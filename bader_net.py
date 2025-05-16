# -*- coding: utf-8 -*-
"""
Created on Wed May 14 15:49:13 2025

@author: Benjamin Kafin
"""

import os
import sys

class BaderParser:
    """
    A class to parse and combine Bader analysis data from a given directory.
    
    The class looks for the following files:
      - ACF.dat: Contains integrated atomic data (e.g., coordinates, charge, minimum distance, and atomic volume).
      - AVF.dat: Provides a mapping from atom numbers to BCF record indices.
      - BCF.dat: Contains detailed data for each individual Bader region (coordinates, charge, and distance).
      - POSCAR (optional): Contains the structural information including element symbols and the number of atoms per type.
    
    The POSCAR is parsed (if available) to automatically build a mapping from atom numbers to element types.
    Then, using an expected-valence dictionary (which you may extend), the net electron (valence) count is computed by subtracting the expected number from the ACF charge.
    
    The final combined record for each atom will include:
      - ACF data: coordinates, integrated charge, minimum distance, and atomic volume.
      - BCF data: coordinates, charge, and the distance (taken from the BCF file).
      - Atom type and net charge (derived from the POSCAR mapping and expected valence).
    """

    # Define expected valence electrons for select elements.
    EXPECTED_VALENCE = {"Au": 11.0, "N": 5.0, "C": 4.0, "H": 1.0}

    def __init__(self, directory, atom_types=None):
        """
        Initialize with a directory and, optionally, a manual mapping from atom number to element type.
        If a POSCAR file is found in the directory, it will be parsed to build this mapping automatically.
        """
        if not os.path.isdir(directory):
            raise ValueError(f"Error: {directory} is not a valid directory.")
        
        self.directory = directory
        self.acf_file = self.find_file("acf")
        self.avf_file = self.find_file("avf")
        self.bcf_file = self.find_file("bcf")
        
        if not self.acf_file or not self.avf_file or not self.bcf_file:
            raise FileNotFoundError("Could not locate one or more required files (ACF.dat, AVF.dat, BCF.dat) in the directory.")
        
        # These attributes will hold parsed Bader data.
        self.acf_data = None
        self.avf_map = None    # List of mappings: [{'atom': int, 'bcf_index': int}, ...]
        self.bcf_data = None
        
        # If a POSCAR is found, parse it; otherwise use the provided atom_types mapping.
        self.atom_types = atom_types if atom_types is not None else {}
        poscar_mapping = self.parse_poscar()
        if poscar_mapping:
            self.atom_types = poscar_mapping

    def find_file(self, keyword):
        """
        Searches the directory for a file whose name contains the given keyword (case insensitive)
        and ends with '.dat'.
        
        Returns:
            The full path to the first matching file, or None if not found.
        """
        keyword = keyword.lower()
        for fname in os.listdir(self.directory):
            if fname.lower().endswith(".dat") and keyword in fname.lower():
                return os.path.join(self.directory, fname)
        return None

    def parse_poscar(self):
        """
        Parse the POSCAR file (if found in the directory) to build a mapping from atom number to element symbol.
        
        POSCAR format (VASP 5+):
          - Line 1: Comment/title.
          - Line 2: Scaling factor.
          - Lines 3-5: Lattice vectors.
          - Line 6: Element symbols (space separated).
          - Line 7: Number of atoms for each element (space separated).
          - (Then optional lines and coordinates.)
        
        Returns:
            A dictionary mapping atom numbers (starting at 1) to element symbols.
            For example: {1: "Au", 2: "N", 3: "N", 4: "C", ...}
            If the POSCAR file is not found or not in the expected format, returns {}.
        """
        poscar_path = None
        for fname in os.listdir(self.directory):
            if fname.lower() == "poscar":
                poscar_path = os.path.join(self.directory, fname)
                break
        if not poscar_path:
            return {}
        
        try:
            with open(poscar_path, 'r') as f:
                lines = f.readlines()
            if len(lines) < 7:
                return {}
            # Assume line 6 contains element symbols and line 7 contains counts.
            symbols = lines[5].strip().split()
            counts = [int(x) for x in lines[6].split()]
            mapping = {}
            atom_number = 1
            for sym, count in zip(symbols, counts):
                for i in range(count):
                    mapping[atom_number] = sym
                    atom_number += 1
            return mapping
        except Exception:
            # If something goes wrong, return an empty mapping.
            return {}

    def parse_acf(self):
        """
        Parse the ACF file.
        
        Expected columns (in order):
            Atom, X, Y, Z, CHARGE, MIN DIST, ATOMIC VOL
        
        Returns:
            A list of dictionaries, each containing:
              {
                 'atom': int,
                 'acf_coords': (x, y, z),
                 'acf_charge': float,
                 'min_dist': float,
                 'acf_volume': float
              }
        """
        records = []
        with open(self.acf_file, 'r') as f:
            lines = f.readlines()
        for line in lines:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#") or set(line) == {"-"}:
                continue
            if line.startswith("VACUUM"):
                break
            parts = line.split()
            if len(parts) < 7:
                continue
            try:
                atom = int(parts[0])
                x = float(parts[1])
                y = float(parts[2])
                z = float(parts[3])
                charge = float(parts[4])
                min_dist = float(parts[5])
                acf_vol = float(parts[6])
            except ValueError:
                continue
            records.append({
                "atom": atom,
                "acf_coords": (x, y, z),
                "acf_charge": charge,
                "min_dist": min_dist,
                "acf_volume": acf_vol
            })
        self.acf_data = records
        return records

    def parse_avf(self):
        """
        Parse the AVF file, which contains a mapping from atom numbers to BCF record indices.
        
        Expected format:
          - A header line (e.g., containing "Atom" and "Volume(s)") and a dashed separator.
          - Data lines where the first token is the atom number and the last token is interpreted
            as the BCF record index for that atom.
        
        Returns:
            A list of dictionaries, each containing:
              {
                 'atom': int,
                 'bcf_index': int
              }
        """
        mapping = []
        with open(self.avf_file, 'r') as f:
            lines = f.readlines()
        for line in lines:
            line = line.strip()
            if not line:
                continue
            if "Atom" in line or "Volume" in line or line.startswith("-"):
                continue
            tokens = line.split()
            if len(tokens) < 2:
                continue
            try:
                atom = int(tokens[0])
                bcf_index = int(tokens[-1])
            except ValueError:
                continue
            mapping.append({"atom": atom, "bcf_index": bcf_index})
        self.avf_map = mapping
        return mapping

    def parse_bcf(self):
        """
        Parse the BCF file, which contains detailed information for each individual Bader region.
        
        Expected columns:
            Index, X, Y, Z, CHARGE, ATOM, DISTANCE
        
        Returns:
            A list of dictionaries, each containing:
              {
                 'index': int,
                 'bcf_coords': (x, y, z),
                 'bcf_charge': float,
                 'atom': int,           # Not used for mapping; we use the record's index.
                 'bcf_distance': float  # Directly from the distance column.
              }
        """
        records = []
        with open(self.bcf_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith("#"):
                    continue
                if set(line) == {"-"}:
                    continue
                tokens = line.split()
                if len(tokens) < 7:
                    continue
                try:
                    rec_index = int(tokens[0])
                    x = float(tokens[1])
                    y = float(tokens[2])
                    z = float(tokens[3])
                    bcf_charge = float(tokens[4])
                    # The BCF file has an "ATOM" column, but for our mapping we use the unique record index.
                    bcf_atom = int(tokens[5])
                    bcf_distance = float(tokens[6])
                except ValueError:
                    continue
                records.append({
                    "index": rec_index,
                    "bcf_coords": (x, y, z),
                    "bcf_charge": bcf_charge,
                    "atom": bcf_atom,
                    "bcf_distance": bcf_distance
                })
        self.bcf_data = records
        return records

    def combine_data(self):
        """
        Combine the parsed ACF data with the corresponding BCF data using the AVF mapping.
        
        For each atom number:
          - From ACF: acf_coords, acf_charge, min_dist, acf_volume.
          - From AVF: We obtain the corresponding BCF record index.
          - From BCF: bcf_coords, bcf_charge, bcf_distance.
          
        Additionally, if an element type is determined (from the POSCAR or provided externally),
        the net charge is computed as:
             net_charge = acf_charge - expected_valence_electrons
        using the EXPECTED_VALENCE dictionary.
        
        Returns:
            A dictionary keyed by atom number.
        """
        if self.acf_data is None:
            self.parse_acf()
        if self.avf_map is None:
            self.parse_avf()
        if self.bcf_data is None:
            self.parse_bcf()
        
        combined = {}
        # Build a BCF lookup keyed by the record's "index".
        bcf_lookup = {rec["index"]: rec for rec in self.bcf_data}
        
        # Start by incorporating ACF data.
        for rec in self.acf_data:
            atom = rec["atom"]
            combined[atom] = {
                "acf_coords": rec["acf_coords"],
                "acf_charge": rec["acf_charge"],
                "min_dist": rec["min_dist"],
                "acf_volume": rec["acf_volume"],
                "atom_type": self.atom_types.get(atom, "N/A"),
                "net_charge": None,  # Will be computed if expected valence is known.
                "bcf": None        # To be filled in using the mapping.
            }
            elem = self.atom_types.get(atom)
            if elem and elem in self.EXPECTED_VALENCE:
                expected = self.EXPECTED_VALENCE[elem]
                combined[atom]["net_charge"] = rec["acf_charge"] - expected
        
        # Use the AVF mapping to add the corresponding BCF record.
        for mapping in self.avf_map:
            atom = mapping["atom"]
            bcf_index = mapping["bcf_index"]
            bcf_rec = bcf_lookup.get(bcf_index)
            if bcf_rec is not None:
                combined_bcf = {
                    "bcf_coords": bcf_rec["bcf_coords"],
                    "bcf_charge": bcf_rec["bcf_charge"],
                    "bcf_distance": bcf_rec["bcf_distance"]
                }
                if atom in combined:
                    combined[atom]["bcf"] = combined_bcf
                else:
                    combined[atom] = {"bcf": combined_bcf}
        return combined

def main():


    directory = 'C:/Users/Benjamin Kafin/Documents/VASP/NHC/IPR/lone/NHC/NHC_iPr/4layers/freegold1/freegold2/kpoints551/dipole_correction/'

    try:
        # We do not necessarily provide a manual mapping since the POSCAR (if present) will supply it.
        parser = BaderParser(directory)
    except Exception as e:
        print(e)
        sys.exit(1)
    
    print("Found files:")
    print("ACF file:", parser.acf_file)
    print("AVF file:", parser.avf_file)
    print("BCF file:", parser.bcf_file)
    
    # Optionally, show the POSCAR-derived mapping.
    if parser.atom_types:
        print("\nAtom type mapping (from POSCAR):")
        for atom in sorted(parser.atom_types.keys()):
            print(f"  Atom {atom:3d}: {parser.atom_types[atom]}")
    
    print("\nParsing data...")
    parser.parse_acf()
    parser.parse_avf()
    parser.parse_bcf()
    
    combined = parser.combine_data()
    
    print("\nCombined Data:")
    for atom in sorted(combined.keys()):
        data = combined[atom]
        print(f"Atom {atom:3d}:")
        print(f"  Atom Type: {data.get('atom_type', 'N/A')}")
        print(f"  ACF: Coordinates {data.get('acf_coords', 'N/A')}, Charge = {data.get('acf_charge', 'N/A')},")
        print(f"       Min Dist = {data.get('min_dist', 'N/A')}, Atomic Volume = {data.get('acf_volume', 'N/A')}")
        if data.get("net_charge") is not None:
            print(f"       Net Charge (ACF charge - expected valence) = {data['net_charge']:.6f}")
        if data.get("bcf") is not None:
            bcf = data["bcf"]
            print(f"  BCF: Coordinates {bcf.get('bcf_coords', 'N/A')}, Charge = {bcf.get('bcf_charge', 'N/A')},")
            print(f"       Distance from BCF file = {bcf.get('bcf_distance', 'N/A')}")
        else:
            print("  BCF: No matching BCF record found via AVF mapping.")
        print()

if __name__ == '__main__':
    main()