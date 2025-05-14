import os
import sys

class BaderParser:
    """
    A class to parse and combine Bader analysis data from a given directory.
    
    When an instance is created with a directory, it locates three required files:
      - ACF.dat: Contains integrated atomic data (e.g., coordinates, charge, 
                minimum distance, and atomic volume).
      - AVF.dat: Provides a mapping from each atom number to the corresponding 
                BCF record index.
      - BCF.dat: Contains detailed data for each individual Bader region including
                its coordinates, charge, and the distance value (from the file’s
                distance column).
                
    The final combined dataset is arranged by atom number and includes:
      - ACF data: (coordinates, charge, min distance, atomic volume)
      - BCF data: (coordinates, charge, distance) taken directly from the BCF file.
    """

    def __init__(self, directory):
        if not os.path.isdir(directory):
            raise ValueError(f"Error: {directory} is not a valid directory.")
        
        self.directory = directory
        self.acf_file = self.find_file("acf")
        self.avf_file = self.find_file("avf")
        self.bcf_file = self.find_file("bcf")
        
        if not self.acf_file or not self.avf_file or not self.bcf_file:
            raise FileNotFoundError("Could not locate one or more required files (ACF.dat, AVF.dat, BCF.dat) in the directory.")
        
        # These attributes will hold parsed data.
        self.acf_data = None
        self.avf_map = None  # Mapping: a list of dicts with keys 'atom' and 'bcf_index'
        self.bcf_data = None

    def find_file(self, keyword):
        """
        Searches the directory for a file whose name contains the given keyword (case insensitive)
        and ends with '.dat'.
        
        Returns:
            The full path to the first matching file, or None if no match is found.
        """
        keyword = keyword.lower()
        for fname in os.listdir(self.directory):
            if fname.lower().endswith(".dat") and keyword in fname.lower():
                return os.path.join(self.directory, fname)
        return None

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
            # Skip header lines (such as those starting with '#' or lines of dashes).
            if line.startswith("#") or set(line) == {"-"}:
                continue
            if line.startswith("VACUUM"):  # Stop parsing at footer
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
                acf_volume = float(parts[6])
            except ValueError:
                continue

            records.append({
                "atom": atom,
                "acf_coords": (x, y, z),
                "acf_charge": charge,
                "min_dist": min_dist,
                "acf_volume": acf_volume
            })
        self.acf_data = records
        return records

    def parse_avf(self):
        """
        Parse the AVF file, which provides a mapping from atom numbers to BCF record indices.
        
        Expected format:
          - A header line with column titles (often "Atom   Volume(s)") and a dashed separator.
          - Data lines where the first token is the atom number and the last token (despite the header)
            is treated as the BCF record index corresponding to that atom.
        
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
                # Here the header is "Volume(s)" but we use the last token as the BCF index.
                bcf_index = int(tokens[-1])
            except ValueError:
                continue

            mapping.append({"atom": atom, "bcf_index": bcf_index})
        self.avf_map = mapping
        return mapping

    def parse_bcf(self):
        """
        Parse the BCF file which contains detailed data for each individual Bader region.
        
        Expected columns:
          Index, X, Y, Z, CHARGE, ATOM, DISTANCE
        
        Returns:
            A list of dictionaries, each containing:
              {
                'index': int,          # This is the unique identifier for the BCF record.
                'bcf_coords': (x, y, z),
                'bcf_charge': float,
                'atom': int,           # The atom field in the BCF file (may not be the same as the ACF atom number).
                'bcf_distance': float  # Taken directly from the BCF file's distance column.
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
                    # Although the BCF file has an "ATOM" column, for our mapping we rely on the record's "index".
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
        Combine the parsed ACF data with the corresponding BCF data by using the AVF mapping.
        
        For each atom:
          - The ACF record provides the atom number, its coordinates, integrated charge,
            minimum distance, and atomic volume.
          - The AVF mapping tells us which BCF record (by its index) maps to this atom.
          - The corresponding BCF record provides its own coordinates, charge, and the 
            value from its distance column.
        
        The final combined record (keyed by the atom number) contains:
          - 'acf_coords', 'acf_charge', 'min_dist', 'acf_volume' (from the ACF file)
          - 'bcf_coords', 'bcf_charge', 'bcf_distance' (from the BCF record as indicated by the AVF mapping)
        
        Returns:
            A dictionary keyed by atom number.
        """
        # Ensure that all data are parsed.
        if self.acf_data is None:
            self.parse_acf()
        if self.avf_map is None:
            self.parse_avf()
        if self.bcf_data is None:
            self.parse_bcf()

        combined = {}

        # Create a lookup for BCF records keyed by their "index".
        bcf_lookup = {rec["index"]: rec for rec in self.bcf_data}

        # First, fill in the combined record with ACF data.
        for rec in self.acf_data:
            atom = rec["atom"]
            combined[atom] = {
                "acf_coords": rec["acf_coords"],
                "acf_charge": rec["acf_charge"],
                "min_dist": rec["min_dist"],
                "acf_volume": rec["acf_volume"],
                "bcf": None  # This will hold the corresponding BCF record details if mapping is found.
            }

        # Use the AVF mapping to attach the corresponding BCF record.
        for mapping in self.avf_map:
            atom = mapping["atom"]
            bcf_index = mapping["bcf_index"]
            # Look up the BCF record directly using its "index" (which is provided by the AVF file).
            bcf_rec = bcf_lookup.get(bcf_index)
            if bcf_rec is not None:
                # We include the BCF coordinates, charge, and the distance value from the BCF file.
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
        parser = BaderParser(directory)
    except Exception as e:
        print(e)
        sys.exit(1)

    print("Found files:")
    print("ACF file:", parser.acf_file)
    print("AVF file:", parser.avf_file)
    print("BCF file:", parser.bcf_file)

    print("\nParsing data...")
    parser.parse_acf()
    parser.parse_avf()
    parser.parse_bcf()

    combined = parser.combine_data()

    print("\nCombined Data:")
    for atom in sorted(combined.keys()):
        data = combined[atom]
        print(f"Atom {atom:3d}:")
        print(f"  ACF:  Coordinates {data.get('acf_coords', 'N/A')}, Charge = {data.get('acf_charge', 'N/A')},")
        print(f"        Min Dist = {data.get('min_dist', 'N/A')}, Atomic Vol = {data.get('acf_volume', 'N/A')}")
        if data.get("bcf") is not None:
            bcf = data["bcf"]
            print(f"  BCF:  Coordinates {bcf.get('bcf_coords', 'N/A')}, Charge = {bcf.get('bcf_charge', 'N/A')},")
            print(f"        Distance = {bcf.get('bcf_distance', 'N/A')}")
        else:
            print("  BCF:  No matching BCF record found via AVF mapping.")
        print()

if __name__ == '__main__':
    main()