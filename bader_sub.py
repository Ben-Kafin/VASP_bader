#!/usr/bin/env python
import os
from bader_net import BaderParser  # Assumes bader_net.py defines BaderParser

class BaderSubtractor:
    """
    This class uses BaderParser to parse Bader output (ACF/AVF/BCF and POSCAR mapping)
    from a full system directory and one or more component directories.
    
    It then—by matching atoms based on their ACF coordinates (within a tolerance)—subtracts
    the component values from the full system values.
    """
    
    def __init__(self, full_system_dir, component_dirs, tol=1e-4):
        self.full_system_dir = full_system_dir
        # If a single directory is given as a string, convert it into a one-element list.
        if isinstance(component_dirs, str):
            self.component_dirs = [component_dirs]
        else:
            self.component_dirs = component_dirs
        self.tol = tol

    def get_combined_data(self, directory):
        """
        Create a BaderParser instance for the given directory, call its parsing methods, 
        and return the combined data.
        """
        parser = BaderParser(directory)
        parser.parse_acf()
        parser.parse_avf()
        parser.parse_bcf()
        return parser.combine_data()

    def coords_equal(self, coord1, coord2):
        """Compare two coordinate tuples element-wise within a tolerance."""
        return all(abs(a - b) < self.tol for a, b in zip(coord1, coord2))
    
    def match_by_coord(self, target_coord, dataset):
        """
        In a given dataset (a dict keyed by atom numbers containing an 'acf_coords' field),
        return the record whose coordinates match target_coord (within tolerance).
        Returns None if no record is found.
        """
        for atom, record in dataset.items():
            coord = record.get("acf_coords")
            if coord and self.coords_equal(coord, target_coord):
                return record
        return None

    def subtract_components(self):
        """
        For each atom in the full system data, look for matching atoms (by ACF coordinates)
        in each component's combined data, sum the component values for additive quantities,
        and subtract these sums from the full system values.
        
        Returns a dictionary keyed by the atom number with the differences.
        """
        # Get full system combined data.
        full_data = self.get_combined_data(self.full_system_dir)
        # Get component combined data for each component directory.
        comp_data_list = [self.get_combined_data(comp_dir) for comp_dir in self.component_dirs]
        
        differences = {}
        for atom, full_rec in full_data.items():
            target_coord = full_rec.get("acf_coords")
            # Initialize sums.
            sum_acf_charge = 0.0
            sum_acf_volume = 0.0
            sum_net_charge = 0.0
            comp_count_with_net = 0
            
            # Loop over each component's data and look for a matching atom.
            for comp_data in comp_data_list:
                comp_rec = self.match_by_coord(target_coord, comp_data)
                if comp_rec is not None:
                    sum_acf_charge += comp_rec.get("acf_charge", 0.0)
                    sum_acf_volume += comp_rec.get("acf_volume", 0.0)
                    if comp_rec.get("net_charge") is not None:
                        sum_net_charge += comp_rec.get("net_charge")
                        comp_count_with_net += 1
            
            diff_acf_charge = full_rec.get("acf_charge", 0.0) - sum_acf_charge
            diff_acf_volume = full_rec.get("acf_volume", 0.0) - sum_acf_volume
            diff_net_charge = None
            if full_rec.get("net_charge") is not None and comp_count_with_net > 0:
                diff_net_charge = full_rec.get("net_charge") - sum_net_charge
            
            differences[atom] = {
                "full_acf_charge": full_rec.get("acf_charge"),
                "component_acf_charge": sum_acf_charge,
                "acf_charge_diff": diff_acf_charge,
                "full_acf_volume": full_rec.get("acf_volume"),
                "component_acf_volume": sum_acf_volume,
                "acf_volume_diff": diff_acf_volume,
                "full_net_charge": full_rec.get("net_charge"),
                "component_net_charge": sum_net_charge if comp_count_with_net > 0 else None,
                "net_charge_diff": diff_net_charge,
                "acf_coords": target_coord,
                "atom_type": full_rec.get("atom_type", "N/A")
            }
        return differences

def write_output_file(differences, out_filename="bader_output.txt"):
    """
    Write the differences to an output file in tab-delimited format.
    
    The file includes a header and one line per atom, with columns:
      Atom, Atom_Type, X, Y, Z, Full_ACF_Charge, Component_ACF_Charge, ACF_Charge_Diff,
      Full_ACF_Volume, Component_ACF_Volume, ACF_Volume_Diff, Full_Net_Charge,
      Component_Net_Charge, Net_Charge_Diff
    """
    header = (
        "Atom\tAtom_Type\tX\tY\tZ\tFull_ACF_Charge\tComponent_ACF_Charge\tACF_Charge_Diff\t"
        "Full_ACF_Volume\tComponent_ACF_Volume\tACF_Volume_Diff\tFull_Net_Charge\t"
        "Component_Net_Charge\tNet_Charge_Diff"
    )
    with open(out_filename, "w") as f:
        f.write(header + "\n")
        for atom in sorted(differences.keys()):
            diff = differences[atom]
            # Unpack coordinates.
            coords = diff.get("acf_coords", (None, None, None))
            x, y, z = coords if all(v is not None for v in coords) else ("N/A", "N/A", "N/A")
            # Format numbers; if value is None, we output "N/A".
            def fmt(val):
                return f"{val:.4f}" if isinstance(val, (int, float)) else str(val)
            full_acf_charge = fmt(diff.get("full_acf_charge", "N/A"))
            component_acf_charge = fmt(diff.get("component_acf_charge", "N/A"))
            acf_charge_diff = fmt(diff.get("acf_charge_diff", "N/A"))
            full_acf_volume = fmt(diff.get("full_acf_volume", "N/A"))
            component_acf_volume = fmt(diff.get("component_acf_volume", "N/A"))
            acf_volume_diff = fmt(diff.get("acf_volume_diff", "N/A"))
            full_net_charge = fmt(diff.get("full_net_charge", "N/A"))
            component_net_charge = fmt(diff.get("component_net_charge", "N/A"))
            net_charge_diff = fmt(diff.get("net_charge_diff", "N/A"))
            atom_type = diff.get("atom_type", "N/A")
            
            line = (
                f"{atom}\t{atom_type}\t{x}\t{y}\t{z}\t"
                f"{full_acf_charge}\t{component_acf_charge}\t{acf_charge_diff}\t"
                f"{full_acf_volume}\t{component_acf_volume}\t{acf_volume_diff}\t"
                f"{full_net_charge}\t{component_net_charge}\t{net_charge_diff}"
            )
            f.write(line + "\n")

if __name__ == '__main__':
    # Hard-coded directories (edit these paths as needed)
    full_system_dir = "C:/Users/Benjamin Kafin/Documents/VASP/NHC/IPR/lone/NHC/NHC_iPr/4layers/freegold1/freegold2/kpoints551/dipole_correction/"
    component_dirs = [
        "C:/Users/Benjamin Kafin/Documents/VASP/NHC/IPR/lone/NHC/NHC_iPr/4layers/freegold1/freegold2/kpoints551/NHC/dipole_correction/",
        "C:/Users/Benjamin Kafin/Documents/VASP/NHC/IPR/lone/NHC/NHC_iPr/4layers/pristine_gold/kpoints551/"
    ]
    
    subtractor = BaderSubtractor(full_system_dir, component_dirs)
    differences = subtractor.subtract_components()

    # Write the combined differences to a tab-delimited output file.
    output_filename = "C:/Users/Benjamin Kafin/Documents/VASP/NHC/IPR/lone/NHC/NHC_iPr/4layers/freegold1/freegold2/kpoints551/dipole_correction/bader_output_1.txt"
    write_output_file(differences, out_filename=output_filename)
    
    print(f"Bader Data Differences have been written to '{output_filename}'")