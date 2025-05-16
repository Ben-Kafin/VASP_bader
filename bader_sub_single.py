from bader_net import BaderParser  # Assumes "bader_net.py" defines BaderParser

class BaderSubtractorSingle:
    """
    This class uses BaderParser to parse the Bader output from a full system directory
    and one component directory. It assumes that the atom numbering is consistent between
    the full system and the component, so matching is done solely by atom number.
    
    The component values (for example, ACF charge, volume, net charge) are subtracted from
    the full system values for each atom.
    """
    
    def __init__(self, full_system_dir, component_dir):
        self.full_system_dir = full_system_dir
        # If a list or tuple is provided, take only the first element.
        if isinstance(component_dir, (list, tuple)):
            self.component_dir = component_dir[0]
        else:
            self.component_dir = component_dir
            
    def get_combined_data(self, directory):
        """
        Creates a BaderParser instance for the given directory, calls its parsing methods,
        and returns the combined data as a dictionary keyed by atom number.
        """
        parser = BaderParser(directory)
        parser.parse_acf()
        parser.parse_avf()
        parser.parse_bcf()
        return parser.combine_data()
    
    def subtract_component(self):
        """
        Subtracts the component values from the full system values using atom number matching.
        For each atom in the full system data, if there is a record with the same atom number in
        the component data, its ACf charge, volume, and net charge (if available) are used in the
        subtraction. If not, zero is assumed.
        
        Returns a dictionary keyed by atom number with the following keys:
          - full_acf_charge, component_acf_charge, acf_charge_diff
          - full_acf_volume, component_acf_volume, acf_volume_diff
          - full_net_charge, component_net_charge, net_charge_diff
          - acf_coords (the full system coordinates) and atom_type.
        """
        # Get combined data from the full system.
        full_data = self.get_combined_data(self.full_system_dir)
        # Get combined data from the component.
        comp_data = self.get_combined_data(self.component_dir)
        
        differences = {}
        for atom, full_rec in full_data.items():
            # Use atom number for matching:
            comp_rec = comp_data.get(atom)
            
            comp_acf_charge = comp_rec.get("acf_charge", 0.0) if comp_rec is not None else 0.0
            comp_acf_volume = comp_rec.get("acf_volume", 0.0) if comp_rec is not None else 0.0
            comp_net_charge = comp_rec.get("net_charge", 0.0) if (comp_rec is not None and "net_charge" in comp_rec) else None
            
            diff_acf_charge = full_rec.get("acf_charge", 0.0) - comp_acf_charge
            diff_acf_volume = full_rec.get("acf_volume", 0.0) - comp_acf_volume
            diff_net_charge = None
            if full_rec.get("net_charge") is not None and comp_net_charge is not None:
                diff_net_charge = full_rec.get("net_charge") - comp_net_charge
            
            differences[atom] = {
                "full_acf_charge": full_rec.get("acf_charge"),
                "component_acf_charge": comp_acf_charge,
                "acf_charge_diff": diff_acf_charge,
                "full_acf_volume": full_rec.get("acf_volume"),
                "component_acf_volume": comp_acf_volume,
                "acf_volume_diff": diff_acf_volume,
                "full_net_charge": full_rec.get("net_charge"),
                "component_net_charge": comp_net_charge,
                "net_charge_diff": diff_net_charge,
                "acf_coords": full_rec.get("acf_coords"),
                "atom_type": full_rec.get("atom_type", "N/A")
            }
            
        return differences

def write_output_file(differences, out_filename="bader_output.txt"):
    """
    Writes the differences to an output file in tab-delimited format.
    
    The output file includes a header and one line per atom with the following columns:
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
            def fmt(val):
                return f"{val:.4f}" if isinstance(val, (int, float)) else str(val)
            full_acf_charge = fmt(diff.get("full_acf_charge", "N/A"))
            comp_acf_charge = fmt(diff.get("component_acf_charge", "N/A"))
            acf_charge_diff = fmt(diff.get("acf_charge_diff", "N/A"))
            full_acf_volume = fmt(diff.get("full_acf_volume", "N/A"))
            comp_acf_volume = fmt(diff.get("component_acf_volume", "N/A"))
            acf_volume_diff = fmt(diff.get("acf_volume_diff", "N/A"))
            full_net_charge = fmt(diff.get("full_net_charge", "N/A"))
            comp_net_charge = fmt(diff.get("component_net_charge", "N/A"))
            net_charge_diff = fmt(diff.get("net_charge_diff", "N/A"))
            atom_type = diff.get("atom_type", "N/A")
            
            line = (
                f"{atom}\t{atom_type}\t{x}\t{y}\t{z}\t"
                f"{full_acf_charge}\t{comp_acf_charge}\t{acf_charge_diff}\t"
                f"{full_acf_volume}\t{comp_acf_volume}\t{acf_volume_diff}\t"
                f"{full_net_charge}\t{comp_net_charge}\t{net_charge_diff}"
            )
            f.write(line + "\n")

if __name__ == '__main__':
    # Hard-coded directories (edit these paths as needed)
    full_system_dir = "C:/Users/Benjamin Kafin/Documents/VASP/NHC/IPR/lone/NHC/NHC_iPr/4layers/freegold1/freegold2/kpoints551/dipole_correction/"
    component_dir = "C:/Users/Benjamin Kafin/Documents/VASP/NHC/IPR/lone/NHC/NHC_iPr/4layers/pristine_gold/kpoints551/"
    
    subtractor = BaderSubtractorSingle(full_system_dir, component_dir)
    differences = subtractor.subtract_component()
    
    output_filename = "C:/Users/Benjamin Kafin/Documents/VASP/NHC/IPR/lone/NHC/NHC_iPr/4layers/freegold1/freegold2/kpoints551/dipole_correction/bader_output.txt"
    write_output_file(differences, out_filename=output_filename)
    
    print(f"Bader Data Differences have been written to '{output_filename}'")