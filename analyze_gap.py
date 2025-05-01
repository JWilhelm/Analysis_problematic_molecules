import os

# Define basis sets
basis_sets = [
    "aug-cc-pVDZ",
    "aug-cc-pVTZ",
    "aug-cc-pVQZ",
    "aug-cc-pV5Z",
    "aug-SZVSR",
    "aug-SZV",
    "aug-DZVP",
    "aug-TZVP"
]

# Define atom names up to Krypton
atom_names = [
    'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
    'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar',
    'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni',
    'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr'
]

# Prepare result containers
GW_HOMO_LUMO_gap = {}
vir_4th_value = {}
atom_counts_all = {}

# Get all top-level molecule directories
molecules = [d for d in os.listdir('.') if os.path.isdir(d)]
molecules.sort()

for molecule in molecules:
    mol_path = os.path.join('.', molecule)
    if not os.path.isdir(mol_path):
        continue

    GW_HOMO_LUMO_gap[molecule] = {}
    vir_4th_value[molecule] = {}
    atom_counts_all[molecule] = {atom: 0 for atom in atom_names}

    for basis in basis_sets:
        file_path = os.path.join(mol_path, molecule, basis, 'GW+BSE', 'cp2k.out')
        if not os.path.isfile(file_path):
            file_path = os.path.join(mol_path, molecule, basis, 'GW+BSE', 'RILIM', 'cp2k.out')
            if not os.path.isfile(file_path):
                continue

        with open(file_path, 'r') as f:
            reading_atoms = False
            for line in f:
                if 'SCF PARAMETERS' in line:
                    break
                if 'ATOMIC COORDINATES' in line:
                    reading_atoms = True
                    continue
                if reading_atoms:
                    tokens = line.strip().split()
                    if len(tokens) >= 3 and tokens[2] in atom_names:
                        atom_counts_all[molecule][tokens[2]] += 1

        found_gap = False
        found_vir = False
        try:
            with open(file_path, 'r') as f:
                for line in f:
                    if not found_gap and "G0W0 HOMO-LUMO gap (eV)" in line:
                        tokens = line.strip().split()
                        if len(tokens) >= 5:
                            try:
                                gap_value = float(tokens[4])
                                GW_HOMO_LUMO_gap[molecule][basis] = gap_value
                                found_gap = True
                            except ValueError:
                                pass
                    if not found_vir and "( vir )" in line:
                        tokens = line.strip().split()
                        if len(tokens) >= 4:
                            try:
                                vir_value = float(tokens[4])
                                vir_4th_value[molecule][basis] = vir_value
                                found_vir = True
                            except ValueError:
                                pass
                    if found_gap and found_vir:
                        break
        except Exception as e:
            print(f"Error reading {file_path}: {e}")

# Print GW_HOMO_LUMO_gap and vir_4th_value
for molecule in molecules:
    counts = atom_counts_all.get(molecule, {})
    formula_parts = [f"{atom}{count // 8}" for atom, count in counts.items() if count > 0]
    formula = ''.join(formula_parts) if formula_parts else 'N/A'

    print(f"Molecule: {molecule}  Formula: {formula}")
    GW_data = GW_HOMO_LUMO_gap.get(molecule, {})
    vir_data = vir_4th_value.get(molecule, {})

    for basis in basis_sets:
        gw_value = GW_data.get(basis, 'N/A')
        vir_value = vir_data.get(basis, 'N/A')
        if isinstance(gw_value, float):
            gw_str = f"{gw_value:7.3f}"
        else:
            gw_str = "   N/A"
        if isinstance(vir_value, float):
            vir_str = f"{vir_value:7.3f}"
        else:
            vir_str = "   N/A"
        print(f" {basis:<12}: GW gap {gw_str}   DFT LUMO {vir_str}")
    print()

