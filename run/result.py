# Result file for the pure python plasma code. Should use a provided bondingarray file to know which bonds to investigate. Then it can look through xyz.all and 
# then make a dissociation.out as well as the C-C.out order files both in the output folder and in the total folder of the a whole set of trajectories. 
import os 
import json
from init import Molecule
import networkx as nx
import shutil

def process_results(rep):
    # 1.read in the bonding array
    bondarr = read_bondarr()
    # # 2. detect dissociation
    bond = detect_dissociation(bondarr,rep)
    # 3. Find fragments.
    frag = fragments(rep)

    return bond, frag

def compile(bonds,frags,bond_groups=None):
    # 4. compile dissociation data
    compile_results(bonds,bond_groups=None)
    # 5. Compile fragment data
    combine_fragments(frags)

def read_bondarr():
    bondarr = {}
    with open('results/bondarr.txt', 'r') as file:
        for line in file:
            line = line.strip()
            if line:  # Skip empty lines
                atoms, bond_type = line.split(':')
                atom1, atom2 = atoms.split('-')
                atom1 = int(atom1)
                atom2 = int(atom2)
                bondarr[(atom1, atom2)] = bond_type
    
    return bondarr

def detect_dissociation(bondarr,rep):

    # Read the input data from a file
    with open('rep-'+str(rep)+'output/xyz.all', 'r') as f:
        lines = f.readlines()

    # Initialize variables to store the current timestep, atom data, and dissociation flag
    timestep = None
    atoms = {}
    dissociated_bonds = set()  # Keep track of dissociated bond pairs
    store=[]
    # Open the output file for writing
    with open('rep-'+str(rep)+'output/dissociation.out', 'w') as output:
        # Iterate through the lines and process the data
        for line in lines:
            parts = line.split()
            if len(parts) == 2 and parts[1].isdigit():
                if atoms:
                    for bonded_pair, bond_type in bondarr.items(): 
                        if bonded_pair not in dissociated_bonds: 
                            i, j = bonded_pair
                            atom1 = atoms[i]        
                            atom2 = atoms[j]
                            distance = ((atom1[1] - atom2[1])**2 + (atom1[2] - atom2[2])**2 + (atom1[3] - atom2[3])**2)**0.5
                            if distance > 4.76:  # Adjust this threshold as needed
                                broken_bond_str = f"{i}-{j}:{bond_type}"  # Include bond type in the output
                                if broken_bond_str not in dissociated_bonds:
                                    output.write(f"Dissociation detected at timestep {timestep}, Broken bond: {broken_bond_str}\n")
                                    dissociated_bonds.add(bonded_pair)
                                    store.append([timestep,bond_type,[i,j]])
                    atoms = {}
                timestep = float(parts[1])
            elif len(parts) == 5 and parts[0].isdigit():
                atom_num = int(parts[0])
                atom_info = [parts[1]] + [float(x) for x in parts[2:]]
                atoms[atom_num] = atom_info

    return store

def compile_results(bonds):
    bonds.sort(key=lambda x:x[0])
    for i in range(len(bonds)):
        bond=bonds[i]
        timestep=bond[0]
        bond_type=bond[1]
        bond_number=f"{bond[2][0]}-{bond[2][1]}"
        bond_number_file = os.path.join("results", f"{bond_number}.out")
        with open(bond_number_file, "a") as f_out:
            f_out.write(f"{timestep}\n")
        
        # Write to bond-type-specific output file
        bond_type_file = os.path.join("results", f"{bond_type}.out")
        with open(bond_type_file, "a") as f_out:
            f_out.write(f"{timestep}\n")

       


def fragments(rep):
    molecule= Molecule.from_json('rep-'+str(rep)+'output/molecule.json')
    coordinates = molecule.coordinates
    fragment_formulas = {} 
    atom_coordinates = [(i, molecule.symbols[i], molecule.coordinates[i, 0], molecule.coordinates[i, 1], molecule.coordinates[i, 2]) for i in range(len(molecule.symbols))]
    atom_distances = {}
    for i in range(len(atom_coordinates)):
        for j in range(i + 1, len(atom_coordinates)):
            dist = distance(atom_coordinates[i], atom_coordinates[j])
            atom_distances[(i, j)] = dist

    G = nx.Graph()
    for (i, j), dist in atom_distances.items():
        if dist < 5:
            G.add_edge(i, j)

    # Initialize a set to keep track of used atoms
    used_atoms = set()
    connected_fragments = list(nx.connected_components(G))
    fragments = []

    for idx, fragment in enumerate(connected_fragments):
        fragment_atoms = []
        for atom_idx in fragment:
            atom_number, atom_name, x, y, z = atom_coordinates[atom_idx]
            fragment_atoms.append((atom_number, atom_name, x, y, z))
            used_atoms.add(atom_idx)
        fragments.append(fragment_atoms)

        # Calculate the molecular formula for the fragment
        formula = {}
        for atom in fragment_atoms:
            atom_name = atom[1]
            if atom_name in formula:
                formula[atom_name] += 1
            else:
                formula[atom_name] = 1

        # Organize formula by CFHO (include 'O' in the order)
        cfho_formula = ""
        if 'C' in formula:
            cfho_formula += f"C{formula['C']}"
        if 'F' in formula:
            cfho_formula += f"F{formula['F']}"
        if 'H' in formula:
            cfho_formula += f"H{formula['H']}"
        if 'O' in formula:
            cfho_formula += f"O{formula['O']}"
        if 'Si' in formula:
            cfho_formula += f"Si{formula['Si']}"

        if cfho_formula in fragment_formulas:
            fragment_formulas[cfho_formula] += 1
        else:
            fragment_formulas[cfho_formula] = 1

    # Process isolated atoms as separate fragments
    isolated_atoms = set(range(len(atom_coordinates))) - used_atoms
    for atom_idx in isolated_atoms:
        atom_number, atom_name, x, y, z = atom_coordinates[atom_idx]
        fragments.append([(atom_number, atom_name, x, y, z)])

        # Calculate the molecular formula for the isolated atom
        formula = {atom_name: 1}

        # Organize formula by CFH
        cfho_formula = ""
        if 'C' in formula:
            cfho_formula += f"C{formula['C']}"
        if 'F' in formula:
            cfho_formula += f"F{formula['F']}"
        if 'H' in formula:
            cfho_formula += f"H{formula['H']}"
        if 'O' in formula:
            cfho_formula += f"O{formula['O']}"
        if 'Si' in formula:
            cfho_formula += f"Si{formula['Si']}"

        if cfho_formula in fragment_formulas:
            fragment_formulas[cfho_formula] += 1
        else:
            fragment_formulas[cfho_formula] = 1
    store=[]
    with open('rep-'+str(rep)+'output/fragments.out', "w") as count_file:
        count_file.write("Fragment Formulas:\n")
        for formula, count in fragment_formulas.items():
            count_file.write(f"{formula}: {count}\n")
            store.append([formula,count])
    
    return store

def distance(point1, point2):
    x1, y1, z1 = point1[2:]
    x2, y2, z2 = point2[2:]
    return ((x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2) ** 0.5

def combine_fragments(frags):
    combo=[]
    for i in range(len(frags)):
        frag=frags[i][0]
        cnt=frags[i][1]
        if frag not in combo:
            combo.append(frags[i])
        else: 
            index=combo.index(frag)
            combo[index[0],1]=combo[index[0],1]+cnt
    
    combo.sort(key=lambda x:x[1], reverse=True)
    with open("results/fragments.out", 'w') as file:
        file.write("Fragment Formulas (Most to Least Common):\n")
        for i in range(len(combo)):
            file.write(f"{combo[i,0]}: {combo[i,1]}\n")


def combine_fragment_counts(file1, file2, output_file):
    # Read fragment counts from the first file
    fragment_formulas1 = read_fragment_counts(file1)

    # Read fragment counts from the second file
    fragment_formulas2 = read_fragment_counts(file2)

    # Combine fragment counts from both files
    combined_fragment_formulas = combine_counts(fragment_formulas1, fragment_formulas2)

    # Write the combined fragment counts to the output file
    write_fragment_counts(output_file, combined_fragment_formulas)

def read_fragment_counts(file_path):
    fragment_formulas = {}
    with open(file_path, 'r') as file:
        lines = file.readlines()

        # Skip the first line (title or header)
        lines = lines[1:]

        for line in lines:
            parts = line.strip().split(':')
            if len(parts) == 2:
                formula, count = parts[0], int(parts[1])
                fragment_formulas[formula] = count
    return fragment_formulas

def combine_counts(counts1, counts2):
    combined_counts = counts1.copy()
    for formula, count in counts2.items():
        if formula in combined_counts:
            combined_counts[formula] += count
        else:
            combined_counts[formula] = count
    return combined_counts

def write_fragment_counts(file_path, fragment_formulas):
    # Sort fragment formulas based on counts in descending order
    sorted_formulas = sorted(fragment_formulas.items(), key=lambda x: x[1], reverse=True)

    with open(file_path, 'w') as file:
        file.write("Fragment Formulas (Most to Least Common):\n")
        for formula, count in sorted_formulas:
            file.write(f"{formula}: {count}\n")

if __name__ == "__main__":
    with open('inputs.json') as f:
        inputs=json.load(f)
    reps=inputs["run"]["repeats"]
    bonds=[]
    frags=[]
    for i in range(1,reps+1):
        bond, frag = process_results(i)
        bonds.append(bond)
        frags.append(frag)

    
    compile(bonds,frags)
