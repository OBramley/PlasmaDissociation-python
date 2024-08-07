import os

def parse_bond_info(dissociation_info):
    start_index = dissociation_info.find('[')
    end_index = dissociation_info.find(']')
    if start_index != -1 and end_index != -1:
        bonds_text = dissociation_info[start_index + 1 : end_index]
        pairs = bonds_text.split(', ')
        broken_bonds = [tuple(map(int, pair.strip('()').split(','))) for pair in pairs]
        return broken_bonds
    else:
        return []

def order_and_number_rows(input_file):
    # Read the content of the input file
    with open(input_file, "r") as f:
        lines = f.readlines()

    # Remove newline characters and convert strings to integers
    timesteps = [int(line.strip()) for line in lines]

    # Sort timesteps in ascending order
    sorted_timesteps = sorted(timesteps)

    # Create a new content with row numbers
    new_content = [f"{row_number}\t{timestep}\n" for row_number, timestep in enumerate(sorted_timesteps, start=1)]

    # Write the sorted and numbered content back to the input file
    with open(input_file, "w") as f:
        f.writelines(new_content)


def collate_dissociation_data(folder_path, bonded_array):
    collated_file = 'collated_dissociation.txt'
    
    with open(collated_file, 'w') as output:
        for i in range(1, 101):  # Assuming you have 100 folders t1, t2, ..., t100
            trajectory_folder = f't{i}'
            dissociation_file = os.path.join(trajectory_folder, 'dissociation.out')

            if os.path.exists(dissociation_file):
                with open(dissociation_file, 'r') as input_file:
                    dissociation_info = input_file.read().strip()
                    if dissociation_info:
                        output.write(f"Trajectory{i}:\n")
                        output.write(f"{dissociation_info}\n")
                    
                    else:
                        output.write(f"Trajectory{i}: No dissociation\n")
                    
                    output.write('-' * 40 + '\n')  # Add a separator line
        
        print("Data collation completed.")

def process_trajectories(input_file):
    # Read the input file and extract relevant information
    with open(input_file, "r") as f:
        lines = f.readlines()

    bonded_array = {
        '1-2': 'C-H',
        '1-3': 'C-H',
        '1-4': 'C-H',
        '1-5': 'C-C',
        '5-6': 'C-H',
        '5-7': 'C=C',
        '7-8': 'C-H',
        '7-9': 'C-H'
        # Add more bonded pairs as needed
    }

    broken_bonds_info = {}

    for line in lines:
        if line.startswith("Dissociation detected"):
            parts = line.split(", ")
            timestep = int(float(parts[0].split(" ")[-1]))
            bond_info = parts[1].split(": ")[1]
            bond = bond_info.replace("[", "").replace("]", "").strip()  # Remove spaces
            
            if bond not in broken_bonds_info:
                broken_bonds_info[bond] = []
            
            broken_bonds_info[bond].append(timestep)

    # Create a "Compiled results" directory if it doesn't exist
    compiled_results_dir = "Compiled results"
    if not os.path.exists(compiled_results_dir):
        os.makedirs(compiled_results_dir)

    # Write the extracted information to separate output files for each bond
    for bond, timesteps in broken_bonds_info.items():
        bond_type = bonded_array.get(bond, "Unknown")

        # Save to the bond-type-specific output file
        bond_type_output_file = os.path.join(compiled_results_dir, f"{bond_type}.out")
        with open(bond_type_output_file, "a") as f:
            for timestep in timesteps:
                f.write(f"{timestep}\n")

        # Save to the old output file format
        old_output_file = os.path.join(compiled_results_dir, f"{bond}.out")
        with open(old_output_file, "a") as f:
            for timestep in timesteps:
                f.write(f"{timestep}\n")

def read_data(file_path):
    with open(file_path, "r") as f:
        lines = f.readlines()
    return [line.strip().split("\t") for line in lines]

def combine_files(input_files,output):
    combined_data = []

    # Read data from each file and combine into a single list
    for file_path in input_files:
        data = read_data(file_path)
        combined_data.extend(data)

    # Sort the combined data based on timestep
    sorted_data = sorted(combined_data, key=lambda x: int(x[1]))

    # Update row numbers based on sorted order
    for new_row_number, row in enumerate(sorted_data, start=1):
        row[0] = str(new_row_number)

    # Write the sorted and renumbered data to a new output file

    with open(output, "w") as f:
        for row in sorted_data:
            f.write("\t".join(row) + "\n")

# Provide the path to the folder containing t1, t2, t3, ... folders
folder = os.getcwd()
# collate_dissociation_data(folder, bonded_array)

input_files = ["Compiled results/7-8.out", "Compiled results/7-9.out"]
combine_files(input_files, 'Compiled results/thirdcarbon.out')

# Call the function with the input file name

# process_trajectories('collated_dissociation.txt')

# order_and_number_rows("Compiled results/7-8.out")
