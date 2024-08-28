import matplotlib.pyplot as plt
import matplotlib as mpl
import json

def plot_data(filename, label, no_bonds, color, total_bonds,ax):
    # Read data from file
    with open(filename, 'r') as file:
        lines = file.readlines()

    # Extract timesteps and sort them
    timesteps = [int(float(line.strip())) * 0.02418884254 for line in lines]
    timesteps.sort()

    # Calculate cumulative number of bonds broken
    cumulative_bonds = [0] * len(timesteps)
    for i, timestep in enumerate(timesteps):
        cumulative_bonds[i] = cumulative_bonds[i - 1] + 1 if i > 0 else 1

    # Calculate average number of bonds broken
    avg_bonds = [100 * count / (no_bonds * total_bonds) for count in cumulative_bonds]

    # Plot average number of bonds broken
    if timesteps[-1] < 500:
        timesteps.append(500)
        avg_bonds.append(avg_bonds[-1])

    ax.plot(timesteps, avg_bonds, label=label, color=color)

def generate_graphs_from_json(json_file,reps):
    # Load the JSON file
    with open(json_file, 'r') as file:
        graphs = json.load(file)
    
    # Read the total number of bonds from allbonds.out
    # with open('../results/bonds/allbonds.out', 'r') as all_file:
    #     lines = all_file.readlines() 

    # Iterate through each graph configuration in the JSON
    for graph_name, graph_data in graphs.items():
        mpl.rcParams['font.family']='DejaVu Sans'
        plt.rcParams['font.size']=18
        plt.rcParams['axes.linewidth']=2
        fig=plt.figure(figsize=(3.37,5.055))
        ax=fig.add_axes([0,0,2,1])
       

        # Iterate through each file info in the current graph data
        for file_info in graph_data['files']:
            # Adjust the path to look for the files in ../results/bonds/
            bond_file_path = file_info['filename']
            plot_data(bond_file_path, file_info['label'], file_info['no_bonds'], 
                      file_info['color'], reps,ax)

        # Set the graph labels and title
        ax.set_xlabel('Femtoseconds (fs)')
        ax.set_ylabel('Percentage of bonds broken (%)')
        ax.set_title(f"Bonds broken for {graph_name}")
        ax.legend()
        ax.set_ylim(0)  # Set the y-axis limits
        ax.set_xlim(0,500)  # Set the x-axis lower limit to 0
        # Save the plot to ../results/graphs/
        output_path = graph_data['output_file']
        # plt.savefig(output_path)
        plt.savefig(output_path,dpi=300, transparent=False,bbox_inches='tight')
        plt.close()

def create_graphs(reps):
    json_file = '../results/graphs_config.json'  # Path to your JSON file
    generate_graphs_from_json(json_file,reps)

if __name__ == "__main__":
    reps=int(input("Number of repeats: "))
    create_graphs(reps)