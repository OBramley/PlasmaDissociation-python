import re
import numpy as np 
import math
def convert_to_bohr_units(r0_angstrom):
    # Convert coordinates from Angstrom to Bohr units
    bohr_conversion_factor = 0.529177
    return r0_angstrom * bohr_conversion_factor

def make_geometry_input(inp):
    with open(inp, "r") as file:
        lines = file.readlines()

    natom, nst = map(int, lines[0].split())
    mc, mult = map(int, lines[3].split())
    # Write the geometry data to geom.in
    with open("geom.in", "w") as geom_file:
        # geom_file.write(f"{natom:2d} {nst:2d}\n")
        geom_file.write(f"{mc:2d} {mult:2d}\n")
        current_line = 4  # Start with the first line after the header

        for _ in range(natom):
            line = lines[current_line]
            data = line.split()
            if len(data) >= 4:  # Check if the line has at least four elements (atom and three coordinates)
                atom = data[0]
                r0_angstrom = [float(val) for val in data[1:4]]
                r0_bohr = [convert_to_bohr_units(coord) for coord in r0_angstrom]  # Convert each coordinate separately
                geom_file.write(f"{atom:1s} {r0_bohr[0]:15.8f} {r0_bohr[1]:15.8f} {r0_bohr[2]:15.8f}\n")
            else:
                print(f"WARNING: Skipping improperly formatted line: {line}")
            current_line += 1

def find_line(filename, target_line):
    try:
        with open(filename, 'r') as file:
            for line in file:
                if target_line in line:
                    return line.strip(), line.split()  # Return the found line and its content
    except FileNotFoundError:
        print("File not found:", filename)
        return None, None

def q_to_prop(output_file):
    try:
        with open('t.0', 'r') as file:
            first_line = file.readline().strip()
            natoms = int(first_line.split()[0])
            nst = int(first_line.split()[1])
            # print(natoms,nst)
    except (FileNotFoundError, ValueError, IndexError):
        print("Error reading the number of atoms from the file.")
        return None
    ndim = 3*natoms
    l1t = ' Excited state   1: excitation energy (eV) ='
    l2t = ' Gradient of the state energy (including CIS Excitation Energy)'

    with open('f.out', 'r') as file:
        f = np.zeros((1,nst*ndim))
        e= np.zeros(nst)
        C = np.zeros((1,ndim))
        l1_found = False
        l2_found = False
        found_target = False
        for line in file:
            if found_target:
                # Read the line below the target line
                data_line = line.strip()
                match = re.search(r'-?\d+\.\d+', data_line)
                if match:
                    e[0] = float(match.group())
                else:
                    print("Number not found in the line.")
                break
            if l1t in line:
                found_target = True
        
        found_target = False
        lines_read = 0 
        skip_counter = 0 
        lines_to_read = 4*(math.ceil(natoms/6))-1
        start_index = 0
        for line in file:
            if found_target and lines_read < lines_to_read:
                if skip_counter != 3:  # Skip every fourth line
                    data_line = line.strip()
                    lines_read += 1
                    start_index += 1
                    # print('lines_read', lines_read)
                    parts = data_line.split("  ")
                    # print(parts)
                    m = start_index
                    for j in range(1,len(parts)):
                        f[0,m-1] = parts[j]
                        # print(m)
                        m = m+3
                if skip_counter == 3: 
                    lines_read += 1
                    start_index = 6*start_index
                skip_counter = (skip_counter + 1) % 4
            if l2t in line:
                found_target = True
                # Skip the header line
                next(file)
        f = -f
        f = np.where(f == -0.0, 0.0, f)
        with open(output_file, "a") as out:
            for i in range(nst):
                out.write(f'{e[i]:.16e} {i:8d}\n')
            out.write(" \n")
            for i in range(0,2*ndim):
                if i<ndim:
                    print(i)
                    out.write(f'{f[0, i]:.16e} 1 {i+1}\n')
                if i>=ndim:
                    print(i,'2')
                    out.write(f'{f[0, i]:.16e} 2 {i-ndim+1}\n')
            out.write(" \n")
            for i in range(ndim):
                out.write(f'{C[0, i]:.16e} 1 2 {i+1}\n')
