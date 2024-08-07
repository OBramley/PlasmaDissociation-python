# Qchem.py 15/11/2023

import os
import numpy as np 
import math
from pyqchem import QchemInput, Structure
from pyqchem import get_output_from_qchem

np.set_printoptions(precision =30)

# Electronic structure via QChem.
        
def create_qchem_input(molecule,spin_flip,scf_algorithm="DIIS", Guess=True):
    # Filter indices based on dissociation flag
    molecule = Structure(coordinates=molecule.active_coordinates, symbols=active_atoms, multiplicity=molecule.multiplicity)
    if spin_flip==0:
        if Guess:             
            qc_inp=QchemInput(molecule,
                        scf_guess = molecule.elecinfo,
                        jobtype='force',
                        exchange='BHHLYP',
                        basis='6-31+G*',
                        unrestricted=True,
                        max_scf_cycles=500,
                        sym_ignore=True,
                        scf_algorithm=scf_algorithm,
                        extra_rem_keywords={'input_bohr':'true'}
                        # set_iter=500,
                        )
        else:
            qc_inp=QchemInput(molecule,
                        jobtype='force',
                        exchange='BHHLYP',
                        basis='6-31+G*',
                        unrestricted=True,
                        max_scf_cycles=500,
                        sym_ignore=True,
                        scf_algorithm=scf_algorithm,
                        extra_rem_keywords={'input_bohr':'true'}
                        # set_iter=500,
                        )  
    elif spin_flip==1:   
        if Guess:             
            qc_inp=QchemInput(molecule,
                        scf_guess = molecule.elecinfo,
                        jobtype='force',
                        exchange='BHHLYP',
                        basis='6-31+G*',
                        unrestricted=True,
                        max_scf_cycles=500,
                        sym_ignore=True,
                        scf_algorithm=scf_algorithm,
                        extra_rem_keywords={'input_bohr':'true','spin_flip':'true','set_iter':500},
                        # set_iter=500,
                        cis_n_roots=1,
                        cis_state_deriv=1
                        )
        else:
            qc_inp=QchemInput(molecule,
                        jobtype='force',
                        exchange='BHHLYP',
                        basis='6-31+G*',
                        unrestricted=True,
                        max_scf_cycles=500,
                        sym_ignore=True,
                        scf_algorithm=scf_algorithm,
                        extra_rem_keywords={'input_bohr':'true','spin_flip':'true','set_iter':500},
                        # set_iter=500,
                        cis_n_roots=1,
                        cis_state_deriv=1
                        )
       
    return qc_inp
                      
def run_qchem(ncpu, molecule, spin_flip, Guess=True): 
    qc_inp=create_qchem_input(molecule,spin_flip, scf_algorithm="DIIS", Guess=Guess)
    try:
        output, ee = get_output_from_qchem(qc_inp,processors=ncpu,return_electronic_structure=True)
        molecule.elecinfo=(ee['coefficients'])
    except:
        print('Using DIIS_GDM algorithm')
        # Retry with a different setup
        qc_inp=create_qchem_input(molecule,spin_flip, scf_algorithm="DIIS_GDM", Guess=False)
        try:
            output, ee = get_output_from_qchem(qc_inp,processors=ncpu,return_electronic_structure=True)
            molecule.elecinfo=(ee['coefficients'])
        except:
            with open("ERROR", "w") as file:
                file.write("Error occurred during QChem job. Help.\n" + os.getcwd())
                file.write(output)
            exit()
    # Job completed successfully
    readqchem(output, molecule,spin_flip)
    # Append f.out content to f.all
    with open("f.all", "a") as f_all:
        f_all.write(output)
    return molecule


def readqchem(output, molecule,spin_flip):
    ndim = 3 * molecule.active_no_atoms
    if spin_flip==1:
        enum = output.find('Total energy for state  1:')
        scf_erg = float(output[enum: enum+100].split()[5])
        molecule.update_scf_energy(scf_erg)
        l2t = ' Gradient of the state energy (including CIS Excitation Energy)'
    else:
        enum = output.find('Total energy in the final basis set')
        scf_erg = float(output[enum: enum+100].split()[8])
        molecule.update_scf_energy(scf_erg)
        l2t = ' Gradient of SCF Energy'
    
    output_lines = output.split("\n")
    enum = output.find(l2t)
    output_lines = output[enum:-1].split("\n")
    lines_to_read = 4 * (math.ceil(molecule.active_no_atoms / 6)) +1
    forces = output_lines[1:lines_to_read]
    forces = [line.split() for line in forces]
    f = np.zeros(ndim,dtype = np.float64)
    strt = 0
    for i in range(molecule.active_no_atoms):
        f[i,0] = float(forces[1 + 4 * (i // 6)][i % 6 + 1])
        f[i,1] = float(forces[2 + 4 * (i // 6)][i % 6 + 1])
        f[i,2] = float(forces[3 + 4 * (i // 6)][i % 6 + 1])
       
    f = -f
    f = np.where(f == -0.0, 0.0, f)
    # Update the forces in the Molecule object
    molecule.update_active_forces(f)
    
# -------------------------------------------------------------------------


