import sys
import os
import json
import time 
from mol import Molecule
import mol
# import elec
# import prop
# import output as out
# import result


start_time = time.time()
if __name__ == "__main__":
    # Load inputs
    with open('../inputs.json') as f:
        inputs=json.load(f)
    # Check basic arguements
    # reps=inputs["setup"]["repeats"]
    # ncpu=inputs["setup"]["cores"]
    # nstates=inputs["run"]["States"]
    # nbranch=inputs["run"]["Branches"]
    # increment=inputs["run"]["Timestep"]
    # endstep=inputs["run"]["Tot_timesteps"]
    # geom_start=inputs["run"]["Geom_start"]
    # spin_flip=inputs["run"]["Spin_flip"]
    # method = inputs["run"]["method"]

    nstates=1
    spin_flip=0
    filename = 'molecule.json'
    molecule1 = Molecule.setup_from_json(nstates,spin_flip,filename)
    molecule1 = elec.run_elec_structure(molecule1,ncpu,n,nstates,spin_flip,method,Guess=False)
    if(molecule1.timestep==0):
        startstep = 1
    else: 
        startstep = molecule1.timestep / increment

Guess = True
for i in range(int(startstep), endstep+1):
    molecule2 = prop_prelim(molecule1,increment)
    molecule1 = elec.run_elec_structure(molecule2,ncpu,n,nstates,spin_flip,method,Guess=Guess)
    molecule1.elecinfo = molecule2.elecinfo
    
    molecule1 = prop.prop_2(molecule1, molecule2, n, nstates, increment)
    molecule1, dissociated = prop.fragements(molecule1,spin_flip)
    molecule1 = prop.prop_diss(molecule1,increment)
    out.output_molecule(molecule1)
    if dissociated == 0:
        Guess = True
    else:
        Guess = False
result.process_results()
end_time=time.time()
print("Time taken to run: ", end_time-start_time)