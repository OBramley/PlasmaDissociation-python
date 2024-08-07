import os
import json
import getpass
import subprocess

with open('../inputs.json') as f:
    inputs=json.load(f)
reps=400
for i in range(1,reps+1):
    EXDIR1 = '../rep-'+str(i)
    os.chdir(EXDIR1)
    with open('t1.all','r')as file:
        data= file.read()
    if '25000.000000' not in data:
        print(i) 

        file1="Plasma_"+inputs["setup"]["Runfolder"]+"_"+str(i)+".sh"
        f=open(file1,"w")
        f.write("#$ -cwd -V \n")
        f.write("#$ -l h_vmem=1G,h_rt=10:00:00 \n")
        f.write("#$ -N Plasma_"+inputs["setup"]["Runfolder"]+"_"+str(i)+"\n")
        f.write("#$ -pe smp "+str(inputs["setup"]["cores"])+" \n") #Use shared memory parallel environemnt
        # if(inputs["run"]["GPU"]==1):
            # f.write('#$ -l coproc_p100=1 \n')
        # if(inputs["run"]["method"]=="QChem"):
        f.write("module load test qchem \n")
        f.write("module load qchem \n")
        f.write("mkdir $TMPDIR/qchemlocal\n")
        f.write('tar -xzvf /nobackup/'+getpass.getuser()+'/scatter/qchem.tar.gz -C $TMPDIR/qchemlocal\n')
        f.write('qchemlocal=$TMPDIR/qchemlocal/apps/applications/qchem/6.0.1/1/default\n')
        f.write('export QCHEM_HOME="$qchemlocal"\n')
        f.write('export QC="$qchemlocal"\n')
        f.write('export QCAUX="$QC/qcaux"\n')
        f.write('export QCPROG="$QC/exe/qcprog.exe"\n')
        f.write('export QCPROG_S="$QC/exe/qcprog.exe_s"\n')
        f.write('export PATH="$PATH:$QC/exe:$QC/bin"\n')
        f.write("export QCSCRATCH="+EXDIR1+"/tmp \n")
        f.write("python ./../code/main.py")
        f.close()
        command = ['qsub','-N','Plasma_'+inputs["setup"]["Runfolder"]+str(i), file1]
        subprocess.call(command)