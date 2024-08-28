import os
import json
import getpass
import subprocess
from random import randrange

with open('inputs.json') as f:
    inputs=json.load(f)
reps=inputs["setup"]["repeats"]
rep_flg=0
hold_ids=''
EXDIR=os.getcwd()
for i in range(1,reps+1):
    EXDIR1=EXDIR+'/rep-'+str(i)
    filename = EXDIR1+"/output/complete.all"
    if(i%50==0):
        print(i,' repetitions checked')
        
    if not os.path.isfile(filename):
        name='Plasma_'+inputs["setup"]["Runfolder"]+"_"+str(i)+"_"+str(randrange(1000))
        if(rep_flg==0):
            rep_flg=1
            hold_ids=name
        else: 
            hold_ids=hold_ids+' '+name
        os.chdir(EXDIR1)
        print('Repetition ',i,' incomplete') 
        file1="Plasma_"+inputs["setup"]["Runfolder"]+"_"+str(i)+".sh"
        f=open(file1,"w")
        f.write("#$ -cwd -V \n")
        f.write("#$ -l h_vmem=1G,h_rt=48:00:00 \n")
        f.write("#$ -N "+name+" \n")
        f.write("#$ -pe smp "+str(inputs["setup"]["cores"])+" \n") #Use shared memory parallel environemnt
        f.write("module swith gnu intel \n")
        f.write("module add intel \n") 
        f.write("module switch intel/18.0.2 \n")
        f.write("module load cuda/11.1.1 \n")
        f.write("module add mkl \n")
        f.write("module add anaconda \n")
        f.write("source activate scatter \n")
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
        command = ['qsub','-N',name, file1]
        subprocess.call(command)

if(rep_flg==1):
    os.chdir(EXDIR)
    command = ['qsub','-N','Check_'+str(randrange(1000)), '-hold_jid', hold_ids, 'recheck.sh'] 
    subprocess.call(command)
    print(command)
else:
    file2="complete.sh"
    f=open(file2,"w")
    f.write("#$ -cwd -V \n")
    f.write("#$ -l h_vmem=1G,h_rt=00:30:00 \n")
    f.write("find -name "*.o*" -delete")
    f.write("find -name "*.e*" -delete")
    f.write("module add anaconda \n")
    f.write("source activate scatter \n")
    f.write("cd code \n")
    f.write("python result.py \n")
    f.write("python graphs.py \n")
    f.close()
    command = ['qsub','-N','Complete_'+EXDIR, file2]
    subprocess.call(command) 