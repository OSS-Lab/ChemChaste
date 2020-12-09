# This script generates and submits jobs for a parameter sweep of our new force
import multiprocessing
import subprocess
from ChemChasteDefinitions import *


# generate a list of bash commands
command_list = []

# config files for each of the 
config = "/home/chaste/projects/ChemChaste/DataInput/ChemChasteConfig.txt"
config0 = "/home/chaste/projects/ChemChaste/DataInput/ChemChasteConfigConstDiff.txt"
config1 = "/home/chaste/projects/ChemChaste/DataInput/ChemChasteConfigStepDiff.txt"
config2 = "/home/chaste/projects/ChemChaste/DataInput/ChemChasteConfigER.txt"
config3 = "/home/chaste/projects/ChemChaste/DataInput/ChemChasteConfigMP.txt"



#print(determineExecutable("/home/chaste/projects/ChemChaste/DataInput/ChemChasteConfig.txt"))

# ... and populate it with calls to a Chaste executable
for simulation_id in range(2):
    command = str(determineExecutable(config0)) + str(simulation_id)+" --config="+config0
    command_list.append(command)

for simulation_id in range(2):
    command = str(determineExecutable(config1)) + str(simulation_id)+" --config="+config1
    command_list.append(command.strip())

for simulation_id in range(2):
    command = str(determineExecutable(config2)) + str(simulation_id)+" --config="+config2
    command_list.append(command.strip())

for simulation_id in range(2):
    command = str(determineExecutable(config3)) + str(simulation_id)+" --config="+config3
    command_list.append(command.strip())
    

# use `count' no of processes 
count = 4 # use multiprocessing.cpu_count() for the number of cores on your machine

# generate a pool of workers
pool = multiprocessing.Pool(processes=count)

# ... and pass the list of bash commands to the pool
pool.map(execute_command, command_list)


