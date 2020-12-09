# This script generates and submits jobs for a parameter sweep of our new force
import multiprocessing
import subprocess
from ChemChasteDefinitions import *

# generate a list of bash commands
command_list = []

# config files for each of the simulations
config0 = "/home/chaste/projects/ChemChaste/DataInput/ChemChasteConfigConstDiff.txt"
config1 = "/home/chaste/projects/ChemChaste/DataInput/ChemChasteConfigStepDiff.txt"
config2 = "/home/chaste/projects/ChemChaste/DataInput/ChemChasteConfigER.txt"
config3 = "/home/chaste/projects/ChemChaste/DataInput/ChemChasteConfigMP.txt"
config4 = "/home/chaste/projects/ChemChaste/DataInput/ChemChasteConfigCrossFeeding.txt"
config5 = "/home/chaste/projects/ChemChaste/DataInput/ChemChasteConfigFisherWave.txt"
config6 = "/home/chaste/projects/ChemChaste/DataInput/ChemChasteConfigSchnackenbergCaseA.txt"
config7 = "/home/chaste/projects/ChemChaste/DataInput/ChemChasteConfigSchnackenbergCaseB.txt"
config8 = "/home/chaste/projects/ChemChaste/DataInput/ChemChasteConfigSchnackenbergCaseC.txt"


timeStepList = [1,0.1,0.01,0.001,0.0001]

# populate command_list with ChemChaste simulation calls
for simulation_id in range(1):
    simulationExecutable = str(determineExecutable(config0))
    command = simulationExecutable + str(simulation_id)
    # add config
    command +=  " --config="+config0
    # add simulation type (default "coupled_cell")
    command += " --simulation_type=coupled_cell"
    # add additional commands to override config
    command += " --simulation_timestep=1e-3"

    # add simulation to the list
    command_list.append(command)

# add aditional simulations to the list
for simulation_id in range(1):
    simulationExecutable = str(determineExecutable(config1))
    command = simulationExecutable + str(simulation_id)
    # add config
    command +=  " --config="+config1
    # add simulation type (default "coupled_cell")
    command += " --simulation_type=coupled_cell"
    # add additional commands to override config
    command += " --simulation_timestep=1e-3"

    # add simulation to the list
    command_list.append(command)

for simulation_id in range(1):
    simulationExecutable = str(determineExecutable(config2))
    command = simulationExecutable + str(simulation_id)
    # add config
    command +=  " --config="+config2
    # add simulation type (default "coupled_cell")
    command += " --simulation_type=coupled_cell"
    # add additional commands to override config
    command += " --simulation_timestep=1e-3"

    # add simulation to the list
    command_list.append(command)

for simulation_id in range(1):
    simulationExecutable = str(determineExecutable(config3))
    command = simulationExecutable + str(simulation_id)
    # add config
    command +=  " --config="+config3
    # add simulation type (default "coupled_cell")
    command += " --simulation_type=coupled_cell"
    # add additional commands to override config
    command += " --simulation_timestep=1e-3"

    # add simulation to the list
    command_list.append(command)

for simulation_id in range(1):
    simulationExecutable = str(determineExecutable(config4))
    command = simulationExecutable + str(simulation_id)
    # add config
    command +=  " --config="+config4
    # add simulation type (default "coupled_cell")
    command += " --simulation_type=coupled_cell"
    # add additional commands to override config
    command += " --simulation_timestep=1e-3"

    # add simulation to the list
    command_list.append(command)


for timestep in timeStepList:
    simulationExecutable = str(determineExecutable(config5))
    command = simulationExecutable + str(0)
    # add config
    command +=  " --config="+config5
    # add simulation type (default "coupled_cell")
    command += " --simulation_type=fisher_equation"
    # add additional commands to override config
    command += " --simulation_timestep="+str(timestep)

    # add simulation to the list
    command_list.append(command)

for simulation_id in range(1):
    simulationExecutable = str(determineExecutable(config6))
    command = simulationExecutable + str(simulation_id)
    # add config
    command +=  " --config="+config6
    # add simulation type (default "coupled_cell")
    command += " --simulation_type=domain_only"
    # add additional commands to override config
    command += " --simulation_timestep=1e-3"

    # add simulation to the list
    command_list.append(command)

for timestep in timeStepList:
    simulationExecutable = str(determineExecutable(config7))
    command = simulationExecutable + str(0)
    # add config
    command +=  " --config="+config7
    # add simulation type (default "coupled_cell")
    command += " --simulation_type=domain_only"
    # add additional commands to override config
    command += " --simulation_timestep="+str(timestep)

    # add simulation to the list
    command_list.append(command)

for simulation_id in range(1):
    simulationExecutable = str(determineExecutable(config8))
    command = simulationExecutable + str(simulation_id)
    # add config
    command +=  " --config="+config8
    # add simulation type (default "coupled_cell")
    command += " --simulation_type=domain_only"
    # add additional commands to override config
    command += " --simulation_timestep=1e-3"

    # add simulation to the list
    command_list.append(command)



# use `count' no of processes 
count = 4 # use multiprocessing.cpu_count() for the number of cores on your machine

# generate a pool of workers
pool = multiprocessing.Pool(processes=count)

# ... and pass the list of bash commands to the pool
pool.map(execute_command, command_list)


