# This script generates and submits jobs for a parameter sweep of our new force
import multiprocessing
import subprocess
from ChemChasteDefinitions import *

# generate a list of bash commands
command_list = []

# config files for each of the simulations

# Figure 2
configFisher1 = "/home/chaste/projects/ChemChaste/DataInput/ChemChasteConfigFisherWaveConvergence_1.txt"
configFisher01 = "/home/chaste/projects/ChemChaste/DataInput/ChemChasteConfigFisherWaveConvergence_01.txt"
configFisher001 = "/home/chaste/projects/ChemChaste/DataInput/ChemChasteConfigFisherWaveConvergence_001.txt"
configFisher0001 = "/home/chaste/projects/ChemChaste/DataInput/ChemChasteConfigFisherWaveConvergence_0001.txt"
configFisher00001 = "/home/chaste/projects/ChemChaste/DataInput/ChemChasteConfigFisherWaveConvergence_00001.txt"
configFisher02 = "/home/chaste/projects/ChemChaste/DataInput/ChemChasteConfigFisherWaveConvergence_02.txt"
configFisher002 = "/home/chaste/projects/ChemChaste/DataInput/ChemChasteConfigFisherWaveConvergence_002.txt"
configFisher0002 = "/home/chaste/projects/ChemChaste/DataInput/ChemChasteConfigFisherWaveConvergence_0002.txt"
configFisher00002 = "/home/chaste/projects/ChemChaste/DataInput/ChemChasteConfigFisherWaveConvergence_00002.txt"
configFisher04 = "/home/chaste/projects/ChemChaste/DataInput/ChemChasteConfigFisherWaveConvergence_04.txt"
configFisher004 = "/home/chaste/projects/ChemChaste/DataInput/ChemChasteConfigFisherWaveConvergence_004.txt"
configFisher0004 = "/home/chaste/projects/ChemChaste/DataInput/ChemChasteConfigFisherWaveConvergence_0004.txt"
configFisher00004 = "/home/chaste/projects/ChemChaste/DataInput/ChemChasteConfigFisherWaveConvergence_00004.txt"
configFisher06 = "/home/chaste/projects/ChemChaste/DataInput/ChemChasteConfigFisherWaveConvergence_06.txt"
configFisher006 = "/home/chaste/projects/ChemChaste/DataInput/ChemChasteConfigFisherWaveConvergence_006.txt"
configFisher0006 = "/home/chaste/projects/ChemChaste/DataInput/ChemChasteConfigFisherWaveConvergence_0006.txt"
configFisher00006 = "/home/chaste/projects/ChemChaste/DataInput/ChemChasteConfigFisherWaveConvergence_00006.txt"
configFisher08 = "/home/chaste/projects/ChemChaste/DataInput/ChemChasteConfigFisherWaveConvergence_08.txt"
configFisher008 = "/home/chaste/projects/ChemChaste/DataInput/ChemChasteConfigFisherWaveConvergence_008.txt"
configFisher0008 = "/home/chaste/projects/ChemChaste/DataInput/ChemChasteConfigFisherWaveConvergence_0008.txt"
configFisher00008 = "/home/chaste/projects/ChemChaste/DataInput/ChemChasteConfigFisherWaveConvergence_00008.txt"

# Figure 3
configPatterning = "/home/chaste/projects/ChemChaste/DataInput/ChemChasteConfigSchnackenbergCaseCCFL_domain100.txt"

# Figure 4-5
configCrossfeedingPaper = "/home/chaste/projects/ChemChaste/DataInput/ChemChasteConfigCrossFeedingEnzymePaper.txt"

# timesteps for time parameter sweeping, Figure 2
timesteps = [1,0.1,0.01,0.001,0.0001,0.08,0.06,0.04,0.02,0.008,0.006,0.004,0.002,0.0008,0.0006,0.0004,0.0002]

# Figure 2
for dt in timesteps:
    simulationExecutable = str(determineExecutable(configFisher1))
    command = simulationExecutable + str(1)
    # add config
    command +=  " --config="+configFisher1
    # add simulation type (default "coupled_cell")
    command += " --simulation_type=fisher_equation"
    # add additional commands to override config
    command += " --simulation_timestep="+str(dt)
    command += " --sampling_timestep=1e-1"

    # add simulation to the list
    command_list.append(command)


for dt in timesteps:
    simulationExecutable = str(determineExecutable(configFisher1))
    command = simulationExecutable + str(1)
    # add config
    command +=  " --config="+configFisher1
    # add simulation type (default "coupled_cell")
    command += " --simulation_type=fisher_equation"
    # add additional commands to override config
    command += " --simulation_timestep="+str(dt)
    command += " --sampling_timestep=1e-1"

    # add simulation to the list
    command_list.append(command)


for dt in timesteps:
    simulationExecutable = str(determineExecutable(configFisher01))
    command = simulationExecutable + str(1)
    # add config
    command +=  " --config="+configFisher01
    # add simulation type (default "coupled_cell")
    command += " --simulation_type=fisher_equation"
    # add additional commands to override config
    command += " --simulation_timestep="+str(dt)
    command += " --sampling_timestep=1e-1"

    # add simulation to the list
    command_list.append(command)

for dt in timesteps:
    simulationExecutable = str(determineExecutable(configFisher001))
    command = simulationExecutable + str(1)
    # add config
    command +=  " --config="+configFisher001
    # add simulation type (default "coupled_cell")
    command += " --simulation_type=fisher_equation"
    # add additional commands to override config
    command += " --simulation_timestep="+str(dt)
    command += " --sampling_timestep=1e-1"

    # add simulation to the list
    command_list.append(command)

for dt in timesteps:
    simulationExecutable = str(determineExecutable(configFisher0001))
    command = simulationExecutable + str(1)
    # add config
    command +=  " --config="+configFisher0001
    # add simulation type (default "coupled_cell")
    command += " --simulation_type=fisher_equation"
    # add additional commands to override config
    command += " --simulation_timestep="+str(dt)
    command += " --sampling_timestep=1e-1"

    # add simulation to the list
    command_list.append(command)

for dt in timesteps:
    simulationExecutable = str(determineExecutable(configFisher00001))
    command = simulationExecutable + str(1)
    # add config
    command +=  " --config="+configFisher00001
    # add simulation type (default "coupled_cell")
    command += " --simulation_type=fisher_equation"
    # add additional commands to override config
    command += " --simulation_timestep="+str(dt)
    command += " --sampling_timestep=1e-1"

    # add simulation to the list
    command_list.append(command)

for dt in timesteps:
    simulationExecutable = str(determineExecutable(configFisher00001))
    command = simulationExecutable + str(1)
    # add config
    command +=  " --config="+configFisher00001
    # add simulation type (default "coupled_cell")
    command += " --simulation_type=fisher_equation"
    # add additional commands to override config
    command += " --simulation_timestep="+str(dt)
    command += " --sampling_timestep=1e-1"

    # add simulation to the list
    command_list.append(command)

for dt in timesteps:
    simulationExecutable = str(determineExecutable(configFisher02))
    command = simulationExecutable + str(1)
    # add config
    command +=  " --config="+configFisher02
    # add simulation type (default "coupled_cell")
    command += " --simulation_type=fisher_equation"
    # add additional commands to override config
    command += " --simulation_timestep="+str(dt)
    command += " --sampling_timestep=1e-1"

    # add simulation to the list
    command_list.append(command)

for dt in timesteps:
    simulationExecutable = str(determineExecutable(configFisher002))
    command = simulationExecutable + str(1)
    # add config
    command +=  " --config="+configFisher002
    # add simulation type (default "coupled_cell")
    command += " --simulation_type=fisher_equation"
    # add additional commands to override config
    command += " --simulation_timestep="+str(dt)
    command += " --sampling_timestep=1e-1"

    # add simulation to the list
    command_list.append(command)

for dt in timesteps:
    simulationExecutable = str(determineExecutable(configFisher0002))
    command = simulationExecutable + str(1)
    # add config
    command +=  " --config="+configFisher0002
    # add simulation type (default "coupled_cell")
    command += " --simulation_type=fisher_equation"
    # add additional commands to override config
    command += " --simulation_timestep="+str(dt)
    command += " --sampling_timestep=1e-1"

    # add simulation to the list
    command_list.append(command)

for dt in timesteps:
    simulationExecutable = str(determineExecutable(configFisher00002))
    command = simulationExecutable + str(1)
    # add config
    command +=  " --config="+configFisher00002
    # add simulation type (default "coupled_cell")
    command += " --simulation_type=fisher_equation"
    # add additional commands to override config
    command += " --simulation_timestep="+str(dt)
    command += " --sampling_timestep=1e-1"

    # add simulation to the list
    command_list.append(command)

for dt in timesteps:
    simulationExecutable = str(determineExecutable(configFisher00002))
    command = simulationExecutable + str(1)
    # add config
    command +=  " --config="+configFisher00002
    # add simulation type (default "coupled_cell")
    command += " --simulation_type=fisher_equation"
    # add additional commands to override config
    command += " --simulation_timestep="+str(dt)
    command += " --sampling_timestep=1e-1"

    # add simulation to the list
    command_list.append(command)

for dt in timesteps:
    simulationExecutable = str(determineExecutable(configFisher04))
    command = simulationExecutable + str(1)
    # add config
    command +=  " --config="+configFisher04
    # add simulation type (default "coupled_cell")
    command += " --simulation_type=fisher_equation"
    # add additional commands to override config
    command += " --simulation_timestep="+str(dt)
    command += " --sampling_timestep=1e-1"

    # add simulation to the list
    command_list.append(command)

for dt in timesteps:
    simulationExecutable = str(determineExecutable(configFisher004))
    command = simulationExecutable + str(1)
    # add config
    command +=  " --config="+configFisher004
    # add simulation type (default "coupled_cell")
    command += " --simulation_type=fisher_equation"
    # add additional commands to override config
    command += " --simulation_timestep="+str(dt)
    command += " --sampling_timestep=1e-1"

    # add simulation to the list
    command_list.append(command)

for dt in timesteps:
    simulationExecutable = str(determineExecutable(configFisher0004))
    command = simulationExecutable + str(1)
    # add config
    command +=  " --config="+configFisher0004
    # add simulation type (default "coupled_cell")
    command += " --simulation_type=fisher_equation"
    # add additional commands to override config
    command += " --simulation_timestep="+str(dt)
    command += " --sampling_timestep=1e-1"

    # add simulation to the list
    command_list.append(command)

for dt in timesteps:
    simulationExecutable = str(determineExecutable(configFisher00004))
    command = simulationExecutable + str(1)
    # add config
    command +=  " --config="+configFisher00004
    # add simulation type (default "coupled_cell")
    command += " --simulation_type=fisher_equation"
    # add additional commands to override config
    command += " --simulation_timestep="+str(dt)
    command += " --sampling_timestep=1e-1"

    # add simulation to the list
    command_list.append(command)

for dt in timesteps:
    simulationExecutable = str(determineExecutable(configFisher00004))
    command = simulationExecutable + str(1)
    # add config
    command +=  " --config="+configFisher00004
    # add simulation type (default "coupled_cell")
    command += " --simulation_type=fisher_equation"
    # add additional commands to override config
    command += " --simulation_timestep="+str(dt)
    command += " --sampling_timestep=1e-1"

    # add simulation to the list
    command_list.append(command)

for dt in timesteps:
    simulationExecutable = str(determineExecutable(configFisher06))
    command = simulationExecutable + str(1)
    # add config
    command +=  " --config="+configFisher06
    # add simulation type (default "coupled_cell")
    command += " --simulation_type=fisher_equation"
    # add additional commands to override config
    command += " --simulation_timestep="+str(dt)
    command += " --sampling_timestep=1e-1"

    # add simulation to the list
    command_list.append(command)

for dt in timesteps:
    simulationExecutable = str(determineExecutable(configFisher006))
    command = simulationExecutable + str(1)
    # add config
    command +=  " --config="+configFisher006
    # add simulation type (default "coupled_cell")
    command += " --simulation_type=fisher_equation"
    # add additional commands to override config
    command += " --simulation_timestep="+str(dt)
    command += " --sampling_timestep=1e-1"

    # add simulation to the list
    command_list.append(command)

for dt in timesteps:
    simulationExecutable = str(determineExecutable(configFisher0006))
    command = simulationExecutable + str(1)
    # add config
    command +=  " --config="+configFisher0006
    # add simulation type (default "coupled_cell")
    command += " --simulation_type=fisher_equation"
    # add additional commands to override config
    command += " --simulation_timestep="+str(dt)
    command += " --sampling_timestep=1e-1"

    # add simulation to the list
    command_list.append(command)

for dt in timesteps:
    simulationExecutable = str(determineExecutable(configFisher00006))
    command = simulationExecutable + str(1)
    # add config
    command +=  " --config="+configFisher00006
    # add simulation type (default "coupled_cell")
    command += " --simulation_type=fisher_equation"
    # add additional commands to override config
    command += " --simulation_timestep="+str(dt)
    command += " --sampling_timestep=1e-1"

    # add simulation to the list
    command_list.append(command)

for dt in timesteps:
    simulationExecutable = str(determineExecutable(configFisher00006))
    command = simulationExecutable + str(1)
    # add config
    command +=  " --config="+configFisher00006
    # add simulation type (default "coupled_cell")
    command += " --simulation_type=fisher_equation"
    # add additional commands to override config
    command += " --simulation_timestep="+str(dt)
    command += " --sampling_timestep=1e-1"

    # add simulation to the list
    command_list.append(command)


for dt in timesteps:
    simulationExecutable = str(determineExecutable(configFisher08))
    command = simulationExecutable + str(1)
    # add config
    command +=  " --config="+configFisher08
    # add simulation type (default "coupled_cell")
    command += " --simulation_type=fisher_equation"
    # add additional commands to override config
    command += " --simulation_timestep="+str(dt)
    command += " --sampling_timestep=1e-1"

    # add simulation to the list
    command_list.append(command)

for dt in timesteps:
    simulationExecutable = str(determineExecutable(configFisher008))
    command = simulationExecutable + str(1)
    # add config
    command +=  " --config="+configFisher008
    # add simulation type (default "coupled_cell")
    command += " --simulation_type=fisher_equation"
    # add additional commands to override config
    command += " --simulation_timestep="+str(dt)
    command += " --sampling_timestep=1e-1"

    # add simulation to the list
    command_list.append(command)

for dt in timesteps:
    simulationExecutable = str(determineExecutable(configFisher0008))
    command = simulationExecutable + str(1)
    # add config
    command +=  " --config="+configFisher0008
    # add simulation type (default "coupled_cell")
    command += " --simulation_type=fisher_equation"
    # add additional commands to override config
    command += " --simulation_timestep="+str(dt)
    command += " --sampling_timestep=1e-1"

    # add simulation to the list
    command_list.append(command)

for dt in timesteps:
    simulationExecutable = str(determineExecutable(configFisher00008))
    command = simulationExecutable + str(1)
    # add config
    command +=  " --config="+configFisher00008
    # add simulation type (default "coupled_cell")
    command += " --simulation_type=fisher_equation"
    # add additional commands to override config
    command += " --simulation_timestep="+str(dt)
    command += " --sampling_timestep=1e-1"

    # add simulation to the list
    command_list.append(command)

for dt in timesteps:
    simulationExecutable = str(determineExecutable(configFisher00008))
    command = simulationExecutable + str(1)
    # add config
    command +=  " --config="+configFisher00008
    # add simulation type (default "coupled_cell")
    command += " --simulation_type=fisher_equation"
    # add additional commands to override config
    command += " --simulation_timestep="+str(dt)
    command += " --sampling_timestep=1e-1"

    # add simulation to the list
    command_list.append(command)


# Figure 3
for simulation_id in range(1):
    simulationExecutable = str(determineExecutable(configPatterning))
    command = simulationExecutable + str(simulation_id)
    # add config
    command +=  " --config="+configPatterning
    # add simulation type (default "coupled_cell")
    command += " --simulation_type=domain_only"
    # add additional commands to override config
    command += " --simulation_timestep=1e-3"
    command += " --sampling_timestep=1e-1"

    # add simulation to the list
    command_list.append(command)

# Figure 4-5
for simulation_id in range(1):
    simulationExecutable = str(determineExecutable(configCrossfeedingPaper))
    command = simulationExecutable + str(simulation_id)
    # add config
    command +=  " --config="+configCrossfeedingPaper
    # add simulation type (default "coupled_cell")
    command += " --simulation_type=complex_cell"
    # add additional commands to override config
    command += " --simulation_timestep=1e-2"
    command += " --sampling_timestep=1e-1"

    # add simulation to the list
    command_list.append(command)


# use `count' no of processes 
count = 3 # use multiprocessing.cpu_count() for the number of cores on your machine

# generate a pool of workers
pool = multiprocessing.Pool(processes=count)

# ... and pass the list of bash commands to the pool
pool.map(execute_command, command_list)


