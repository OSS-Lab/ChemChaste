# This script generates and submits jobs for a parameter sweep of our new force
import multiprocessing
import subprocess
from ChemChasteDefinitions import *

# generate a list of bash commands
command_list = []

# config files for each of the simulations

configCrossfeedingPaper = "/home/chaste/projects/ChemChaste/DataInput/ChemChasteConfigCrossFeedingEnzymePaper.txt"
configCrossfeedingPaperDirichlet = "/home/chaste/projects/ChemChaste/DataInput/ChemChasteConfigCrossFeedingEnzymePaperDirichlet.txt"

configCrossfeeding1 = "/home/chaste/projects/ChemChaste/DataInput/ChemChasteConfigCrossFeedingEnzyme1.txt"
configCrossfeeding2 = "/home/chaste/projects/ChemChaste/DataInput/ChemChasteConfigCrossFeedingEnzyme2.txt"
configCrossfeeding3 = "/home/chaste/projects/ChemChaste/DataInput/ChemChasteConfigCrossFeedingEnzyme3.txt"
configCrossfeeding4 = "/home/chaste/projects/ChemChaste/DataInput/ChemChasteConfigCrossFeedingEnzyme4.txt"
configCrossfeeding5 = "/home/chaste/projects/ChemChaste/DataInput/ChemChasteConfigCrossFeedingEnzyme5.txt"
configCrossfeeding7 = "/home/chaste/projects/ChemChaste/DataInput/ChemChasteConfigCrossFeedingEnzyme7.txt"

configCrossfeedingDirichlet1 = "/home/chaste/projects/ChemChaste/DataInput/ChemChasteConfigCrossFeedingEnzymeDirichlet1.txt"
configCrossfeedingDirichlet2 = "/home/chaste/projects/ChemChaste/DataInput/ChemChasteConfigCrossFeedingEnzymeDirichlet2.txt"
configCrossfeedingDirichlet3 = "/home/chaste/projects/ChemChaste/DataInput/ChemChasteConfigCrossFeedingEnzymeDirichlet3.txt"
configCrossfeedingDirichlet4 = "/home/chaste/projects/ChemChaste/DataInput/ChemChasteConfigCrossFeedingEnzymeDirichlet4.txt"
configCrossfeedingDirichlet5 = "/home/chaste/projects/ChemChaste/DataInput/ChemChasteConfigCrossFeedingEnzymeDirichlet5.txt"
configCrossfeedingDirichlet7 = "/home/chaste/projects/ChemChaste/DataInput/ChemChasteConfigCrossFeedingEnzymeDirichlet7.txt"

configEnvironmentProject = "/home/chaste/projects/ChemChaste/DataInput/ChemChasteConfigEnvironmentProject.txt"
configEnvironmentProject2 = "/home/chaste/projects/ChemChaste/DataInput/ChemChasteConfigEnvironmentProject2.txt"
configEnvironmentProject10 = "/home/chaste/projects/ChemChaste/DataInput/ChemChasteConfigEnvironmentProject10.txt"

configEnvironmentProjectReverse = "/home/chaste/projects/ChemChaste/DataInput/ChemChasteConfigEnvironmentProjectReverse.txt"
configEnvironmentProjectReverse2 = "/home/chaste/projects/ChemChaste/DataInput/ChemChasteConfigEnvironmentProject2Reverse.txt"
configEnvironmentProjectReverse10 = "/home/chaste/projects/ChemChaste/DataInput/ChemChasteConfigEnvironmentProject10Reverse.txt"

configEnvironmentParameterSweep = "/home/chaste/projects/ChemChaste/DataInput/ChemChasteConfigEnvironmentParameterSweep.txt"

configCrossfeedingPaper = "/home/chaste/projects/ChemChaste/DataInput/ChemChasteConfigCrossFeedingEnzymePaper.txt"

configsweep1 = "/home/chaste/projects/ChemChaste/DataInput/ChemChasteConfigAlternateDiffusionGore1BiStable.txt"
configsweep2 = "/home/chaste/projects/ChemChaste/DataInput/ChemChasteConfigAlternateDiffusionGore2BiStable.txt"
configsweep3 = "/home/chaste/projects/ChemChaste/DataInput/ChemChasteConfigAlternateDiffusionGore3BiStable.txt"

configsweep4 = "/home/chaste/projects/ChemChaste/DataInput/ChemChasteConfigCompetitionAlternateDiffusionGore1BiStable.txt"
configsweep5 = "/home/chaste/projects/ChemChaste/DataInput/ChemChasteConfigCompetitionAlternateDiffusionGore2BiStable.txt"
configsweep6 = "/home/chaste/projects/ChemChaste/DataInput/ChemChasteConfigCompetitionAlternateDiffusionGore3BiStable.txt"

configsweep7 = "/home/chaste/projects/ChemChaste/DataInput/ChemChasteConfigCompetitionDiffusionGore1BiStable.txt"
configsweep8 = "/home/chaste/projects/ChemChaste/DataInput/ChemChasteConfigCompetitionDiffusionGore2BiStable.txt"
configsweep9 = "/home/chaste/projects/ChemChaste/DataInput/ChemChasteConfigCompetitionDiffusionGore3BiStable.txt"

configsweep10 = "/home/chaste/projects/ChemChaste/DataInput/ChemChasteConfigLayerDiffusionGore1BiStable.txt"
configsweep11 = "/home/chaste/projects/ChemChaste/DataInput/ChemChasteConfigLayerDiffusionGore2BiStable.txt"
configsweep12 = "/home/chaste/projects/ChemChaste/DataInput/ChemChasteConfigLayerDiffusionGore3BiStable.txt"

configsweep13 = "/home/chaste/projects/ChemChaste/DataInput/ChemChasteConfigDecayDiffusionGore1BiStable.txt"
configsweep14 = "/home/chaste/projects/ChemChaste/DataInput/ChemChasteConfigDecayDiffusionGore3BiStable.txt"


sweep1 = ParameterSweeping(configsweep1)
sweep2 = ParameterSweeping(configsweep2)
sweep3 = ParameterSweeping(configsweep3)
sweep4 = ParameterSweeping(configsweep4)
sweep5 = ParameterSweeping(configsweep5)
sweep6 = ParameterSweeping(configsweep6)
sweep7 = ParameterSweeping(configsweep7)
sweep8 = ParameterSweeping(configsweep8)
sweep9 = ParameterSweeping(configsweep9)
sweep10 = ParameterSweeping(configsweep10)
sweep11 = ParameterSweeping(configsweep11)
sweep12 = ParameterSweeping(configsweep12)




timesteps = [1,0.1,0.01,0.001,0.0001,0.08,0.06,0.04,0.02,0.008,0.006,0.004,0.002,0.0008,0.0006,0.0004,0.0002]

sweepConfigList = ParameterSweeping(configEnvironmentParameterSweep)

for paramSimConfig in sweep1:
    simulationExecutable = str(determineExecutable(paramSimConfig))
    command = simulationExecutable + str(1)
    # add config
    command +=  " --config="+paramSimConfig
    # add simulation type (default "coupled_cell")
    command += " --simulation_type=environment_cell"
    # add additional commands to override config
    command += " --simulation_timestep=1e-3"
    command += " --sampling_timestep=1e-1"

    # add simulation to the list
    command_list.append(command)

for paramSimConfig in sweep1:
    simulationExecutable = str(determineExecutable(paramSimConfig))
    command = simulationExecutable + str(1)
    # add config
    command +=  " --config="+paramSimConfig
    # add simulation type (default "coupled_cell")
    command += " --simulation_type=environment_cell"
    # add additional commands to override config
    command += " --simulation_timestep=1e-3"
    command += " --sampling_timestep=1e-1"

    # add simulation to the list
    command_list.append(command)

for paramSimConfig in sweep2:
    simulationExecutable = str(determineExecutable(paramSimConfig))
    command = simulationExecutable + str(1)
    # add config
    command +=  " --config="+paramSimConfig
    # add simulation type (default "coupled_cell")
    command += " --simulation_type=environment_cell"
    # add additional commands to override config
    command += " --simulation_timestep=1e-3"
    command += " --sampling_timestep=1e-1"

    # add simulation to the list
    command_list.append(command)

for paramSimConfig in sweep3:
    simulationExecutable = str(determineExecutable(paramSimConfig))
    command = simulationExecutable + str(1)
    # add config
    command +=  " --config="+paramSimConfig
    # add simulation type (default "coupled_cell")
    command += " --simulation_type=environment_cell"
    # add additional commands to override config
    command += " --simulation_timestep=1e-3"
    command += " --sampling_timestep=1e-1"

    # add simulation to the list
    command_list.append(command)

for paramSimConfig in sweep4:
    simulationExecutable = str(determineExecutable(paramSimConfig))
    command = simulationExecutable + str(1)
    # add config
    command +=  " --config="+paramSimConfig
    # add simulation type (default "coupled_cell")
    command += " --simulation_type=environment_cell"
    # add additional commands to override config
    command += " --simulation_timestep=1e-3"
    command += " --sampling_timestep=1e-1"

    # add simulation to the list
    command_list.append(command)

for paramSimConfig in sweep5:
    simulationExecutable = str(determineExecutable(paramSimConfig))
    command = simulationExecutable + str(1)
    # add config
    command +=  " --config="+paramSimConfig
    # add simulation type (default "coupled_cell")
    command += " --simulation_type=environment_cell"
    # add additional commands to override config
    command += " --simulation_timestep=1e-3"
    command += " --sampling_timestep=1e-1"

    # add simulation to the list
    command_list.append(command)

for paramSimConfig in sweep6:
    simulationExecutable = str(determineExecutable(paramSimConfig))
    command = simulationExecutable + str(1)
    # add config
    command +=  " --config="+paramSimConfig
    # add simulation type (default "coupled_cell")
    command += " --simulation_type=environment_cell"
    # add additional commands to override config
    command += " --simulation_timestep=1e-3"
    command += " --sampling_timestep=1e-1"

    # add simulation to the list
    command_list.append(command)

for paramSimConfig in sweep7:
    simulationExecutable = str(determineExecutable(paramSimConfig))
    command = simulationExecutable + str(1)
    # add config
    command +=  " --config="+paramSimConfig
    # add simulation type (default "coupled_cell")
    command += " --simulation_type=environment_cell"
    # add additional commands to override config
    command += " --simulation_timestep=1e-3"
    command += " --sampling_timestep=1e-1"

    # add simulation to the list
    command_list.append(command)

for paramSimConfig in sweep8:
    simulationExecutable = str(determineExecutable(paramSimConfig))
    command = simulationExecutable + str(1)
    # add config
    command +=  " --config="+paramSimConfig
    # add simulation type (default "coupled_cell")
    command += " --simulation_type=environment_cell"
    # add additional commands to override config
    command += " --simulation_timestep=1e-3"
    command += " --sampling_timestep=1e-1"

    # add simulation to the list
    command_list.append(command)

for paramSimConfig in sweep9:
    simulationExecutable = str(determineExecutable(paramSimConfig))
    command = simulationExecutable + str(1)
    # add config
    command +=  " --config="+paramSimConfig
    # add simulation type (default "coupled_cell")
    command += " --simulation_type=environment_cell"
    # add additional commands to override config
    command += " --simulation_timestep=1e-3"
    command += " --sampling_timestep=1e-1"

    # add simulation to the list
    command_list.append(command)

for paramSimConfig in sweep10:
    simulationExecutable = str(determineExecutable(paramSimConfig))
    command = simulationExecutable + str(1)
    # add config
    command +=  " --config="+paramSimConfig
    # add simulation type (default "coupled_cell")
    command += " --simulation_type=environment_cell"
    # add additional commands to override config
    command += " --simulation_timestep=1e-3"
    command += " --sampling_timestep=1e-1"

    # add simulation to the list
    command_list.append(command)

for paramSimConfig in sweep11:
    simulationExecutable = str(determineExecutable(paramSimConfig))
    command = simulationExecutable + str(1)
    # add config
    command +=  " --config="+paramSimConfig
    # add simulation type (default "coupled_cell")
    command += " --simulation_type=environment_cell"
    # add additional commands to override config
    command += " --simulation_timestep=1e-3"
    command += " --sampling_timestep=1e-1"

    # add simulation to the list
    command_list.append(command)


for paramSimConfig in sweep12:
    simulationExecutable = str(determineExecutable(paramSimConfig))
    command = simulationExecutable + str(1)
    # add config
    command +=  " --config="+paramSimConfig
    # add simulation type (default "coupled_cell")
    command += " --simulation_type=environment_cell"
    # add additional commands to override config
    command += " --simulation_timestep=1e-3"
    command += " --sampling_timestep=1e-1"

    # add simulation to the list
    command_list.append(command)

# use `count' no of processes 
count = multiprocessing.cpu_count()# for the number of cores on your machine

# generate a pool of workers
pool = multiprocessing.Pool(processes=count)

# ... and pass the list of bash commands to the pool
pool.map(execute_command, command_list)


