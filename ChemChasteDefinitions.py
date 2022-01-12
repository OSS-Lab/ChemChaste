import multiprocessing
import os
import subprocess
import csv
from distutils.dir_util import copy_tree


# This is a helper function that allows to run a bash command in a separate process
def execute_command(cmd):
    return subprocess.call(cmd, shell=True)

def TestDimensions(problem_dim,space_dim=2,element_dim=2):

    # cast dimensions if not already strings
    problem_dim = str(problem_dim)
    space_dim = str(space_dim)
    element_dim = str(element_dim)
    
    config_file = open("/home/chaste/projects/ChemChaste/DataInput/ChemChaste_configuration.txt", 'r') 
    pde_substring = "domain_problem_dimensions ="
    space_substring = "domain_space_dimensions ="
    element_substring = "domain_FE_element_dimensions ="

    while True: 

        # Get next line from file 
        line = config_file.readline() 
        # if line is empty 
        # end of file is reached 
        if not line: 
            break
        line=line.strip()

        if pde_substring in line:
            pde_number_container = line.split("=")[1].strip()
            pde_number_list = pde_number_container.split(',')
            if problem_dim not in pde_number_list:
                print("not found "+ problem_dim + " in "+ pde_number_container)
                return False
        elif space_substring in line:
            space_number_container = line.split("=")[1].strip()
            space_number_list = space_number_container.split(',')
            if space_dim not in space_number_list:
                print("found "+ space_dim + " in "+ space_number_container)
                return False
        elif element_substring in line:
            element_number_container = line.split("=")[1].strip()
            element_number_list = element_number_container.split(',')
            if element_dim not in element_number_list:
                print("found "+ element_dim + " in "+ element_number_container)
                return False

    return True

def UpdateConfigDimensions(problem_dim,space_dim=2,element_dim=2):

    # cast dimensions if not already strings
    problem_dim = str(problem_dim)
    space_dim = str(space_dim)
    element_dim = str(element_dim)
    
    config_filename = "/home/chaste/projects/ChemChaste/DataInput/ChemChaste_configuration.txt"

    config_file = open(config_filename, 'r') 
    pde_substring = "domain_problem_dimensions ="
    space_substring = "domain_space_dimensions ="
    element_substring = "domain_FE_element_dimensions ="

    mod_config_file_content = ""
    for line in config_file:

        line=line.strip()

        if pde_substring in line:
            pde_number_container = line.split("=")[1].strip()
            pde_number_list = pde_number_container.split(',')
            if problem_dim not in pde_number_list:
                newLine = line+","+problem_dim
            else:
                newLine = line   
        elif space_substring in line:
            space_number_container = line.split("=")[1].strip()
            space_number_list = space_number_container.split(',')
            if space_dim not in space_number_list:
                newLine = line+","+space_dim
            else:
                newLine = line
        elif element_substring in line:
            element_number_container = line.split("=")[1].strip()
            element_number_list = element_number_container.split(',')
            if element_dim not in element_number_list:
                newLine = line+","+element_dim
            else:
                newLine = line
        else:
            newLine = line

        mod_config_file_content += newLine +"\n"


    config_file.close()


    config_file = open(config_filename, "w")

    config_file.write(mod_config_file_content)

    config_file.close()

    return True

def TestCellVirtual(config_file):

    original_cell_str = "~Cell();"
    virtual_cell_str  = "virtual ~Cell();"

    cell_file_loc = "/home/chaste/src/cell_based/src/cell/Cell.hpp"

    cell_file = open(cell_file_loc, 'r') 

    for line in cell_file:

        line=line.strip()

        if original_cell_str in line:
            if virtual_cell_str in line:
                return True
            else:
                return False

    return False


def EditCppVirtual():

    original_cell_str = "~Cell();"
    virtual_cell_str  = "virtual ~Cell();"

    original_divide_str = "CellPtr Divide();"
    virtual_divide_str = "virtual CellPtr Divide();"

    original_protected_str = "protected:"
    mCanDivide_str  = "bool mCanDivide;"
    placeNextLine = False
    
    cell_file_loc = "/home/chaste/src/cell_based/src/cell/Cell.hpp"
    temp_file_loc = "/home/chaste/src/cell_based/src/cell/CellTemp.hpp"

    # do we need to also move mCanDivide from private to protected?

    cell_file = open(cell_file_loc, 'r')

    temp_file_content = ""
    for line in cell_file:

        line=line.strip()

        if original_cell_str in line:
            newLine = virtual_cell_str
        elif original_divide_str in line:
            newLine = virtual_divide_str
        elif mCanDivide_str in line:
            # assume mCanDivide_str comes before original_protected_str
            newLine = ""
        elif original_protected_str in line:
            # place mCanDivide_str on next line
            newLine = line
            placeNextLine = True
        else:
            if(placeNextLine):
                newLine = mCanDivide_str
                placeNextLine = False
            newLine = line

        temp_file_content += newLine +"\n"

    cell_file.close()


    temp_file = open(temp_file_loc, "w")

    temp_file.write(temp_file_content)

    temp_file.close()

    # delete old cell file
    os.remove(cell_file_loc)

    # rename temp file
    os.rename(temp_file_loc,cell_file_loc)


    return True

def EditCppDimensions(problem_dim,space_dim=2,element_dim=2):

    # cast dimensions if not already a string
    problem_dim = str(problem_dim).strip()
    space_dim = str(space_dim).strip()
    element_dim = str(element_dim).strip()

    config_file = "/home/chaste/projects/ChemChaste/DataInput/ChemChaste_configuration.txt"

    substring = "ChemChaste_cpp_location ="
    app_location = ""
    config = open(config_file, 'r') 

    while True: 
        # Get next line from file 
        line = config.readline() 
        # if line is empty 
        # end of file is reached 
        if not line: 
            break
        line=line.strip()

        if substring in line:
            print()
            app_location = line.split("=")[1].strip()
            print("cpp loc ="+app_location)

    config.close() 

    exe_filename = app_location
    substring = "ChemChaste_cpp_location ="
    exe_file = open(exe_filename, 'r')

    mod_exe_filename = exe_filename.split(".cpp")[0].strip()+"_"+element_dim+"_"+space_dim+"_"+problem_dim+".cpp"

    problem_dim_string = "const unsigned probDim ="
    space_dim_string = "const unsigned spaceDim="
    element_dim_string = "const unsigned elementDim="

    mod_file_content = ""
    for line in exe_file:

        line=line.strip()

        if problem_dim_string in line:
            newLine = problem_dim_string+problem_dim+";"
        elif space_dim_string in line:
            newLine = space_dim_string+space_dim+";"
        elif element_dim_string in line:
            newLine = element_dim_string+element_dim+";"
        else:
            newLine = line

        mod_file_content += newLine +"\n"


    exe_file.close()


    mod_exe_file = open(mod_exe_filename, "w")

    mod_exe_file.write(mod_file_content)

    mod_exe_file.close()

    return True

def CompileCPP(problem_dim,space_dim=2,element_dim=2):
    # new .cpp formed, make 
    mod_exe = ("ChemChaste"+"_"+element_dim+"_"+space_dim+"_"+problem_dim).strip()
    exe_string = "cd /home/chaste/lib && sudo make -j4 "+mod_exe
    os.system("cd /home/chaste/lib && sudo cmake /home/chaste/src")
    os.system(exe_string)
    os.system("cd /home/chaste/projects/ChemChaste")

    return True

def RetriveConfigLocation():

    return location

def determineExecutable(config_file):

    #with open(config_file) as file:
    #    file_contents = file.read()
    #    print(file_contents)
    pde_substring = "number_of_reaction_pdes ="
    space_substring = "spatial_dimensions ="
    element_substring = "FE_element_dimension ="

    toCompile=False

    pde_number = 0
    space_dim = 2
    element_dim = 2
    config = open(config_file, 'r') 
  
    while True: 
        # Get next line from file 
        line = config.readline() 
        # if line is empty 
        # end of file is reached 
        if not line: 
            break
        line=line.strip()

        if pde_substring in line:
            pde_number = line.split("=")[1].strip()
            #print("found number_of_reaction_pdes ="+str(pde_number))
        elif space_substring in line:
            space_dim = line.split("=")[1].strip()
            #print("found spatial_dimensions ="+str(space_dim))
        elif element_substring in line:
            element_dim = line.split("=")[1].strip()
            #print("found FE_element_dimension ="+str(element_dim))

    config.close() 

    executableName = 'hello'

    if TestDebugMode(config_file):
        # recompile anyway
        print("Debugging")
        toCompile=True

    else:
        if TestDimensions(pde_number,space_dim,element_dim):
            print("Template dimensions known")
        else:
            print("Makeing new executable app")
            EditCppDimensions(pde_number,space_dim,element_dim)
            UpdateConfigDimensions(pde_number,space_dim,element_dim)
            toCompile=True

        if not TestCellVirtual(config_file):
            EditCppVirtual()
            toCompile=True

    if toCompile:
        CompileCPP(pde_number,space_dim,element_dim)


    executableName = "cd /home/chaste/lib/projects/ChemChaste/apps/ && ./ChemChaste"+"_"+str(element_dim)+"_"+str(space_dim)+"_"+str(pde_number)+" --ID "
    
    return executableName


def TestDebugMode(config_filename):
    config_filename = "/home/chaste/projects/ChemChaste/DataInput/ChemChaste_configuration.txt"
    config_file = open(config_filename, 'r') 
    debug_substring = "debug_mode ="
    
    while True: 
        # Get next line from file 
        line = config_file.readline() 
        # if line is empty 
        # end of file is reached 
        if not line: 
            break
        line=line.strip()
        print(line)
        if debug_substring in line:
            debug_value = str(line.split("=")[1].strip())
            print(debug_value)
            if debug_value ==str(1):
                return True
            else:
                return False
    return False


def DetermineConfigFileName(simulation_config,clean_parameter_names,parameterSet):
    
    configFilename = simulation_config.split('.')[0]
    fileCode=""
    for i in range(len(clean_parameter_names)):
        fileCode=fileCode+"_"+clean_parameter_names[i]+"_"+parameterSet[i]
    
    
    configFilename=configFilename+fileCode+'.'+simulation_config.split('.')[1]
    
    fileCode = fileCode[1:] 
    
    return [configFilename, fileCode]

def CreateNewConfigFile(simulation_config,clean_parameter_names,parameterSet,configFilename,fileCode):
    
    simulationConfig= open(simulation_config,'r')
    newConfig = open(configFilename,'w')

    for line in simulationConfig:
        line=line.strip()

        if 'output_filename' in line:
            newLine = line + '/'+ fileCode
        elif 'domain_file_root' in line:
            splitLine1 = line.rsplit('/',3)
            splitLine2 = line.rsplit('/',2)
            newLine=splitLine1[0]+'/'+fileCode+'/'+splitLine2[1]+'/'
        elif 'cell_file_root' in line:
            splitLine1 = line.rsplit('/',3)
            splitLine2 = line.rsplit('/',2)
            newLine=splitLine1[0]+'/'+fileCode+'/'+splitLine2[1]+'/'
        else:
            newLine=line
        
        newConfig.write(newLine+'\n')
        newLine=""
        
    simulationConfig.close()
    newConfig.close()


    return configFilename


def UpdateFile(filePath,parameterNamesList,parameterSet):
    
    newFileString=""
    
    newFile = open(filePath,'r')
    
    for line in newFile:
        line=line.strip()
        
        for i in range(len(parameterNamesList)):

            para_name= parameterNamesList[i]
            para_value=parameterSet[i]

            if para_name in line:
                line = line.replace(para_name, para_value)
        
        newFileString += line +"\n"
        
    newFile.close()
    
    newFile = open(filePath,'w')
    newFile.write(newFileString)
    newFile.close()
    
    return 

def CreateNewSimulationFiles(parameterNamesList,parameterSet,simulationDirectory, fileCode):
    
    newDirectoryRoot = simulationDirectory.rsplit('/',1)[0] + '/' + fileCode
    
    print("SimulationDirectory: ",simulationDirectory)
    print("newDirectoryRoot: ",newDirectoryRoot)

    copy_tree(simulationDirectory, newDirectoryRoot)

    for root, dirs, files in os.walk(newDirectoryRoot):
        for file in files:

            UpdateFile(os.path.join(root, file),parameterNamesList,parameterSet)
    
    return 

def ParseSweepingFile(sweeping_filename):

    parameterNames = []
    parameterValues =[]

    sweep_file = open(sweeping_filename, 'r') 
    csv_reader = csv.reader(sweep_file, delimiter=',')

    for lineList in csv_reader:
        #print(lineList)
        if lineList:
            if not lineList[0].count('#'): 
                parameterNames.append(lineList[0].strip())
                del lineList[0]

                lineValuesList=[]

                for entry in lineList:
                    entry=entry.strip()
                    testRangeCharacters = entry.count(':')

                    if testRangeCharacters==0:
                        # single value given
                        lineValuesList.append(entry.strip())

                    elif testRangeCharacters==2:
                        # first entry, stepsize, last entry
                        entryList=entry.split(':')
                        firstValue=float(entryList[0].strip())
                        stepSize=float(entryList[1].strip())
                        lastValue=float(entryList[2].strip())

                        lineValuesList.append(str(firstValue))
                        numericParameter=firstValue+stepSize
                        while numericParameter<=lastValue:
                            lineValuesList.append(str(numericParameter))
                            numericParameter=numericParameter+stepSize

                    elif(testRangeCharacters==1):
                        # first entry, last entry, provide an intermediate value as a third parameter
                        entryList=entry.split(':')
                        firstValue=float(entryList[0].strip())
                        lastValue=float(entryList[1].strip())
                        intermediateValue=firstValue + 0.5*(lastValue-firstValue)

                        lineValuesList.append(str(firstValue))
                        lineValuesList.append(str(intermediateValue))
                        lineValuesList.append(str(lastValue))



                    else:
                        # error
                        print('ParseSweepingFile: Error parsing'+entry )
                parameterValues.append(lineValuesList)


    sweep_file.close()


    return [parameterNames,parameterValues]

def ConstructParameterSets(parameterNamesList,valuesListOfLists):

    listOfParameterSets=[]

    depth = len(valuesListOfLists)-1
    startingDepth=depth-1
    container = []

    while startingDepth >=0:
        #print(valuesListOfLists[startingDepth])
        for i in range(len(valuesListOfLists[startingDepth])):
            entry = valuesListOfLists[startingDepth][i]
            for entryLower in valuesListOfLists[startingDepth+1]:
                container.append(entry+','+entryLower)

        valuesListOfLists[startingDepth]= container

        if startingDepth != 0:
            container =[]

        startingDepth=startingDepth-1

    #print("Container: ",container)
 
    listOfParameterSets=container

    return listOfParameterSets


def ParameterSweeping(simulation_config):

    parametersFilename = ""

    simulationDirectory= ""

    configFile = open(simulation_config, 'r')
    for line in configFile:
            line=line.strip()

            if 'parameter_filename' in line:
                parametersFilename= line.split("=")[1].strip()

    configFile.close()

    simulationDirectory = parametersFilename.rsplit('/',1)[0]

    print("parametersFilename: ",parametersFilename)
    print("simulationDirectory: ",simulationDirectory)
    
    #parametersFilename  = simulationDirectory + 'Parameters.csv'
    # read parameter sweeping file and parse parameters into a names 
    # list and values lists of lists

    [parameterNamesList,valuesListOfLists] = ParseSweepingFile(parametersFilename)

    clean_parameter_names = []
    
    for entry in parameterNamesList:
        
        clean_parameter_names.append(entry.replace('&',''))
    
    
    ListOfSimulationParameters = ConstructParameterSets(parameterNamesList,valuesListOfLists)

    print(ListOfSimulationParameters)

    config_list = []
    
    for parameterSet in ListOfSimulationParameters:

        parameterSetSplit = parameterSet.split(',')

        print("parameter names: ",clean_parameter_names)
        print("parameter set: ",parameterSetSplit)
        [configFilename, fileCode] = DetermineConfigFileName(simulation_config,clean_parameter_names,parameterSetSplit)

        newConfigFilename = CreateNewConfigFile(simulation_config,clean_parameter_names,parameterSetSplit,configFilename, fileCode)

        print("New config name: ",newConfigFilename)
        config_list.append(newConfigFilename)
        CreateNewSimulationFiles(parameterNamesList,parameterSetSplit,simulationDirectory, fileCode)

    

    return config_list
