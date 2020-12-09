import multiprocessing
import os
import subprocess

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

def UpdateConfig(problem_dim,space_dim=2,element_dim=2):

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

def EditCppAndCompile(problem_dim,space_dim=2,element_dim=2):

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


    # new .cpp formed, make 
    mod_exe = ("ChemChaste"+"_"+element_dim+"_"+space_dim+"_"+problem_dim).strip()
    exe_string = "cd ~/lib && sudo make -j4 "+mod_exe
    os.system("cd ~/lib && sudo cmake ~/src")
    os.system(exe_string)
    os.system("cd ~/projects/ChemChaste")

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
            print("found number_of_reaction_pdes ="+str(pde_number))
        elif space_substring in line:
            space_dim = line.split("=")[1].strip()
            print("found spatial_dimensions ="+str(space_dim))
        elif element_substring in line:
            element_dim = line.split("=")[1].strip()
            print("found FE_element_dimension ="+str(element_dim))

    config.close() 

    executableName = 'hello'
    if TestDimensions(pde_number,space_dim,element_dim):
        print("Template dimensions known")
    else:
        print("Makeing new executable app")
        EditCppAndCompile(pde_number,space_dim,element_dim)
        UpdateConfig(pde_number,space_dim,element_dim)

    executableName = "cd ~/lib/projects/ChemChaste/apps/ && ./ChemChaste"+"_"+str(element_dim)+"_"+str(space_dim)+"_"+str(pde_number)+" --ID "
    
    return executableName


