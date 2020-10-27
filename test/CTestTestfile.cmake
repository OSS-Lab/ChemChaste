# CMake generated Testfile for 
# Source directory: /home/chaste/src/projects/ChemicalChaste/test
# Build directory: /home/chaste/projects/ChemicalChaste/test
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(TestHello "/home/chaste/projects/ChemicalChaste/test/TestHello")
set_tests_properties(TestHello PROPERTIES  LABELS "Continuous_project_ChemicalChaste" PROCESSORS "1" WORKING_DIRECTORY "/home/chaste/src/")
add_test(TestChemicalChaste "/home/chaste/projects/ChemicalChaste/test/TestChemicalChaste")
set_tests_properties(TestChemicalChaste PROPERTIES  LABELS "Continuous_project_ChemicalChaste" PROCESSORS "1" WORKING_DIRECTORY "/home/chaste/src/")
