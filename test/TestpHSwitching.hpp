#ifndef TESTPHSWITCHING_HPP_
#define TESTPHSWITCHING_HPP_

// chaste includes
#include "ComplexCell.hpp"
#include "ChasteHeaders.hpp"

// general includes
#include "GeneralHeaders.hpp"

//ChemChaste includes
#include "ChemChasteHeaders.hpp"

// test specific include
#include "ChemicalStructuresForTests.hpp"

//#include "Cell_virtual.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "BoundaryConditionsContainer_extended.hpp"
#include "ChemChasteFeAssemblerCommon.hpp"
#include "ChemChasteVolumeAssembler.hpp"
#include "ChemChasteSurfaceAssembler.hpp"


#include <cxxtest/TestSuite.h>

//#include <cxxtest/GlobalFixture.h>
#include "CheckpointArchiveTypes.hpp"
#include "SmartPointers.hpp"
#include <string>
#include <fstream>
#include <iostream>
#include <chrono>

//#include "StateSwitchingCellProperty.hpp"

//#include "StateSwitchingCellPropertyFromFile.hpp"


//ChemChaste includes
#include "ChemChasteExecutableHeaders.hpp"

#include "ChemicalCellFromFile.hpp"
#include "CustomCellFromFile.hpp"

struct CellSwitchingStruct {
    std::string dataFileRoot = "/home/chaste/projects/ChemChaste/src/Data/StateSwitchingModel/Cell/";
    std::string cellConfig = "cell_configuration.txt";
    std::string cellEnvironment = "Environmnet.csv";
    std::string initialConditionFilename = "InitialCellConcentrations.csv";
    std::string membranePropertyFilename = "MembraneReactions.txt";
    std::string transportPropertyFilename = "TransportReactions.txt";
    std::string environmentPropertyFilename = "Environment.csv";
    std::string srnFilename = "Srn.txt";
    std::string cellCycleFilename = "SpeciesThreshold.csv";
    std::string cellDivisionRules = "SpeciesDivisionRules.csv";

    std::string cellStatesCSV = "CellStates.csv";
    std::string stateSwitchFileName = "StateSwitchingConditions.csv";

    std::string cellName = "CellA";

    std::vector<std::string> cellTypeNames = {"state_high_pH","state_low_pH"};

} stateSwitch;

class TestpHSwitching : public AbstractCellBasedTestSuite
{
public:

    void TestStateSwitchingPropertyFromFile()
    {
/*        
        // make cell
        std::string srnState = "state_high_pH/"+stateSwitch.srnFilename;

        std::string particularCellDirectory = stateSwitch.dataFileRoot + stateSwitch.cellName + "/";

        ComplexCellFromFile* p_cell_reader = new ComplexCellFromFile(
                                    particularCellDirectory,
                                    stateSwitch.cellCycleFilename, 
                                    stateSwitch.cellDivisionRules, 
                                    srnState, 
                                    stateSwitch.initialConditionFilename, 
                                    stateSwitch.transportPropertyFilename, 
                                    stateSwitch.membranePropertyFilename, 
                                    stateSwitch.environmentPropertyFilename,
                                    1,
                                    true
                                    );

        p_cell_reader->SetUp();

        std::cout<<"has membrane: "<<p_cell_reader -> GetCellPtr() -> rGetCellPropertyCollection().HasProperty<MembraneCellProperty>()<<std::endl;
        std::cout<<"has transport: "<<p_cell_reader -> GetCellPtr() -> rGetCellPropertyCollection().HasProperty<TransportCellProperty>()<<std::endl;


        StateSwitchingCellPropertyFromFile* pStateSwitchingFromFile = new StateSwitchingCellPropertyFromFile(
                                    particularCellDirectory,
                                    stateSwitch.cellStatesCSV,
                                    stateSwitch.stateSwitchFileName
        );

        //StateSwitchingCellProperty* pStateProperty = new StateSwitchingCellProperty();
        boost::shared_ptr<StateSwitchingCellProperty> pStateProperty( new StateSwitchingCellProperty());
    
        std::cout<<"TestPh - here"<<std::endl;
        CellPropertyCollection p_collection = p_cell_reader -> GetCellPtr() -> rGetCellPropertyCollection();
    
        p_cell_reader -> GetCellPtr() -> AddCellProperty(pStateProperty);

        std::cout<<"TestPh - here2"<<std::endl;

        pStateSwitchingFromFile -> SetUpStateSwitchingProperty(p_cell_reader -> GetCellPtr());

        std::cout<<"TestPh - here3"<<std::endl;
        pStateProperty->DetermineState();

        std::cout<<"TestPh - here4"<<std::endl;

        pStateProperty->SetUpSwitchedProperties();

        std::cout<<"TestPh - change environment"<<std::endl;

        boost::shared_ptr<EnvironmentCellProperty> pEnvironmentCellProperty = boost::static_pointer_cast<EnvironmentCellProperty>(p_cell_reader -> GetCellPtr()->rGetCellPropertyCollection().GetPropertiesType<EnvironmentCellProperty>().GetProperty());

        pEnvironmentCellProperty -> SetEnvironmentValueByName("pH", 9.0);

        pStateProperty->DetermineState();

        pStateProperty->SwitchState();

        std::cout<<"#####################################################################"<<std::endl;
*/
    }

    void TestCustomCellFromFile()
    {
/*
    std::string dataFileRoot = "/home/chaste/projects/ChemChaste/src/Data/StateSwitchingModel/Cell/";
    std::string cellConfig = "cell_configuration.txt";
    std::string cellEnvironment = "Environmnet.csv";
    std::string initialConditionFilename = "InitialCellConcentrations.csv";
    std::string membranePropertyFilename = "MembraneReactions.txt";
    std::string transportPropertyFilename = "TransportReactions.txt";
    std::string environmentPropertyFilename = "Environment.csv";
    std::string srnFilename = "Srn.txt";
    std::string cellCycleFilename = "SpeciesThreshold.csv";
    std::string cellDivisionRules = "SpeciesDivisionRules.csv";

    std::string cellStatesCSV = "CellStates.csv";
    std::string stateSwitchFileName = "StateSwitchingConditions.csv";

    std::string cellName = "CellA";
        
        CustomCellFromFile* p_cell_reader = new CustomCellFromFile(
            stateSwitch.dataFileRoot,
            stateSwitch.cellName,
            stateSwitch.cellCycleFilename, 
            stateSwitch.cellDivisionRules, 
            "Srn.txt",
            stateSwitch.initialConditionFilename, 
            stateSwitch.transportPropertyFilename, 
            stateSwitch.membranePropertyFilename, 
            stateSwitch.environmentPropertyFilename,
            stateSwitch.cellStatesCSV,
            "cell_configuration.txt",
            1,
            true                                                                    
            );

        p_cell_reader->SetUp();

        std::cout<<"has membrane property: "<<p_cell_reader -> GetCellPtr() -> rGetCellPropertyCollection().HasProperty<MembraneCellProperty>()<<std::endl;
        std::cout<<"has transport property: "<<p_cell_reader -> GetCellPtr() -> rGetCellPropertyCollection().HasProperty<TransportCellProperty>()<<std::endl;
        std::cout<<"has environment property: "<<p_cell_reader -> GetCellPtr() -> rGetCellPropertyCollection().HasProperty<EnvironmentCellProperty>()<<std::endl;
        std::cout<<"has state switch property: "<<p_cell_reader -> GetCellPtr() -> rGetCellPropertyCollection().HasProperty<StateSwitchingCellProperty>()<<std::endl;

        boost::shared_ptr<StateSwitchingCellProperty> pStateProperty = boost::static_pointer_cast<StateSwitchingCellProperty>(p_cell_reader -> GetCellPtr()->rGetCellPropertyCollection().GetPropertiesType<StateSwitchingCellProperty>().GetProperty());


        std::vector<std::string> stateNames = pStateProperty -> GetStateNames();

        std::cout<<"State names:"<<std::endl;
        for(unsigned i=0; i<stateNames.size(); i++)
        {
            std::cout<<stateNames[i]<<std::endl;
        }

        std::cout<<"#################### Cell division ########################"<<std::endl;

        CellPtr p_new_cell = p_cell_reader -> GetCellPtr() -> Divide();

*/
        
    }

};

#endif