

#ifndef STATESWITCHINGCELLPROPERTY_HPP
#define STATESWITCHINGCELLPROPERTY_HPP

//class ChemicalSrnModel;

//general includes
#include <vector>
#include <boost/shared_ptr.hpp>

#include "ChemicalSRNFromFile.hpp"
#include "AbstractTransportReactionSystemFromFile.hpp"
#include "AbstractMembraneReactionSystemFromFile.hpp"
#include "PreferredCellEnvironmentFromFile.hpp"
#include "SimpleChemicalThresholdCellCycleFromFile.hpp"

#include "ChemicalCellProperty.hpp"

// chaste includes
#include "Cell.hpp"
//#include "AbstractCellProperty.hpp"
//#include "StateVariableRegister.hpp"
//#include "RandomNumberGenerator.hpp"
//#include "AbstractChemistry.hpp"
//#include "ChastePoint.hpp"


//class ComplexCell;


class StateSwitchingCellProperty : public AbstractCellProperty
{
protected:

    // CellPtr to a given cell for access to cell data and properties
    CellPtr mThis_cellPtr;

    // Cell state
    unsigned mCellState;

    std::string mBaseDirectoryName;

    std::vector<bool> mStateActive;

    std::vector<std::string> mStateNames;

    std::vector<std::vector<std::string>> mStateFiles;


    // conditions that determine what the state is
    std::vector<std::vector<std::string>> mStateConditionsChemicalNames;

    std::vector<std::vector<std::string>> mStateConditionTypes;

    // whether there has been a maximum value set, minimum will default to 0.0
    std::vector<std::vector<bool>> mIsStateConditionsMaxValueSet;

    std::vector<std::vector<bool>> mIsStateConditionsMinValueSet;

    std::vector<std::vector<double>> mStateConditionsMaxValue;

    std::vector<std::vector<double>> mStateConditionsMinValue;

    SimpleChemicalThresholdCellCycleModel*  mpSimpleChemicalThresholdCellCycleModel;

    ChemicalSrnModel*   mpChemicalSrnModel;



public:

    StateSwitchingCellProperty();

    virtual ~StateSwitchingCellProperty(){};

    StateSwitchingCellProperty(const StateSwitchingCellProperty&);

    // virtual methods

    virtual void SetUp(CellPtr);

    virtual void PreparePostDivisionParent(double);
    
    virtual void PreparePostDivisionDaughter(const StateSwitchingCellProperty&, double);

    virtual void SetUpSwitchedProperties();

    // concrete methods

    void DetermineState();

    void SwitchState();

    std::string RetrieveStateName();

    std::vector<std::string> RetrieveStateFiles();

    std::vector<std::vector<std::string>> ReadMatrix(std::string);

    std::vector<std::string> parseMatrixLineString(std::string);

    // set methods

    void SetCellState(unsigned);

    void SetCellPtr(CellPtr);

    void SetBaseDirectoryName(std::string);

    void SetStateNames(std::vector<std::string>);

    void SetStateFiles(std::vector<std::vector<std::string>>);

    

    void SetConditionChemicalNames(std::vector<std::vector<std::string>>);

    void SetConditionTypes(std::vector<std::vector<std::string>>);

    void SetIsStateConditionsMaxSet(std::vector<std::vector<bool>>);

    void SetIsStateConditionsMinSet(std::vector<std::vector<bool>>);

    void SetStateConditionsMax(std::vector<std::vector<double>>);

    void SetStateConditionsMin(std::vector<std::vector<double>>);
    
    void SetChemicalSrnModel(ChemicalSrnModel*);

    void SetChemicalCellCycleModel(SimpleChemicalThresholdCellCycleModel*);

    // get methods

    unsigned GetCellState();

    CellPtr GetCellPtr();

    std::string GetBaseDirectoryName();

    std::vector<std::string> GetStateNames();

    std::vector<std::vector<std::string>> GetStateFiles();

    std::vector<std::vector<std::string>> GetConditionChemicalNames();

    std::vector<std::vector<std::string>> GetConditionTypes();

    std::vector<std::vector<bool>> GetIsStateConditionsMaxSet();

    std::vector<std::vector<bool>> GetIsStateConditionsMinSet();

    std::vector<std::vector<double>> GetStateConditionsMax();

    std::vector<std::vector<double>> GetStateConditionsMin();

    std::string GetStateNamesByIndex(unsigned);

    std::vector<std::string> GetStateFilesByIndex(unsigned);

    std::vector<std::string> GetConditionChemicalNamesByIndex(unsigned);

    std::vector<std::string> GetConditionTypesByIndex(unsigned);

    std::vector<bool> GetIsStateConditionsMaxSetByIndex(unsigned);

    std::vector<bool> GetIsStateConditionsMinSetByIndex(unsigned);

    std::vector<double> GetStateConditionsMaxByIndex(unsigned);

    std::vector<double> GetStateConditionsMinByIndex(unsigned);

    ChemicalSrnModel* GetChemicalSrnModel();

    SimpleChemicalThresholdCellCycleModel* GetChemicalCellCycleModel();
};

StateSwitchingCellProperty::StateSwitchingCellProperty()
    : AbstractCellProperty()
{
    std::cout<<"StateSwitchingCellProperty::StateSwitchingCellProperty() - start"<<std::endl;
    std::cout<<"StateSwitchingCellProperty::StateSwitchingCellProperty() - end"<<std::endl;
}

StateSwitchingCellProperty::StateSwitchingCellProperty(const StateSwitchingCellProperty& existingProperty)
{
    mCellState = existingProperty.mCellState;

    mBaseDirectoryName = existingProperty.mBaseDirectoryName;
    mStateActive = existingProperty.mStateActive;
    mStateNames = existingProperty.mStateNames;

    mStateFiles = existingProperty.mStateFiles;

    mStateConditionsChemicalNames = existingProperty.mStateConditionsChemicalNames;

    mStateConditionTypes = existingProperty.mStateConditionTypes;

    mIsStateConditionsMaxValueSet = existingProperty.mIsStateConditionsMaxValueSet;

    mIsStateConditionsMinValueSet = existingProperty.mIsStateConditionsMinValueSet;;

    mStateConditionsMaxValue = existingProperty.mStateConditionsMaxValue;

    mStateConditionsMinValue = existingProperty.mStateConditionsMinValue;
}

void StateSwitchingCellProperty::SetUp(CellPtr)
{
    std::cout<<"StateSwitchingCellProperty::SetUp() - start"<<std::endl;
    this -> SetUpSwitchedProperties();
    std::cout<<"StateSwitchingCellProperty::SetUp() - end"<<std::endl;
}

void StateSwitchingCellProperty::PreparePostDivisionParent(double splitRatio)
{
    // split any properties that are shared
}
    
void  StateSwitchingCellProperty::PreparePostDivisionDaughter(const StateSwitchingCellProperty& parentProperty, double splitRatio)
{
    // split any properties that are shared
}

void StateSwitchingCellProperty::SetUpSwitchedProperties()
{
    std::cout<<"StateSwitchingCellProperty::SetUpSwitchedProperties() - start"<<std::endl;
    std::vector<std::string> stateFiles = GetStateFilesByIndex(mCellState);

    std::string stateName = GetStateNamesByIndex(mCellState);

    std::string directoryName = mBaseDirectoryName+"/"+stateName+"/";

    std::string fileName="";
    std::string fileNameFull="";

    
    // properties should only be changed if new stage wants the change
    for(unsigned fileIndex=0; fileIndex<stateFiles.size(); fileIndex++)
    {
        fileName = stateFiles[fileIndex];

        std::cout<< "For cell state: "<<mCellState<<" mutable file: "<<fileName<<std::endl;
        if(fileName=="Srn.txt")
        {
            fileNameFull=directoryName+fileName;
            std::cout<<"Srn file: "<<fileNameFull<<std::endl;

            ChemicalSRNFromFile* p_srn_reaction_system_from_file = new ChemicalSRNFromFile(fileNameFull);

            ChemicalSrnModel *p_chemical_srn = static_cast<ChemicalSrnModel*>(p_srn_reaction_system_from_file -> GetChemicalSrnModel());
      
            SetChemicalSrnModel(p_chemical_srn);
            p_chemical_srn->GetReactionSystem()->SetCell(mThis_cellPtr); // should work or else type cast to ChemicalSrnModel*
            p_chemical_srn->GetReactionSystem()->DistributeCellPtr();

        }
        else if(fileName=="TransportReactions.txt")
        {
            fileNameFull=directoryName+fileName;
            boost::shared_ptr<TransportCellProperty> p_transport = boost::static_pointer_cast<TransportCellProperty>(mThis_cellPtr->rGetCellPropertyCollection().GetPropertiesType<TransportCellProperty>().GetProperty());
            
            AbstractTransportReactionSystemFromFile* p_transport_system_from_file = new AbstractTransportReactionSystemFromFile(fileNameFull);
        
            p_transport -> SetUp(p_transport_system_from_file, mThis_cellPtr);
        
        }
        else if(fileName=="MembraneReactions.txt")
        {
            fileNameFull=directoryName+fileName;
            boost::shared_ptr<MembraneCellProperty> p_membrane = boost::static_pointer_cast<MembraneCellProperty>(mThis_cellPtr->rGetCellPropertyCollection().GetPropertiesType<MembraneCellProperty>().GetProperty());
            AbstractMembraneReactionSystemFromFile* p_membrane_system_from_file = new AbstractMembraneReactionSystemFromFile(fileNameFull);
        
            p_membrane -> SetUp(p_membrane_system_from_file, mThis_cellPtr);

        }
        else if(fileName=="SpeciesDivisionRules.csv")
        {
            fileNameFull=directoryName+fileName;

            std::vector<std::string> chemicalNames;
            std::vector<std::string> chemicalDivisionRules;

            std::vector<std::vector<std::string>> divisionDatabase;

            divisionDatabase = ReadMatrix(fileNameFull);

            unsigned numberOfChemicals = divisionDatabase.size();

            for(unsigned row=0; row<numberOfChemicals;row++)
            {
                chemicalNames.push_back(divisionDatabase[row][0]);
                chemicalDivisionRules.push_back(divisionDatabase[row][1]);
            }
           
            boost::shared_ptr<ComplexCell> p_cell_new = boost::static_pointer_cast<ComplexCell>(mThis_cellPtr);
            p_cell_new->SetChemicalNames(chemicalNames);

            p_cell_new->SetChemicalDivsionRules(chemicalDivisionRules);


        }
        else if(fileName=="Environment.csv")
        {
            fileNameFull=directoryName+fileName;
            boost::shared_ptr<EnvironmentCellProperty> p_environment = boost::static_pointer_cast<EnvironmentCellProperty>(mThis_cellPtr->rGetCellPropertyCollection().GetPropertiesType<EnvironmentCellProperty>().GetProperty());
            
            PreferredCellEnvironmentFromFile* p_preferred_environment_from_file = new PreferredCellEnvironmentFromFile(fileNameFull);

            p_environment -> SetUp(mThis_cellPtr, p_preferred_environment_from_file->GetEnvironmentNamesVector());
            p_environment -> SetPreferredEnvironmentVector(p_preferred_environment_from_file->GetEnvironmentValueVector());

        }
        else if(fileName=="SpeciesThreshold.csv")
        {
            fileNameFull=directoryName+fileName;
            SimpleChemicalThresholdCellCycleFromFile* p_cell_cycle = static_cast<SimpleChemicalThresholdCellCycleFromFile*>(const_cast <AbstractCellCycleModel*>(mThis_cellPtr->GetCellCycleModel()));
         
            p_cell_cycle->SetCellCycleFilename(fileNameFull);
            p_cell_cycle->SetUp();

        }
        else
        {
            std::cout<<"Error: StateSwitchingCellProperty::SetUpSwitchedProperties() unexpected filename "<<std::endl;
        }
    }
    std::cout<<"StateSwitchingCellProperty::SetUpSwitchedProperties() - end"<<std::endl;
}


// concrete methods

void StateSwitchingCellProperty::DetermineState()
{
    std::cout<<"StateSwitchingCellProperty::DetermineState() - start"<<std::endl;
    // current cell state
    unsigned cellState = GetCellState();

    unsigned numberOfCellStates = mStateNames.size();
std::cout<<"here0"<<std::endl;
    // use ChemicalCellProperty to check cell data for 
    assert(pCell->rGetCellPropertyCollection().HasProperty<ChemicalCellProperty>());
    boost::shared_ptr<ChemicalCellProperty> pChemicalCellProperty = boost::static_pointer_cast<ChemicalCellProperty>(mThis_cellPtr->rGetCellPropertyCollection().GetPropertiesType<ChemicalCellProperty>().GetProperty());
std::cout<<"here1"<<std::endl;
    StateVariableRegister* pStateRegister =  pChemicalCellProperty->GetStateVariableRegister();


    std::vector<bool> mStateActive(numberOfCellStates, true);
    bool isStateActive=false; 
    std::string chemicalName="";
    double this_concentration=0.0;
std::cout<<"here2"<<std::endl;
    for(unsigned stateIndex =0; stateIndex<numberOfCellStates; stateIndex++)
    {
        // run through each of the states to check whether they are potentially active
        isStateActive = true;
        std::cout<<"stateIndex: "<<stateIndex<<std::endl;
        for(unsigned conditionIndex=0; conditionIndex<mStateConditionTypes[stateIndex].size(); conditionIndex++)
        {
            // run through each condition and determine whether the condition is met
            std::cout<<"conditionIndex: "<<conditionIndex<<std::endl;
            chemicalName = GetConditionChemicalNamesByIndex(stateIndex)[conditionIndex];
            std::cout<<"Chemical name: "<<chemicalName<<" state condition: "<<mStateConditionTypes[stateIndex][conditionIndex]<<std::endl;
            if(mStateConditionTypes[stateIndex][conditionIndex]=="CellConcentration")
            {
                std::cout<<"CellConcentration"<<std::endl;
                // get concentration from cellData
                if(pStateRegister->IsStateVariablePresent(chemicalName))
                {
                    this_concentration = mThis_cellPtr->GetCellData()->GetItem(chemicalName);

                    if(GetIsStateConditionsMaxSetByIndex(stateIndex)[conditionIndex])
                    {
                        if(this_concentration>GetStateConditionsMaxByIndex(stateIndex)[conditionIndex])
                        {
                            isStateActive = false;
                            break;
                        }
                    }
                    if(GetIsStateConditionsMinSetByIndex(stateIndex)[conditionIndex])
                    {
                        if(this_concentration<GetStateConditionsMinByIndex(stateIndex)[conditionIndex])
                        {
                            isStateActive = false;
                            break;
                        } 
                    }
                }
            }
            else if(mStateConditionTypes[stateIndex][conditionIndex]=="EnvironmentConcentration")
            {
                std::cout<<"EnvironmentConcentration"<<std::endl;
                if(mThis_cellPtr->rGetCellPropertyCollection().HasProperty<EnvironmentCellProperty>())
                {
                    std::cout<<"cell has environment"<<std::endl;
                    boost::shared_ptr<EnvironmentCellProperty> pEnvironmentCellProperty = boost::static_pointer_cast<EnvironmentCellProperty>(mThis_cellPtr->rGetCellPropertyCollection().GetPropertiesType<EnvironmentCellProperty>().GetProperty());

                    this_concentration = pEnvironmentCellProperty->GetEnvironmentValueByName(chemicalName);

                    // check the environment concentration is within the bounds
                    if(GetIsStateConditionsMaxSetByIndex(stateIndex)[conditionIndex])
                    {
                        std::cout<<"env conc. "<<this_concentration<<" max condition: "<<GetStateConditionsMaxByIndex(stateIndex)[conditionIndex]<<std::endl;
                        if(this_concentration>GetStateConditionsMaxByIndex(stateIndex)[conditionIndex])
                        {
                            isStateActive = false;
                            break;
                        } 
                    }
                    if(GetIsStateConditionsMinSetByIndex(stateIndex)[conditionIndex])
                    {
                        std::cout<<"env conc. "<<this_concentration<<" min condition: "<<GetStateConditionsMinByIndex(stateIndex)[conditionIndex]<<std::endl;
                        if(this_concentration<GetStateConditionsMinByIndex(stateIndex)[conditionIndex])
                        {
                            isStateActive = false;
                            break;
                        } 
                    }
                }
            }
            else
            {
                std::cout<<"StateSwitchingCellProperty::DetermineState() state condition not recognised"<<std::endl;
            }
            
        }
        std::cout<<"here2"<<std::endl;
        mStateActive[stateIndex] = isStateActive;
        std::cout<<"Is state active: "<<isStateActive<<std::endl;
    }
    std::cout<<"select a state"<<std::endl;
    // for each active state determine the selected state (random if more than one)
    unsigned numberOfActiveStates=0;
    unsigned lastActiveState=0;
    for(unsigned stateIndex=0; stateIndex<mStateActive.size(); stateIndex++)
    {
        if(mStateActive[stateIndex])
        {
            lastActiveState = stateIndex;
            numberOfActiveStates++;
        }
    }

    if(numberOfActiveStates==0)
    {
        // No state active, retain previous state
        std::cout<<"no active state, retain previous state"<<std::endl;
        SetCellState(cellState);
    }
    else if(numberOfActiveStates==1)
    {
        // switch to sole active state
        std::cout<<"one active state, switch state"<<std::endl;
        SetCellState(lastActiveState);
    }
    else
    {
        std::cout<<"many active states, random select"<<std::endl;
        // more than one state is possible, select state randomly
        double statePortion = 1.0/numberOfActiveStates;
        double uniformRandomNumber = RandomNumberGenerator::Instance()->ranf();

        for(unsigned stateIndex=0; stateIndex<mStateActive.size(); stateIndex++)
        {
            if(fabs(uniformRandomNumber - stateIndex*statePortion)<1e-6)
            {
                SetCellState(stateIndex);
                break;
            }
        }
    }

    std::cout<<"Active state: "<< mCellState << " "<<mStateNames[mCellState]<<std::endl;

    std::cout<<"StateSwitchingCellProperty::DetermineState() - end"<<std::endl;
}

void StateSwitchingCellProperty::SwitchState()
{
    std::cout<<"StateSwitchingCellProperty::SwitchState() - start"<<std::endl;
    this -> SetUpSwitchedProperties();
    // SetUpSwitchedProperties();

    /*
    cell_label = p_Pde_field->GetCellLabelByIndex(i);
    cell_key = p_Pde_field->ReturnCellKeyFromCellLabel(cell_label);
    numericalCellID = p_Pde_field->ReturnUnsignedIDFromCellKeyString(cell_key);
    given_cell_root = variables_map["cell_file_root"].as<std::string>()+cell_key+"/";

    CustomCellFromFile* p_cell_reader = new CustomCellFromFile(
                        given_cell_root,
                        "SpeciesThreshold.csv", 
                        "SpeciesDivisionRules.csv", 
                        "Srn.txt",
                        "InitialCellConcentrations.csv",
                        "TransportReactions.txt",
                        "MembraneReactions.txt",
                        "Environment.csv",
                        "StateSwitches.csv",
                        "cell_configuration.txt",
                        numericalCellID,
                        true                                                                    
                        );

    p_cell_reader->SetUp();

    cells.push_back(p_cell_reader -> GetCellPtr());
    */
   
    std::cout<<"StateSwitchingCellProperty::SwitchState() - end"<<std::endl;
}

std::string StateSwitchingCellProperty::RetrieveStateName()
{
    std::cout<<"here"<<std::endl;
    return mStateNames[mCellState];
}

std::vector<std::string> StateSwitchingCellProperty::RetrieveStateFiles()
{
    return mStateFiles[mCellState];
}


std::vector<std::vector<std::string>> StateSwitchingCellProperty::ReadMatrix(std::string filename)
{
    // parse a matrix file (.csv) line by line, ignore escape line,s containing file information
    // that is lines starting with '#' 
    
    std::string line;
    std::ifstream inputFile(filename);
    
    // read all data types as std::string therefore return the matrix of strings for personalised
    // methods down stream
    std::vector<std::vector<std::string>> outputMatrix = std::vector<std::vector<std::string>>();

    // check file exists and is openable
    if(inputFile.is_open()){
        // open the matrix file
        while (getline(inputFile,line)){
            // while the file still has lines not read.
            // read line left to right, top to bottom.
            if(!line.empty()){
                if(line.at(0)=='#')
                {
                    //std::cout<<"Escape line: "<<line<<std::endl;
                }
                else
                {
                    outputMatrix.push_back(parseMatrixLineString(line));
                }   
            }
        }
        inputFile.close();

        return outputMatrix;
    }else{
        std::cout<<"Error: Unable to open file: "<<filename<<std::endl;
        return outputMatrix;
    }
}


std::vector<std::string> StateSwitchingCellProperty::parseMatrixLineString(std::string line)
{
    // for a line string in the matrix read, parse into vector data entries based on delimiters ','
    std::vector<std::string> rowVector = std::vector<std::string>();

    // delimiter, may be modified by further methods
    std::string delim = ",";
    std::string matrixCell;

    // determine the position of the delimiter
    size_t posSnew=line.find(delim);

    bool IsEndOfLine = false;
    
    while(!IsEndOfLine)
    {
        // while not at the end of the file, sample sub strings from the posiiton of the delimiter
        if(posSnew == std::string::npos)
        {
            IsEndOfLine = true;
        }
        
        // sample substring from begining of the string to the delimiter position, store as data entry
        matrixCell = line.substr(0,posSnew);

        // remove the sampled entry from the string
        line = line.substr(posSnew+1,std::string::npos);

        rowVector.push_back(matrixCell);

        // update delimiter position
        posSnew=line.find(delim);
    }
    return rowVector;
}



// set methods

void StateSwitchingCellProperty::SetCellState(unsigned cellState)
{
    mCellState = cellState;
}

void StateSwitchingCellProperty::SetCellPtr(CellPtr this_cellPtr)
{
    mThis_cellPtr = this_cellPtr;
}

void StateSwitchingCellProperty::SetBaseDirectoryName(std::string name)
{
    mBaseDirectoryName = name;
}

void StateSwitchingCellProperty::SetStateNames(std::vector<std::string> names)
{
    std::cout<<"StateSwitchingCellProperty::SetStateNames(std::vector<std::string> names)"<<std::endl;
    mStateNames = names;
}

void StateSwitchingCellProperty::SetStateFiles(std::vector<std::vector<std::string>> stateFiles)
{
    mStateFiles = stateFiles;
}


void StateSwitchingCellProperty::SetConditionChemicalNames(std::vector<std::vector<std::string>> names)
{
    mStateConditionsChemicalNames = names;
}

void StateSwitchingCellProperty::SetConditionTypes(std::vector<std::vector<std::string>> types)
{
    mStateConditionTypes = types;
}

void StateSwitchingCellProperty::SetIsStateConditionsMaxSet(std::vector<std::vector<bool>>maxValuesSet)
{
    mIsStateConditionsMaxValueSet = maxValuesSet;
}

void StateSwitchingCellProperty::SetIsStateConditionsMinSet(std::vector<std::vector<bool>>minValuesSet)
{
    mIsStateConditionsMinValueSet = minValuesSet;
}

void StateSwitchingCellProperty::SetStateConditionsMax(std::vector<std::vector<double>>maxValues)
{
    mStateConditionsMaxValue = maxValues;
}

void StateSwitchingCellProperty::SetStateConditionsMin(std::vector<std::vector<double>> minValues)
{
    mStateConditionsMinValue = minValues;
}

void StateSwitchingCellProperty::SetChemicalSrnModel(ChemicalSrnModel* p_chemicalSrn)
{
    mpChemicalSrnModel = p_chemicalSrn;
}

void StateSwitchingCellProperty::SetChemicalCellCycleModel(SimpleChemicalThresholdCellCycleModel* p_chemicalCellCycleModel)
{
    mpSimpleChemicalThresholdCellCycleModel = p_chemicalCellCycleModel;
}

// get methods

unsigned StateSwitchingCellProperty::GetCellState()
{
    return mCellState;
}

CellPtr StateSwitchingCellProperty::GetCellPtr()
{
    return mThis_cellPtr;
}

std::string StateSwitchingCellProperty::GetBaseDirectoryName()
{
    return mBaseDirectoryName;
}

std::vector<std::vector<std::string>> StateSwitchingCellProperty::GetStateFiles()
{
    return mStateFiles;
}

std::vector<std::string> StateSwitchingCellProperty::GetStateNames()
{
    return mStateNames;
}

std::vector<std::vector<std::string>> StateSwitchingCellProperty::GetConditionChemicalNames()
{
    return mStateConditionsChemicalNames;
}

std::vector<std::vector<std::string>> StateSwitchingCellProperty::GetConditionTypes()
{
    return mStateConditionTypes;
}

std::vector<std::vector<bool>> StateSwitchingCellProperty::GetIsStateConditionsMaxSet()
{
    return mIsStateConditionsMaxValueSet;
}

std::vector<std::vector<double>> StateSwitchingCellProperty::GetStateConditionsMax()
{
    return mStateConditionsMaxValue;
}

std::vector<std::vector<double>> StateSwitchingCellProperty::GetStateConditionsMin()
{
    return mStateConditionsMinValue;
}

std::string StateSwitchingCellProperty::GetStateNamesByIndex(unsigned index)
{
    return mStateNames[index];
}

std::vector<std::string> StateSwitchingCellProperty::GetStateFilesByIndex(unsigned index)
{
    return mStateFiles[index];
}

std::vector<std::string> StateSwitchingCellProperty::GetConditionChemicalNamesByIndex(unsigned index)
{
    return mStateConditionsChemicalNames[index];
}

std::vector<std::string> StateSwitchingCellProperty::GetConditionTypesByIndex(unsigned index)
{
    return mStateConditionTypes[index];
}

std::vector<bool> StateSwitchingCellProperty::GetIsStateConditionsMaxSetByIndex(unsigned index)
{
    return mIsStateConditionsMaxValueSet[index];
}

std::vector<bool> StateSwitchingCellProperty::GetIsStateConditionsMinSetByIndex(unsigned index)
{
    return mIsStateConditionsMinValueSet[index];
}

std::vector<double> StateSwitchingCellProperty::GetStateConditionsMaxByIndex(unsigned index)
{
    return mStateConditionsMaxValue[index];
}

std::vector<double> StateSwitchingCellProperty::GetStateConditionsMinByIndex(unsigned index)
{
    return mStateConditionsMinValue[index];
}

ChemicalSrnModel* StateSwitchingCellProperty::GetChemicalSrnModel()
{
    return mpChemicalSrnModel;
}

SimpleChemicalThresholdCellCycleModel* StateSwitchingCellProperty::GetChemicalCellCycleModel()
{
    return mpSimpleChemicalThresholdCellCycleModel;
}

#endif