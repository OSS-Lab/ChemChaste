/*
#ifndef CUSTOMCELLFROMFILE_HPP_
#define CUSTOMCELLFROMFILE_HPP_

#include "ComplexCell.hpp"
#include "ComplexCellFromFile.hpp"
#include <vector>
#include <string>
#include <boost/utility.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/enable_shared_from_this.hpp>
#include "StateSwitchingCellPropertyFromFile.hpp"
#include "StateSwitchingCellProperty.hpp"

//class ComplexCell;

class CustomCellFromFile : public ComplexCellFromFile
{
protected:
    using ComplexCellFromFile::SetUp;
    using ComplexCellFromFile::SetUpSRNandCellCycle;
    using ComplexCellFromFile::BuildCellPropertyCollection;

    std::string mConfigurationFilename = "";

    std::string mConfigurationDelimiter = " = ";
    std::string mConfigurationDelimiterShort="=";

    std::string mCellModelType = "";

    std::string mInitialCellStateName = "";

    std::vector<std::vector<std::string>> mConfigurationParameters;

    boost::shared_ptr<StateSwitchingCellProperty> mp_cell_state_switching;

    std::string mStateSwitchingPropertyFilename;

    bool mIsStateSwitchingPropertySet = false;

    std::vector<std::string> mStateNames;

    unsigned mInitialState = 0;
    


public:

    CustomCellFromFile(     std::string cellFileRoot="",
                            std::string cellModelType="",
                            std::string cellCycleFilename="",
                            std::string divisionRulesFilename="", 
                            std::string srnFilename="",
                            std::string initialConditionFilename="",
                            std::string transportPropertyFilename="",
                            std::string membranePropertyFilename="",
                            std::string environmentPropertyFilename="",
                            std::string stateSwitchingPropertyFilename="",
                            std::string configurationFilename="",
                            unsigned cellID =0,
                            bool isCellIDSet = false);

    virtual ~CustomCellFromFile()
    {
    };

    virtual void SetUp();

    virtual void BuildCellPropertyCollection(CellPropertyCollection& collection);

    void ReadStateSwitchingReplaceFileNames();

    virtual void SetUpSRNandCellCycle();

    void ParseConfigurationFile();

    void SetUpCustomCellProperties();

    bool IsParameterInConfiguration(std::string);

    std::string ReturnFromConfiguration(std::string);

    bool ReturnBoolFromConfiguration(std::string);

    double ReturnDoubleFromConfiguration(std::string);

    int ReturnIntFromConfiguration(std::string);

    std::string ReturnStringFromConfiguration(std::string);


    // set/get methods

    void SetCellModelType(std::string);

    std::string GetCellModelType();

    void SetConfigurationParameters(std::vector<std::vector<std::string>>);

    std::vector<std::vector<std::string>> GetConfigurationParameters();

    void SetStateSwitchingPropertyFilename(std::string);

    std::string GetStateSwitchingPropertyFilename();

    void GetIsStateSwitchingPropertySet(bool);

    bool GetIsStateSwitchingPropertySet();

    void SetStateSwitchingCellProperty(boost::shared_ptr<StateSwitchingCellProperty>);

    boost::shared_ptr<StateSwitchingCellProperty> GetStateSwitchingCellProperty();

    void SetCellStateNames(std::vector<std::string>);

    std::vector<std::string> GetCellStateNames();

    void SetInitialCellState(unsigned);

    unsigned GetInitialCellState();
    
    void SetInitialCellStateName(std::string);

    std::string GetInitialCellStateName();


    

};

CustomCellFromFile::CustomCellFromFile(
                        std::string cellFileRoot,
                        std::string cellModelType,
                        std::string cellCycleFilename,
                        std::string divisionRulesFilename,
                        std::string srnFilename,
                        std::string initialConditionFilename,
                        std::string transportPropertyFilename,
                        std::string membranePropertyFilename,
                        std::string environmentPropertyFilename,
                        std::string stateSwitchingPropertyFilename,
                        std::string configurationFilename,
                        unsigned cellID,
                        bool isCellIDSet)
        :   ComplexCellFromFile(cellFileRoot+cellModelType+"/",
                                cellCycleFilename,
                                divisionRulesFilename,
                                srnFilename,
                                initialConditionFilename,
                                transportPropertyFilename,
                                membranePropertyFilename,
                                environmentPropertyFilename,
                                cellID,
                                isCellIDSet),
            mCellModelType(cellModelType),
            mConfigurationFilename(configurationFilename),
            mStateSwitchingPropertyFilename(stateSwitchingPropertyFilename)
{
    std::cout<<"CustomCellFromFile::CustomCellFromFile::CustomCellFromFile::CustomCellFromFile() - start"<<std::endl;
    ParseConfigurationFile();

    if(IsParameterInConfiguration("split_ratio"))
    {
        
        this->mSplitRatio = ReturnDoubleFromConfiguration("split_ratio");
      
        //boost::shared_ptr<Cell> 
        //mpCell->SetSplitRatio(ReturnDoubleFromConfiguration("split_ratio"));
    }

    // if there is a cell state file set up state switching property
    if(mStateSwitchingPropertyFilename != "")
    {
  //      std::cout<<"################test if there is and set transprot"<<std::endl;
        mIsStateSwitchingPropertySet = true;
    }
    std::cout<<"CustomCellFromFile::CustomCellFromFile::CustomCellFromFile::CustomCellFromFile() - end"<<std::endl;
}

void CustomCellFromFile::SetUp()
{
    std::cout<<"CustomCellFromFile::SetUp() - start"<<std::endl;
    // for properties that are switched, overide the property filename paths while setting the property to be set 

    if(mIsStateSwitchingPropertySet)
    {
        ReadStateSwitchingReplaceFileNames();
    }

    ComplexCellFromFile::SetUp();

    SetUpCustomCellProperties();
    
    std::cout<<"CustomCellFromFile::SetUp() - end"<<std::endl;
}

void CustomCellFromFile::ReadStateSwitchingReplaceFileNames()
{
    std::cout<<"CustomCellFromFile::ReadStateSwitchingReplaceFileNames() - start"<<std::endl;
    std::cout<<"here state switching property is set: "<<mStateSwitchingPropertyFilename<<std::endl;
    StateSwitchingCellPropertyFromFile* p_state_switching_property_from_file = new StateSwitchingCellPropertyFromFile(this->mCellFileRoot,mStateSwitchingPropertyFilename,"StateSwitchingConditions.csv");
    std::cout<<"here0"<<std::endl;
    std::vector<std::string> fileNamesVector = p_state_switching_property_from_file -> GetUniquePropertyNamesVector();
    std::cout<<p_state_switching_property_from_file->GetCellStatesVector()[GetInitialCellState()]<<std::endl;
    SetInitialCellStateName(p_state_switching_property_from_file->GetCellStatesVector()[GetInitialCellState()]); //
    //SetStateSwitchingCellProperty(p_state_switching_property_from_file ->GetStateSwitchingCellProperty());
    SetCellStateNames(p_state_switching_property_from_file->GetCellStatesVector());

    boost::shared_ptr<StateSwitchingCellProperty> p_cell_switch(new StateSwitchingCellProperty());
    SetStateSwitchingCellProperty(p_cell_switch);


    std::cout<<"here01"<<std::endl;
    //p_state_switching_property_from_file ->GetStateSwitchingCellProperty() -> SetStateNames(p_state_switching_property_from_file->GetCellStatesVector());

    std::string fileName = "";
    std::cout<<"here1"<<std::endl;
    // search the state switching property for files that may be switched
    if(fileNamesVector.size()>0)
    {
        // replace swithced properties with null equivalents
        for(unsigned i=0; i<fileNamesVector.size(); i++)
        {
            fileName = fileNamesVector[i];
            // if a file is defined within the switching property then don't read in complex cell
            if(fileName=="Srn.txt")
            {
                this->mIsSRNSet = true;
                this->mSrnFilename = "";
                
                //ChemicalSRNFromFile* p_srn_reaction_system_from_file = new ChemicalSRNFromFile(p_state_switching_property_from_file->GetFileDirectory() + fileName);
                //SetChemicalSrnModel(p_srn_reaction_system_from_file -> GetChemicalSrnModel());

            }
            else if(fileName=="TransportReactions.txt")
            {
                this->mTransportPropertyFilename = "";
                this->mIsTransportPropertySet = true;
                boost::shared_ptr<TransportCellProperty> p_cell_transport(new TransportCellProperty());

                this->SetTransportCellProperty(p_cell_transport);
            }
            else if(fileName=="MembraneReactions.txt")
            {
                this->mMembranePropertyFilename = "";
                this->mIsMembranePropertySet = true;

                boost::shared_ptr<MembraneCellProperty> p_cell_membrane(new MembraneCellProperty());

                this->SetMembraneCellProperty(p_cell_membrane);
            }
            else if(fileName=="SpeciesDivisionRules.csv")
            {
                this->mIsDivisionRulesSet = true;
                this->mDivisionRulesFilename = "";
            }
            else if(fileName=="Environment.csv")
            {
                this->mEnvironmentPropertyFilename = "";
                this->mIsEnvironmentPropertySet = true;

                boost::shared_ptr<EnvironmentCellProperty> p_cell_environment(new EnvironmentCellProperty());

                this->SetEnvironmentCellProperty(p_cell_environment);

            }
            else if(fileName=="SpeciesThreshold.csv")
            {
                this->mCellCycleFilename = "";
                this->mIsCellCycleSet = true;
                //SimpleChemicalThresholdCellCycleFromFile* pCellCycle = new SimpleChemicalThresholdCellCycleFromFile(p_state_switching_property_from_file->GetFileDirectory() + fileName);

                //pCellCycle -> SetUp();
                //this->SetChemicalCellCycleModel(pCellCycle);
            }
        }
    }
    std::cout<<"CustomCellFromFile::ReadStateSwitchingReplaceFileNames() - end"<<std::endl;
}



void CustomCellFromFile::SetUpSRNandCellCycle()
{
    std::cout<<"CustomCellFromFile::SetUpSRNandCellCycle() - start"<<std::endl;
    if(mIsSRNSet)
    {
        std::cout<<"set SRN"<<std::endl;
        std::string srnFilename="";

        if(mSrnFilename != "")
        {
            // as normal
            //mSrnFilename = mCellFileRoot+ mSrnFilename;
            srnFilename = mSrnFilename;
        }
        else if(mIsStateSwitchingPropertySet)
        {
            std::cout<<"switching"<<std::endl;
            // pull the srn name as a cell switched property
            srnFilename = mCellFileRoot+mInitialCellStateName + "/Srn.txt";
            
            //mp_cell_state_switching->RetrieveStateName();// + "/srn.txt";

        }
        else
        {
            srnFilename="";
            //ccFilename = ccFilename;
        }
        
        std::cout<<"setting file: "<<srnFilename<<std::endl;
        ChemicalSRNFromFile* p_srn_reaction_system_from_file = new ChemicalSRNFromFile(srnFilename);
        SetChemicalSrnModel(p_srn_reaction_system_from_file -> GetChemicalSrnModel());  
    }


    if(mIsCellCycleSet)
    {
        std::cout<<"set cc"<<std::endl;
        std::string ccFilename="";

        if(mCellCycleFilename != "")
        {
            // as normal
            //mCellCycleFilename  =mCellFileRoot+mCellCycleFilename;
            ccFilename = mCellCycleFilename;
        }
        else if(mIsStateSwitchingPropertySet)
        {
            std::cout<<"switching"<<std::endl;
            // pull the srn name as a cell switched property
            ccFilename = mCellFileRoot+mInitialCellStateName + "/SpeciesThreshold.csv";//mp_cell_state_switching->RetrieveStateName() + "/SpeciesThreshold.csv";

        }
        else
        {
            ccFilename="";
            //ccFilename = ccFilename;
        }
        std::cout<<"setting file: "<<ccFilename<<std::endl;
        SimpleChemicalThresholdCellCycleFromFile* pCellCycle = new SimpleChemicalThresholdCellCycleFromFile(mCellCycleFilename);

        pCellCycle -> SetUp();
        SetChemicalCellCycleModel(pCellCycle);
    }
    std::cout<<"CustomCellFromFile::SetUpSRNandCellCycle() - end"<<std::endl;

}

void CustomCellFromFile::BuildCellPropertyCollection(CellPropertyCollection& collection)
{
    std::cout<<"CustomCellFromFile::BuildCellPropertyCollection(CellPropertyCollection& collection) - start"<<std::endl;
    if(mIsStateSwitchingPropertySet)
    {
        std::cout<<"add switching property"<<std::endl;
        collection.AddProperty(mp_cell_state_switching);
        std::cout<<"-------has switch: "<<collection.HasProperty<StateSwitchingCellProperty>()<<std::endl;
    }
    ComplexCellFromFile::BuildCellPropertyCollection(collection);
    std::cout<<"CustomCellFromFile::BuildCellPropertyCollection(CellPropertyCollection& collection) - end"<<std::endl;
}

void CustomCellFromFile::ParseConfigurationFile()
{
    std::cout<<"CustomCellFromFile::ParseConfigurationFile() - start"<<std::endl;
    std::ifstream inputFile(mCellFileRoot+mConfigurationFilename);

    std::string line;

    std::vector<std::string> parameterString;

    std::vector<std::vector<std::string>> configParameters;

    if(inputFile.is_open())
    {   
        // open the reaction file
        while (getline(inputFile,line))
        {

            // for each non-empty reation file line, parse the reactions into
            // data structures on line by line basis
            if(!line.empty())
            {
                // look for the delimiter in the string
                size_t configuration_delimiter_position = line.find(mConfigurationDelimiter);
                size_t configuration_delimiter_short_position = line.find(mConfigurationDelimiterShort);
                std::string parameterInfo;
                std::string parameterValue;
                if(configuration_delimiter_position != std::string::npos)
                {

                    // take substring from the delimiter to end
                    parameterValue = line.substr(configuration_delimiter_position+3,std::string::npos);

                    line.erase(configuration_delimiter_position, std::string::npos);

                    parameterInfo = line;

                    parameterString = {parameterInfo,parameterValue};

                    configParameters.push_back(parameterString);

                }
                else if(configuration_delimiter_short_position != std::string::npos)
                {
                    // take substring from the delimiter to end
                    parameterValue = line.substr(configuration_delimiter_short_position+1,std::string::npos);

                    line.erase(configuration_delimiter_short_position, std::string::npos);

                    parameterInfo = line;

                    parameterString = {parameterInfo,parameterValue};

                    configParameters.push_back(parameterString);

                }
                // line doesn't have the delimiter so skip
            }
        }

        SetConfigurationParameters(configParameters);


    inputFile.close();
    }
    else
    {
        std::cout<<"Error: Filename not found: "<<mConfigurationFilename<<std::endl;
    }
    std::cout<<"CustomCellFromFile::ParseConfigurationFile() - end"<<std::endl;
}

void CustomCellFromFile::SetUpCustomCellProperties()
{
    std::cout<<"CustomCellFromFile::SetUpCustomCellProperties() - start"<<std::endl;
    StateVariableRegister* pFullStateVariableRegister = GetFullChemicalStateRegister();

    if(mIsStateSwitchingPropertySet)
    {
        std::string stateSwitchFileName = "StateSwitchingConditions.csv";
        

        StateSwitchingCellPropertyFromFile* pStateSwitchingFromFile = new StateSwitchingCellPropertyFromFile(
                                    mCellFileRoot,
                                    mStateSwitchingPropertyFilename,
                                    stateSwitchFileName
        );

        //boost::shared_ptr<StateSwitchingCellProperty> pStateProperty(new StateSwitchingCellProperty());

        pStateSwitchingFromFile -> SetUpStateSwitchingProperty(mpCell);

        boost::shared_ptr<StateSwitchingCellProperty> pStateProperty = boost::static_pointer_cast<StateSwitchingCellProperty>(mpCell->rGetCellPropertyCollection().GetPropertiesType<StateSwitchingCellProperty>().GetProperty());

        pStateProperty = pStateSwitchingFromFile -> GetStateSwitchingCellProperty();
        pStateProperty->SetCellPtr(mpCell);

        pStateProperty->DetermineState();
        pStateProperty->SetUpSwitchedProperties();

        SetStateSwitchingCellProperty(pStateProperty);


    //    std::vector<std::string> cell_state_switching_chemical_names = mp_cell_state_switching -> GetCellStateVariableRegister() -> GetStateVariableRegisterVector();
    //    pFullStateVariableRegister -> AddStateVariableVector(cell_state_switching_chemical_names);
        std::cout<<"set state switching property"<<std::endl;
    } 

    SetFullChemicalStateRegister(pFullStateVariableRegister);
    std::cout<<"CustomCellFromFile::SetUpCustomCellProperties() - end"<<std::endl;
}


bool CustomCellFromFile::IsParameterInConfiguration(std::string query)
{
    bool value = false;
    for(unsigned i=0;i<mConfigurationParameters.size();i++)
    {
        if(mConfigurationParameters[i][0] == query)
        {
            value = true;
            break;
        }
    }

    return value;
}

std::string CustomCellFromFile::ReturnFromConfiguration(std::string query)
{
    std::cout<<"CustomCellFromFile::ReturnFromConfiguration(std::string query) - start"<<std::endl;
    // assume IsParameterInConfiguration already true   
    std::string value;

    for(unsigned i=0;i<mConfigurationParameters.size();i++)
    {
        if(mConfigurationParameters[i][0] == query)
        {
            value = mConfigurationParameters[i][1];
            break;
        }
    }
    std::cout<<value<<std::endl;
    std::cout<<"CustomCellFromFile::ReturnFromConfiguration(std::string query) - end"<<std::endl;
    return value;
}

bool CustomCellFromFile::ReturnBoolFromConfiguration(std::string query)
{
    // assume IsParameterInConfiguration already true

    bool value=false;

    if(ReturnFromConfiguration(query)=="true")
    {
        value = true;
    }

    return value;
}

double CustomCellFromFile::ReturnDoubleFromConfiguration(std::string query)
{
    // assume IsParameterInConfiguration already true
    double value;

    value = std::stod(ReturnFromConfiguration(query));
    std::cout<<"return double: "<<value<<std::endl;
    return value;
}

int CustomCellFromFile::ReturnIntFromConfiguration(std::string query)
{
    // assume IsParameterInConfiguration already true   
    int value;

    value = std::stoi(ReturnFromConfiguration(query));

    return value;
}

std::string CustomCellFromFile::ReturnStringFromConfiguration(std::string query)
{
    // assume IsParameterInConfiguration already true   
    std::string value;

    value = ReturnFromConfiguration(query);

    return value;
}


void CustomCellFromFile::SetConfigurationParameters(std::vector<std::vector<std::string>> configParameters)
{
    mConfigurationParameters = configParameters;
}

std::vector<std::vector<std::string>> CustomCellFromFile::GetConfigurationParameters()
{
    return mConfigurationParameters;
}

void CustomCellFromFile::SetCellModelType(std::string modelType)
{
    mCellModelType = modelType;
}

std::string CustomCellFromFile::GetCellModelType()
{
    return mCellModelType;
}


void CustomCellFromFile::SetStateSwitchingPropertyFilename(std::string filename)
{
    mStateSwitchingPropertyFilename = filename;
}

std::string CustomCellFromFile::GetStateSwitchingPropertyFilename()
{
    return mStateSwitchingPropertyFilename;
}

void CustomCellFromFile::GetIsStateSwitchingPropertySet(bool isStateSwitchingSet)
{
    mIsStateSwitchingPropertySet = isStateSwitchingSet;
}


bool CustomCellFromFile::GetIsStateSwitchingPropertySet()
{
    return mIsStateSwitchingPropertySet;
}


void CustomCellFromFile::SetStateSwitchingCellProperty(boost::shared_ptr<StateSwitchingCellProperty> cell_state_switching)
{
    mp_cell_state_switching = cell_state_switching;
}

boost::shared_ptr<StateSwitchingCellProperty> CustomCellFromFile::GetStateSwitchingCellProperty()
{
    return mp_cell_state_switching;
}

void CustomCellFromFile::SetInitialCellStateName(std::string cellState)
{
    mInitialCellStateName = cellState;
}

std::string CustomCellFromFile::GetInitialCellStateName()
{
    return mInitialCellStateName;
}


void CustomCellFromFile::SetCellStateNames(std::vector<std::string> names)
{
    mStateNames = names;
}

std::vector<std::string> CustomCellFromFile::GetCellStateNames()
{
    return mStateNames;
}

void CustomCellFromFile::SetInitialCellState(unsigned state)
{
    mInitialState = state;
}

unsigned CustomCellFromFile::GetInitialCellState()
{
    return mInitialState;
}
    


#endif

*/