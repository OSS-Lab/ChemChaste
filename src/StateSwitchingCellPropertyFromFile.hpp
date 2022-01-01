#ifndef STATESWITCHINGCELLPROPERTYFROMFILE_HPP
#define STATESWITCHINGCELLPROPERTYFROMFILE_HPP

//general includes
#include <string>
#include <tuple>
#include <vector>
#include <iostream>
#include <stdlib.h> 
#include <boost/shared_ptr.hpp>

// reaction includes
#include "AbstractChemistry.hpp"
#include "StateSwitchingCellProperty.hpp"
#include "ChemicalCell.hpp"
#include "PreferredCellEnvironmentFromFile.hpp"



class StateSwitchingCellPropertyFromFile
{
protected:

    std::string mFileDirectory;

    std::string mCellStatesFilename;

    std::string mStateSwitchingFilename;

    boost::shared_ptr<StateSwitchingCellProperty> mpStateSwitchingProperty;

    std::vector<std::string> mStateNamesVector;

    std::vector<std::vector<std::string>> mPropertyNamesVector;

    std::vector<std::string> mUniquePropertyNamesVector;

    std::vector<std::vector<std::vector<std::string>>> mSwitchConditions;


public:

    StateSwitchingCellPropertyFromFile(std::string fileDirectory ="", std::string statesFilename ="", std::string stateSwitchingFilename="");

    ~StateSwitchingCellPropertyFromFile()
    {
    };

    void SetUpStateSwitchingProperty(CellPtr);

    std::vector<std::vector<std::string>> ReadMatrix(std::string);

    std::vector<std::string> parseMatrixLineString(std::string);

    std::vector<std::string> ReturnUnique(std::vector<std::vector<std::string>> );


    void SetUpStateSwitchinProperty(CellPtr);

    void SetFileDirectory(std::string);

    void SetCellStatesFilename(std::string);

    void SetStateSwitchinProperty(boost::shared_ptr<StateSwitchingCellProperty>);

    void SetStateSwitchingProperty(boost::shared_ptr<StateSwitchingCellProperty>);

    void SetCellStatesVector(std::vector<std::string>);

    void SetPropertyNamesVector(std::vector<std::vector<std::string>>);

    void SetUniquePropertyNamesVector(std::vector<std::string>);

    void SetSwitchConditions(std::vector<std::vector<std::vector<std::string>>>);

    std::string GetFileDirectory();

    std::string GetCellStatesFilename();

    boost::shared_ptr<StateSwitchingCellProperty> GetStateSwitchingCellProperty();

    std::vector<std::string> GetCellStatesVector();

    std::vector<std::vector<std::string>> GetPropertyNamesVector();

    std::vector<std::string> GetUniquePropertyNamesVector();

    std::vector<std::vector<std::vector<std::string>>> GetSwitchConditions();

};

StateSwitchingCellPropertyFromFile::StateSwitchingCellPropertyFromFile(std::string fileDirectory, std::string statesFilename, std::string stateSwitchingFilename)
    :   mFileDirectory(fileDirectory),
        mCellStatesFilename(statesFilename),
        mStateSwitchingFilename(stateSwitchingFilename)
{
    // CellStates.csv
    
    // this constructor contents can be deleted?
    if(statesFilename == "")
    {
        std::cout<<"Error: StateSwitchingCellPropertyFromFile::StateSwitchingCellPropertyFromFile No state file passed"<<std::endl;
    }

    if(stateSwitchingFilename == "")
    {
        std::cout<<"Error: StateSwitchingCellPropertyFromFile::StateSwitchingCellPropertyFromFile No state switching file name passed"<<std::endl;
    }

    std::vector<std::string> stateNamesVector;
    std::vector<std::vector<std::string>> propertyNamesVector;
  
    std::vector<std::vector<std::string>> cellStateDatabase;
    // "CellStates.csv"
    cellStateDatabase = ReadMatrix(fileDirectory+statesFilename);

    // fileDirectory + "CELL_STATE/stateSwitchingFilename"
    std::string stateSwitchingFile="";
    std::vector<std::vector<std::string>> switchConditionDatabase;

    std::vector<std::vector<std::vector<std::string>>> switchConditions;

    // parse the matrix
    for(unsigned i=0; i<cellStateDatabase.size(); i++)
    {
        // the first element of the read in database is the identifying names, linking to cell subdirectories
        stateNamesVector.push_back(cellStateDatabase[i][0]);

        // vector of the property file names changed during switching, located within the state directory 
        std::vector<std::string> cellPropertyFiles = std::vector<std::string>();
        if(cellStateDatabase[i].size()>1)
        {
            for(unsigned j=0; j<cellStateDatabase[i].size()-1; j++)
            {
                cellPropertyFiles.push_back(cellStateDatabase[i][j+1]);
            }
        }
        propertyNamesVector.push_back(cellPropertyFiles);

        // for each of the states read the particualr switching conditions
        
        stateSwitchingFile = fileDirectory+stateNamesVector[i] + "/" + stateSwitchingFilename;

        switchConditionDatabase  = ReadMatrix(stateSwitchingFile);

        switchConditions.push_back(switchConditionDatabase);

    }


    // store the information
    SetCellStatesVector(stateNamesVector);
    SetPropertyNamesVector(propertyNamesVector);
    SetUniquePropertyNamesVector(ReturnUnique(propertyNamesVector));

    SetSwitchConditions(switchConditions);

}

void StateSwitchingCellPropertyFromFile::SetUpStateSwitchingProperty(CellPtr p_cell)
{
    std::cout<<" StateSwitchingCellPropertyFromFile::SetUpStateSwitchingProperty - start"<<std::endl;
    // set properties depending on directory names
    boost::shared_ptr<StateSwitchingCellProperty> p_state_switching = boost::static_pointer_cast<StateSwitchingCellProperty>(p_cell->rGetCellPropertyCollection().GetPropertiesType<StateSwitchingCellProperty>().GetProperty());
    
    std::cout<<"here1"<<std::endl;
    p_state_switching -> SetCellPtr(p_cell);

    p_state_switching -> SetBaseDirectoryName(mFileDirectory);

    p_state_switching -> SetStateNames(mStateNamesVector);

    p_state_switching -> SetStateFiles(mPropertyNamesVector);

    // for each state parse the conditions, assume each condition is of the form:
    // TYPE CHEMICAL_NAME MAX_CONCENTRATION MIN_CONCENTRATION
    // N.B. a state may have multiple conditons 
    std::vector<std::vector<std::string>> stateConditionTypes;

    std::vector<std::vector<std::string>> stateConditionsChemicalNames;

    // whether there has been a maximum value set
    std::vector<std::vector<bool>> isStateConditionsMaxValueSet;

    // whether there has been a minimum value set
    std::vector<std::vector<bool>> isStateConditionsMinValueSet;

    std::vector<std::vector<double>> stateConditionsMaxValue;

    std::vector<std::vector<double>> stateConditionsMinValue;

    for(unsigned state_index=0; state_index<mSwitchConditions.size(); state_index++)
    {

        std::vector<std::string> conditionTypes;
        std::vector<std::string> conditionChemicalNames;
        std::vector<bool> conditionMaxSet;
        std::vector<bool> conditionMinSet;
        std::vector<double> conditionMax;
        std::vector<double> conditionMin;

        std::cout<<"state number: "<<state_index<<std::endl;

        // parse switch condition database
        for(unsigned condition_index=0; condition_index<mSwitchConditions[state_index].size(); condition_index++)
        {
            std::cout<<"condition number: "<<state_index<<std::endl;
            // for state state_index 
            bool isMinSet=false;
            bool isMaxSet=false;
            double minVal = 0.0;
            double maxVal = 0.0;

            if(mSwitchConditions[state_index][condition_index].size()>=4)
            {
                // full condition specified
                conditionTypes.push_back(mSwitchConditions[state_index][condition_index][0]);
                conditionChemicalNames.push_back(mSwitchConditions[state_index][condition_index][1]);
                std::cout<<"Type: "<<mSwitchConditions[state_index][condition_index][0]<<std::endl;
                std::cout<<"Chemical: "<<mSwitchConditions[state_index][condition_index][1]<<std::endl;
                // if value == "0.0" then not set 
                
                if(mSwitchConditions[state_index][condition_index][2] != "0.0")
                {
                    isMaxSet = true;
                }
                if(mSwitchConditions[state_index][condition_index][3] != "0.0")
                {
                    isMinSet = true;
                }
                if(isMaxSet && isMinSet)
                {
                    maxVal = std::stod (mSwitchConditions[state_index][condition_index][2]);
                    minVal = std::stod (mSwitchConditions[state_index][condition_index][3]);

                    // do not allow for negative values
                    if(maxVal<=0.0)
                    {
                        maxVal = 0.0;
                        isMaxSet = false;
                    }
                    if(minVal<=0.0)
                    {
                        minVal = 0.0;
                        isMinSet = false;
                    }

                    // if maximum < minimum then maximum is not set
                    if(maxVal<=minVal && minVal)
                    {
                        maxVal = 0.0;
                        isMaxSet = false;
                    }
                }
                std::cout<<"isMaxSet: "<<isMaxSet<<std::endl;
                std::cout<<"isMinSet: "<<isMinSet<<std::endl;
                std::cout<<"maxVal: "<<maxVal<<std::endl;
                std::cout<<"minVal: "<<minVal<<std::endl;
                conditionMaxSet.push_back(isMaxSet);
                conditionMinSet.push_back(isMinSet);
                conditionMax.push_back(maxVal);
                conditionMin.push_back(minVal);

            }
            else if(mSwitchConditions[state_index][condition_index].size()==3)
            {
                // assume maximum is specified 
                conditionTypes.push_back(mSwitchConditions[state_index][condition_index][0]);
                conditionChemicalNames.push_back(mSwitchConditions[state_index][condition_index][1]);

                if(mSwitchConditions[state_index][condition_index][2] != "0.0")
                {
                    isMaxSet = true;
                    maxVal = std::stod (mSwitchConditions[state_index][condition_index][2]);

                    // do not allow for negative values
                    if(maxVal<=0.0)
                    {
                        maxVal = 0.0;
                        isMaxSet = false;
                    }
                }

                conditionMaxSet.push_back(isMaxSet);
                conditionMinSet.push_back(isMinSet);
                conditionMax.push_back(maxVal);
                conditionMin.push_back(minVal);

            }
            else if(mSwitchConditions[state_index][condition_index].size()==3)
            {
                // no max/min conditions set
                conditionTypes.push_back(mSwitchConditions[state_index][condition_index][0]);
                conditionChemicalNames.push_back(mSwitchConditions[state_index][condition_index][1]);

                // push default
                conditionMaxSet.push_back(isMaxSet);
                conditionMinSet.push_back(isMinSet);
                conditionMax.push_back(maxVal);
                conditionMin.push_back(minVal);
            }
            else
            {
                // no parameters and/or type specified
                std::cout<<"Error: StateSwitchingCellPropertyFromFile::SetUpStateSwitchingProperty No conditions passed"<<std::endl;
            }
        }
        
        stateConditionTypes.push_back(conditionTypes);

        stateConditionsChemicalNames.push_back(conditionChemicalNames);

        isStateConditionsMaxValueSet.push_back(conditionMaxSet);

        isStateConditionsMinValueSet.push_back(conditionMinSet);

        stateConditionsMaxValue.push_back(conditionMax);

        stateConditionsMinValue.push_back(conditionMin);

    }


    p_state_switching -> SetConditionTypes(stateConditionTypes);
    p_state_switching -> SetConditionChemicalNames(stateConditionsChemicalNames);
    p_state_switching -> SetIsStateConditionsMaxSet(isStateConditionsMaxValueSet);
    p_state_switching -> SetIsStateConditionsMinSet(isStateConditionsMinValueSet);
    p_state_switching -> SetStateConditionsMax(stateConditionsMaxValue);
    p_state_switching -> SetStateConditionsMin(stateConditionsMinValue);

    p_state_switching -> SetCellState(0); // initial state
    std::cout<<"here2"<<std::endl;
    p_state_switching -> SetUp(p_cell);
    std::cout<<"here3"<<std::endl;
    SetStateSwitchingProperty(p_state_switching);
    std::cout<<" StateSwitchingCellPropertyFromFile::SetUpStateSwitchingProperty - end"<<std::endl;
}


std::vector<std::vector<std::string>> StateSwitchingCellPropertyFromFile::ReadMatrix(std::string filename)
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

std::vector<std::string> StateSwitchingCellPropertyFromFile::parseMatrixLineString(std::string line)
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
        
        // sample substring from begining of the string to the delimiter positioon, store as data entry
        matrixCell = line.substr(0,posSnew);

        // remove the sampled entry from the string
        line = line.substr(posSnew+1,std::string::npos);

        rowVector.push_back(matrixCell);

        // update delimiter position
        posSnew=line.find(delim);
    }
    return rowVector;
}

std::vector<std::string> StateSwitchingCellPropertyFromFile::ReturnUnique(std::vector<std::vector<std::string>> candidateMatrix)
{
    // serach through a matrix of strings and store the first occurance of each unique string
    std::vector<std::string> resultUnique = std::vector<std::string>();
    resultUnique.push_back(candidateMatrix[0][0]);
    unsigned uniqueCount =1;
    bool IsFound = false;

    for(unsigned j=0; j<candidateMatrix.size(); j++)
    {
        for(unsigned i=0; i<candidateMatrix[j].size(); i++)
        {
            IsFound = false;
            for(unsigned k=0; k<uniqueCount; k++)
            {
                if(candidateMatrix[j][i] == resultUnique[k])
                {
                    IsFound =true;
                    break;
                }
            }
            if(!IsFound)
            {
                resultUnique.push_back(candidateMatrix[j][i]);
                uniqueCount++;
            }
        }
    }
    return resultUnique;
}

void StateSwitchingCellPropertyFromFile::SetFileDirectory(std::string fileDirectory)
{
    mFileDirectory = fileDirectory;
}


void StateSwitchingCellPropertyFromFile::SetCellStatesFilename(std::string filename)
{
    mCellStatesFilename = filename;
}

void StateSwitchingCellPropertyFromFile::SetStateSwitchingProperty(boost::shared_ptr<StateSwitchingCellProperty> pStateSwitchingProperty)
{
    mpStateSwitchingProperty = pStateSwitchingProperty;
}

void StateSwitchingCellPropertyFromFile::SetCellStatesVector(std::vector<std::string> namesVector)
{
    mStateNamesVector = namesVector;
}

void StateSwitchingCellPropertyFromFile::SetPropertyNamesVector(std::vector<std::vector<std::string>> namesVector)
{
    mPropertyNamesVector = namesVector;
}

void StateSwitchingCellPropertyFromFile::SetUniquePropertyNamesVector(std::vector<std::string> namesVector)
{
    mUniquePropertyNamesVector = namesVector;
}

void StateSwitchingCellPropertyFromFile::SetSwitchConditions(std::vector<std::vector<std::vector<std::string>>> switchConditions)
{
    mSwitchConditions = switchConditions;
}

std::string StateSwitchingCellPropertyFromFile::GetFileDirectory()
{
    return mFileDirectory;
}

std::string StateSwitchingCellPropertyFromFile::GetCellStatesFilename()
{
    return mCellStatesFilename;
}

boost::shared_ptr<StateSwitchingCellProperty> StateSwitchingCellPropertyFromFile::GetStateSwitchingCellProperty()
{
    return mpStateSwitchingProperty;

}


std::vector<std::string> StateSwitchingCellPropertyFromFile::GetCellStatesVector()
{
    std::cout<<"StateSwitchingCellPropertyFromFile::GetCellStatesVector()"<<std::endl;
    return mStateNamesVector;
}

std::vector<std::vector<std::string>> StateSwitchingCellPropertyFromFile::GetPropertyNamesVector()
{
    return mPropertyNamesVector;
}

std::vector<std::string> StateSwitchingCellPropertyFromFile::GetUniquePropertyNamesVector()
{
    return mUniquePropertyNamesVector;
}

std::vector<std::vector<std::vector<std::string>>> StateSwitchingCellPropertyFromFile::GetSwitchConditions()
{
    return mSwitchConditions;
}

#endif
