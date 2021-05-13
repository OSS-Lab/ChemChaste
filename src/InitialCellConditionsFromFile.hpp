#ifndef INITIALCELLCONDITIONSFROMFILE_HPP_
#define INITIALCELLCONDITIONSFROMFILE_HPP_

// general includes
#include <string>
#include <tuple>
#include <vector>
#include <iostream>
#include <stdlib.h> 

#include "RandomNumberGenerator.hpp"

// reaction includes
#include "AbstractChemistry.hpp"

 
class InitialCellConditionsFromFile
{
protected:

    std::vector<std::string> mChemicalNamesVector;

    std::vector<double> mConcentrationVector;

    std::vector<bool> mPerturbConcentrationVector;

    std::vector<std::string> mBoolTestDictionary = {"true","True","TRUE","1"};

public:

    InitialCellConditionsFromFile(std::string);

    virtual ~InitialCellConditionsFromFile()
    {
    };

    std::vector<std::vector<std::string>> ReadMatrix(std::string);

    std::vector<std::string> parseMatrixLineString(std::string);

    bool StringToBool(std::string);

    void SetChemicalNamesVector(std::vector<std::string>);

    void SetConcentrationVector(std::vector<double>);

    void SetPerturbConcentrationVector(std::vector<bool>);

    std::vector<std::string> GetChemicalNamesVector();

    std::vector<double> GetConcentrationVector();

    std::vector<bool> GetPerturbConcentrationVector();

};

InitialCellConditionsFromFile::InitialCellConditionsFromFile(std::string conditionFilename) 
{
    std::vector<std::string> chemicalNamesVector;
    std::vector<double> concentrationVector;
    std::vector<bool> perturbConcentrationVector;

    // read and store the diffusion data as a matrix of strings 
    std::vector<std::vector<std::string>> conditionsDatabase;
    conditionsDatabase = ReadMatrix(conditionFilename);

    // parse the matrix
    for(unsigned i=0; i<conditionsDatabase.size(); i++)
    {
        // the first element of the read in database is the identifying names
        chemicalNamesVector.push_back(conditionsDatabase[i][0]);

        if(conditionsDatabase[0].size()==2)
        {
            // data for species concentration
            concentrationVector.push_back(atof(conditionsDatabase[i][1].c_str()));
        }
        else if(conditionsDatabase[0].size()==3)
        {
            // data for species concentration and whether to perturb the initial value
            
            perturbConcentrationVector.push_back(StringToBool(conditionsDatabase[i][2]));
            if(perturbConcentrationVector[i]==true)
            {
                concentrationVector.push_back(fabs(atof(conditionsDatabase[i][1].c_str()) + RandomNumberGenerator::Instance()->ranf()));
            }
            else
            {
                concentrationVector.push_back(atof(conditionsDatabase[i][1].c_str()));
            }
        }
    }


    // store the information
    SetChemicalNamesVector(chemicalNamesVector);
    SetConcentrationVector(concentrationVector);
    SetPerturbConcentrationVector(perturbConcentrationVector);
}

std::vector<std::vector<std::string>> InitialCellConditionsFromFile::ReadMatrix(std::string filename)
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

std::vector<std::string> InitialCellConditionsFromFile::parseMatrixLineString(std::string line)
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

bool InitialCellConditionsFromFile::StringToBool(std::string test_string)
{
    // take in a test string and determine if a true character string is present

    for(unsigned i=0; i<mBoolTestDictionary.size(); i++)
    {
        if(test_string == mBoolTestDictionary[i])
        {
            return true;
        }
    }

    return false;
}



void InitialCellConditionsFromFile::SetChemicalNamesVector(std::vector<std::string> namesVector)
{
    mChemicalNamesVector = namesVector;
}

void InitialCellConditionsFromFile::SetConcentrationVector(std::vector<double> concentrationVector)
{
    mConcentrationVector = concentrationVector;
}

void InitialCellConditionsFromFile::SetPerturbConcentrationVector(std::vector<bool> perturbConcentrationVector)
{
    mPerturbConcentrationVector = perturbConcentrationVector;
}


std::vector<std::string> InitialCellConditionsFromFile::GetChemicalNamesVector()
{
    return mChemicalNamesVector;
}

std::vector<double> InitialCellConditionsFromFile::GetConcentrationVector()
{
    return mConcentrationVector;
}

std::vector<bool> InitialCellConditionsFromFile::GetPerturbConcentrationVector()
{
    return mPerturbConcentrationVector;
}


#endif