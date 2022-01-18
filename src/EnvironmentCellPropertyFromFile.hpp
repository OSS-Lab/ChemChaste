#ifndef ENVIRONMENTCELLPROPERTYFROMFILE_HPP
#define ENVIRONMENTCELLPROPERTYFROMFILE_HPP

//general includes
#include <string>
#include <tuple>
#include <vector>
#include <iostream>
#include <stdlib.h> 
#include <boost/shared_ptr.hpp>

// reaction includes
#include "AbstractChemistry.hpp"
#include "EnvironmentCellProperty.hpp"
#include "ChemicalCell.hpp"
#include "PreferredCellEnvironmentFromFile.hpp"



class EnvironmentCellPropertyFromFile
{
protected:

    std::string mEnvironmentFilename;

    boost::shared_ptr<EnvironmentCellProperty> mpEnvironmentProperty;

    std::vector<std::string> mChemicalNamesVector;

    std::vector<double> mConcentrationVector;

    std::vector<bool> mPerturbConcentrationVector;

    std::vector<std::string> mBoolTestDictionary = {"true","True","TRUE","1"};


public:

    EnvironmentCellPropertyFromFile(std::string filename ="");

    ~EnvironmentCellPropertyFromFile()
    {
    };


    std::vector<std::vector<std::string>> ReadMatrix(std::string);

    std::vector<std::string> parseMatrixLineString(std::string);

    bool StringToBool(std::string);

    void SetChemicalNamesVector(std::vector<std::string>);

    void SetConcentrationVector(std::vector<double>);

    void SetPerturbConcentrationVector(std::vector<bool>);

    void SetUpEnvironmentProperty(CellPtr);

    void SetEnvironmentFilename(std::string);

    void SetEnvironmentProperty(boost::shared_ptr<EnvironmentCellProperty>);

    std::vector<std::string> GetChemicalNamesVector();

    std::vector<double> GetConcentrationVector();

    std::vector<bool> GetPerturbConcentrationVector();


    std::string GetEnvironmentFilename();

    boost::shared_ptr<EnvironmentCellProperty> GetEnvironmentProperty();

};

EnvironmentCellPropertyFromFile::EnvironmentCellPropertyFromFile(std::string filename)
    : mEnvironmentFilename(filename)
{
    // this constructor contents can be deleted?
    if(filename == "")
    {
        std::cout<<"Error: No Environment file passed"<<std::endl;
    }

    std::vector<std::string> chemicalNamesVector;
    std::vector<double> concentrationVector;
    std::vector<bool> perturbConcentrationVector;

    // read and store the diffusion data as a matrix of strings 
    std::vector<std::vector<std::string>> environmentDatabase;
    environmentDatabase = ReadMatrix(filename);

    // parse the matrix
    for(unsigned i=0; i<environmentDatabase.size(); i++)
    {
        // the first element of the read in database is the identifying names
        chemicalNamesVector.push_back(environmentDatabase[i][0]);

        if(environmentDatabase[0].size()==2)
        {
            // data for species concentration
            concentrationVector.push_back(atof(environmentDatabase[i][1].c_str()));
        }
        else if(environmentDatabase[0].size()==3)
        {
            // data for species concentration and whether to perturb the initial value
            
            perturbConcentrationVector.push_back(StringToBool(environmentDatabase[i][2]));
            if(perturbConcentrationVector[i]==true)
            {
                concentrationVector.push_back(fabs(atof(environmentDatabase[i][1].c_str()) + RandomNumberGenerator::Instance()->ranf()));
            }
            else
            {
                concentrationVector.push_back(atof(environmentDatabase[i][1].c_str()));
            }
        }
    }


    // store the information
    SetChemicalNamesVector(chemicalNamesVector);
    SetConcentrationVector(concentrationVector);
    SetPerturbConcentrationVector(perturbConcentrationVector);


}

void EnvironmentCellPropertyFromFile::SetUpEnvironmentProperty(CellPtr p_cell)
{
    boost::shared_ptr<EnvironmentCellProperty> p_environment = boost::static_pointer_cast<EnvironmentCellProperty>(p_cell->rGetCellPropertyCollection().GetPropertiesType<EnvironmentCellProperty>().GetProperty());
    
    
    PreferredCellEnvironmentFromFile* p_preferred_environment_from_file = new PreferredCellEnvironmentFromFile(mEnvironmentFilename);

    p_environment -> SetUp(p_cell, p_preferred_environment_from_file->GetEnvironmentNamesVector());
    p_environment -> SetPreferredEnvironmentVector(p_preferred_environment_from_file->GetEnvironmentValueVector());

    SetEnvironmentProperty(p_environment);

}


std::vector<std::vector<std::string>> EnvironmentCellPropertyFromFile::ReadMatrix(std::string filename)
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

std::vector<std::string> EnvironmentCellPropertyFromFile::parseMatrixLineString(std::string line)
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

bool EnvironmentCellPropertyFromFile::StringToBool(std::string test_string)
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

void EnvironmentCellPropertyFromFile::SetEnvironmentFilename(std::string filename)
{
    mEnvironmentFilename = filename;
}

void EnvironmentCellPropertyFromFile::SetEnvironmentProperty(boost::shared_ptr<EnvironmentCellProperty> pEnvironmentProperty)
{
    mpEnvironmentProperty = pEnvironmentProperty;
}

void EnvironmentCellPropertyFromFile::SetChemicalNamesVector(std::vector<std::string> namesVector)
{
    mChemicalNamesVector = namesVector;
}

void EnvironmentCellPropertyFromFile::SetConcentrationVector(std::vector<double> concentrationVector)
{
    mConcentrationVector = concentrationVector;
}

void EnvironmentCellPropertyFromFile::SetPerturbConcentrationVector(std::vector<bool> perturbConcentrationVector)
{
    mPerturbConcentrationVector = perturbConcentrationVector;
}


std::vector<std::string> EnvironmentCellPropertyFromFile::GetChemicalNamesVector()
{
    return mChemicalNamesVector;
}

std::vector<double> EnvironmentCellPropertyFromFile::GetConcentrationVector()
{
    return mConcentrationVector;
}

std::vector<bool> EnvironmentCellPropertyFromFile::GetPerturbConcentrationVector()
{
    return mPerturbConcentrationVector;
}


std::string EnvironmentCellPropertyFromFile::GetEnvironmentFilename()
{
    return mEnvironmentFilename;
}

boost::shared_ptr<EnvironmentCellProperty> EnvironmentCellPropertyFromFile::GetEnvironmentProperty()
{
    return mpEnvironmentProperty;

}

#endif