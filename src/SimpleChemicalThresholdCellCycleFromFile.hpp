#ifndef SIMPLECHEMICALTHRESHOLDCELLCYCLEFROMFILE_HPP_
#define SIMPLECHEMICALTHRESHOLDCELLCYCLEFROMFILE_HPP_

// general includes
#include <string>
#include <tuple>
#include <vector>
#include <iostream>

// reaction includes
#include "AbstractChemistry.hpp"
#include "SimpleChemicalThresholdCellCycleModel.hpp"
 
class SimpleChemicalThresholdCellCycleFromFile : public SimpleChemicalThresholdCellCycleModel
{
protected:

    using SimpleChemicalThresholdCellCycleModel::SetUp;
    
    std::string mCellCycleFilename;

    double mSmallestMaximumThreshold = 1e-6;

public:

    SimpleChemicalThresholdCellCycleFromFile(std::string);

    virtual ~SimpleChemicalThresholdCellCycleFromFile()
    {
    };

    virtual void SetUp();

    std::vector<std::vector<std::string>> ReadMatrix(std::string);

    std::vector<std::string> parseMatrixLineString(std::string);
};

SimpleChemicalThresholdCellCycleFromFile::SimpleChemicalThresholdCellCycleFromFile(std::string thresholdFilename) 
    :   SimpleChemicalThresholdCellCycleModel(),
        mCellCycleFilename(thresholdFilename)

{
}

void SimpleChemicalThresholdCellCycleFromFile::SetUp()
{
    std::vector<std::string> thresholdChemicalVector;
    std::vector<double> maximumThresholdVector;
    std::vector<double> minimumThresholdVector;
    std::vector<bool> maximumThresholdCheck;
    std::vector<bool> minimumThresholdCheck;

    // read and store the diffusion data as a matrix of strings 
    std::vector<std::vector<std::string>> thresholdDatabase;
    thresholdDatabase = ReadMatrix(mCellCycleFilename);

    // parse the matrix
    for(unsigned i=0; i<thresholdDatabase.size(); i++)
    {
        // the first element of the read in database is the identifying names
        thresholdChemicalVector.push_back(thresholdDatabase[i][0]);

        if(thresholdDatabase[i].size()==3)
        {
            // data for both the maximum and minimum thresholds

            if(std::stoul(thresholdDatabase[i][1].c_str())<std::stoul(thresholdDatabase[i][2].c_str()))
            {
                // then can't possibly divide as will be apoptotic,so don't check
                maximumThresholdVector.push_back(std::stoul(thresholdDatabase[i][1].c_str()));
                minimumThresholdVector.push_back(std::stoul(thresholdDatabase[i][2].c_str()));
                maximumThresholdCheck.push_back(false);
                minimumThresholdCheck.push_back(true);
            }
            else
            {
                maximumThresholdVector.push_back(std::stoul(thresholdDatabase[i][1].c_str()));
                minimumThresholdVector.push_back(std::stoul(thresholdDatabase[i][2].c_str()));
                maximumThresholdCheck.push_back(true);
                minimumThresholdCheck.push_back(true);
            }
        }
        else if(thresholdDatabase[i].size()==2)
        {
            // data for only the maximum the thresholds
            maximumThresholdVector.push_back(std::stoul(thresholdDatabase[i][1].c_str()));
            minimumThresholdVector.push_back(0.0);
            maximumThresholdCheck.push_back(true);
            minimumThresholdCheck.push_back(false);
        }
        else
        {
            // then no data for either the maximum and minimum thresholds
            maximumThresholdVector.push_back(0.0);
            minimumThresholdVector.push_back(0.0);
            maximumThresholdCheck.push_back(false);
            minimumThresholdCheck.push_back(false);
        }
        // error checking
        if(maximumThresholdVector[i] < mSmallestMaximumThreshold)
        {
            maximumThresholdVector[i] = 0.0;
            maximumThresholdCheck[i] = false;
        }
    }



    // store the information
    AbstractChemistry* pThresholdChemistry = new AbstractChemistry();

    for(unsigned i=0; i<thresholdChemicalVector.size(); i++)
    {
        pThresholdChemistry -> AddChemical(new AbstractChemical(thresholdChemicalVector[i]));
    }
    

    SimpleChemicalThresholdCellCycleModel::SetUp(pThresholdChemistry);


    SetUp(pThresholdChemistry);
    SetMaximumSpeciesThreshold(maximumThresholdVector);
    SetMinimumSpeciesThreshold(minimumThresholdVector);
    SetNumberThresholdSpecies(thresholdChemicalVector.size());
    SetMaximumThresholdCheck(maximumThresholdCheck);
    SetMinimumThresholdCheck(minimumThresholdCheck);

}

std::vector<std::vector<std::string>> SimpleChemicalThresholdCellCycleFromFile::ReadMatrix(std::string filename)
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

std::vector<std::string> SimpleChemicalThresholdCellCycleFromFile::parseMatrixLineString(std::string line)
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


#endif