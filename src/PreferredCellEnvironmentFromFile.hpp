#ifndef PREFERREDCELLENVRIONMENTFROMFILE_HPP_
#define PREFERREDCELLENVRIONMENTFROMFILE_HPP_

// general includes
#include <string>
#include <tuple>
#include <vector>
#include <iostream>
#include <stdlib.h> 

class PreferredCellEnvironmentFromFile
{
protected:

    std::vector<std::string> mEnvironmentNamesVector;

    std::vector<double> mValueVector;

public:

    PreferredCellEnvironmentFromFile(std::string);

    virtual ~PreferredCellEnvironmentFromFile()
    {
    };

    std::vector<std::vector<std::string>> ReadMatrix(std::string);

    std::vector<std::string> parseMatrixLineString(std::string);

    void SetEnvironmentNamesVector(std::vector<std::string>);

    void SetEnvironmentValueVector(std::vector<double> );

    std::vector<std::string> GetEnvironmentNamesVector();

    std::vector<double> GetEnvironmentValueVector();
};

PreferredCellEnvironmentFromFile::PreferredCellEnvironmentFromFile(std::string environmentFilename) 
{
    std::vector<std::string> environmentNamesVector;
    std::vector<double> valueVector;

    // read and store the diffusion data as a matrix of strings 
    std::vector<std::vector<std::string>> environmentDatabase;
    environmentDatabase = ReadMatrix(environmentFilename);

    // parse the matrix
    for(unsigned i=0; i<environmentDatabase.size(); i++)
    {
        // the first element of the read in database is the identifying names
        environmentNamesVector.push_back(environmentDatabase[i][0]);
        // data for species concentration
        valueVector.push_back(atof(environmentDatabase[i][1].c_str()));
    }


    // store the information
    SetEnvironmentNamesVector(environmentNamesVector);
    SetEnvironmentValueVector(valueVector);
}

std::vector<std::vector<std::string>> PreferredCellEnvironmentFromFile::ReadMatrix(std::string filename)
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

std::vector<std::string> PreferredCellEnvironmentFromFile::parseMatrixLineString(std::string line)
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


void PreferredCellEnvironmentFromFile::SetEnvironmentNamesVector(std::vector<std::string> namesVector)
{
    mEnvironmentNamesVector = namesVector;
}

void PreferredCellEnvironmentFromFile::SetEnvironmentValueVector(std::vector<double> valueVector)
{
    mValueVector = valueVector;
}

std::vector<std::string> PreferredCellEnvironmentFromFile::GetEnvironmentNamesVector()
{
    return mEnvironmentNamesVector;
}

std::vector<double> PreferredCellEnvironmentFromFile::GetEnvironmentValueVector()
{
    return mValueVector;
}


#endif