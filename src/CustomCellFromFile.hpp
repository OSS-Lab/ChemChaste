#ifndef CUSTOMCELLFROMFILE_HPP_
#define CUSTOMCELLFROMFILE_HPP_

#include "ComplexCellFromFile.hpp"
#include <vector>
#include <string>
#include <boost/utility.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/enable_shared_from_this.hpp>

class CustomCellFromFile : public ComplexCellFromFile
{
protected:
    //using ComplexCellFromFile::Divide;

    std::string mConfigurationFilename = "";

    std::string mConfigurationDelimiter = " = ";
    std::string mConfigurationDelimiterShort="=";

    std::string mCellModelType = "";

    std::vector<std::vector<std::string>> mConfigurationParameters;




public:

    CustomCellFromFile(     std::string cellCycleFilename="",
                            std::string divisionRulesFilename="", 
                            std::string srnFilename="",
                            std::string initialConditionFilename="",
                            std::string transportPropertyFilename="",
                            std::string membranePropertyFilename="",
                            std::string configurationFilename="",
                            unsigned cellID =0,
                            bool isCellIDSet = false);

    virtual ~CustomCellFromFile()
    {
    };

    void ParseConfigurationFile();

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

    

};

CustomCellFromFile::CustomCellFromFile(
                        std::string cellCycleFilename,
                        std::string divisionRulesFilename,
                        std::string srnFilename,
                        std::string initialConditionFilename,
                        std::string transportPropertyFilename,
                        std::string membranePropertyFilename,
                        std::string configurationFilename,
                        unsigned cellID,
                        bool isCellIDSet)
        :   ComplexCellFromFile(cellCycleFilename,
                                divisionRulesFilename,
                                srnFilename,
                                initialConditionFilename,
                                transportPropertyFilename,
                                membranePropertyFilename,
                                cellID,
                                isCellIDSet),
            mConfigurationFilename(configurationFilename)
{
    ParseConfigurationFile();

    if(IsParameterInConfiguration("split_ratio"))
    {
        boost::shared_ptr<ComplexCell> cellPtr = boost::dynamic_pointer_cast<ComplexCell>(mpCell);
        cellPtr->SetSplitRatio(ReturnDoubleFromConfiguration("split_ratio"));
        //boost::shared_ptr<Cell> 
        //mpCell->SetSplitRatio(ReturnDoubleFromConfiguration("split_ratio"));
    }

    if(IsParameterInConfiguration("cell_model"))
    {
        SetCellModelType(ReturnStringFromConfiguration("cell_model"));
    }

}

void CustomCellFromFile::ParseConfigurationFile()
{


    std::ifstream inputFile(mConfigurationFilename);

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

//std::vector<std::vector<std::string>> CustomCellFromFile::GetConfigurationParameters()
//{
//    return mConfigurationParameters;
//}

void CustomCellFromFile::SetCellModelType(std::string modelType)
{
    mCellModelType = modelType;
}

std::string CustomCellFromFile::GetCellModelType()
{
    return mCellModelType;
}

#endif