#ifndef ENVIRONMENTCELLPROPERTYFROMFILE_HPP
#define ENVIRONMENTCELLPROPERTYFROMFILE_HPP

//general includes
#include <vector>
#include <boost/shared_ptr.hpp>

#include "EnvironmentCellProperty.hpp"
#include "ChemicalCell.hpp"


class EnvironmentCellPropertyFromFile
{
protected:

    std::string mEnvironmentFilename;

    boost::shared_ptr<EnvironmentCellProperty> mpEnvironmentProperty;


public:

    EnvironmentCellPropertyFromFile(std::string filename ="");

    ~EnvironmentCellPropertyFromFile()
    {
    };


    void SetUpEnvironmentProperty(CellPtr);


    void SetEnvironmentFilename(std::string);

    void SetEnvironmentProperty(boost::shared_ptr<EnvironmentCellProperty>);


    std::string GetEnvironmentFilename();

    boost::shared_ptr<EnvironmentCellProperty> GetEnvironmentProperty();

};

    EnvironmentCellPropertyFromFile::EnvironmentCellPropertyFromFile(std::string filename)
        : mEnvironmentFilename(filename)
    {

        if(filename != "")
        {
            SetEnvironmentFilename(filename);
        }

    }

    void EnvironmentCellPropertyFromFile::SetUpEnvironmentProperty(CellPtr p_cell)
    {
        boost::shared_ptr<EnvironmentCellProperty> p_environment = boost::static_pointer_cast<EnvironmentCellProperty>(p_cell->rGetCellPropertyCollection().GetPropertiesType<EnvironmentCellProperty>().GetProperty());
        
        
        PreferredCellEnvironmentFromFile* p_preferred_environment_from_file = new PreferredCellEnvironmentFromFile(mEnvironmentFilename);

        p_environment -> SetUp(p_cell, p_preferred_environment_from_file->GetEnvironmentStateVariableRegister()->GetStateVariableRegisterVector());
        p_environment -> SetEnvironmentVector(p_preferred_environment_from_file->GetEnvironmentValueVector());
        std::cout<<"Environment property from file"<<std::endl;
        SetEnvironmentProperty(p_environment);
    }

    void EnvironmentCellPropertyFromFile::SetEnvironmentFilename(std::string filename)
    {
        mEnvironmentFilename = filename;
    }

    void EnvironmentCellPropertyFromFile::SetEnvironmentProperty(boost::shared_ptr<EnvironmentCellProperty> pEnvironmentProperty)
    {
        mpEnvironmentProperty = pEnvironmentProperty;
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