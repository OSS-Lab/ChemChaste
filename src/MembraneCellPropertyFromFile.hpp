#ifndef MEMBRANECELLPROPERTYFROMFILE_HPP
#define MEMBRANECELLPROPERTYFROMFILE_HPP

//general includes
#include <vector>
#include <boost/shared_ptr.hpp>

#include "MembraneCellProperty.hpp"
#include "ChemicalCell.hpp"
#include "AbstractMembraneReactionSystemFromFile.hpp"

class MembraneCellPropertyFromFile
{
protected:

    std::string mMembraneFilename;

    boost::shared_ptr<MembraneCellProperty> mpMembraneProperty;

public:

    MembraneCellPropertyFromFile(std::string filename ="");

    ~MembraneCellPropertyFromFile()
    {
    };


    void SetUpMembraneProperty(CellPtr);


    void SetMembraneFilename(std::string);

    void SetMembraneProperty(boost::shared_ptr<MembraneCellProperty>);


    std::string GetMembraneFilename();

    boost::shared_ptr<MembraneCellProperty> GetMembraneProperty();

};

    MembraneCellPropertyFromFile::MembraneCellPropertyFromFile(std::string filename)
        : mMembraneFilename(filename)
    {

        if(filename != "")
        {
            SetMembraneFilename(filename);
        }

    }

    void MembraneCellPropertyFromFile::SetUpMembraneProperty(CellPtr p_cell)
    {
        boost::shared_ptr<MembraneCellProperty> p_membrane = boost::static_pointer_cast<MembraneCellProperty>(p_cell->rGetCellPropertyCollection().GetPropertiesType<MembraneCellProperty>().GetProperty());
        AbstractMembraneReactionSystemFromFile* p_membrane_system_from_file = new AbstractMembraneReactionSystemFromFile(mMembraneFilename);

  
        p_membrane -> SetUp(p_membrane_system_from_file, p_cell);
        SetMembraneProperty(p_membrane);
    }

    void MembraneCellPropertyFromFile::SetMembraneFilename(std::string filename)
    {
        mMembraneFilename = filename;
    }

    void MembraneCellPropertyFromFile::SetMembraneProperty(boost::shared_ptr<MembraneCellProperty> pMembraneProperty)
    {
        mpMembraneProperty = pMembraneProperty;
    }

    std::string MembraneCellPropertyFromFile::GetMembraneFilename()
    {
        return mMembraneFilename;
    }

    boost::shared_ptr<MembraneCellProperty> MembraneCellPropertyFromFile::GetMembraneProperty()
    {
        return mpMembraneProperty;
    }

#endif