#ifndef TRANSPORTCELLPROPERTYFROMFILE_HPP
#define TRANSPORTCELLPROPERTYFROMFILE_HPP

//general includes
#include <vector>
#include <boost/shared_ptr.hpp>

#include "TransportCellProperty.hpp"
#include "ChemicalCell.hpp"
#include "AbstractTransportReactionSystemFromFile.hpp"

class TransportCellPropertyFromFile
{
protected:

    std::string mTransportFilename;

    boost::shared_ptr<TransportCellProperty> mpTransportProperty;


public:

    TransportCellPropertyFromFile(std::string filename ="");

    ~TransportCellPropertyFromFile()
    {
    };


    void SetUpTransportProperty(CellPtr);


    void SetTransportFilename(std::string);

    void SetTransportProperty(boost::shared_ptr<TransportCellProperty>);


    std::string GetTransportFilename();

    boost::shared_ptr<TransportCellProperty> GetTransportProperty();

};

    TransportCellPropertyFromFile::TransportCellPropertyFromFile(std::string filename)
        : mTransportFilename(filename)
    {

        if(filename != "")
        {
            SetTransportFilename(filename);
        }

    }

    void TransportCellPropertyFromFile::SetUpTransportProperty(CellPtr p_cell)
    {
        boost::shared_ptr<TransportCellProperty> p_transport = boost::static_pointer_cast<TransportCellProperty>(p_cell->rGetCellPropertyCollection().GetPropertiesType<TransportCellProperty>().GetProperty());
        AbstractTransportReactionSystemFromFile* p_transport_system_from_file = new AbstractTransportReactionSystemFromFile(mTransportFilename);
        
        //boost::shared_ptr<TransportCellProperty> pTransportProperty(new TransportCellProperty());
        p_transport -> SetUp(p_transport_system_from_file, p_cell);
        std::cout<<"Transport property from file"<<std::endl;
        SetTransportProperty(p_transport);
    }

    void TransportCellPropertyFromFile::SetTransportFilename(std::string filename)
    {
        mTransportFilename = filename;
    }

    void TransportCellPropertyFromFile::SetTransportProperty(boost::shared_ptr<TransportCellProperty> pTransportProperty)
    {
        mpTransportProperty = pTransportProperty;
    }

    std::string TransportCellPropertyFromFile::GetTransportFilename()
    {
        return mTransportFilename;
    }

    boost::shared_ptr<TransportCellProperty> TransportCellPropertyFromFile::GetTransportProperty()
    {
        return mpTransportProperty;
    }


#endif