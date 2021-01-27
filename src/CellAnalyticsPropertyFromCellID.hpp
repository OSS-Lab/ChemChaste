#ifndef CELLANALYTICSPROPERTYFROMCELLID_HPP
#define CELLANALYTICSPROPERTYFROMCELLID_HPP

//general includes
#include <vector>
#include <boost/shared_ptr.hpp>

#include "CellAnalyticsProperty.hpp"
#include "ChemicalCell.hpp"
#include "AbstractTransportReactionSystemFromFile.hpp"

class CellAnalyticsPropertyFromCellID
{
protected:


    unsigned mCellID;

    boost::shared_ptr<CellAnalyticsProperty> mpCellAnalyticsProperty;


public:

    CellAnalyticsPropertyFromCellID(unsigned);

    ~CellAnalyticsPropertyFromCellID()
    {
    };


    void SetUpCellAnalyticsProperty(CellPtr);


    void SetCellID(unsigned);

    void SetCellAnalyticsProperty(boost::shared_ptr<CellAnalyticsProperty>);


    unsigned GetCellID();

    boost::shared_ptr<CellAnalyticsProperty> GetCellAnalyticsProperty();

};

    CellAnalyticsPropertyFromCellID::CellAnalyticsPropertyFromCellID(unsigned cellID)
        : mCellID(cellID)
    {

        SetCellID(cellID);

    }

    void CellAnalyticsPropertyFromCellID::SetUpCellAnalyticsProperty(CellPtr p_cell)
    {
        boost::shared_ptr<CellAnalyticsProperty> p_analytics = boost::static_pointer_cast<CellAnalyticsProperty>(p_cell->rGetCellPropertyCollection().GetPropertiesType<CellAnalyticsProperty>().GetProperty());
       
        p_analytics -> SetUp(p_cell,mCellID);

        SetCellAnalyticsProperty(p_analytics);
    }

    void CellAnalyticsPropertyFromCellID::SetCellID(unsigned cellID)
    {
        mCellID = cellID;
    }

    void CellAnalyticsPropertyFromCellID::SetCellAnalyticsProperty(boost::shared_ptr<CellAnalyticsProperty> pCellAnalyticsProperty)
    {
        mpCellAnalyticsProperty = pCellAnalyticsProperty;
    }

    unsigned CellAnalyticsPropertyFromCellID::GetCellID()
    {
        return mCellID;
    }

    boost::shared_ptr<CellAnalyticsProperty> CellAnalyticsPropertyFromCellID::GetCellAnalyticsProperty()
    {
        return mpCellAnalyticsProperty;
    }


#endif