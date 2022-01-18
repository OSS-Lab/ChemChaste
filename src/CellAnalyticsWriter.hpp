#ifndef CELLANALYTICSWRITER_HPP_
#define CELLANALYTICSWRITER_HPP_

#include "AbstractCellWriter.hpp"
#include "CellAnalyticsProperty.hpp"
#include "AbstractCellPopulation.hpp"
#include <boost/serialization/base_object.hpp>


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class CellAnalyticsWriter : public AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>
{

public:

    CellAnalyticsWriter();
    

    double GetCellDataForVtkOutput(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation);
    

    virtual void VisitCell(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation);
    
};

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
CellAnalyticsWriter<ELEMENT_DIM, SPACE_DIM>::CellAnalyticsWriter()
    : AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>("cellanalyticsresults.dat")
{
    this->mVtkCellDataName = "Strain type";
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double CellAnalyticsWriter<ELEMENT_DIM, SPACE_DIM>::GetCellDataForVtkOutput(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    if(pCell->rGetCellPropertyCollection().HasProperty<CellAnalyticsProperty>())
    {

        boost::shared_ptr<CellAnalyticsProperty> cellAnalyticsProperty = boost::static_pointer_cast<CellAnalyticsProperty>(pCell->rGetCellPropertyCollection().GetPropertiesType<CellAnalyticsProperty>().GetProperty());

        unsigned typeID = cellAnalyticsProperty -> GetCellID();

        return typeID;
    }

    return 0;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellAnalyticsWriter<ELEMENT_DIM, SPACE_DIM>::VisitCell(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{

    if(pCell->rGetCellPropertyCollection().HasProperty<CellAnalyticsProperty>())
    {
        boost::shared_ptr<CellAnalyticsProperty> cellAnalyticsProperty = boost::static_pointer_cast<CellAnalyticsProperty>(pCell->rGetCellPropertyCollection().GetPropertiesType<CellAnalyticsProperty>().GetProperty());

        unsigned typeID = cellAnalyticsProperty -> GetCellID();

        *this->mpOutStream << typeID << " ";
    }else
    {
        *this->mpOutStream << 0 << " ";
    }

    
}






#endif
