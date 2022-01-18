#ifndef CELLSTATEWRITER_HPP_
#define CELLSTATEWRITER_HPP_

#include "AbstractCellWriter.hpp"
//#include "CellStateSwitchingProperty.hpp"
#include "AbstractCellPopulation.hpp"
#include <boost/serialization/base_object.hpp>


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class CellStateWriter : public AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>
{

public:

    CellStateWriter();
    

    double GetCellDataForVtkOutput(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation);
    

    virtual void VisitCell(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation);
    
};

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
CellStateWriter<ELEMENT_DIM, SPACE_DIM>::CellStateWriter()
    : AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>("cellstateresults.dat")
{
    this->mVtkCellDataName = "Cell State";
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double CellStateWriter<ELEMENT_DIM, SPACE_DIM>::GetCellDataForVtkOutput(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    /*
    if(pCell->rGetCellPropertyCollection().HasProperty<CellStateSwitchingProperty>())
    {

        boost::shared_ptr<CellStateSwitchingProperty> cellStateSwitchingProperty = boost::static_pointer_cast<CellStateSwitchingProperty>(pCell->rGetCellPropertyCollection().GetPropertiesType<CellStateSwitchingProperty>().GetProperty());

        unsigned stateID = cellStateSwitchingProperty -> GetCellStateID();

        return stateID;
    }
*/
    return 0;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellStateWriter<ELEMENT_DIM, SPACE_DIM>::VisitCell(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
/*
    if(pCell->rGetCellPropertyCollection().HasProperty<CellStateSwitchingProperty>())
    {

        //boost::shared_ptr<CellStateSwitchingProperty> CellStateSwitchingProperty = boost::static_pointer_cast<CellStateSwitchingProperty>(pCell->rGetCellPropertyCollection().GetPropertiesType<CellStateSwitchingProperty>().GetProperty());
        boost::shared_ptr<CellStateSwitchingProperty> cellStateSwitchingProperty = boost::static_pointer_cast<CellStateSwitchingProperty>(pCell->rGetCellPropertyCollection().GetPropertiesType<CellStateSwitchingProperty>().GetProperty());


        unsigned stateID = cellStateSwitchingProperty -> GetCellStateID();

        *this->mpOutStream << stateID << " ";
    }else
    {
*/
        *this->mpOutStream << 0 << " ";
//    }

    
}






#endif