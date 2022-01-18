#ifndef CELLSTATETRACKINGMODIFIER_HPP_
#define CELLSTATETRACKINGMODIFIER_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/shared_ptr.hpp>

#include "AbstractCellBasedSimulationModifier.hpp"
#include "ChemicalSrnModel.hpp"
#include "AbstractCellPopulation.hpp"
#include "AbstractCellProperty.hpp"
#include "TransportCellProperty.hpp"
#include "MembraneCellProperty.hpp"
//#include "StateSwitchingCellProperty.hpp"

#include "ComplexCellFromFile.hpp"
#include "CustomCellFromFile.hpp"



template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class CellStateTrackingModifier : public AbstractCellBasedSimulationModifier<ELEMENT_DIM,SPACE_DIM>
{
private:
    double mDelta_error;

    bool mIsCheckConcentration;

public:

    CellStateTrackingModifier();

    virtual ~CellStateTrackingModifier();

    virtual void UpdateAtEndOfTimeStep(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation);

    virtual void SetupSolve(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation, std::string outputDirectory);

    void UpdateCellData(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation);


};

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
CellStateTrackingModifier<ELEMENT_DIM,SPACE_DIM>::CellStateTrackingModifier()
    : AbstractCellBasedSimulationModifier<SPACE_DIM>()
{

}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
CellStateTrackingModifier<ELEMENT_DIM,SPACE_DIM>::~CellStateTrackingModifier()
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellStateTrackingModifier<ELEMENT_DIM, SPACE_DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation)
{
    UpdateCellData(rCellPopulation);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellStateTrackingModifier<ELEMENT_DIM,SPACE_DIM>::UpdateCellData(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation)
{
    // run through each cell, if switch type copy cell data and reinitlialise the cell with the switched state
/*
    for (typename AbstractCellPopulation<SPACE_DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {   
        if (cell_iter-> rGetCellPropertyCollection(). template HasProperty<StateSwitchingCellProperty>())
        {
            boost::shared_ptr<StateSwitchingCellProperty> stateSwitching_cell_property = boost::static_pointer_cast<StateSwitchingCellProperty>(cell_iter-> template rGetCellPropertyCollection(). template GetPropertiesType<StateSwitchingCellProperty>().GetProperty());


            unsigned previous_cell_state = stateSwitching_cell_property -> GetCellState();

            stateSwitching_cell_property -> DetermineState();

            unsigned new_cell_state = stateSwitching_cell_property -> GetCellState();

            if(new_cell_state != previous_cell_state)
            {
                // renew cell state
                stateSwitching_cell_property -> SwitchState();

            }
        }
    }
    */
}


#endif 