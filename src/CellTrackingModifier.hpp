#ifndef CELLTRACKINGMODIFIER_HPP_
#define CELLTRACKINGMODIFIER_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractCellBasedSimulationModifier.hpp"
#include "CellSrnModel.hpp"
#include "AbstractCellPopulation.hpp"
#include "AbstractCellProperty.hpp"
#include "TransportCellProperty.hpp"


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class CellTrackingModifier : public AbstractCellBasedSimulationModifier<ELEMENT_DIM,SPACE_DIM>
{
public:

    CellTrackingModifier();

    virtual ~CellTrackingModifier();

    virtual void UpdateAtEndOfTimeStep(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation);

    virtual void SetupSolve(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation, std::string outputDirectory);


    void UpdateCellData(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation);

    void OutputSimulationModifierParameters(out_stream& rParamsFile);
};

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
CellTrackingModifier<ELEMENT_DIM,SPACE_DIM>::CellTrackingModifier()
    : AbstractCellBasedSimulationModifier<SPACE_DIM>()
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
CellTrackingModifier<ELEMENT_DIM,SPACE_DIM>::~CellTrackingModifier()
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellTrackingModifier<ELEMENT_DIM, SPACE_DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation)
{
    UpdateCellData(rCellPopulation);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellTrackingModifier<ELEMENT_DIM,SPACE_DIM>::SetupSolve(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation, std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    //std::cout<<"CellTrackingModifier - SetupSolve"<<std::endl;
    UpdateCellData(rCellPopulation);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellTrackingModifier<ELEMENT_DIM,SPACE_DIM>::UpdateCellData(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation)
{
   // std::cout<<"CellTrackingModifier - UpdateCellData"<<std::endl;
    // Make sure the cell population is updated
    rCellPopulation.Update();
    //std::cout<<"CellTrackingModifier - UpdateCellData - access transport property"<<std::endl;
    for (typename AbstractCellPopulation<SPACE_DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {

        CellPropertyCollection collection = cell_iter->rGetCellPropertyCollection();

        bool transport_active = collection.HasProperty<TransportCellProperty>();
        if (transport_active)
        {

            CellPropertyCollection transport_collection = collection.GetPropertiesType<TransportCellProperty>();
            boost::shared_ptr<TransportCellProperty> transportProp = boost::static_pointer_cast<TransportCellProperty>(transport_collection.GetProperty());
    
            std::vector<double> transportIn =  transportProp->GetTransportIn();
            std::vector<double> transportOut = transportProp->GetTransportOut();
            std::vector<std::string> speciesNameVector = transportProp->GetSpeciesNameVector();

            unsigned numberSpecies = transportIn.size();
            double tempCell=0;
            double tempNode=0;
            double cellVal=0;
            double nodeVal=0;

            for(unsigned species=0; species<numberSpecies; species++){
      
                cellVal = cell_iter->GetCellData()->GetItem("var_cell_"+speciesNameVector[species]);

                nodeVal = cell_iter->GetCellData()->GetItem("var_node_"+speciesNameVector[species]);
            
                // subtract transport out
                tempCell = -transportOut[species]*cellVal + transportIn[species]*nodeVal;
                tempNode = transportOut[species]*cellVal - transportIn[species]*nodeVal;
        
                // update
                cell_iter->GetCellData()->SetItem("var_cell_"+speciesNameVector[species], tempCell);
                cell_iter->GetCellData()->SetItem("var_node_"+speciesNameVector[species], tempNode);
            }
            
        }
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellTrackingModifier<ELEMENT_DIM,SPACE_DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<SPACE_DIM>::OutputSimulationModifierParameters(rParamsFile);
}

#endif 
