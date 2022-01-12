#ifndef CELLNEIGHBOURTRACKINGMODIFIER_HPP_
#define CELLNEIGHBOURTRACKINGMODIFIER_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/shared_ptr.hpp>

#include "AbstractCellBasedSimulationModifier.hpp"
#include "AbstractCellPopulation.hpp"
#include "AbstractCellProperty.hpp"
#include <vector>
#include <string>


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class CellNeighbourTrackingModifier : public AbstractCellBasedSimulationModifier<ELEMENT_DIM,SPACE_DIM>
{
private:
    
    std::vector<std::string> mCellTypeNameDictionary;

    unsigned mNumberOfCellTypes;

public:

    CellNeighbourTrackingModifier(std::vector<std::string>);

    virtual ~CellNeighbourTrackingModifier();

    virtual void UpdateAtEndOfTimeStep(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation);

    virtual void SetupSolve(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation, std::string outputDirectory);

    void UpdateCellData(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation);

    void OutputSimulationModifierParameters(out_stream& rParamsFile);

    void SetCellTypeNameDictionary(std::vector<std::string>);

    std::vector<std::string> GetCellTypeNameDictionary();

};

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
CellNeighbourTrackingModifier<ELEMENT_DIM,SPACE_DIM>::CellNeighbourTrackingModifier(std::vector<std::string> cellTypeNameDictionary)
    : AbstractCellBasedSimulationModifier<SPACE_DIM>(),
    mCellTypeNameDictionary(cellTypeNameDictionary),
    mNumberOfCellTypes(cellTypeNameDictionary.size())
{
    std::cout<<"CellNeighbourTrackingModifier<ELEMENT_DIM,SPACE_DIM>::CellNeighbourTrackingModifier - start"<<std::endl;
    std::cout<<"CellNeighbourTrackingModifier<ELEMENT_DIM,SPACE_DIM>::CellNeighbourTrackingModifier - end"<<std::endl;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
CellNeighbourTrackingModifier<ELEMENT_DIM,SPACE_DIM>::~CellNeighbourTrackingModifier()
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellNeighbourTrackingModifier<ELEMENT_DIM, SPACE_DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation)
{
    std::cout<<"CellNeighbourTrackingModifier<ELEMENT_DIM,SPACE_DIM>::UpdateAtEndOfTimeStep - start"<<std::endl;
    UpdateCellData(rCellPopulation);
    std::cout<<"CellNeighbourTrackingModifier<ELEMENT_DIM,SPACE_DIM>::UpdateAtEndOfTimeStep - end"<<std::endl;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellNeighbourTrackingModifier<ELEMENT_DIM,SPACE_DIM>::SetupSolve(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation, std::string outputDirectory)
{
    std::cout<<"CellNeighbourTrackingModifier<ELEMENT_DIM,SPACE_DIM>::SetupSolve - start"<<std::endl;
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    UpdateCellData(rCellPopulation);
     std::cout<<"CellNeighbourTrackingModifier<ELEMENT_DIM,SPACE_DIM>::SetupSolve - end"<<std::endl;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellNeighbourTrackingModifier<ELEMENT_DIM,SPACE_DIM>::UpdateCellData(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation)
{
    //std::cout<<"CellNeighbourTrackingModifier<ELEMENT_DIM,SPACE_DIM>::UpdateCellData - start"<<std::endl;
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    rCellPopulation.Update();

    /// run through each cell in the population and update the internal chemical concentrations in cellData
    for (typename AbstractCellPopulation<SPACE_DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {   
        //std::cout<<"here0"<<std::endl;
        // update the cell data based on any transport properties
        if (cell_iter-> rGetCellPropertyCollection(). template HasProperty<CellAnalyticsProperty>())
        {
            boost::shared_ptr<CellAnalyticsProperty> analytics_cell_property = boost::static_pointer_cast<CellAnalyticsProperty>(cell_iter-> template rGetCellPropertyCollection(). template GetPropertiesType<CellAnalyticsProperty>().GetProperty());
            //std::cout<<"here1"<<std::endl;
             // Get the set of neighbouring location indices
            std::set<unsigned> neighbour_indices = rCellPopulation.GetNeighbouringLocationIndices(*cell_iter);

            std::string cellTypeName = "";
            std::vector<unsigned> neighbourCellCount(mNumberOfCellTypes,0);
            //std::cout<<"here2"<<std::endl;
            if (!neighbour_indices.empty())
            {
                //std::cout<<"here3"<<std::endl;
                for (std::set<unsigned>::iterator iter = neighbour_indices.begin();
                    iter != neighbour_indices.end();
                    ++iter)
                {
                    CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(*iter);
                    //std::cout<<"here4"<<std::endl;
                    if (p_cell-> rGetCellPropertyCollection(). HasProperty<CellAnalyticsProperty>())
                    {
                        boost::shared_ptr<CellAnalyticsProperty> neighbour_analytics_cell_property = boost::static_pointer_cast<CellAnalyticsProperty>(p_cell-> rGetCellPropertyCollection(). GetPropertiesType<CellAnalyticsProperty>().GetProperty());

                        cellTypeName = neighbour_analytics_cell_property -> GetCellTypeName();
                        //std::cout<<"here5"<<std::endl;
                        for(unsigned i=0; i<mNumberOfCellTypes; i++)
                        {
                            //std::cout<<"here5.1"<<std::endl;
                            //std::cout<<mNumberOfCellTypes<<std::endl;
                            //std::cout<<cellTypeName<<std::endl;
                            //std::cout<<mCellTypeNameDictionary[i]<<std::endl;
                            if(cellTypeName == mCellTypeNameDictionary[i])
                            {
                                //std::cout<<"here5.2"<<std::endl;
                                neighbourCellCount[i] +=1;
                            }
                        }
                        //std::cout<<"here6"<<std::endl;
                    }
                }
            }
            //std::cout<<"here7"<<std::endl;
            analytics_cell_property -> SetNeighbourTypes(neighbourCellCount);
            //std::cout<<"here8"<<std::endl;
        }
    }
    //std::cout<<"CellNeighbourTrackingModifier<ELEMENT_DIM,SPACE_DIM>::UpdateCellData - start"<<std::endl;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellNeighbourTrackingModifier<ELEMENT_DIM,SPACE_DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    std::cout<<"CellNeighbourTrackingModifier<ELEMENT_DIM,SPACE_DIM>::void CellNeighbourTrackingModifier<ELEMENT_DIM,SPACE_DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)- start"<<std::endl;
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<SPACE_DIM>::OutputSimulationModifierParameters(rParamsFile);
    std::cout<<"CellNeighbourTrackingModifier<ELEMENT_DIM,SPACE_DIM>::void CellNeighbourTrackingModifier<ELEMENT_DIM,SPACE_DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)- end"<<std::endl;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellNeighbourTrackingModifier<ELEMENT_DIM,SPACE_DIM>::SetCellTypeNameDictionary(std::vector<std::string> cellTypeNameDictionary)
{
    mCellTypeNameDictionary = cellTypeNameDictionary;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<std::string> CellNeighbourTrackingModifier<ELEMENT_DIM,SPACE_DIM>::GetCellTypeNameDictionary()
{
    return mCellTypeNameDictionary;
}


#endif 