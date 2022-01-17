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
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
CellNeighbourTrackingModifier<ELEMENT_DIM,SPACE_DIM>::~CellNeighbourTrackingModifier()
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellNeighbourTrackingModifier<ELEMENT_DIM, SPACE_DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation)
{
    //std::cout<<"CellNeighbourTrackingModifier<ELEMENT_DIM,SPACE_DIM>::UpdateAtEndOfTimeStep - start"<<std::endl;
    UpdateCellData(rCellPopulation);
    //std::cout<<"CellNeighbourTrackingModifier<ELEMENT_DIM,SPACE_DIM>::UpdateAtEndOfTimeStep - end"<<std::endl;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellNeighbourTrackingModifier<ELEMENT_DIM,SPACE_DIM>::SetupSolve(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation, std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    UpdateCellData(rCellPopulation);
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

    double meanCellDiversityValue=0.0;
    unsigned numberOfCells=0;

    /// run through each cell in the population and update the internal chemical concentrations in cellData
    for (typename AbstractCellPopulation<SPACE_DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {   
        // update the cell data based on any transport properties
        if (cell_iter-> rGetCellPropertyCollection(). template HasProperty<CellAnalyticsProperty>())
        {
            boost::shared_ptr<CellAnalyticsProperty> analytics_cell_property = boost::static_pointer_cast<CellAnalyticsProperty>(cell_iter-> template rGetCellPropertyCollection(). template GetPropertiesType<CellAnalyticsProperty>().GetProperty());

             // Get the set of neighbouring location indices
            std::set<unsigned> neighbour_indices = rCellPopulation.GetNeighbouringLocationIndices(*cell_iter);

            std::string cellTypeName = "";
            std::string thisCellTypeName = analytics_cell_property -> GetCellTypeName();
            double cellDiversityValue = 0.0;
            std::vector<unsigned> neighbourCellCount(mNumberOfCellTypes,0);
            std::vector<double> neighbourTypeProbability(mNumberOfCellTypes,0.0);
            unsigned numberOfNeighbours=0;

            if (!neighbour_indices.empty())
            {
                for (std::set<unsigned>::iterator iter = neighbour_indices.begin();
                    iter != neighbour_indices.end();
                    ++iter)
                {
                    CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(*iter);
                    if (p_cell-> rGetCellPropertyCollection(). HasProperty<CellAnalyticsProperty>())
                    {
                        boost::shared_ptr<CellAnalyticsProperty> neighbour_analytics_cell_property = boost::static_pointer_cast<CellAnalyticsProperty>(p_cell-> rGetCellPropertyCollection(). GetPropertiesType<CellAnalyticsProperty>().GetProperty());

                        numberOfNeighbours +=1;

                        cellTypeName = neighbour_analytics_cell_property -> GetCellTypeName();

                        for(unsigned i=0; i<mNumberOfCellTypes; i++)
                        {
                            if(cellTypeName == mCellTypeNameDictionary[i])
                            {
                                neighbourCellCount[i] +=1;
                            }
                        }
                    }
                }
            }
            
            if(numberOfNeighbours==0)
            {
                for(unsigned i=0; i<mNumberOfCellTypes; i++)
                {
                    neighbourTypeProbability[i] = 0.0;
                }
            }
            else
            {
                for(unsigned i=0; i<mNumberOfCellTypes; i++)
                {
                    neighbourTypeProbability[i] = (double) neighbourCellCount[i]/ (double) numberOfNeighbours;
                }
            }
            
            for(unsigned i=0; i<mNumberOfCellTypes; i++)
            {
                if(thisCellTypeName == mCellTypeNameDictionary[i])
                {
                    cellDiversityValue = neighbourTypeProbability[i];
                }
            }

            analytics_cell_property -> SetNeighbourTypes(neighbourCellCount);
            analytics_cell_property -> SetNeighbourTypeProbability(neighbourTypeProbability);
            analytics_cell_property -> SetCellDiversity(cellDiversityValue);
            numberOfCells +=1;
            meanCellDiversityValue += cellDiversityValue;
        }
    }
    for (typename AbstractCellPopulation<SPACE_DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {   
        // update the cell data based on any transport properties
        if (cell_iter-> rGetCellPropertyCollection(). template HasProperty<CellAnalyticsProperty>())
        {
            boost::shared_ptr<CellAnalyticsProperty> analytics_cell_property = boost::static_pointer_cast<CellAnalyticsProperty>(cell_iter-> template rGetCellPropertyCollection(). template GetPropertiesType<CellAnalyticsProperty>().GetProperty());

            if(numberOfCells == 0)
            {
                analytics_cell_property -> SetMeanCellDiversity(0.0);
            }
            else
            {
            analytics_cell_property -> SetMeanCellDiversity(meanCellDiversityValue/ (double) numberOfCells);
            }
        }
    }
    //std::cout<<"CellNeighbourTrackingModifier<ELEMENT_DIM,SPACE_DIM>::UpdateCellData - start"<<std::endl;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellNeighbourTrackingModifier<ELEMENT_DIM,SPACE_DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<SPACE_DIM>::OutputSimulationModifierParameters(rParamsFile);
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
