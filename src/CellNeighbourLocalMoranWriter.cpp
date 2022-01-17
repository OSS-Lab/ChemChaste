#include "CellNeighbourLocalMoranWriter.hpp"
#include "AbstractCellPopulation.hpp"
#include "CellAnalyticsProperty.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
CellNeighbourLocalMoranWriter<ELEMENT_DIM, SPACE_DIM>::CellNeighbourLocalMoranWriter()
    : AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>("cellNeighbourLocalMoran.dat")
{
    this->mVtkCellDataName = "LocalMoranIndex";
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double CellNeighbourLocalMoranWriter<ELEMENT_DIM, SPACE_DIM>::GetCellDataForVtkOutput(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    return CalculateNeighbourLocalMoranIndex(pCell, pCellPopulation);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellNeighbourLocalMoranWriter<ELEMENT_DIM, SPACE_DIM>::VisitCell(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    // Write location index corresponding to cell
    *this->mpOutStream << pCellPopulation->GetLocationIndexUsingCell(pCell) << " ";

    // Write cell location
    c_vector<double, SPACE_DIM> cell_location = pCellPopulation->GetLocationOfCellCentre(pCell);
    for (unsigned i=0; i<SPACE_DIM; i++)
    {
        *this->mpOutStream << cell_location[i] << " ";
    }

    // Write cell age
    *this->mpOutStream << CalculateNeighbourLocalMoranIndex(pCell, pCellPopulation) << " ";
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double CellNeighbourLocalMoranWriter<ELEMENT_DIM, SPACE_DIM>::CalculateNeighbourLocalMoranIndex(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    double localMoranIndex = 0.0;
    unsigned numberOfCells=0;
    double weightingFactor=1.0;

    if (pCell-> rGetCellPropertyCollection(). template HasProperty<CellAnalyticsProperty>())
    {
        boost::shared_ptr<CellAnalyticsProperty> analytics_cell_property = boost::static_pointer_cast<CellAnalyticsProperty>(pCell-> rGetCellPropertyCollection(). GetPropertiesType<CellAnalyticsProperty>().GetProperty());
    
        std::vector<unsigned> cellTypeCounts = analytics_cell_property-> GetNeighbourTypes();
    
        double meanCellDiversity = analytics_cell_property-> GetMeanCellDiversity();
        
        double neighbourCellDiversity = 0.0;

        std::set<unsigned> neighbour_indices = pCellPopulation -> GetNeighbouringLocationIndices(pCell);

        if (!neighbour_indices.empty())
        {
            for (std::set<unsigned>::iterator iter = neighbour_indices.begin();
                iter != neighbour_indices.end();
                ++iter)
            {
                CellPtr p_cell = pCellPopulation ->GetCellUsingLocationIndex(*iter);

                if (p_cell-> rGetCellPropertyCollection(). HasProperty<CellAnalyticsProperty>())
                {
                    boost::shared_ptr<CellAnalyticsProperty> neighbour_analytics_cell_property = boost::static_pointer_cast<CellAnalyticsProperty>(p_cell-> rGetCellPropertyCollection(). GetPropertiesType<CellAnalyticsProperty>().GetProperty());

                    neighbourCellDiversity = neighbour_analytics_cell_property -> GetCellDiversity();

            
                    localMoranIndex += weightingFactor*(neighbourCellDiversity - meanCellDiversity);

                
                }
            }
        }

        localMoranIndex = localMoranIndex*(analytics_cell_property -> GetCellDiversity() - meanCellDiversity);
    }

    return localMoranIndex;
}

// Explicit instantiation
template class CellNeighbourLocalMoranWriter<1,1>;
template class CellNeighbourLocalMoranWriter<1,2>;
template class CellNeighbourLocalMoranWriter<2,2>;
template class CellNeighbourLocalMoranWriter<1,3>;
template class CellNeighbourLocalMoranWriter<2,3>;
template class CellNeighbourLocalMoranWriter<3,3>;
