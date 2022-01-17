#include "CellNeighbourGetisOrdWriter.hpp"
#include "AbstractCellPopulation.hpp"
#include "CellAnalyticsProperty.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
CellNeighbourGetisOrdWriter<ELEMENT_DIM, SPACE_DIM>::CellNeighbourGetisOrdWriter()
    : AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>("cellNeighbourGetisOrd.dat")
{
    this->mVtkCellDataName = "GetisOrdIndex";
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double CellNeighbourGetisOrdWriter<ELEMENT_DIM, SPACE_DIM>::GetCellDataForVtkOutput(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    return CalculateNeighbourGetisOrdIndex(pCell, pCellPopulation);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellNeighbourGetisOrdWriter<ELEMENT_DIM, SPACE_DIM>::VisitCell(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    //std::cout<<"CellNeighbourGetisOrdWriter<ELEMENT_DIM, SPACE_DIM>::VisitCell - start"<<std::endl;
    // Write location index corresponding to cell
    *this->mpOutStream << pCellPopulation->GetLocationIndexUsingCell(pCell) << " ";

    // Write cell location
    c_vector<double, SPACE_DIM> cell_location = pCellPopulation->GetLocationOfCellCentre(pCell);
    for (unsigned i=0; i<SPACE_DIM; i++)
    {
        *this->mpOutStream << cell_location[i] << " ";
    }

    // Write cell age
    *this->mpOutStream << CalculateNeighbourGetisOrdIndex(pCell, pCellPopulation) << " ";
    //std::cout<<"CellNeighbourGetisOrdWriter<ELEMENT_DIM, SPACE_DIM>::VisitCell - end"<<std::endl;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double CellNeighbourGetisOrdWriter<ELEMENT_DIM, SPACE_DIM>::CalculateNeighbourGetisOrdIndex(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    std::cout<<"CellNeighbourGetisOrdWriter<ELEMENT_DIM, SPACE_DIM>::CalculateNeighbourGetisOrdIndex - start"<<std::endl;
    double getisOrdIndex = 0.0;
    unsigned numberOfCells=0;
    double weightingFactor=1.0;
    double sumWeightingFactor=0.0;

    if (pCell-> rGetCellPropertyCollection(). template HasProperty<CellAnalyticsProperty>())
    {
        boost::shared_ptr<CellAnalyticsProperty> analytics_cell_property = boost::static_pointer_cast<CellAnalyticsProperty>(pCell-> rGetCellPropertyCollection(). GetPropertiesType<CellAnalyticsProperty>().GetProperty());
    
        std::set<unsigned> neighbour_indices = pCellPopulation -> GetNeighbouringLocationIndices(pCell);

        double neighbourCellDiversity = 0;
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

            
                    getisOrdIndex += weightingFactor*neighbourCellDiversity;

                    sumWeightingFactor +=weightingFactor;
                }
            }
        }
    }

    std::cout<<"CellNeighbourGetisOrdWriter<ELEMENT_DIM, SPACE_DIM>::CalculateNeighbourGetisOrdIndex - start"<<std::endl;
    return getisOrdIndex/sumWeightingFactor;
}

// Explicit instantiation
template class CellNeighbourGetisOrdWriter<1,1>;
template class CellNeighbourGetisOrdWriter<1,2>;
template class CellNeighbourGetisOrdWriter<2,2>;
template class CellNeighbourGetisOrdWriter<1,3>;
template class CellNeighbourGetisOrdWriter<2,3>;
template class CellNeighbourGetisOrdWriter<3,3>;