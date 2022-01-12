#include "CellNeighbourDiversityWriter.hpp"
#include "AbstractCellPopulation.hpp"
#include "CellAnalyticsProperty.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
CellNeighbourDiversityWriter<ELEMENT_DIM, SPACE_DIM>::CellNeighbourDiversityWriter()
    : AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>("cellNeighbourDiversity.dat")
{
    this->mVtkCellDataName = "DiversityIndex";
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double CellNeighbourDiversityWriter<ELEMENT_DIM, SPACE_DIM>::GetCellDataForVtkOutput(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    return CalculateNeighbourDiversityIndex(pCell, pCellPopulation);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellNeighbourDiversityWriter<ELEMENT_DIM, SPACE_DIM>::VisitCell(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    //std::cout<<"CellNeighbourDiversityWriter<ELEMENT_DIM, SPACE_DIM>::VisitCell - start"<<std::endl;
    // Write location index corresponding to cell
    *this->mpOutStream << pCellPopulation->GetLocationIndexUsingCell(pCell) << " ";

    // Write cell location
    c_vector<double, SPACE_DIM> cell_location = pCellPopulation->GetLocationOfCellCentre(pCell);
    for (unsigned i=0; i<SPACE_DIM; i++)
    {
        *this->mpOutStream << cell_location[i] << " ";
    }

    // Write cell age
    *this->mpOutStream << CalculateNeighbourDiversityIndex(pCell, pCellPopulation) << " ";
    //std::cout<<"CellNeighbourDiversityWriter<ELEMENT_DIM, SPACE_DIM>::VisitCell - end"<<std::endl;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double CellNeighbourDiversityWriter<ELEMENT_DIM, SPACE_DIM>::CalculateNeighbourDiversityIndex(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    //std::cout<<"CellNeighbourDiversityWriter<ELEMENT_DIM, SPACE_DIM>::CalculateNeighbourDiversityIndex - start"<<std::endl;
    double diversityIndex = 0.0;
    unsigned numberOfCells=0;

    if (pCell-> rGetCellPropertyCollection(). template HasProperty<CellAnalyticsProperty>())
    {
        boost::shared_ptr<CellAnalyticsProperty> analytics_cell_property = boost::static_pointer_cast<CellAnalyticsProperty>(pCell-> rGetCellPropertyCollection(). GetPropertiesType<CellAnalyticsProperty>().GetProperty());
    
        std::vector<unsigned> cellTypeCounts = analytics_cell_property-> GetNeighbourTypes();
    

        unsigned sumNumerator=0;

        for(unsigned i=0; i<cellTypeCounts.size(); i++)
        {
            sumNumerator += cellTypeCounts[i]*(cellTypeCounts[i]-1); 
            numberOfCells += cellTypeCounts[i];
        }
    //    std::cout<<"sumNumerator: "<<sumNumerator<<std::endl;
    //    std::cout<<"numberOfCells: "<<numberOfCells<<std::endl;
        if(numberOfCells>0)
        {
            diversityIndex = 1 - ( (double) sumNumerator /((double)numberOfCells*(double) (numberOfCells-1)));
        }
        else
        {
            std::cout<<"Error: CellNeighbourDiversityWriter<ELEMENT_DIM, SPACE_DIM>::CalculateNeighbourDiversityIndex - number of neighbours ==0"<<std::endl;
        }
    }

    //std::cout<<"CellNeighbourDiversityWriter<ELEMENT_DIM, SPACE_DIM>::CalculateNeighbourDiversityIndex - start"<<std::endl;
    return diversityIndex;
}

// Explicit instantiation
template class CellNeighbourDiversityWriter<1,1>;
template class CellNeighbourDiversityWriter<1,2>;
template class CellNeighbourDiversityWriter<2,2>;
template class CellNeighbourDiversityWriter<1,3>;
template class CellNeighbourDiversityWriter<2,3>;
template class CellNeighbourDiversityWriter<3,3>;