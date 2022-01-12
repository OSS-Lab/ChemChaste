#include "PopulationDiversityWriter.hpp"

#include "AbstractCellPopulation.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "CaBasedCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "CellAnalyticsProperty.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
PopulationDiversityWriter<ELEMENT_DIM, SPACE_DIM>::PopulationDiversityWriter()
    : AbstractCellPopulationCountWriter<ELEMENT_DIM, SPACE_DIM>("populationDiversity.dat")
{
}

 
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double PopulationDiversityWriter<ELEMENT_DIM, SPACE_DIM>::CalculateDiversityIndex(unsigned numberOfCells)
{
    //std::cout<<"PopulationDiversityWriter<ELEMENT_DIM, SPACE_DIM>::CalculateDiversityIndex - start"<<std::endl;
    unsigned sumNumerator=0;
    double diversityIndex=0.0;

    for(unsigned i=0; i<mCellTypesString.size(); i++)
    {
        sumNumerator += mNumberOfCellsOfType[i]*(mNumberOfCellsOfType[i]-1); 

    }

    diversityIndex = 1 - ( (double) sumNumerator /((double)numberOfCells*(double) (numberOfCells-1)));
    //std::cout<<"sumNumerator: "<<sumNumerator<<std::endl;
    //std::cout<<"diversityIndex: "<<diversityIndex<<std::endl;
    //std::cout<<"PopulationDiversityWriter<ELEMENT_DIM, SPACE_DIM>::CalculateDiversityIndex - end"<<std::endl;
    return diversityIndex;
    
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void PopulationDiversityWriter<ELEMENT_DIM, SPACE_DIM>::WriteHeader(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{

    //initialise from the first cell in the population to have an analytic property
    for (typename AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>::Iterator cell_iter = pCellPopulation -> Begin();
         cell_iter != pCellPopulation -> End();
         ++cell_iter)
    {
        if (cell_iter-> rGetCellPropertyCollection(). template HasProperty<CellAnalyticsProperty>())
        {
            boost::shared_ptr<CellAnalyticsProperty> analytics_cell_property = boost::static_pointer_cast<CellAnalyticsProperty>(cell_iter-> template rGetCellPropertyCollection(). template GetPropertiesType<CellAnalyticsProperty>().GetProperty());

            SetCellTypesString(analytics_cell_property -> GetPopulationCellTypeNames());
            std::vector<unsigned> numberOfCellsOfType(analytics_cell_property -> GetPopulationCellTypeNames().size(),0);

            SetNumberOfCellsOfType(numberOfCellsOfType);
            
            break;
        }    
    }

    if (PetscTools::AmMaster())
    {
        //pCellPopulation->SetDefaultCellMutationStateAndProliferativeTypeOrdering();

        *this->mpOutStream << "Time\t ";

        for (unsigned i=0; i<mCellTypesString.size(); i++)
        {
            *this->mpOutStream << mCellTypesString[i] << "\t ";
        }
        *this->mpOutStream << "Diversity_Index" << "\t ";
        this->WriteNewline();
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void PopulationDiversityWriter<ELEMENT_DIM, SPACE_DIM>::Visit(NodeBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    //std::cout<<"PopulationDiversityWriter<ELEMENT_DIM, SPACE_DIM>::Visit - start"<<std::endl;
  
    std::string cellTypeName = "";
    
    unsigned numberOfCells =0;

    double diversityIndex=1.0;

    std::fill(mNumberOfCellsOfType.begin(), mNumberOfCellsOfType.end(), 0);

    //std::cout<<"here0"<<std::endl;
    if (PetscTools::AmMaster())
    {
        for (typename AbstractCellPopulation<SPACE_DIM>::Iterator cell_iter = pCellPopulation -> Begin();
         cell_iter != pCellPopulation -> End();
         ++cell_iter)
        {
            if (cell_iter-> rGetCellPropertyCollection(). template HasProperty<CellAnalyticsProperty>())
            {
                boost::shared_ptr<CellAnalyticsProperty> analytics_cell_property = boost::static_pointer_cast<CellAnalyticsProperty>(cell_iter-> template rGetCellPropertyCollection(). template GetPropertiesType<CellAnalyticsProperty>().GetProperty());

                cellTypeName = analytics_cell_property-> GetCellTypeName();
            
                for(unsigned i=0; i<mCellTypesString.size(); i++)
                {
                    //std::cout<<"CellName: "<<cellTypeName<<" test: "<<mCellTypesString[i]<<std::endl;
                    if(cellTypeName == mCellTypesString[i])
                    {
                        mNumberOfCellsOfType[i] += 1;
                        break;
                    }

                }
                numberOfCells +=1;
            }
        }
        diversityIndex = CalculateDiversityIndex(numberOfCells);
        //std::cout<<"DiversityIndex: "<<diversityIndex<<std::endl;
        for (unsigned i=0; i<mCellTypesString.size(); i++)
        {
            *this->mpOutStream << mNumberOfCellsOfType[i] << "\t ";
        }
        *this->mpOutStream << diversityIndex << "\t ";
        
    }
    //std::cout<<"PopulationDiversityWriter<ELEMENT_DIM, SPACE_DIM>::Visit - end"<<std::endl;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void PopulationDiversityWriter<ELEMENT_DIM, SPACE_DIM>::Visit(MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void PopulationDiversityWriter<ELEMENT_DIM, SPACE_DIM>::Visit(CaBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void PopulationDiversityWriter<ELEMENT_DIM, SPACE_DIM>::Visit(PottsBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void PopulationDiversityWriter<ELEMENT_DIM, SPACE_DIM>::Visit(VertexBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void PopulationDiversityWriter<ELEMENT_DIM, SPACE_DIM>::SetNumberOfCellsOfType(std::vector<unsigned> numberOfCellsOfType)
{
    mNumberOfCellsOfType = numberOfCellsOfType;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void PopulationDiversityWriter<ELEMENT_DIM, SPACE_DIM>::SetCellTypesString(std::vector<std::string> cellTypes)
{
    mCellTypesString = cellTypes;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<unsigned>  PopulationDiversityWriter<ELEMENT_DIM, SPACE_DIM>::GetNumberOfCellsOfType()
{
    return mNumberOfCellsOfType;
}



// Explicit instantiation
template class PopulationDiversityWriter<1,1>;
template class PopulationDiversityWriter<1,2>;
template class PopulationDiversityWriter<2,2>;
template class PopulationDiversityWriter<1,3>;
template class PopulationDiversityWriter<2,3>;
template class PopulationDiversityWriter<3,3>;