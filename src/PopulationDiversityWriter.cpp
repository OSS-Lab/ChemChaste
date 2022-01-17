#include "PopulationDiversityWriter.hpp"

#include "AbstractCellPopulation.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "CaBasedCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "CellAnalyticsProperty.hpp"
#include <math.h> 

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
PopulationDiversityWriter<ELEMENT_DIM, SPACE_DIM>::PopulationDiversityWriter()
    : AbstractCellPopulationCountWriter<ELEMENT_DIM, SPACE_DIM>("populationHeterogeneity.dat")
{
}

 
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double PopulationDiversityWriter<ELEMENT_DIM, SPACE_DIM>::CalculateDiversityIndex(unsigned numberOfCells)
{
    unsigned sumNumerator=0;
    double diversityIndex=0.0;

    for(unsigned i=0; i<mCellTypesString.size(); i++)
    {
        sumNumerator += mNumberOfCellsOfType[i]*(mNumberOfCellsOfType[i]-1); 

    }

    diversityIndex = 1 - ( (double) sumNumerator /((double)numberOfCells*(double) (numberOfCells-1)));
    
    return diversityIndex;
    
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double PopulationDiversityWriter<ELEMENT_DIM, SPACE_DIM>::CalculateShannonIndex()
{
    double shannonIndex=0.0; 

    for(unsigned i=0; i<mProbabilityOfType.size(); i++)
    {
        if(mProbabilityOfType[i] > 1e-5)
        {
            // in case mProbabilityOfType[i] resolves to zero where log undefined
            shannonIndex -= mProbabilityOfType[i]*log(mProbabilityOfType[i]);
        }
    }

    return shannonIndex;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double PopulationDiversityWriter<ELEMENT_DIM, SPACE_DIM>::CalculateGiniSimpsonIndex()
{
    double giniSimpsonIndex=1.0;

    for(unsigned i=0; i<mProbabilityOfType.size(); i++)
    {
        giniSimpsonIndex -= mProbabilityOfType[i]*mProbabilityOfType[i];
    }
    
    return giniSimpsonIndex;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double PopulationDiversityWriter<ELEMENT_DIM, SPACE_DIM>::CalculateLeeOyburnIndex(NodeBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    // calculate the Lee-Oyburn Index where the connected neighbour cells are given a weighting of 1 and other cells 0 i.e are not considered.
    double leeOyburnIndex=0.0;
    double cellLeeOyburnIndex=0.0;
    double thisCellTypeProbability=0.0;
    double neighbourCellTypeProbability=0.0;

    std::string thisCellTypeName="";
    std::string neighbourCellTypeName="";

    double weightingFactor =1.0;
    double indicatorFunction=0.0;
    double sumWeightingFactor=0.0;
    double thisCellDiversity=0.0;

    for (typename AbstractCellPopulation<SPACE_DIM>::Iterator cell_iter = pCellPopulation -> Begin();
         cell_iter != pCellPopulation -> End();
         ++cell_iter)
    {
        cellLeeOyburnIndex =0.0;
        thisCellTypeName="";
        thisCellTypeProbability=0.0;
        if (cell_iter-> rGetCellPropertyCollection(). template HasProperty<CellAnalyticsProperty>())
        {
            boost::shared_ptr<CellAnalyticsProperty> analytics_cell_property = boost::static_pointer_cast<CellAnalyticsProperty>(cell_iter-> template rGetCellPropertyCollection(). template GetPropertiesType<CellAnalyticsProperty>().GetProperty());

            thisCellDiversity = analytics_cell_property -> GetCellDiversity();
            thisCellTypeName = analytics_cell_property-> GetCellTypeName();

            std::set<unsigned> neighbour_indices = pCellPopulation -> GetNeighbouringLocationIndices(*cell_iter);

            if (!neighbour_indices.empty())
            {
                for (std::set<unsigned>::iterator iter = neighbour_indices.begin();
                    iter != neighbour_indices.end();
                    ++iter)
                {
                    neighbourCellTypeName="";
                    neighbourCellTypeProbability=0.0;
                    CellPtr p_cell = pCellPopulation -> GetCellUsingLocationIndex(*iter);

                    if (p_cell-> rGetCellPropertyCollection(). HasProperty<CellAnalyticsProperty>())
                    {
                        boost::shared_ptr<CellAnalyticsProperty> neighbour_analytics_cell_property = boost::static_pointer_cast<CellAnalyticsProperty>(p_cell-> rGetCellPropertyCollection(). GetPropertiesType<CellAnalyticsProperty>().GetProperty());

                        neighbourCellTypeName = neighbour_analytics_cell_property-> GetCellTypeName();


                        indicatorFunction=0.0;
                        if(thisCellTypeName == neighbourCellTypeName)
                        {
                            indicatorFunction=1.0;
                        }

                        for(unsigned i=0; i<mProbabilityOfType.size();i++)
                        {
                            if(neighbourCellTypeName == mCellTypesString[i])
                            {
                                neighbourCellTypeProbability = mProbabilityOfType[i];
                                break;
                            }
                        }

                        cellLeeOyburnIndex += weightingFactor*(2*indicatorFunction - 1.0)/neighbourCellTypeProbability;

                        sumWeightingFactor +=weightingFactor;
                    }
                }
            }

            for(unsigned i=0; i<mProbabilityOfType.size();i++)
            {
                if(thisCellTypeName == mCellTypesString[i])
                {
                    thisCellTypeProbability = mProbabilityOfType[i];
                    break;
                }
            }

            leeOyburnIndex += cellLeeOyburnIndex/thisCellTypeProbability;
        }
    }
    
    return leeOyburnIndex/sumWeightingFactor;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double PopulationDiversityWriter<ELEMENT_DIM, SPACE_DIM>::CalculateMoranIndex(NodeBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    // calculate the Moran Index where the connected neighbour cells are given a weighting of 1 and other cells 0 i.e are not considered.
    double moranIndex=0.0;
    double cellMoranIndex=0.0;
    double thisCellDiversity=0.0;
    double neighbourCellDiversity=0.0;

    double weightingFactor =1.0;

    double sumWeightingFactor=0.0;

    for (typename AbstractCellPopulation<SPACE_DIM>::Iterator cell_iter = pCellPopulation -> Begin();
         cell_iter != pCellPopulation -> End();
         ++cell_iter)
    {
        cellMoranIndex =0.0;
        if (cell_iter-> rGetCellPropertyCollection(). template HasProperty<CellAnalyticsProperty>())
        {
            boost::shared_ptr<CellAnalyticsProperty> analytics_cell_property = boost::static_pointer_cast<CellAnalyticsProperty>(cell_iter-> template rGetCellPropertyCollection(). template GetPropertiesType<CellAnalyticsProperty>().GetProperty());

            thisCellDiversity = analytics_cell_property -> GetCellDiversity();

            std::set<unsigned> neighbour_indices = pCellPopulation  -> GetNeighbouringLocationIndices(*cell_iter);

            if (!neighbour_indices.empty())
            {
                //std::cout<<"here3"<<std::endl;
                for (std::set<unsigned>::iterator iter = neighbour_indices.begin();
                    iter != neighbour_indices.end();
                    ++iter)
                {
                    CellPtr p_cell = pCellPopulation -> GetCellUsingLocationIndex(*iter);

                    if (p_cell-> rGetCellPropertyCollection(). HasProperty<CellAnalyticsProperty>())
                    {
                        boost::shared_ptr<CellAnalyticsProperty> neighbour_analytics_cell_property = boost::static_pointer_cast<CellAnalyticsProperty>(p_cell-> rGetCellPropertyCollection(). GetPropertiesType<CellAnalyticsProperty>().GetProperty());

                        neighbourCellDiversity = neighbour_analytics_cell_property -> GetCellDiversity();

              
                        cellMoranIndex += weightingFactor*( thisCellDiversity - analytics_cell_property -> GetMeanCellDiversity())*( neighbourCellDiversity - analytics_cell_property -> GetMeanCellDiversity());

                        sumWeightingFactor +=weightingFactor;



                    }
                }
            }
            moranIndex += cellMoranIndex/(( thisCellDiversity - analytics_cell_property -> GetMeanCellDiversity())*( thisCellDiversity - analytics_cell_property -> GetMeanCellDiversity()));

        }
    }
    
    return moranIndex * mNumberOfCells/sumWeightingFactor;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double PopulationDiversityWriter<ELEMENT_DIM, SPACE_DIM>::CalculateGearyIndex(NodeBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    // calculate the Geary Index where the connected neighbour cells are given a weighting of 1 and other cells 0 i.e are not considered.
    double gearyIndex=0.0;
    double cellGearyIndex=0.0;
    double thisCellDiversity=0.0;
    double neighbourCellDiversity=0.0;

    double weightingFactor =1.0;

    double sumWeightingFactor=0.0;

    for (typename AbstractCellPopulation<SPACE_DIM>::Iterator cell_iter = pCellPopulation -> Begin();
         cell_iter != pCellPopulation -> End();
         ++cell_iter)
    {
        cellGearyIndex =0.0;
        if (cell_iter-> rGetCellPropertyCollection(). template HasProperty<CellAnalyticsProperty>())
        {
            boost::shared_ptr<CellAnalyticsProperty> analytics_cell_property = boost::static_pointer_cast<CellAnalyticsProperty>(cell_iter-> template rGetCellPropertyCollection(). template GetPropertiesType<CellAnalyticsProperty>().GetProperty());

            thisCellDiversity = analytics_cell_property -> GetCellDiversity();

            std::set<unsigned> neighbour_indices = pCellPopulation -> GetNeighbouringLocationIndices(*cell_iter);

            if (!neighbour_indices.empty())
            {
                //std::cout<<"here3"<<std::endl;
                for (std::set<unsigned>::iterator iter = neighbour_indices.begin();
                    iter != neighbour_indices.end();
                    ++iter)
                {
                    CellPtr p_cell = pCellPopulation -> GetCellUsingLocationIndex(*iter);

                    if (p_cell-> rGetCellPropertyCollection(). HasProperty<CellAnalyticsProperty>())
                    {
                        boost::shared_ptr<CellAnalyticsProperty> neighbour_analytics_cell_property = boost::static_pointer_cast<CellAnalyticsProperty>(p_cell-> rGetCellPropertyCollection(). GetPropertiesType<CellAnalyticsProperty>().GetProperty());

                        neighbourCellDiversity = neighbour_analytics_cell_property -> GetCellDiversity();

              
                        cellGearyIndex += weightingFactor*( thisCellDiversity - neighbourCellDiversity)*( thisCellDiversity - neighbourCellDiversity);

                        sumWeightingFactor +=weightingFactor;
                    }
                }
            }
            gearyIndex += cellGearyIndex/(( thisCellDiversity - analytics_cell_property -> GetMeanCellDiversity())*( thisCellDiversity - analytics_cell_property -> GetMeanCellDiversity()));

        }
    }
    
    return gearyIndex * (mNumberOfCells-1)/(2*sumWeightingFactor);
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
            std::vector<double> probabilityOfType(analytics_cell_property -> GetPopulationCellTypeNames().size(),0.0);

            SetNumberOfCellsOfType(numberOfCellsOfType);
            SetProbabilityOfType(probabilityOfType);
            
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

        *this->mpOutStream << "Shannon_Index" << "\t ";

        *this->mpOutStream << "Gini-Simpson_Index" << "\t ";

        *this->mpOutStream << "Lee-Oyburn_Index" << "\t ";

        *this->mpOutStream << "Moran_Index" << "\t ";

        *this->mpOutStream << "Geary_Index" << "\t ";

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
    double Shannon_Index=1.0;
    double Gini_Simpson_Index=1.0;
    double Lee_Oyburn_Index=1.0;

    double Moran_Index = 1.0;
    double Geary_Index = 1.0;

    std::fill(mNumberOfCellsOfType.begin(), mNumberOfCellsOfType.end(), 0);

    std::fill(mProbabilityOfType.begin(), mProbabilityOfType.end(), 0.0);

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

        mNumberOfCells = numberOfCells;

        for(unsigned i=0; i<mNumberOfCellsOfType.size(); i++)
        {
            mProbabilityOfType[i] = mNumberOfCellsOfType[i]/numberOfCells;
        }

        diversityIndex = CalculateDiversityIndex(numberOfCells);
        //std::cout<<"DiversityIndex: "<<diversityIndex<<std::endl;

        Shannon_Index = CalculateShannonIndex();

        Gini_Simpson_Index = CalculateGiniSimpsonIndex();

        Lee_Oyburn_Index = CalculateLeeOyburnIndex(pCellPopulation);

        Moran_Index = CalculateMoranIndex(pCellPopulation);

        Geary_Index = CalculateGearyIndex(pCellPopulation);

        for (unsigned i=0; i<mCellTypesString.size(); i++)
        {
            *this->mpOutStream << mNumberOfCellsOfType[i] << "\t ";
        }
        *this->mpOutStream << diversityIndex << "\t ";

        *this->mpOutStream << Shannon_Index << "\t ";

        *this->mpOutStream << Gini_Simpson_Index << "\t ";

        *this->mpOutStream << Lee_Oyburn_Index << "\t ";

        *this->mpOutStream << Moran_Index << "\t ";

        *this->mpOutStream << Geary_Index << "\t ";
        
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
void PopulationDiversityWriter<ELEMENT_DIM, SPACE_DIM>::SetProbabilityOfType(std::vector<double> probabilityOfType)
{
    mProbabilityOfType = probabilityOfType;
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

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<double>  PopulationDiversityWriter<ELEMENT_DIM, SPACE_DIM>::GetProbabilityOfType()
{
    return mProbabilityOfType;
}

// Explicit instantiation
template class PopulationDiversityWriter<1,1>;
template class PopulationDiversityWriter<1,2>;
template class PopulationDiversityWriter<2,2>;
template class PopulationDiversityWriter<1,3>;
template class PopulationDiversityWriter<2,3>;
template class PopulationDiversityWriter<3,3>;