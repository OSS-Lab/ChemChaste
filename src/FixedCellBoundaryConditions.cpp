#include "FixedCellBoundaryConditions.hpp"
#include "AbstractCentreBasedCellPopulation.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
FixedCellBoundaryCondition<ELEMENT_DIM,SPACE_DIM>::FixedCellBoundaryCondition(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>* pCellPopulation,
                                                    std::vector<double> boundaryMax,
                                                    std::vector<double> boundaryMin)
        : AbstractCellPopulationBoundaryCondition<ELEMENT_DIM,SPACE_DIM>(pCellPopulation),
          mBoundaryMax(boundaryMax),
          mBoundaryMin(boundaryMin)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void FixedCellBoundaryCondition<ELEMENT_DIM,SPACE_DIM>::ImposeBoundaryCondition(const std::map<Node<SPACE_DIM>*, c_vector<double, SPACE_DIM> >& rOldLocations)
{

    for (typename AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>::Iterator cell_iter = this->mpCellPopulation->Begin();
            cell_iter != this->mpCellPopulation->End();
            ++cell_iter)
    {
        unsigned node_index = this->mpCellPopulation->GetLocationIndexUsingCell(*cell_iter);
        Node<SPACE_DIM>* p_node = this->mpCellPopulation->GetNode(node_index);

        c_vector<double, SPACE_DIM> node_location = p_node->rGetLocation();

        for(unsigned dim=0; dim<SPACE_DIM; dim++)
        {
            if(node_location[dim]>mBoundaryMax[dim])
            {
                p_node->rGetModifiableLocation()[dim] = mBoundaryMax[dim];
            }
            else if(node_location[dim]<mBoundaryMin[dim])
            {
                p_node->rGetModifiableLocation()[dim] = mBoundaryMin[dim];
            }

        }
    }    
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool FixedCellBoundaryCondition<ELEMENT_DIM,SPACE_DIM>::VerifyBoundaryCondition()
{
    //std::cout<<"FixedCellBoundaryCondition<ELEMENT_DIM,SPACE_DIM>::VerifyBoundaryCondition() - start"<<std::endl;
    bool condition_satisfied = true;
    //std::cout<<"mBoundaryMax: x= "<<mBoundaryMax[0]<<" y= "<<mBoundaryMax[1]<<std::endl;
    //std::cout<<"mBoundaryMin: x= "<<mBoundaryMin[0]<<" y= "<<mBoundaryMin[1]<<std::endl;
    if (SPACE_DIM == 1)
    {
        EXCEPTION("FixedCellBoundaryCondition is not implemented in 1D");
    }
    else
    {
        unsigned count=0;
        for (typename AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::Iterator cell_iter = this->mpCellPopulation->Begin();
             cell_iter != this->mpCellPopulation->End();
             ++cell_iter)
        {
            c_vector<double, SPACE_DIM> cell_location = this->mpCellPopulation->GetLocationOfCellCentre(*cell_iter);
            //std::cout<<"Cell number: "<<count<<" location: x= "<<cell_location[0]<<" y= "<<cell_location[1]<<std::endl;
            count +=1;
            for(unsigned dim=0; dim<SPACE_DIM; dim++)
            {
                if(cell_location[dim]>mBoundaryMax[dim])
                {
                    condition_satisfied = false;
                    break;
                }
                else if(cell_location[dim]<mBoundaryMin[dim])
                {
                    condition_satisfied = false;
                    break;
                }   
            }
        }
    }
    //std::cout<<"condition_satisfied = "<<condition_satisfied<<std::endl;
    //std::cout<<"FixedCellBoundaryCondition<ELEMENT_DIM,SPACE_DIM>::VerifyBoundaryCondition() - end"<<std::endl;
    return condition_satisfied;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void FixedCellBoundaryCondition<ELEMENT_DIM,SPACE_DIM>::OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile)
{
   
    // Call method on direct parent class
    AbstractCellPopulationBoundaryCondition<ELEMENT_DIM,SPACE_DIM>::OutputCellPopulationBoundaryConditionParameters(rParamsFile);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void FixedCellBoundaryCondition<ELEMENT_DIM,SPACE_DIM>::SetBoundaryMax(std::vector<double> boundaryMax)
{
    mBoundaryMax = boundaryMax;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void FixedCellBoundaryCondition<ELEMENT_DIM,SPACE_DIM>::SetBoundaryMin(std::vector<double> boundaryMin)
{
    mBoundaryMin = boundaryMin;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<double> FixedCellBoundaryCondition<ELEMENT_DIM,SPACE_DIM>::GetBoundaryMax()
{
    return mBoundaryMax;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<double> FixedCellBoundaryCondition<ELEMENT_DIM,SPACE_DIM>::GetBoundaryMin()
{
    return mBoundaryMin;
}

// Explicit instantiation
template class FixedCellBoundaryCondition<1,1>;
template class FixedCellBoundaryCondition<1,2>;
template class FixedCellBoundaryCondition<2,2>;
template class FixedCellBoundaryCondition<1,3>;
template class FixedCellBoundaryCondition<2,3>;
template class FixedCellBoundaryCondition<3,3>;