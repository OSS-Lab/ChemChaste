#ifndef PDECONSUMERPRODUCER_HPP_
#define PDECONSUMERPRODUCER_HPP_

#include "AbstractLinearParabolicPdeSystemForCoupledOdeSystem.hpp"
#include "ChastePoint.hpp"
#include "Element.hpp"



template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
class PdeConsumerProducer : public AbstractLinearParabolicPdeSystemForCoupledOdeSystem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>
{
private:

    double mD1;     
    double mD2;      


public:

    PdeConsumerProducer(double d1=1.0,
                        double d2=1.0)
        : AbstractLinearParabolicPdeSystemForCoupledOdeSystem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>(),
          mD1(d1),
          mD2(d2)
    {
    }

    double ComputeDuDtCoefficientFunction(const ChastePoint<SPACE_DIM>& rX, unsigned index)
    {
        return 1.0;
    }

    double ComputeSourceTerm(const ChastePoint<SPACE_DIM>& rX, c_vector<double,PROBLEM_DIM>& rU, std::vector<double>& rOdeSolution, unsigned pdeIndex)
    {
        return rOdeSolution[pdeIndex];
    }

    c_matrix<double, SPACE_DIM, SPACE_DIM> ComputeDiffusionTerm(const ChastePoint<SPACE_DIM>& rX, unsigned pdeIndex, Element<ELEMENT_DIM,SPACE_DIM>* pElement=NULL)
    {
        assert(pdeIndex == 0 || pdeIndex == 1);

        c_matrix<double, SPACE_DIM, SPACE_DIM> diffusion_term;
        if (pdeIndex == 0)
        {
            diffusion_term = mD1*identity_matrix<double>(SPACE_DIM);
        }
        else // pdeIndex == 1
        {
            diffusion_term = mD2*identity_matrix<double>(SPACE_DIM);
        }
        return diffusion_term;
    }
};

#endif