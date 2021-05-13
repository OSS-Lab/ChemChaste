#ifndef INHOMOGENOUSSCHNACKENBERGNONDIMPDE_HPP_
#define INHOMOGENOUSSCHNACKENBERGNONDIMPDE_HPP_

#include "AbstractLinearParabolicPdeSystemForCoupledOdeSystem.hpp"
#include "InhomogenousParabolicPdeForCoupledOdeSystem_templated.hpp"
#include "ChastePoint.hpp"
#include "Element.hpp"
#include "AbstractChemicalOdeForCoupledPdeSystem.hpp"

/**
* The non-dimensionalised Schnackenberg reaction diffusion equation
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
class InhomogenousSchnackenbergNonDimPde : public InhomogenousParabolicPdeForCoupledOdeSystemTemplated<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>
{
private:

    // diffusion rate
    std::vector<double> mDiffusionRates;      

    double mK1;

    double mK_1;
    
    double mK2;

    double mK3;


public:

    InhomogenousSchnackenbergNonDimPde(
        std::vector<double> diffusionRates,
        double k1,
        double k_1,
        double k2,
        double k3
                        )
        : InhomogenousParabolicPdeForCoupledOdeSystemTemplated<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>(diffusionRates),
          mDiffusionRates(diffusionRates),
          mK1(k1),
          mK_1(k_1),
          mK2(k2),
          mK3(k3)
    {
    }

    double ComputeDuDtCoefficientFunction(const ChastePoint<SPACE_DIM>& rX, unsigned index)
    {
        return 1.0;
    }

    double ComputeSourceTerm(const ChastePoint<SPACE_DIM>& rX, c_vector<double,PROBLEM_DIM>& rU, std::vector<double>& rOdeSolution, unsigned pdeIndex)
    {
        assert(pdeIndex<PROBLEM_DIM);

        double source_term_this_pde_index;
        if(pdeIndex==0)
        {
            // U
            source_term_this_pde_index = mK1 - mK_1*rU(0) + mK3*rU(0)*rU(0)*rU(1);
        }
        else
        {
            // V
            source_term_this_pde_index = mK2 - mK3*rU(0)*rU(0)*rU(1);
        }

        return source_term_this_pde_index;

    }

    c_matrix<double, SPACE_DIM, SPACE_DIM> ComputeDiffusionTerm(const ChastePoint<SPACE_DIM>& rX, unsigned pdeIndex, Element<ELEMENT_DIM,SPACE_DIM>* pElement=NULL)
    {
        assert(pdeIndex<PROBLEM_DIM);

        c_matrix<double, SPACE_DIM, SPACE_DIM> diffusion_term;

        diffusion_term = mDiffusionRates[pdeIndex]*identity_matrix<double>(SPACE_DIM);
        
        return diffusion_term;
    }
};

#endif
