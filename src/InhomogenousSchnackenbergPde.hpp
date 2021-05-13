#ifndef INHOMOGENOUSSCHNACKENBERGPDE_HPP_
#define INHOMOGENOUSSCHNACKENBERGPDE_HPP_

#include "AbstractLinearParabolicPdeSystemForCoupledOdeSystem.hpp"
#include "InhomogenousParabolicPdeForCoupledOdeSystem_templated.hpp"
#include "ChastePoint.hpp"
#include "Element.hpp"
#include "AbstractChemicalOdeForCoupledPdeSystem.hpp"

/**
* The non-dimensionalised Schnackenberg reaction diffusion equation
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
class InhomogenousSchnackenbergPde : public InhomogenousParabolicPdeForCoupledOdeSystemTemplated<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>
{
private:

    // diffusion rate
    std::vector<double> mDiffusionRates;      

    double mA;

    double mB;
    
    double mGamma;


public:

    InhomogenousSchnackenbergPde(
        std::vector<double> diffusionRates,
        double a,
        double b,
        double gamma
                        )
        : InhomogenousParabolicPdeForCoupledOdeSystemTemplated<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>(diffusionRates),
          mDiffusionRates(diffusionRates),
          mA(a),
          mB(b),
          mGamma(gamma)
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
            source_term_this_pde_index = mGamma*(mA - rU(0) + rU(0)*rU(0)*rU(1));
        }
        else
        {
            // V
            source_term_this_pde_index = mGamma*(mB - rU(0)*rU(0)*rU(1));
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
