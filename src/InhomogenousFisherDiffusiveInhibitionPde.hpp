#ifndef INHOMOGENOUSFISHERDIFFUSIVEINHIBITIONPDE_HPP_
#define INHOMOGENOUSFISHERDIFFUSIVEINHIBITIONPDE_HPP_

#include "AbstractLinearParabolicPdeSystemForCoupledOdeSystem.hpp"
#include "InhomogenousParabolicPdeForCoupledOdeSystem_templated.hpp"
#include "ChastePoint.hpp"
#include "Element.hpp"
#include "AbstractChemicalOdeForCoupledPdeSystem.hpp"

/**
* The Fisher - KPP reaction diffusion equation with areas with inhibited diffusion
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
class InhomogenousFisherDiffusiveInhibitionPde : public InhomogenousParabolicPdeForCoupledOdeSystemTemplated<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>
{
private:

    // diffusion rate
    std::vector<double> mDiffusionRates;      

    // population growth rate
    std::vector<double> mGrowthRates;

    // carrying capacity
    std::vector<double> mCarryingCapacities;


public:

    InhomogenousFisherDiffusiveInhibitionPde(
        std::vector<double> diffusionRates,
        std::vector<double> growthRates,
        std::vector<double> carryingCapacities
                        )
        : InhomogenousParabolicPdeForCoupledOdeSystemTemplated<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>(diffusionRates),
          mDiffusionRates(diffusionRates),
          mGrowthRates(growthRates),
          mCarryingCapacities(carryingCapacities)
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

        source_term_this_pde_index = mGrowthRates[pdeIndex]*rU(pdeIndex)*(1-rU(pdeIndex)/mCarryingCapacities[pdeIndex]);

        return source_term_this_pde_index;

    }

    c_matrix<double, SPACE_DIM, SPACE_DIM> ComputeDiffusionTerm(const ChastePoint<SPACE_DIM>& rX, unsigned pdeIndex, Element<ELEMENT_DIM,SPACE_DIM>* pElement=NULL)
    {
        assert(pdeIndex<PROBLEM_DIM);

        c_matrix<double, SPACE_DIM, SPACE_DIM> diffusion_term;

        if(std::pow((rX[0]-50),2)+std::pow((rX[1]-45),2) <= 49 || std::pow((rX[0]-50),2)+std::pow((rX[1]-55),2) <= 49 )
        {
            std::cout<<"rX: "<<rX[0]<<" rY: "<<rX[1]<<std::endl;
            diffusion_term = (0.1*mDiffusionRates[pdeIndex])*identity_matrix<double>(SPACE_DIM);
        }
        else
        {
            diffusion_term = mDiffusionRates[pdeIndex]*identity_matrix<double>(SPACE_DIM);
        }
        
   
        return diffusion_term;
    }
};

#endif
