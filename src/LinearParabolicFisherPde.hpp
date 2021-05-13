#ifndef LINEARPARABOLICFISHERPDE_HPP_
#define LINEARPARABOLICFISHERPDE_HPP_

#include "AbstractLinearParabolicPdeSystemForCoupledOdeSystem.hpp"
#include "InhomogenousParabolicPdeForCoupledOdeSystem_templated.hpp"
#include "ChastePoint.hpp"
#include "Element.hpp"
#include "AbstractChemicalOdeForCoupledPdeSystem.hpp"

/**
* The Fisher - KPP reaction diffusion equation
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
class LinearParabolicFisherPde : public AbstractLinearParabolicPdeSystemForCoupledOdeSystem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>
{
private:

    // diffusion rate
    std::vector<double> mDiffusionRates;      

    // population growth rate
    std::vector<double> mGrowthRates;

    // carrying capacity
    std::vector<double> mCarryingCapacities;


public:

    LinearParabolicFisherPde(
        std::vector<double> diffusionRates,
        std::vector<double> growthRates,
        std::vector<double> carryingCapacities
                        )
        : AbstractLinearParabolicPdeSystemForCoupledOdeSystem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>(),
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

        diffusion_term = mDiffusionRates[pdeIndex]*identity_matrix<double>(SPACE_DIM);
        
        return diffusion_term;
    }
};

#endif
