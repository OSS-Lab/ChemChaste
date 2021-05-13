#ifndef PDESCHNACKENBERGCOUPLEDPDEODESYSTEM_HPP_
#define PDESCHNACKENBERGCOUPLEDPDEODESYSTEM_HPP_

#include "AbstractLinearParabolicPdeSystemForCoupledOdeSystem.hpp"
#include "ChastePoint.hpp"
#include "Element.hpp"
#include "AbstractChemicalOdeForCoupledPdeSystem.hpp"

/**
 * Two coupled PDEs defining the Schnackenberg reaction-diffusion system
 *
 * u_t = D1*del^2 u + kappa1 - kappa_1*u + kappa3*u^2*v,
 * v_t = D2*del^2 u + kappa2 - kappa3*u^2*v.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
class PdeSchnackenbergCoupledPdeOdeSystem : public AbstractLinearParabolicPdeSystemForCoupledOdeSystem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>
{
private:
    AbstractOdeSystemForCoupledPdeSystem* mp_odeSystem;

    unsigned mCount=0;
    double mD1;      /**< Parameter D1 in the Schnackenberg system. */
    double mD2;      /**< Parameter D2 in the Schnackenberg system. */


public:

    PdeSchnackenbergCoupledPdeOdeSystem(
        AbstractOdeSystemForCoupledPdeSystem* p_odeSystem,
        double d1=1.0,
                                  double d2=1.0
                                  )
        : AbstractLinearParabolicPdeSystemForCoupledOdeSystem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>(),
          mD1(d1),
          mD2(d2),
            mp_odeSystem(p_odeSystem)
    {
    }

    double ComputeDuDtCoefficientFunction(const ChastePoint<SPACE_DIM>& rX, unsigned index)
    {
        return 1.0;
    }

    double ComputeSourceTerm(const ChastePoint<SPACE_DIM>& rX, c_vector<double,PROBLEM_DIM>& rU, std::vector<double>& rOdeSolution, unsigned pdeIndex)
    {
        std::vector<double> rUvec(PROBLEM_DIM,0.0);
        for(unsigned i=0;i<PROBLEM_DIM;i++)
        {
            rUvec[i] = rU(i);
        }
        
        mp_odeSystem->EvaluateYDerivatives(1.0, rUvec, rOdeSolution);

        return rUvec[pdeIndex];
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
