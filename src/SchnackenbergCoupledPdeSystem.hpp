#ifndef SCHNACKENBERGCOUPLEDPDESYSTEM_HPP_
#define SCHNACKENBERGCOUPLEDPDESYSTEM_HPP_

#include "AbstractLinearParabolicPdeSystemForCoupledOdeSystem.hpp"
#include "ChastePoint.hpp"
#include "Element.hpp"

/**
 * Two coupled PDEs defining the Schnackenberg reaction-diffusion system
 *
 * u_t = D1*del^2 u + kappa1 - kappa_1*u + kappa3*u^2*v,
 * v_t = D2*del^2 u + kappa2 - kappa3*u^2*v.
 */
template<unsigned SPACE_DIM>
class SchnackenbergCoupledPdeSystem : public AbstractLinearParabolicPdeSystemForCoupledOdeSystem<SPACE_DIM, SPACE_DIM, 2>
{
private:
    unsigned mCount=0;
    double mD1;      /**< Parameter D1 in the Schnackenberg system. */
    double mD2;      /**< Parameter D2 in the Schnackenberg system. */
    double mKappa1;  /**< Parameter kappa1 in the Schnackenberg system. */
    double mKappa_1; /**< Parameter kappa_1 in the Schnackenberg system. */
    double mKappa2;  /**< Parameter kappa2 in the Schnackenberg system. */
    double mKappa3;  /**< Parameter kappa3 in the Schnackenberg system. */

public:

    SchnackenbergCoupledPdeSystem(double d1=1.0,
                                  double d2=1.0,
                                  double kappa1=1.0,
                                  double kappa_1=1.0,
                                  double kappa2=1.0,
                                  double kappa3=1.0)
        : AbstractLinearParabolicPdeSystemForCoupledOdeSystem<SPACE_DIM, SPACE_DIM, 2>(),
          mD1(d1),
          mD2(d2),
          mKappa1(kappa1),
          mKappa_1(kappa_1),
          mKappa2(kappa2),
          mKappa3(kappa3)
    {
    }

    double ComputeDuDtCoefficientFunction(const ChastePoint<SPACE_DIM>& rX, unsigned index)
    {
        return 1.0;
    }

    double ComputeSourceTerm(const ChastePoint<SPACE_DIM>& rX, c_vector<double,2>& rU, std::vector<double>& rOdeSolution, unsigned pdeIndex)
    {
        assert(pdeIndex == 0 || pdeIndex == 1);

        double source_term;
        if (pdeIndex == 0)
        {
            
            source_term = mKappa1 - mKappa_1*rU(0) + mKappa3*rU(1)*rU(0)*rU(0);

            
        }
        else // pdeIndex == 1
        {
            source_term = mKappa2 - mKappa3*rU(1)*rU(0)*rU(0);
        }
        return source_term;
    }

    c_matrix<double, SPACE_DIM, SPACE_DIM> ComputeDiffusionTerm(const ChastePoint<SPACE_DIM>& rX, unsigned pdeIndex, Element<SPACE_DIM,SPACE_DIM>* pElement=NULL)
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

#endif /*SCHNACKENBERGCOUPLEDPDESYSTEM_HPP_*/ 