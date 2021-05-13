
#ifndef _SIMPLEHEATEQUATION_HPP_
#define _SIMPLEHEATEQUATION_HPP_

#include "AbstractLinearParabolicPde.hpp"
#include "ChastePoint.hpp"
#include "Element.hpp"

/**
 * A simple parabolic PDE used in tests.
 */
template <int SPACE_DIM>
class SimpleHeatEquation : public AbstractLinearParabolicPde<SPACE_DIM>
{

public:
    double ComputeSourceTerm(const ChastePoint<SPACE_DIM>& , double, Element<SPACE_DIM,SPACE_DIM>* pElement=NULL)
    {
        return 0.0;
    }

    c_matrix<double, SPACE_DIM, SPACE_DIM> ComputeDiffusionTerm(const ChastePoint<SPACE_DIM>& , Element<SPACE_DIM,SPACE_DIM>* pElement=NULL)
    {
        return identity_matrix<double>(SPACE_DIM);
    }

    double ComputeDuDtCoefficientFunction(const ChastePoint<SPACE_DIM>& )
    {
        return 1;
    }
};

#endif 