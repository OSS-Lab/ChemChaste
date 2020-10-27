#ifndef ODECONSUMERPRODUCER_HPP_
#define ODECONSUMERPRODUCER_HPP_

#include "AbstractOdeSystemForCoupledPdeSystem.hpp"
#include "OdeSystemInformation.hpp"


class OdeConsumerProducer: public AbstractOdeSystemForCoupledPdeSystem
{
private:

    double mKappa1;  
    double mKappa2;  
    double mKappa3;
    double mKappa_3; 

public:

    OdeConsumerProducer(double kappa1=1.0,
                        double kappa2=1.0,
                        double kappa3=1.0,
                        double kappa_3=1.0)
        : AbstractOdeSystemForCoupledPdeSystem(2,2),
          mKappa1(kappa1),
          mKappa2(kappa2),
          mKappa3(kappa3),
          mKappa_3(kappa_3)
    {
        mpSystemInfo = OdeSystemInformation<OdeConsumerProducer>::Instance();
        SetStateVariables(GetInitialConditions());
    }

    /**
     * Method to evaluate the derivatives of the system.
     *
     * @param time  the current time
     * @param rY  the current values of the state variables
     * @param rDY  storage for the derivatives of the system; will be filled in on return
     */
    void EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
    {
        //rY is the same as the pde solution at the point, so pde is updated first then the ode is followed
        rDY[0] = mKappa1 + mKappa_3*rY[1]*rY[0] - mKappa3*rY[1]*rY[0];
        rDY[1] = -mKappa2 + mKappa3*rY[1]*rY[0] - mKappa_3*rY[1]*rY[0];
    }  
};

template<>
void OdeSystemInformation<OdeConsumerProducer>::Initialise()
{
    this->mVariableNames.push_back("var_cell_0");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(1.0);
    this->mVariableNames.push_back("var_cell_1");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(1.0);
    this->mInitialised = true;
}



#endif 