#ifndef ODESCHNACKENBERGCOUPLEDPDEODESYSTEM_HPP_
#define ODESCHNACKENBERGCOUPLEDPDEODESYSTEM_HPP_

#include "AbstractOdeSystemForCoupledPdeSystem.hpp"
#include "OdeSystemInformation.hpp"


class OdeSchnackenbergCoupledPdeOdeSystem : public AbstractOdeSystemForCoupledPdeSystem
{
private:

    double mKappa1;  /**< Parameter kappa1 in the Schnackenberg system. */
    double mKappa_1; /**< Parameter kappa_1 in the Schnackenberg system. */
    double mKappa2;  /**< Parameter kappa2 in the Schnackenberg system. */
    double mKappa3;  /**< Parameter kappa3 in the Schnackenberg system. */

public:

    /**
     * Constructor.
     *
     * @param a the value of the parameter mA
     */
    OdeSchnackenbergCoupledPdeOdeSystem(double kappa1=1.0,
                                        double kappa_1=1.0,
                                        double kappa2=1.0,
                                        double kappa3=1.0)
        : AbstractOdeSystemForCoupledPdeSystem(2,2),
          mKappa1(kappa1),
          mKappa_1(kappa_1),
          mKappa2(kappa2),
          mKappa3(kappa3)
    {
        mpSystemInfo = OdeSystemInformation<OdeSchnackenbergCoupledPdeOdeSystem>::Instance();
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
        //std::vector<double> pde_sol = this-> rGetPdeSolution();
        //rY is the same as the pde solution at the point, so pde is updated first then the ode is followed
        rDY[0] = mKappa1 - mKappa_1*rY[0] + mKappa3*rY[1]*rY[0]*rY[0];
        rDY[1] = mKappa2 - mKappa3*rY[1]*rY[0]*rY[0];
    }

    /**
     * @return mkappa1
     */
    double GetKappa1()
    {
        return mKappa1;
    }

    /**
     * @return mkappa_1
     */
    double GetKappa_1()
    {
        return mKappa_1;
    }

    /**
     * @return mkappa2
     */
    double GetKappa2()
    {
        return mKappa2;
    }

    /**
     * @return mkappa3
     */
    double GetKapp3()
    {
        return mKappa3;
    }
    
};

template<>
void OdeSystemInformation<OdeSchnackenbergCoupledPdeOdeSystem>::Initialise()
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
