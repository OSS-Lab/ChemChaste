#ifndef _ABSTRACTINHOMOGENOUSODESYSTEM_HPP_
#define _ABSTRACTINHOMOGENOUSODESYSTEM_HPP_

#include "AbstractOdeSystemForCoupledPdeSystem.hpp"
#include "OdeSystemInformation.hpp"
#include "Exception.hpp"


class AbstractInhomogenousOdeSystem : public AbstractOdeSystemForCoupledPdeSystem
{
protected:

    /**
     * Current solution to the PDE problem.
     */
    std::vector<double> mPdeSolution;

    /**
     * The size of the PDE solution at a point in space.
     */
    unsigned mPdeSolutionSize;

    // populate from node chemistry?
    StateVariableRegister* mpStateVariableRegister; 


public:
    /**
     * Constructor.
     *
     * @param numberOfStateVariables  the number of state variables in the ODE system (defaults to 0)
     * @param pdeSolutionSize  the size of the PDE solution at a point in space (defaults to 0)
     */
    AbstractInhomogenousOdeSystem(unsigned numberOfStateVariables = 0,
                                         unsigned pdeSolutionSize = 0)
        : AbstractOdeSystem(numberOfStateVariables),
          mPdeSolutionSize(pdeSolutionSize)
    {
        mPdeSolution.clear();
    }

    /**
     * @return #mPdeSolution.
     */
    std::vector<double>& rGetPdeSolution()
    {
        return mPdeSolution;
    }

    /**
     * Set #mPdeSolution.
     *
     * @param pdeSolution the PDE solution at a point in space
     */
    void SetPdeSolution(std::vector<double> pdeSolution)
    {

        if (pdeSolution.size() != mPdeSolutionSize)
        {
            EXCEPTION("The supplied vector is not the correct size.");
        }
        
        mPdeSolution = pdeSolution;
    }

    /**
     * @return #mPdeSolutionSize.
     */
    unsigned GetPdeSolutionSize()
    {
        return mPdeSolutionSize;
    }

    double rRetrieveStateVariable(std::string variable_name)
    {
        if(IsStateVariablePresent(variable_name))
        {
            return rGetStateVariable()[mpStateVariableRegister -> RetrieveStateVariableIndex(variable_name)];
        }
        else
        {
            return 0.0;
        }
    }

    unsigned GetNumberOfStateVariables()
    {
        return mpStateVariableRegister -> GetNumberOfStateVariables();
    }

};

#endif 
