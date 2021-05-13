#ifndef _ABSTRACTINHOMOGENOUSODESYSTEMFORCOUPLEDPDESYSTEM_HPP_
#define _ABSTRACTINHOMOGENOUSODESYSTEMFORCOUPLEDPDESYSTEM_HPP_

#include "AbstractOdeSystemForCoupledPdeSystem.hpp"
#include "OdeSystemInformation.hpp"
#include "StateVariableRegister.hpp"
#include <string>


class AbstractInhomogenousOdeSystemForCoupledPdeSystem : public AbstractOdeSystemForCoupledPdeSystem
{
protected:

    // state variable register for all the states under the influence of this ode system
    StateVariableRegister* mpStateVariableRegister;

    // state variable register for all the states in the full domain system
    StateVariableRegister* mpDomainStateVariableRegister;

public:


    AbstractInhomogenousOdeSystemForCoupledPdeSystem(   unsigned numberOfStateVariables = 0,
                                                        unsigned pdeSolutionSize = 0
                                                    )
        : AbstractOdeSystemForCoupledPdeSystem(numberOfStateVariables,pdeSolutionSize)
          
    {
        mpStateVariableRegister = new StateVariableRegister(std::vector<std::string>());
        mpDomainStateVariableRegister = new StateVariableRegister(std::vector<std::string>());
    }


    void SetStateVariableRegister(StateVariableRegister* p_stateVariableRegister)
    {
        mpStateVariableRegister = p_stateVariableRegister;
    }

    StateVariableRegister* GetStateVariableRegister()
    {
        return mpStateVariableRegister;
    }

    void SetStateVariableNameVector(std::vector<std::string> nameVector)
    {
        delete mpStateVariableRegister;
        mpStateVariableRegister = new StateVariableRegister(nameVector);
    }

    void SetDomainStateVariableRegister(StateVariableRegister* p_domainStateVaribleRegister)
    {
        mpDomainStateVariableRegister = p_domainStateVaribleRegister;
        mPdeSolutionSize = p_domainStateVaribleRegister-> GetNumberOfStateVariables();
    }

    StateVariableRegister* GetDomainStateVariableRegister()
    {   
        return mpDomainStateVariableRegister;
    }

};

#endif 