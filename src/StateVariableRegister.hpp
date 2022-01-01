#ifndef STATEVARIABLEREGISTER_HPP_
#define STATEVARIABLEREGISTER_HPP_

#include <string>
#include <vector>

/*
    Class to hold the names of the variables in an inhomogenous system, name index pair
*/
class StateVariableRegister
{
protected:

    unsigned mNumberSystemStateVariables;

    std::vector<std::string> mStateVariableRegister;

public:

    StateVariableRegister();
    
    StateVariableRegister(std::vector<std::string> VariableNameVector = std::vector<std::string>() );

    ~StateVariableRegister();

    void UpdateStateVariableRegister(std::vector<std::string>);

    std::vector<std::string> GetStateVariableRegisterVector();

    void AddStateVariableVector(std::vector<std::string>);

    void AddStateVariable(std::string);

    void RemoveStateVariable(std::string);

    bool IsStateVariablePresent(std::string);

    unsigned RetrieveStateVariableIndex(std::string);

    std::string RetrieveStateVariableName(unsigned);

    void SetNumberOfStateVariables(unsigned);

    unsigned GetNumberOfStateVariables();

    std::vector<unsigned> FindIndicesInThisRegister(StateVariableRegister* );

    std::vector<unsigned> FindIndicesInThatRegister(StateVariableRegister* );

    std::vector<std::string> FindCommonNamesInRegisters(StateVariableRegister* );

};

#endif