#include "EnvironmentCellProperty.hpp"

EnvironmentCellProperty::EnvironmentCellProperty()
    : AbstractCellProperty()
{
}


EnvironmentCellProperty::~EnvironmentCellProperty()
{
}

EnvironmentCellProperty::EnvironmentCellProperty(const EnvironmentCellProperty& existingProperty)
{
    mpEnvironmentStateVariableRegister = existingProperty.mpEnvironmentStateVariableRegister;
    mEnvironmentVector = existingProperty.mEnvironmentVector;
}

// virtual methods

void EnvironmentCellProperty::SetUp(CellPtr this_cellPtr, std::vector<std::string> environmentParameterNames)
{
    SetCellPtr(this_cellPtr);

    StateVariableRegister* p_stateRegister = new StateVariableRegister(environmentParameterNames);

    SetEnvironmentStateVariableRegister(p_stateRegister);

    SetUpEnvironmentVector(p_stateRegister -> GetNumberOfStateVariables());

}

void EnvironmentCellProperty::UpdateEnvironmentVector(std::vector<double> environmentVector)
{
    mEnvironmentVector = environmentVector;
}

void EnvironmentCellProperty::PreparePostDivisionParent(double splitRatio)
{
    // split any properties that are shared
}
    
void  EnvironmentCellProperty::PreparePostDivisionDaughter(const EnvironmentCellProperty& parentProperty, double splitRatio)
{
    // split any properties that are shared
}

void EnvironmentCellProperty::SetUpEnvironmentVector(unsigned numberOfStateVariables)
{
    std::vector<double> reset(numberOfStateVariables,0.0);
    mEnvironmentVector = reset;
}

// set methods

void EnvironmentCellProperty::SetEnvironmentStateVariableRegister(StateVariableRegister* p_stateRegister)
{
    mpEnvironmentStateVariableRegister = p_stateRegister;
}

void EnvironmentCellProperty::SetEnvironmentVector(std::vector<double> valueVector)
{
    mEnvironmentVector = valueVector;
}

void EnvironmentCellProperty::SetCellPtr(CellPtr this_cellPtr)
{
    mThis_cellPtr = this_cellPtr;
}

// get methods

StateVariableRegister* EnvironmentCellProperty::GetEnvironmentStateVariableRegister()
{
    return mpEnvironmentStateVariableRegister;
}


std::vector<double> EnvironmentCellProperty::GetEnvironmentVector()
{
    return mEnvironmentVector;
}

double EnvironmentCellProperty::GetEnvironmentValueByIndex(unsigned index)
{
    if(index<mpEnvironmentStateVariableRegister -> GetNumberOfStateVariables())
    {
        return mEnvironmentVector[index];
    }

    return 0.0; 
}

double EnvironmentCellProperty::GetEnvironmentValueByName(std::string stateName)
{
    unsigned index = mpEnvironmentStateVariableRegister->RetrieveStateVariableIndex(stateName);

    if(index<mpEnvironmentStateVariableRegister -> GetNumberOfStateVariables())
    {
        return mEnvironmentVector[index];
    }

    return 0.0; 
}



CellPtr EnvironmentCellProperty::GetCellPtr()
{
    return mThis_cellPtr;
}