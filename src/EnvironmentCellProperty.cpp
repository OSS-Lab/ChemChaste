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
    mPreferredEnvironmentVector = existingProperty.mPreferredEnvironmentVector;
}

// virtual methods

void EnvironmentCellProperty::SetUp(CellPtr this_cellPtr, std::vector<std::string> environmentParameterNames)
{
    SetCellPtr(this_cellPtr);

    StateVariableRegister* p_stateRegister = new StateVariableRegister(environmentParameterNames);

    SetEnvironmentStateVariableRegister(p_stateRegister);

    SetUpEnvironmentVectors(p_stateRegister -> GetNumberOfStateVariables());

}

void EnvironmentCellProperty::UpdateEnvironmentVector(std::vector<double> environmentVector)
{
    mEnvironmentVector = environmentVector;
}

void EnvironmentCellProperty::UpdatePreferredEnvironmentVector(std::vector<double> environmentVector)
{
    mPreferredEnvironmentVector = environmentVector;
}


void EnvironmentCellProperty::PreparePostDivisionParent(double splitRatio)
{
    // split any properties that are shared
}
    
void  EnvironmentCellProperty::PreparePostDivisionDaughter(const EnvironmentCellProperty& parentProperty, double splitRatio)
{
    // split any properties that are shared
}

void EnvironmentCellProperty::SetUpEnvironmentVectors(unsigned numberOfStateVariables)
{
    std::vector<double> reset(numberOfStateVariables,0.0);
    mEnvironmentVector = reset;
    mPreferredEnvironmentVector = reset;
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

void EnvironmentCellProperty::SetEnvironmentValueByIndex(unsigned name_index,double value)
{
    mEnvironmentVector[name_index] = value;
}

void EnvironmentCellProperty::SetEnvironmentValueByName(std::string name, double value)
{

    if(mpEnvironmentStateVariableRegister -> IsStateVariablePresent(name))
    {
        mEnvironmentVector[mpEnvironmentStateVariableRegister -> RetrieveStateVariableIndex(name)] = value;
    }
    else
    {
        std::cout<<"Error: EnvironmentCellProperty::SetEnvironmentValueByName - name not found"<<std::endl;
    }

    return;
}



void EnvironmentCellProperty::SetPreferredEnvironmentVector(std::vector<double> valueVector)
{
    mPreferredEnvironmentVector = valueVector;
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

std::vector<double> EnvironmentCellProperty::GetPreferredEnvironmentVector()
{
    return mPreferredEnvironmentVector;
}

double EnvironmentCellProperty::GetEnvironmentValueByIndex(unsigned index)
{
    if(index<mpEnvironmentStateVariableRegister -> GetNumberOfStateVariables())
    {
        return mEnvironmentVector[index];
    }

    return 0.0; 
}

double EnvironmentCellProperty::GetPreferredEnvironmentValueByIndex(unsigned index)
{
    if(index<mpEnvironmentStateVariableRegister -> GetNumberOfStateVariables())
    {
        return mPreferredEnvironmentVector[index];
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

double EnvironmentCellProperty::GetPreferredEnvironmentValueByName(std::string stateName)
{
    unsigned index = mpEnvironmentStateVariableRegister->RetrieveStateVariableIndex(stateName);

    if(index<mpEnvironmentStateVariableRegister -> GetNumberOfStateVariables())
    {
        return mPreferredEnvironmentVector[index];
    }

    return 0.0; 
}


CellPtr EnvironmentCellProperty::GetCellPtr()
{
    return mThis_cellPtr;
}