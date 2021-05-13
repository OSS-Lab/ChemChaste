#include "ChemicalCellProperty.hpp"

ChemicalCellProperty::ChemicalCellProperty()
: AbstractCellProperty()
{

}

ChemicalCellProperty::~ChemicalCellProperty()
{
}

void ChemicalCellProperty::InitialiseCell(std::vector<std::string> stateNameVector, std::vector<double> concentrationVector)
{
    SetStateVariableRegister(new StateVariableRegister(stateNameVector));

    UpdateCellConcentrationVector(concentrationVector);
}

void ChemicalCellProperty::InitialiseCell(StateVariableRegister* p_StateVariableRegister, std::vector<double> concentrationVector)
{
    SetStateVariableRegister(p_StateVariableRegister);

    UpdateCellConcentrationVector(concentrationVector);
}

void ChemicalCellProperty::UpdateCellConcentrationVector(std::vector<double>& concentrationVector)
{
    mConcentrationVector = concentrationVector;
}


void ChemicalCellProperty::SetStateVariableRegister(StateVariableRegister* p_register)
{
    mpStateVariableRegister = p_register;
}


StateVariableRegister* ChemicalCellProperty::GetStateVariableRegister()
{
    return mpStateVariableRegister;
}

std::vector<double> ChemicalCellProperty::GetCellConcentrationVector()
{
    return mConcentrationVector;
}

double ChemicalCellProperty::GetCellConcentrationByIndex(unsigned index)
{
    if(index < mConcentrationVector.size())
    {
        return mConcentrationVector[index];
    }
    
    return 0.0;
}

double ChemicalCellProperty::GetCellConcentrationByName(std::string name)
{
    if(mpStateVariableRegister->IsStateVariablePresent(name))
    {
        return mConcentrationVector[mpStateVariableRegister->RetrieveStateVariableIndex(name)];
    }
    
    return 0.0;
}
