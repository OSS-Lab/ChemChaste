#include "SchnackenbergOdeSystem.hpp"
#include "CellwiseOdeSystemInformation.hpp"

SchnackenbergOdeSystem::SchnackenbergOdeSystem(std::vector<double> stateVariables)
    : AbstractOdeSystem(2)
{
    mpSystemInfo.reset(new CellwiseOdeSystemInformation<SchnackenbergOdeSystem>);

    /**
     * The state variables are as follows:
     *
     * 0 - Notch concentration for this cell
     * 1 - Delta concentration for this cell
     *
     * We store the last state variable so that it can be written
     * to file at each time step alongside the others, and visualized.
     */

    SetDefaultInitialCondition(0, 1.0); // soon overwritten
    SetDefaultInitialCondition(1, 1.0); // soon overwritten

    this->mParameters.push_back(0.5);

    if (stateVariables != std::vector<double>())
    {
        SetStateVariables(stateVariables);
    }
}

SchnackenbergOdeSystem::~SchnackenbergOdeSystem()
{
}

void SchnackenbergOdeSystem::EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
{

    rDY[0] =   1.0*rY[1]*rY[0]*rY[0];
    rDY[1] = - 1.0*rY[1]*rY[0]*rY[0];
}

template<>
void CellwiseOdeSystemInformation<SchnackenbergOdeSystem>::Initialise()
{
    this->mVariableNames.push_back("U");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0); // will be filled in later

    this->mVariableNames.push_back("V");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0); // will be filled in later

    this->mInitialised = true;
}

