#include "EulerIvpOdeSolver_test.hpp"


void EulerIvpOdeSolver_test::CalculateNextYValue(AbstractOdeSystem* pAbstractOdeSystem,
                                            double timeStep,
                                            double time,
                                            std::vector<double>& rCurrentYValues,
                                            std::vector<double>& rNextYValues)
{
    // For each timestep in AbstractOneStepIvpSolver calculates a vector containing
    // the next Y value from the current one for each equation in the system.

    const unsigned num_equations = pAbstractOdeSystem->GetNumberOfStateVariables();

    // Yes, this looks weird, but it makes good use of memory!
    pAbstractOdeSystem->EvaluateYDerivatives(time, rCurrentYValues, rNextYValues /*dydt is stored here*/);
    for (unsigned i=0; i<num_equations; i++)
    {
        // rNextYValues contains dY/dt until here
        rNextYValues[i] = rCurrentYValues[i] + timeStep*rNextYValues[i];
    }
}
