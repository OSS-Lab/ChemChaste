#include "SchnackenbergSrnModel.hpp"
#include "CellCycleModelOdeSolver.hpp"

SchnackenbergSrnModel::SchnackenbergSrnModel(boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver)
    : AbstractOdeSrnModel(2, pOdeSolver)
{
    if (mpOdeSolver == boost::shared_ptr<AbstractCellCycleModelOdeSolver>())
    {
#ifdef CHASTE_CVODE
        mpOdeSolver = CellCycleModelOdeSolver<SchnackenbergSrnModel, CvodeAdaptor>::Instance();
        mpOdeSolver->Initialise();
        mpOdeSolver->SetMaxSteps(10000);
#else
        mpOdeSolver = CellCycleModelOdeSolver<SchnackenbergSrnModel, RungeKutta4IvpOdeSolver>::Instance();
        mpOdeSolver->Initialise();
        SetDt(0.001);
#endif //CHASTE_CVODE
    }
    assert(mpOdeSolver->IsSetUp());
}

SchnackenbergSrnModel::SchnackenbergSrnModel(const SchnackenbergSrnModel& rModel)
    : AbstractOdeSrnModel(rModel)
{

    assert(rModel.GetOdeSystem());
    SetOdeSystem(new SchnackenbergOdeSystem(rModel.GetOdeSystem()->rGetStateVariables()));
}

AbstractSrnModel* SchnackenbergSrnModel::CreateSrnModel()
{
    return new SchnackenbergSrnModel(*this);
}

void SchnackenbergSrnModel::SimulateToCurrentTime()
{

    // Run the ODE simulation as needed
    AbstractOdeSrnModel::SimulateToCurrentTime();
}

void SchnackenbergSrnModel::Initialise()
{
    AbstractOdeSrnModel::Initialise(new SchnackenbergOdeSystem);
}

double SchnackenbergSrnModel::GetU()
{
    assert(mpOdeSystem != nullptr);
    double U = mpOdeSystem->rGetStateVariables()[0];
    return U;
}

double SchnackenbergSrnModel::GetV()
{
    assert(mpOdeSystem != nullptr);
    double V = mpOdeSystem->rGetStateVariables()[1];
    return V;
}


void SchnackenbergSrnModel::OutputSrnModelParameters(out_stream& rParamsFile)
{
    // No new parameters to output, so just call method on direct parent class
    AbstractOdeSrnModel::OutputSrnModelParameters(rParamsFile);
}
