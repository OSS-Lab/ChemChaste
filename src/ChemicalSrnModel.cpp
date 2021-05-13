#include "ChemicalSrnModel.hpp"

ChemicalSrnModel::ChemicalSrnModel(AbstractReactionSystem* pReactionSystem,boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver)
    : AbstractOdeSrnModel(pReactionSystem -> GetSystemChemistry() -> GetNumberChemicals(), pOdeSolver),
    mpReactionSystem(pReactionSystem),
    //mpChemicalOdeSystem(new AbstractChemicalOdeSystem(pReactionSystem)),
    mpCellChemistry(pReactionSystem -> GetSystemChemistry())
{
    //std::cout<<"ChemicalSrnModel::ChemicalSrnModel - start"<<std::endl;
    if (mpOdeSolver == boost::shared_ptr<AbstractCellCycleModelOdeSolver>())
    {
#ifdef CHASTE_CVODE
        mpOdeSolver = CellCycleModelOdeSolver<ChemicalSrnModel, CvodeAdaptor>::Instance();
        mpOdeSolver->Initialise();
        mpOdeSolver->SetMaxSteps(10000);
#else
        mpOdeSolver = CellCycleModelOdeSolver<ChemicalSrnModel, RungeKutta4IvpOdeSolver>::Instance();
        mpOdeSolver->Initialise();
        SetDt(0.001);
#endif //CHASTE_CVODE
    }
    assert(mpOdeSolver->IsSetUp());
    //std::cout<<"ChemicalSrnModel::ChemicalSrnModel - end"<<std::endl;
}

ChemicalSrnModel::ChemicalSrnModel(const ChemicalSrnModel& rModel)
    : AbstractOdeSrnModel(rModel)
{//std::cout<<"ChemicalSrnModel::ChemicalSrnModel- copy - start"<<std::endl;
    assert(rModel.GetOdeSystem());

    mpReactionSystem = rModel.mpReactionSystem;

    mpCellChemistry = rModel.mpCellChemistry;

    SetOdeSystem(new AbstractChemicalOdeSystem(mpReactionSystem));

    std::vector<double> stateVector = rModel.GetOdeSystem()->rGetStateVariables();

    mpOdeSystem->SetStateVariables(stateVector);
    //std::cout<<"ChemicalSrnModel::ChemicalSrnModel- copy - end"<<std::endl;
}

AbstractSrnModel* ChemicalSrnModel::CreateSrnModel()
{
    return new ChemicalSrnModel(*this);
}

void ChemicalSrnModel::Initialise()
{
    AbstractOdeSrnModel::Initialise( new AbstractChemicalOdeSystem(mpReactionSystem));
}

void ChemicalSrnModel::SimulateToCurrentTime()
{
    //std::cout<<"ChemicalSrnModel::SimulateToCurrentTime() - start"<<std::endl;
    // Custom behaviour
    UpdateOdeParameters();

    // Run the ODE simulation as needed
    AbstractOdeSrnModel::SimulateToCurrentTime();
    //std::cout<<"ChemicalSrnModel::SimulateToCurrentTime() - end"<<std::endl;
}

void ChemicalSrnModel::SetReactionSystem(AbstractReactionSystem* reactionSystem)
{
    mpReactionSystem = reactionSystem;
}

AbstractReactionSystem* ChemicalSrnModel::GetReactionSystem()
{
    return mpReactionSystem;
}

void ChemicalSrnModel::SetCellChemistry(AbstractChemistry* chemistry)
{
    mpCellChemistry = chemistry;
}

AbstractChemistry* ChemicalSrnModel::GetCellChemistry()
{
    return mpCellChemistry;
}

void ChemicalSrnModel::UpdateOdeStatesFromCellData()
{
    //std::cout<<"ChemicalSrnModel::UpdateOdeStatesFromCellData() - start"<<std::endl;
    unsigned numberOfSpecies = mpCellChemistry -> GetNumberChemicals();

    double species_value=0.0;
    std::string species_name="";
    std::vector<double> current_state_values(numberOfSpecies,0.0);
    for(unsigned i=0; i<numberOfSpecies; i++)
    {

        species_name = mpCellChemistry -> GetChemicalNamesByIndex(i);

        species_value=mpCell->GetCellData()->GetItem(species_name);

        current_state_values[i] = species_value;
    }

    mpOdeSystem->SetStateVariables(current_state_values);
    //std::cout<<"ChemicalSrnModel::UpdateOdeStatesFromCellData() - end"<<std::endl;
}

void ChemicalSrnModel::UpdateOdeParameters()
{
}

double ChemicalSrnModel::GetStateValueByIndex(unsigned index)
{
    return mpOdeSystem->rGetStateVariables()[index];
}

double ChemicalSrnModel::GetStateValueByName(std::string name)
{
    unsigned index = mpCellChemistry -> GetChemicalIndexByName(name);
    return GetStateValueByIndex(index);
}

void ChemicalSrnModel::OutputSrnModelParameters(out_stream& rParamsFile)
{
    // No new parameters to output, so just call method on direct parent class
    AbstractOdeSrnModel::OutputSrnModelParameters(rParamsFile);
}

