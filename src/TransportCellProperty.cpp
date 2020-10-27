#include "TransportCellProperty.hpp"


TransportCellProperty::TransportCellProperty()
    : AbstractCellProperty()
{
}


TransportCellProperty::~TransportCellProperty()
{
}

TransportCellProperty::TransportCellProperty(const TransportCellProperty& existingProperty)
{
    mpBulkStateVariableRegister = existingProperty.mpBulkStateVariableRegister;
    mpCellStateVariableRegister = existingProperty.mpCellStateVariableRegister;
    mpTransportReactionSystem = existingProperty.mpTransportReactionSystem;
    mpTransportOdeSystem = existingProperty.mpTransportOdeSystem;
    mpSolver = existingProperty.mpSolver;
    mIsConstantCellConcentration = existingProperty.mIsConstantCellConcentration;
    mBulkBoundaryConcentrationVector = existingProperty.mBulkBoundaryConcentrationVector;
    mChangeBulkBoundaryConcentrationVector = existingProperty.mChangeBulkBoundaryConcentrationVector;
    mCellBoundaryConcentrationVector = existingProperty.mCellBoundaryConcentrationVector;
    mChangeCellBoundaryConcentrationVector = existingProperty.mChangeCellBoundaryConcentrationVector;
    mIncludeOdeInterpolationOnBoundary = existingProperty.mIncludeOdeInterpolationOnBoundary;
}


// virtual methods

void TransportCellProperty::SetUp(AbstractTransportReactionSystem* transportReactionSystem, CellPtr this_cellPtr)
{
    UpdateTransportReactionSystem(transportReactionSystem);

    SetCellPtr(this_cellPtr);

    AbstractTransportOdeSystem* p_ode_system = new AbstractTransportOdeSystem(transportReactionSystem);

    UpdateTransportOdeSystem(p_ode_system);

    StateVariableRegister* p_bulk_stateRegister = new StateVariableRegister(transportReactionSystem -> GetBulkChemistry() -> GetChemicalNames());

    StateVariableRegister* p_cell_stateRegister = new StateVariableRegister(transportReactionSystem -> GetCellChemistry() -> GetChemicalNames());
   
    SetBulkStateVariableRegister(p_bulk_stateRegister);
    SetCellStateVariableRegister(p_cell_stateRegister);

    SetUpCellConcentrationVector(p_cell_stateRegister -> GetNumberOfStateVariables());
    SetUpChangeCellConcentrationVector(p_cell_stateRegister -> GetNumberOfStateVariables());
    SetUpChangeBulkConcentrationVector(p_bulk_stateRegister -> GetNumberOfStateVariables());
}

void TransportCellProperty::UpdateCellConcentrationVector(std::vector<double> cellBoundaryConcentrationVector)
{
    mCellBoundaryConcentrationVector = cellBoundaryConcentrationVector;
}


void TransportCellProperty::SetUpCellConcentrationVector(unsigned numberOfStateVariables)
{
    std::vector<double> reset(numberOfStateVariables,0.0);
    mCellBoundaryConcentrationVector = reset;
}

void TransportCellProperty::SetUpChangeCellConcentrationVector(unsigned numberOfStateVariables)
{
    std::vector<double> reset(numberOfStateVariables,0.0);
    mChangeCellBoundaryConcentrationVector = reset;
}

void TransportCellProperty::SetUpChangeBulkConcentrationVector(unsigned numberOfStateVariables)
{
    std::vector<double> reset(numberOfStateVariables,0.0);
    mChangeBulkBoundaryConcentrationVector = reset;
}

void TransportCellProperty::PerformTransportSystem(const std::vector<double>& currentBulkConcentration, const std::vector<double>& currentCellConcentration, std::vector<double>& changeBulkConc, std::vector<double>& changeCellConc)
{
    mpTransportReactionSystem->ReactSystem(currentBulkConcentration, currentCellConcentration, changeBulkConc, changeCellConc);
}

void TransportCellProperty::UpdateTransportReactionSystem(AbstractTransportReactionSystem* p_system)
{
    mpTransportReactionSystem = p_system;
    boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);
    SetTransportOdeSolver(p_solver);
}   

void TransportCellProperty::UpdateTransportOdeSystem(AbstractTransportOdeSystem* p_ode_system)
{
    mpTransportOdeSystem = p_ode_system;
}

void TransportCellProperty::PreparePostDivisionParent(double splitRatio)
{
    // split any properties that are shared
}
    
void  TransportCellProperty::PreparePostDivisionDaughter(const TransportCellProperty& parentProperty, double splitRatio)
{
    // split any properties that are shared
}

// concrete methods

double TransportCellProperty::RetrieveBoundarySourceByStateName(std::string stateName)
{
    // if cell is on boundary add the result of the transport Ode system
    unsigned index = mpBulkStateVariableRegister -> RetrieveStateVariableIndex(stateName);
    //15/10/2020
    return GetExternalCellBoundaryConcentrationByIndex(index);// - mInitBulkBoundaryConcentrationVector[index];
}

double TransportCellProperty::RetrieveChangeBoundarySourceByStateName(std::string stateName)
{
    // if cell is on boundary add the result of the transport Ode system
    unsigned index = mpBulkStateVariableRegister -> RetrieveStateVariableIndex(stateName);
    //15/10/2020
    return GetChangeExternalCellBoundaryConcentrationByIndex(index);// - mInitBulkBoundaryConcentrationVector[index];
}



void TransportCellProperty::AppendInternalCellBoundaryConcentrations(std::vector<double>& rY)
{
    rY.insert(rY.end(), mCellBoundaryConcentrationVector.begin(), mCellBoundaryConcentrationVector.end());
}

void TransportCellProperty::ReplaceBoundaryStateVariables(std::vector<double>& rY)
{
    //  partition the ODE state variable vector into the constitient internal and bulk state vectors
    unsigned number_of_bulk_states = mBulkBoundaryConcentrationVector.size();
    unsigned number_of_cell_states = mCellBoundaryConcentrationVector.size();

    for(unsigned i=0; i<number_of_bulk_states; i++)
    {
        mBulkBoundaryConcentrationVector[i] = rY[i];
    }

    for(unsigned i=0; i<number_of_cell_states; i++)
    {
        mCellBoundaryConcentrationVector[i] = rY[i + number_of_bulk_states];
    }

}

void TransportCellProperty::ReplaceChangeBoundaryStateVariables(std::vector<double>& rDY)
{
    //  partition the ODE state variable vector into the constitient internal and bulk state vectors
    unsigned number_of_bulk_states = mBulkBoundaryConcentrationVector.size();
    unsigned number_of_cell_states = mCellBoundaryConcentrationVector.size();

    for(unsigned i=0; i<number_of_bulk_states; i++)
    {

        mChangeBulkBoundaryConcentrationVector[i] = rDY[i];
    }

    for(unsigned i=0; i<number_of_cell_states; i++)
    {

        mChangeCellBoundaryConcentrationVector[i] = rDY[i + number_of_bulk_states];
    }

}


// set methods

void TransportCellProperty::SetBulkStateVariableRegister(StateVariableRegister* p_stateRegister)
{
    mpBulkStateVariableRegister = p_stateRegister;
}

void TransportCellProperty::SetCellStateVariableRegister(StateVariableRegister* p_stateRegister)
{
    mpCellStateVariableRegister = p_stateRegister;
}

void TransportCellProperty::SetConstantCellConcentrationBool(bool isConstantCellConcentration)
{
    mIsConstantCellConcentration = isConstantCellConcentration;
}

void TransportCellProperty::SetCellBoundaryConcentrationVector(std::vector<double> concentrationVector)
{
    mCellBoundaryConcentrationVector = concentrationVector;
}

void TransportCellProperty::SetBulkBoundaryConcentrationVector(std::vector<double> concentrationVector)
{
    mBulkBoundaryConcentrationVector = concentrationVector;
}

//15/10/2020
void TransportCellProperty::SetInitBulkBoundaryConcentrationVector(std::vector<double> concentrationVector)
{
    mInitBulkBoundaryConcentrationVector = concentrationVector;
}

void TransportCellProperty::SetTransportOdeSolver(boost::shared_ptr<AbstractIvpOdeSolver> p_solver)
{
    mpSolver = p_solver;
}

void TransportCellProperty::SetIncludeOdeInterpolationOnBoundary(bool includeOdeInterpolationOnBoundary)
{
    mIncludeOdeInterpolationOnBoundary = includeOdeInterpolationOnBoundary;
}

void TransportCellProperty::SetCellPtr(CellPtr this_cellPtr)
{
    mThis_cellPtr = this_cellPtr;
}

// get methods

StateVariableRegister* TransportCellProperty::GetBulkStateVariableRegister()
{
    return mpBulkStateVariableRegister;
}

StateVariableRegister* TransportCellProperty::GetCellStateVariableRegister()
{
    return mpCellStateVariableRegister;
}

bool TransportCellProperty::GetConstantCellConcentrationBool()
{
    return mIsConstantCellConcentration;
}

std::vector<double> TransportCellProperty::GetInternalCellBoundaryConcentrationVector()
{
    return mCellBoundaryConcentrationVector;
}

double TransportCellProperty::GetInternalCellBoundaryConcentrationByIndex(unsigned index)
{
    if(index<mpCellStateVariableRegister -> GetNumberOfStateVariables())
    {
        return mCellBoundaryConcentrationVector[index];
    }

    return 0.0; 
}

double TransportCellProperty::GetInternalCellBoundaryConcentrationByName(std::string stateName)
{
    unsigned index = mpCellStateVariableRegister->RetrieveStateVariableIndex(stateName);

    if(index<mpCellStateVariableRegister -> GetNumberOfStateVariables())
    {
        return mCellBoundaryConcentrationVector[index];
    }

    return 0.0; 
}


double TransportCellProperty::GetChangeInternalCellBoundaryConcentrationByName(std::string stateName)
{
    unsigned index = mpCellStateVariableRegister->RetrieveStateVariableIndex(stateName);

    if(index<mpCellStateVariableRegister -> GetNumberOfStateVariables())
    {
        return mChangeCellBoundaryConcentrationVector[index];
    }

    return 0.0; 
}


std::vector<double> TransportCellProperty::GetExternalCellBoundaryConcentrationVector()
{
    return mBulkBoundaryConcentrationVector;
}

double TransportCellProperty::GetExternalCellBoundaryConcentrationByIndex(unsigned index)
{
    if(index<mBulkBoundaryConcentrationVector.size())
    {
        return mBulkBoundaryConcentrationVector[index];
    }

    return 0.0; 
}

double TransportCellProperty::GetChangeExternalCellBoundaryConcentrationByIndex(unsigned index)
{
    if(index<mChangeBulkBoundaryConcentrationVector.size())
    {
        return mChangeBulkBoundaryConcentrationVector[index];
    }

    return 0.0; 
}

double TransportCellProperty::GetExternalCellBoundaryConcentrationByName(std::string stateName)
{
    unsigned index = mpBulkStateVariableRegister->RetrieveStateVariableIndex(stateName);

    if(index<mBulkBoundaryConcentrationVector.size())
    {
        return mBulkBoundaryConcentrationVector[index];
    }

    return 0.0; 
}

AbstractTransportReactionSystem* TransportCellProperty::GetTransportReactionSystem()
{
    return mpTransportReactionSystem;
}

AbstractTransportOdeSystem* TransportCellProperty::GetTransportOdeSystem()
{
    return mpTransportOdeSystem;
}

boost::shared_ptr<AbstractIvpOdeSolver> TransportCellProperty::GetTransportOdeSolver()
{
    return mpSolver;
}

bool TransportCellProperty::GetIncludeOdeInterpolationOnBoundary()
{
    return mIncludeOdeInterpolationOnBoundary;
}

CellPtr TransportCellProperty::GetCellPtr()
{
    return mThis_cellPtr;
}