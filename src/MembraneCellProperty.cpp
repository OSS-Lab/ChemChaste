#include "MembraneCellProperty.hpp"

MembraneCellProperty::MembraneCellProperty()
    : AbstractCellProperty()
{

}

MembraneCellProperty::~MembraneCellProperty()
{

}

MembraneCellProperty::MembraneCellProperty(const MembraneCellProperty& existingProperty)
{
    mIsExtent = existingProperty.mIsExtent;
    mMembraneThickness = existingProperty.mMembraneThickness;
    mpMembraneStateVariableRegister = existingProperty.mpMembraneStateVariableRegister;
    mMembraneConcentrationVector = existingProperty.mMembraneConcentrationVector;
    mpBulkStateVariableRegister = existingProperty.mpBulkStateVariableRegister;
    mpCellStateVariableRegister = existingProperty.mpCellStateVariableRegister;
    mpMembraneReactionSystem = existingProperty.mpMembraneReactionSystem;
    mpMembraneOdeSystem = existingProperty.mpMembraneOdeSystem;
    mpMembraneOdeSolver = existingProperty.mpMembraneOdeSolver;
    mIsConstantCellConcentration = existingProperty.mIsConstantCellConcentration;
    mBulkBoundaryConcentrationVector = existingProperty.mBulkBoundaryConcentrationVector;
    mChangeBulkBoundaryConcentrationVector = existingProperty.mChangeBulkBoundaryConcentrationVector;
    mCellBoundaryConcentrationVector = existingProperty.mCellBoundaryConcentrationVector;
    mChangeCellBoundaryConcentrationVector = existingProperty.mChangeCellBoundaryConcentrationVector;
    mIncludeMembraneOdeInterpolationOnBoundary = existingProperty.mIncludeMembraneOdeInterpolationOnBoundary;
}

// virtual methods

void MembraneCellProperty::SetUp(AbstractMembraneReactionSystem* membraneReactionSystem, CellPtr this_cellPtr)
{
    UpdateMembraneReactionSystem(membraneReactionSystem);

    SetCellPtr(this_cellPtr);

    AbstractMembraneOdeSystem* p_ode_system = new AbstractMembraneOdeSystem(membraneReactionSystem);

    UpdateMembraneOdeSystem(p_ode_system);

    StateVariableRegister* p_bulk_stateRegister = new StateVariableRegister(membraneReactionSystem -> GetBulkChemistry() -> GetChemicalNames());

    StateVariableRegister* p_cell_stateRegister = new StateVariableRegister(membraneReactionSystem -> GetCellChemistry() -> GetChemicalNames());
    
    SetBulkStateVariableRegister(p_bulk_stateRegister);
    SetCellStateVariableRegister(p_cell_stateRegister);

    SetUpCellConcentrationVector(p_cell_stateRegister -> GetNumberOfStateVariables());
    SetUpChangeCellConcentrationVector(p_cell_stateRegister -> GetNumberOfStateVariables());
    SetUpChangeBulkConcentrationVector(p_bulk_stateRegister -> GetNumberOfStateVariables());
}

void MembraneCellProperty::UpdateMembraneConcentrationVector(std::vector<double> cellBoundaryConcentrationVector)
{
    mMembraneConcentrationVector = cellBoundaryConcentrationVector;
}

void MembraneCellProperty::UpdateCellConcentrationVector(std::vector<double> cellBoundaryConcentrationVector)
{
    mCellBoundaryConcentrationVector = cellBoundaryConcentrationVector;
}

void MembraneCellProperty::SetUpCellConcentrationVector(unsigned numberOfStateVariables)
{
    std::vector<double> reset(numberOfStateVariables,0.0);
    mCellBoundaryConcentrationVector = reset;
}

void MembraneCellProperty::SetUpChangeCellConcentrationVector(unsigned numberOfStateVariables)
{
    std::vector<double> reset(numberOfStateVariables,0.0);
    mChangeCellBoundaryConcentrationVector = reset;
}

void MembraneCellProperty::SetUpChangeBulkConcentrationVector(unsigned numberOfStateVariables)
{
    std::vector<double> reset(numberOfStateVariables,0.0);
    mChangeBulkBoundaryConcentrationVector = reset;
}

void MembraneCellProperty::PerformMembraneSystem(const std::vector<double>& currentBulkConcentration, const std::vector<double>& currentCellConcentration, std::vector<double>& changeBulkConc, std::vector<double>& changeCellConc)
{
    mpMembraneReactionSystem->ReactSystem(currentBulkConcentration, currentCellConcentration, changeBulkConc, changeCellConc);    
}

void MembraneCellProperty::UpdateMembraneReactionSystem(AbstractMembraneReactionSystem* p_system)
{
    mpMembraneReactionSystem = p_system;
    boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);
    SetMembraneOdeSolver(p_solver);
}

void MembraneCellProperty::UpdateMembraneOdeSystem(AbstractMembraneOdeSystem* p_ode_system)
{
    mpMembraneOdeSystem = p_ode_system;
}

void MembraneCellProperty::PreparePostDivisionParent(double splitRatio)
{
    // split any properties that are shared
    for(unsigned i=0; i<mMembraneConcentrationVector.size();i++)
    {
        mMembraneConcentrationVector[i] =splitRatio*mMembraneConcentrationVector[i];
    }
}

void MembraneCellProperty::PreparePostDivisionDaughter(const MembraneCellProperty& parentProperty ,double splitRatio)
{
    // split any properties that are shared
    std::vector<double> parentConcentrationVector = parentProperty.mMembraneConcentrationVector;

    for(unsigned i=0; i<parentConcentrationVector.size();i++)
    {
        mMembraneConcentrationVector[i] =(1-splitRatio)*parentConcentrationVector[i];
    }
}

void MembraneCellProperty::InitialiseMembrane(std::vector<std::string> stateNameVector, std::vector<double> concentrationVector)
{
    SetMembraneStateVariableRegister(new StateVariableRegister(stateNameVector));

    SetMembraneConcentrationVector(concentrationVector);
}

void MembraneCellProperty::InitialiseMembrane(StateVariableRegister* p_StateVariableRegister, std::vector<double> concentrationVector)
{
    SetMembraneStateVariableRegister(p_StateVariableRegister);

    SetMembraneConcentrationVector(concentrationVector);
}

// concrete methods

double MembraneCellProperty::RetrieveBoundarySourceByStateName(std::string stateName)
{
    // if cell is on boundary add the result of the transport Ode system
    unsigned index = mpBulkStateVariableRegister -> RetrieveStateVariableIndex(stateName);
    //15/10/2020
    return GetExternalCellBoundaryConcentrationByIndex(index);// - mInitBulkBoundaryConcentrationVector[index];
}

double MembraneCellProperty::RetrieveChangeBoundarySourceByStateName(std::string stateName)
{
    // if cell is on boundary add the result of the transport Ode system
    unsigned index = mpBulkStateVariableRegister -> RetrieveStateVariableIndex(stateName);
    //15/10/2020
    return GetChangeExternalCellBoundaryConcentrationByIndex(index);// - mInitBulkBoundaryConcentrationVector[index];
}

void MembraneCellProperty::AppendInternalCellBoundaryConcentrations(std::vector<double>& rY)
{
    //SetBulkBoundaryConcentrationVector(rY);

    rY.insert(rY.end(), mCellBoundaryConcentrationVector.begin(), mCellBoundaryConcentrationVector.end());
}

void MembraneCellProperty::ReplaceBoundaryStateVariables(std::vector<double>& rY)
{
    /*
    //  partition the ODE state variable vector into the constitient internal and bulk state vectors
    unsigned number_of_bulk_states = mBulkBoundaryConcentrationVector.size();
    unsigned number_of_cell_states = mCellBoundaryConcentrationVector.size();
    std::cout<<"ReplaceBoundaryStateVariables - numberBulk:"<<number_of_bulk_states<<std::endl;
    AbstractChemistry* p_bulk_chemistry = mpMembraneReactionSystem->GetBulkChemistry();
    std::string this_state="";
    unsigned number_bulk_reaction_states = p_bulk_chemistry -> GetNumberChemicals();;
    unsigned bulk_chemical_index = 0;
    for(unsigned i=0; i<number_bulk_reaction_states; i++)
    {
        this_state = p_bulk_chemistry ->GetChemicalNamesByIndex(i);
  
        
        if(mpMembraneOdeSystem->GetPdeStateVariableRegister()->IsStateVariablePresent(this_state))
        {
            bulk_chemical_index = mpMembraneOdeSystem->GetPdeStateVariableRegister()->RetrieveStateVariableIndex(this_state);
            mBulkBoundaryConcentrationVector[bulk_chemical_index] = rY.at(bulk_chemical_index);
        }
        std::cout<<"Here1"<<std::endl;
    }

    for(unsigned i=0; i<number_of_cell_states; i++)
    {
        mCellBoundaryConcentrationVector[i] = rY[i + mpMembraneOdeSystem->GetPdeStateVariableRegister()->GetNumberOfStateVariables()];
    }
    */

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

void MembraneCellProperty::ReplaceChangeBoundaryStateVariables(std::vector<double>& rDY)
{
    //  partition the ODE state variable vector into the constitient internal and bulk state vectors
    unsigned number_of_bulk_states = mBulkBoundaryConcentrationVector.size();
    unsigned number_of_cell_states = mCellBoundaryConcentrationVector.size();

 
    AbstractChemistry* p_bulk_chemistry = mpMembraneReactionSystem->GetBulkChemistry();
    std::string this_state="";
    unsigned number_bulk_reaction_states = p_bulk_chemistry -> GetNumberChemicals();


    for(unsigned i=0; i<number_bulk_reaction_states; i++)
    {
        mChangeBulkBoundaryConcentrationVector[i] = rDY[i];
    }

    for(unsigned i=0; i<number_of_cell_states; i++)
    {
        mChangeCellBoundaryConcentrationVector[i] = rDY[i + number_bulk_reaction_states];
    }

}

// set methods

void MembraneCellProperty::SetMembraneExtentBool(bool isExtent)
{
    mIsExtent = isExtent;
}

void MembraneCellProperty::SetMembraneThickness(double thickness)
{
    mMembraneThickness = thickness;
}

void MembraneCellProperty::SetDoubleMembraneBool(bool isDouble)
{
    mIsDoubleMembrane = isDouble;
}

void MembraneCellProperty::SetBulkStateVariableRegister(StateVariableRegister* p_stateRegister)
{
    mpBulkStateVariableRegister = p_stateRegister;
}

void MembraneCellProperty::SetCellStateVariableRegister(StateVariableRegister* p_stateRegister)
{
    mpCellStateVariableRegister = p_stateRegister;
}

void MembraneCellProperty::SetConstantCellConcentrationBool(bool isConstantCellConcentration)
{
    mIsConstantCellConcentration = isConstantCellConcentration;
}

void MembraneCellProperty::SetCellBoundaryConcentrationVector(std::vector<double> concentrationVector)
{
    mCellBoundaryConcentrationVector = concentrationVector;
}

void MembraneCellProperty::SetBulkBoundaryConcentrationVector(std::vector<double> concentrationVector)
{
    mBulkBoundaryConcentrationVector = concentrationVector;
}

//15/10/2020
void MembraneCellProperty::SetInitBulkBoundaryConcentrationVector(std::vector<double> concentrationVector)
{
    mInitBulkBoundaryConcentrationVector = concentrationVector;
}

void MembraneCellProperty::SetMembraneOdeSolver(boost::shared_ptr<AbstractIvpOdeSolver> p_solver)
{
    mpMembraneOdeSolver = p_solver;
}
 
void MembraneCellProperty::SetIncludeMembraneOdeInterpolationOnBoundary(bool includeOdeInterpolationOnBoundary)
{
    mIncludeMembraneOdeInterpolationOnBoundary = includeOdeInterpolationOnBoundary;
}

void MembraneCellProperty::SetCellPtr(CellPtr this_cellPtr)
{
    m_membrane_cellPtr = this_cellPtr;
}

// get methods

bool MembraneCellProperty::GetMembraneExtentBool()
{
    return mIsExtent;
}

double MembraneCellProperty::GetMembraneThickness()
{
    return mMembraneThickness;
}

bool MembraneCellProperty::GetDoubleMembraneBool()
{
    return mIsDoubleMembrane;
}

StateVariableRegister* MembraneCellProperty::GetBulkStateVariableRegister()
{
    return mpBulkStateVariableRegister;
}

StateVariableRegister* MembraneCellProperty::GetCellStateVariableRegister()
{
    return mpCellStateVariableRegister;
}

bool MembraneCellProperty::GetConstantCellConcentrationBool()
{
    return mIsConstantCellConcentration;
}

std::vector<double> MembraneCellProperty::GetInternalCellBoundaryConcentrationVector()
{
    return mCellBoundaryConcentrationVector;
}

double MembraneCellProperty::GetInternalCellBoundaryConcentrationByIndex(unsigned index)
{
    if(index<mpCellStateVariableRegister -> GetNumberOfStateVariables())
    {
        return mCellBoundaryConcentrationVector[index];
    }

    return 0.0; 
}

double MembraneCellProperty::GetInternalCellBoundaryConcentrationByName(std::string stateName)
{
    unsigned index = mpCellStateVariableRegister->RetrieveStateVariableIndex(stateName);

    if(index<mpCellStateVariableRegister -> GetNumberOfStateVariables())
    {
        return mCellBoundaryConcentrationVector[index];
    }

    return 0.0; 
}

double MembraneCellProperty::GetChangeInternalCellBoundaryConcentrationByName(std::string stateName)
{
    unsigned index = mpCellStateVariableRegister->RetrieveStateVariableIndex(stateName);

    if(index<mpCellStateVariableRegister -> GetNumberOfStateVariables())
    {
        return mChangeCellBoundaryConcentrationVector[index];
    }

    return 0.0; 
}

std::vector<double> MembraneCellProperty::GetExternalCellBoundaryConcentrationVector()
{
    return mBulkBoundaryConcentrationVector;
}

double MembraneCellProperty::GetExternalCellBoundaryConcentrationByIndex(unsigned index)
{
    if(index<mBulkBoundaryConcentrationVector.size())
    {
        return mBulkBoundaryConcentrationVector[index];
    }

    return 0.0; 
}

double MembraneCellProperty::GetChangeExternalCellBoundaryConcentrationByIndex(unsigned index)
{
    if(index<mChangeBulkBoundaryConcentrationVector.size())
    {
        return mChangeBulkBoundaryConcentrationVector[index];
    }

    return 0.0; 
}

double MembraneCellProperty::GetExternalCellBoundaryConcentrationByName(std::string stateName)
{
    unsigned index = mpBulkStateVariableRegister->RetrieveStateVariableIndex(stateName);

    if(index<mBulkBoundaryConcentrationVector.size())
    {
        return mBulkBoundaryConcentrationVector[index];
    }

    return 0.0; 
}

AbstractMembraneReactionSystem* MembraneCellProperty::GetMembraneReactionSystem()
{
    return mpMembraneReactionSystem;
}

AbstractMembraneOdeSystem* MembraneCellProperty::GetMembraneOdeSystem()
{
    return mpMembraneOdeSystem;
}

boost::shared_ptr<AbstractIvpOdeSolver> MembraneCellProperty::GetMembraneOdeSolver()
{
    return mpMembraneOdeSolver;
}

bool MembraneCellProperty::GetIncludeMembraneOdeInterpolationOnBoundary()
{
    return mIncludeMembraneOdeInterpolationOnBoundary;
}

CellPtr MembraneCellProperty::GetCellPtr()
{
    return m_membrane_cellPtr;
}


// membrane parameters

void MembraneCellProperty::SetMembraneStateVariableRegister(StateVariableRegister* p_register)
{
    mpMembraneStateVariableRegister = p_register;
}

void MembraneCellProperty::SetMembraneConcentrationVector(std::vector<double> concentrationVector)
{
    mMembraneConcentrationVector = concentrationVector;
}

StateVariableRegister* MembraneCellProperty::GetMembraneStateVariableRegister()
{
    return mpMembraneStateVariableRegister;
}

std::vector<double> MembraneCellProperty::GetMembraneConcentrationVector()
{
    return mMembraneConcentrationVector;
}

double MembraneCellProperty::GetMembraneConcentrationByIndex(unsigned index)
{
    if(index < mMembraneConcentrationVector.size())
    {
        return mMembraneConcentrationVector[index];
    }
    
    return 0.0;
}

double MembraneCellProperty::GetMembraneConcentrationByName(std::string name)
{
    if(mpMembraneStateVariableRegister->IsStateVariablePresent(name))
    {
        return mMembraneConcentrationVector[mpMembraneStateVariableRegister->RetrieveStateVariableIndex(name)];
    }
    
    return 0.0;
}
