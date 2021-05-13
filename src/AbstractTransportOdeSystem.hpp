#ifndef ABSTRACTTRANSPORTODESYSTEM_HPP_
#define ABSTRACTTRANSPORTODESYSTEM_HPP_

#include "AbstractOdeSystem.hpp"
#include "AbstractTransportReactionSystem.hpp"
#include "TransportOdeSystemInformation.hpp"
#include "StateVariableRegister.hpp"


// ode system to model transport phenomena across a membrane seaprating two compartments
// the compartments external (bulk) and internal (cell) are recorded separately

class AbstractTransportOdeSystem: public AbstractOdeSystem
{
private:

    AbstractTransportReactionSystem* mpTransportReactionSystem;

    StateVariableRegister* mPdeStateVariableRegister;

    std::vector<double> mReservedCellConcentration;

    std::vector<double> mChangeCellConc;

    unsigned mNumberOfSpecies;

    unsigned mNumberOfReactions;

    double mDelta_error;

    bool mIsCheckConcentration;


public:

    AbstractTransportOdeSystem(AbstractTransportReactionSystem* pReactionSystem = new AbstractTransportReactionSystem());

    virtual ~AbstractTransportOdeSystem()
    {
    };
    
    AbstractTransportOdeSystem(const AbstractTransportOdeSystem&);

    virtual void EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY);

    virtual void UpdateReservedCellConcentration(std::vector<double>);

    virtual void UpdateChangeCellConcentration(std::vector<double>);

    unsigned GetNumberOfSpecies();

    unsigned GetNumberOfReactions();

    AbstractTransportReactionSystem* GetReactionSystem();

    void SetNumberOfSpecies(unsigned);

    void SetNumberOfReactions(unsigned);

    void SetReactionSystem(AbstractTransportReactionSystem* );

    
    void SetReservedCellConcentration(std::vector<double>);
    
    void SetChangeCellConcentration(std::vector<double>);

    void UpdateReservedCellConcentrationByIndex(unsigned, double);

    void UpdateReservedCellConcentrationByName(std::string, double);

    void UpdateChangeCellConcentrationByIndex(unsigned, double);

    void UpdateChangeCellConcentrationByName(std::string, double);

    void SetDeltaError(double);

    double GetDeltaError();

    void SetPdeStateRegister(StateVariableRegister*);

    void CheckConcentration(const std::vector<double>&);

    void SetIsCheckConcentration(bool);

    bool GetIsCheckConcentration();

    std::vector<double> GetReservedCellConcentration();
   
    std::vector<double> GetChangeCellConcentration();

    StateVariableRegister* GetPdeStateVariableRegister();

    double RetrieveReservedCellConcentrationByIndex(unsigned);

    double RetrieveReservedCellConcentrationByName(std::string);

    double RetrieveChangeCellConcentrationByIndex(unsigned);

    double RetrieveChangeCellConcentrationByName(std::string);

};

AbstractTransportOdeSystem::AbstractTransportOdeSystem(AbstractTransportReactionSystem* pReactionSystem)
        :   AbstractOdeSystem(pReactionSystem ->GetNumberOfBulkStates()+pReactionSystem ->GetNumberOfCellStates()),
            mpTransportReactionSystem(pReactionSystem)
{

    SetReactionSystem(mpTransportReactionSystem);
    SetNumberOfSpecies(mpTransportReactionSystem ->GetNumberOfBulkStates()+mpTransportReactionSystem ->GetNumberOfCellStates());
    SetNumberOfReactions(mpTransportReactionSystem -> GetNumberOfReactions());

    // initialsie cell concentration vectors
    std::vector<double> cellConcentration(mpTransportReactionSystem ->GetNumberOfCellStates(),1.0);
    SetReservedCellConcentration(cellConcentration);
    std::vector<double> cellConcentrationChange(mpTransportReactionSystem ->GetNumberOfCellStates(),0.0);
    SetChangeCellConcentration(cellConcentrationChange);

    mpSystemInfo.reset(new TransportOdeSystemInformation<AbstractTransportOdeSystem>(mpTransportReactionSystem));

    SetStateVariables(GetInitialConditions()); // set up in the solver call; size of NumberOfBulkStates + NumberOfCellStates

    SetIsCheckConcentration(true);
    SetDeltaError(1e-6);
}

AbstractTransportOdeSystem::AbstractTransportOdeSystem(const AbstractTransportOdeSystem& existingSystem)
    : AbstractOdeSystem(existingSystem.mNumberOfSpecies)
{
    mpTransportReactionSystem = existingSystem.mpTransportReactionSystem;

    mNumberOfSpecies = existingSystem.mNumberOfSpecies;

    mNumberOfReactions = existingSystem.mNumberOfReactions;

    mDelta_error = existingSystem.mDelta_error;

    mIsCheckConcentration = existingSystem.mIsCheckConcentration;

}


void AbstractTransportOdeSystem::EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
{
   
    unsigned number_of_bulk_states=mpTransportReactionSystem ->GetNumberOfBulkStates();
    unsigned number_of_cell_states=mpTransportReactionSystem ->GetNumberOfCellStates();

    std::vector<double> bulkConcentrations(number_of_bulk_states,0.0);
    std::vector<double> changeBulkConcentrations(number_of_bulk_states,0.0);

    std::vector<double> cellConcentrations(number_of_cell_states,0.0);
    std::vector<double> changeCellConcentrations(number_of_cell_states,0.0);
    
    // parse the different concentration state vectors from the solver output to the different concentration systems
    
    std::string this_state="";
    for(unsigned i=0; i<number_of_bulk_states; i++)
    {
        bulkConcentrations[i]=rY.at(i);
        //changeBulkConcentrations.push_back(rDY.at(i));
        /*
        this_state = mpTransportReactionSystem->GetBulkChemistry()->GetChemicalNamesByIndex(i);
        if(mPdeStateVariableRegister->IsStateVariablePresent(this_state))
        {
            bulkConcentrations[i] = rY.at(mPdeStateVariableRegister->RetrieveStateVariableIndex(this_state));
        }
        else
        {
            bulkConcentrations[i]=0.0;
        }
        */
    }

    for(unsigned i=0; i<number_of_cell_states; i++)
    {
        //mReservedCellConcentration[i]= rY[i+number_of_bulk_states];
        //cellConcentrations.push_back(rY.at(i+number_of_bulk_states));
        //changeCellConcentrations.push_back(rDY.at(i+number_of_bulk_states));

        cellConcentrations[i] = rY.at(i+number_of_bulk_states);
    }
    
    // reset rDY
    for(unsigned i=0; i<rDY.size();i++)
    {
        rDY[i] = 0.0;
    }
    
    if(mIsCheckConcentration)
    {   
        CheckConcentration(bulkConcentrations);
        CheckConcentration(cellConcentrations);
        //CheckConcentration(mReservedCellConcentration);
    }
    
    // reset mChangeCellConc, bulkConcentrations already reset
    //for(unsigned i=0; i<mNumberOfSpecies; i++)
    //{
    //    mChangeCellConc[i] =0.0;
    //}

    mpTransportReactionSystem -> ReactSystem(bulkConcentrations, cellConcentrations, changeBulkConcentrations, changeCellConcentrations);

    // reform rDY for passing to the solver
    for(unsigned i=0; i<number_of_bulk_states; i++)
    {
        rDY.at(i)= changeBulkConcentrations.at(i);
    }   

    for(unsigned i=0; i<number_of_cell_states; i++)
    {
        //rDY[i+number_of_bulk_states] = mChangeCellConc[i];
        rDY.at(i+number_of_bulk_states) = changeCellConcentrations.at(i);
    }
}  

void AbstractTransportOdeSystem::UpdateReservedCellConcentration(std::vector<double> reservedConcentration)
{
    SetReservedCellConcentration(reservedConcentration);
}

void AbstractTransportOdeSystem::UpdateChangeCellConcentration(std::vector<double> changeCellConc)
{
    mChangeCellConc = changeCellConc;
}

unsigned AbstractTransportOdeSystem::GetNumberOfSpecies()
{
    return mNumberOfSpecies;
}

unsigned AbstractTransportOdeSystem::GetNumberOfReactions()
{
    return mNumberOfReactions;
}

AbstractTransportReactionSystem* AbstractTransportOdeSystem::GetReactionSystem()
{
    return mpTransportReactionSystem;
}

void AbstractTransportOdeSystem::SetNumberOfSpecies(unsigned numberOfSpecies)
{
    mNumberOfSpecies = numberOfSpecies;
}

void AbstractTransportOdeSystem::SetNumberOfReactions(unsigned numberOfReactions)
{
    mNumberOfReactions = numberOfReactions;
}

void AbstractTransportOdeSystem::SetReactionSystem(AbstractTransportReactionSystem* p_reactionSystem)
{
    mpTransportReactionSystem = p_reactionSystem;
}

void AbstractTransportOdeSystem::SetReservedCellConcentration(std::vector<double> cellConcentration)
{
    mReservedCellConcentration = cellConcentration;
}
    
void AbstractTransportOdeSystem::SetChangeCellConcentration(std::vector<double> changeCellConcentration)
{
    mChangeCellConc = changeCellConcentration;
}

void AbstractTransportOdeSystem::UpdateReservedCellConcentrationByIndex(unsigned index, double concentration)
{
    mReservedCellConcentration[index] = concentration;
}

void AbstractTransportOdeSystem::UpdateReservedCellConcentrationByName(std::string name, double concentration)
{
    unsigned index = mpTransportReactionSystem -> GetCellChemistry() -> GetChemicalIndexByName(name);

    mReservedCellConcentration[index] = concentration;
}

void AbstractTransportOdeSystem::UpdateChangeCellConcentrationByIndex(unsigned index, double concentration)
{
    mChangeCellConc[index] = concentration;
}

void AbstractTransportOdeSystem::UpdateChangeCellConcentrationByName(std::string name, double concentration)
{
    unsigned index = mpTransportReactionSystem -> GetCellChemistry() -> GetChemicalIndexByName(name);

    mChangeCellConc[index] = concentration;
}


void AbstractTransportOdeSystem::SetDeltaError(double delta_error)
{
    mDelta_error = delta_error;
}

double AbstractTransportOdeSystem::GetDeltaError()
{
    return mDelta_error;
}

void AbstractTransportOdeSystem::SetPdeStateRegister(StateVariableRegister* p_register)
{
    mPdeStateVariableRegister = p_register;
}

void AbstractTransportOdeSystem::CheckConcentration(const std::vector<double>& rY)
{
    // if chemical concentration gets too low then NaNs can occur, concentration must be +ve

    for(unsigned i=0; i<rY.size(); i++)
    {
        // due to the discrete nature occasionally rY can evaluate to <0
        // ensure rY >= 0
        if(rY[i]<mDelta_error)
        {
            const_cast<double&>(rY[i]) = 0;
        }
    }

}

void AbstractTransportOdeSystem::SetIsCheckConcentration(bool IsCheckConcentration)
{
    mIsCheckConcentration = IsCheckConcentration;
}

bool AbstractTransportOdeSystem::GetIsCheckConcentration()
{
    return mIsCheckConcentration;
}

std::vector<double> AbstractTransportOdeSystem::GetReservedCellConcentration()
{
    return mReservedCellConcentration;
}

std::vector<double> AbstractTransportOdeSystem::GetChangeCellConcentration()
{
    return mChangeCellConc;
}

StateVariableRegister* AbstractTransportOdeSystem::GetPdeStateVariableRegister()
{
    return mPdeStateVariableRegister;
}

double AbstractTransportOdeSystem::RetrieveReservedCellConcentrationByIndex(unsigned index)
{
    return mReservedCellConcentration[index];
}

double AbstractTransportOdeSystem::RetrieveReservedCellConcentrationByName(std::string name)
{
    unsigned index = mpTransportReactionSystem -> GetCellChemistry() -> GetChemicalIndexByName(name);

    return mReservedCellConcentration[index];
}

double AbstractTransportOdeSystem::RetrieveChangeCellConcentrationByIndex(unsigned index)
{
    return mChangeCellConc[index];
}

double AbstractTransportOdeSystem::RetrieveChangeCellConcentrationByName(std::string name)
{
    unsigned index = mpTransportReactionSystem -> GetCellChemistry() -> GetChemicalIndexByName(name);

    return mChangeCellConc[index];
}




template<>
void TransportOdeSystemInformation<AbstractTransportOdeSystem>::Initialise()
{
    // initialise the bulk state varibles
    for( unsigned i=0; i<mp_reaction_system -> GetBulkChemistry() -> GetNumberChemicals(); i++)
    {
        this->mVariableNames.push_back(mp_reaction_system -> GetBulkChemistry() -> GetChemicalNamesByIndex(i));
        this->mVariableUnits.push_back(mp_reaction_system -> GetBulkChemistry() -> GetChemicalDimensionsByIndex(i));
        this->mInitialConditions.push_back(1.0); // will be overridden
    }
    
    // initialise the cell state varibles; appended onto the end of the bulk state vectors
    for( unsigned i=0; i<mp_reaction_system -> GetCellChemistry() -> GetNumberChemicals(); i++)
    {
        this->mVariableNames.push_back(mp_reaction_system -> GetCellChemistry() -> GetChemicalNamesByIndex(i));
        this->mVariableUnits.push_back(mp_reaction_system -> GetCellChemistry() -> GetChemicalDimensionsByIndex(i));
        this->mInitialConditions.push_back(1.0); // will be overridden
    }
    this->mInitialised = true;
}

#endif 