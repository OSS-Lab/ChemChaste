#ifndef AbstractTransportOdeForCoupledPdeSystem_HPP_
#define AbstractTransportOdeForCoupledPdeSystem_HPP_

#include "AbstractOdeSystemForCoupledPdeSystem.hpp"
#include "AbstractTransportReactionSystem.hpp"
#include "ChemicalOdeSystemInformation.hpp"


class AbstractTransportOdeForCoupledPdeSystem: public AbstractOdeSystemForCoupledPdeSystem
{
private:

    // reaction system object which contains the data structures and functions regarding the reactions
    AbstractTransportReactionSystem* mpReactionSystem;

    std::vector<double> mReservedCellConcentration;

    std::vector<double> mChangeCellConc;

    unsigned mNumberOfSpecies;

    unsigned mNumberOfReactions;
    // error tolerence for checking whether the concentrations are >0, i.e =0 if <tolerence
    double mDelta_error; 
    // whether to check for negative concentrations
    bool mIsCheckConcentration;


public:

    AbstractTransportOdeForCoupledPdeSystem(AbstractTransportReactionSystem* pReactionSystem = new AbstractTransportReactionSystem());

    virtual ~AbstractTransportOdeForCoupledPdeSystem()
    {
    };
    
    AbstractTransportOdeForCoupledPdeSystem(const AbstractTransportOdeForCoupledPdeSystem&);

    // virtual methods

    virtual void EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY);

    virtual void UpdateReservedCellConcentration();

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

    void CheckConcentration(const std::vector<double>&);

    void SetIsCheckConcentration(bool);

    bool GetIsCheckConcentration();

    std::vector<double> GetReservedCellConcentration();
   
    std::vector<double> GetChangeCellConcentration();

    double RetrieveReservedCellConcentrationByIndex(unsigned);

    double RetrieveReservedCellConcentrationByName(std::string);

    double RetrieveChangeCellConcentrationByIndex(unsigned);

    double RetrieveChangeCellConcentrationByName(std::string);

};

AbstractTransportOdeForCoupledPdeSystem::AbstractTransportOdeForCoupledPdeSystem(AbstractTransportReactionSystem* pReactionSystem)
        :   AbstractOdeSystemForCoupledPdeSystem(pReactionSystem ->GetNumberOfBulkStates()+pReactionSystem ->GetNumberOfCellStates(),pReactionSystem ->GetNumberOfBulkStates())            
{

    SetReactionSystem(pReactionSystem);
    SetNumberOfSpecies(pReactionSystem ->GetNumberOfBulkStates()+pReactionSystem ->GetNumberOfCellStates());
    SetNumberOfReactions(pReactionSystem -> GetNumberOfReactions());

    // initialise cell concentration vectors
    std::vector<double> cellConcentration(pReactionSystem ->GetNumberOfCellStates(),1.0);
    SetReservedCellConcentration(cellConcentration);
    std::vector<double> cellConcentrationChange(pReactionSystem ->GetNumberOfCellStates(),0.0);
    SetChangeCellConcentration(cellConcentrationChange);

    mpSystemInfo.reset(new TransportOdeSystemInformation<AbstractTransportOdeForCoupledPdeSystem>(pReactionSystem));

    SetStateVariables(GetInitialConditions()); // set up in the solver call; size of NumberOfBulkStates + NumberOfCellStates

    SetIsCheckConcentration(true);
    SetDeltaError(1e-6);
}

AbstractTransportOdeForCoupledPdeSystem::AbstractTransportOdeForCoupledPdeSystem(const AbstractTransportOdeForCoupledPdeSystem& existingSystem)
    : AbstractOdeSystemForCoupledPdeSystem(existingSystem.mNumberOfSpecies,existingSystem.mNumberOfSpecies)
{
    mpReactionSystem = existingSystem.mpReactionSystem;

    mNumberOfSpecies = existingSystem.mNumberOfSpecies;

    mNumberOfReactions = existingSystem.mNumberOfReactions;

    mDelta_error = existingSystem.mDelta_error;

    mIsCheckConcentration = existingSystem.mIsCheckConcentration;
}


void AbstractTransportOdeForCoupledPdeSystem::EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
{
    unsigned number_of_bulk_states=pReactionSystem ->GetNumberOfBulkStates();

    std::vector<double> bulkConcentrations(number_of_bulk_states,0.0);
    std::vector<double> changeBulkConcentrations(number_of_bulk_states,0.0);

    // parse the different concentration state vectors from the solver output to the different concentration systems
    for(unsigned i=0; i<number_of_bulk_states; i++)
    {
        bulkConcentrations[i]= rY[i];
    }

    for(unsigned i=0; i<pReactionSystem ->GetNumberOfCellStates(); i++)
    {
        mReservedCellConcentration[i]= rY[i+number_of_bulk_states];
    }

    if(mIsCheckConcentration)
    {
        CheckConcentration(bulkConcentrations);
        CheckConcentration(mReservedCellConcentration);
    }

    // reset mChangeCellConc, bulkConcentrations already reset
    for(unsigned i=0; i<mNumberOfSpecies; i++)
    {
        mChangeCellConc[i] =0.0;
    }
    
    mpReactionSystem -> ReactSystem(bulkConcentrations, mReservedCellConcentration, changeBulkConcentrations, mChangeCellConc);

    // reform rDY for passing to the solver
    for(unsigned i=0; i<number_of_bulk_states; i++)
    {
        rDY[i]= changeBulkConcentrations[i];
    }
    for(unsigned i=0; i<pReactionSystem ->GetNumberOfCellStates(); i++)
    {
        rDY[i+number_of_bulk_states] = mChangeCellConc[i];
    }

}  

void AbstractTransportOdeForCoupledPdeSystem::UpdateReservedCellConcentration(std::vector<double> reservedConcentration)
{
    mReservedCellConcentration = reservedConcentration
}

void AbstractTransportOdeForCoupledPdeSystem::UpdateChangeCellConcentration(std::vector<double> changeCellConc)
{
    mChangeCellConc = changeCellConc;
}

unsigned AbstractTransportOdeForCoupledPdeSystem::GetNumberOfSpecies()
{
    return mNumberOfSpecies;
}

unsigned AbstractTransportOdeForCoupledPdeSystem::GetNumberOfReactions()
{
    return mNumberOfReactions;
}

AbstractTransportReactionSystem* AbstractTransportOdeForCoupledPdeSystem::GetReactionSystem()
{
    return mpReactionSystem;
}

void AbstractTransportOdeForCoupledPdeSystem::SetNumberOfSpecies(unsigned numberOfSpecies)
{
    mNumberOfSpecies = numberOfSpecies;
}

void AbstractTransportOdeForCoupledPdeSystem::SetNumberOfReactions(unsigned numberOfReactions)
{
    mNumberOfReactions = numberOfReactions;
}

void AbstractTransportOdeForCoupledPdeSystem::SetReactionSystem(AbstractTransportReactionSystem* p_reactionSystem)
{
    mpReactionSystem = p_reactionSystem;
}

void AbstractTransportOdeForCoupledPdeSystem::SetReservedCellConcentration(std::vector<double> cellConcentration)
{
    mReservedCellConcentration = cellConcentration;
}
    
void AbstractTransportOdeForCoupledPdeSystem::SetChangeCellConcentration(std::vector<double> changeCellConcentration)
{
    mChangeCellConc = changeCellConcentration;
}

void AbstractTransportOdeForCoupledPdeSystem::UpdateReservedCellConcentrationByIndex(unsigned index, double concentration)
{
    mReservedCellConcentration[index] = concentration;
}

void AbstractTransportOdeForCoupledPdeSystem::UpdateReservedCellConcentrationByName(std::string name, double concentration)
{
    unsigned index = mpTransportReactionSystem -> GetSystemChemistry() -> GetChemicalIndexByName(name);

    mReservedCellConcentration[index] = concentration;
}

void AbstractTransportOdeForCoupledPdeSystem::UpdateChangeCellConcentrationByIndex(unsigned index, double concentration)
{
    mChangeCellConc[index] = concentration;
}

void AbstractTransportOdeForCoupledPdeSystem::UpdateChangeCellConcentrationByName(std::string name, double concentration)
{
    unsigned index = mpTransportReactionSystem -> GetSystemChemistry() -> GetChemicalIndexByName(name);

    mChangeCellConc[index] = concentration;
}


void AbstractTransportOdeForCoupledPdeSystem::SetDeltaError(double delta_error)
{
    mDelta_error = delta_error;
}

double AbstractTransportOdeForCoupledPdeSystem::GetDeltaError()
{
    return mDelta_error;
}

void AbstractTransportOdeForCoupledPdeSystem::CheckConcentration(const std::vector<double>& rY)
{
    // if chemical concentration gets too low then nan can occur, concentration must be +ve

    for(unsigned i=0; i<mNumberOfSpecies; i++)
    {
        // due to the discrete nature occasionally rY can evaluate to <0
        // ensure rY >= 0
        if(rY[i]<mDelta_error)
        {
            const_cast<double&>(rY[i]) = 0;
        }
    }

}

void AbstractTransportOdeForCoupledPdeSystem::SetIsCheckConcentration(bool IsCheckConcentration)
{
    mIsCheckConcentration = IsCheckConcentration;
}

bool AbstractTransportOdeForCoupledPdeSystem::GetIsCheckConcentration()
{
    return mIsCheckConcentration;
}

std::vector<double> AbstractTransportOdeForCoupledPdeSystem::GetReservedCellConcentration()
{
    return mReservedCellConcentration;
}

std::vector<double> AbstractTransportOdeForCoupledPdeSystem::GetChangeCellConcentration()
{
    return mChangeCellConc;
}

double AbstractTransportOdeForCoupledPdeSystem::RetrieveReservedCellConcentrationByIndex(unsigned index)
{
    return mReservedCellConcentration[index];
}

double AbstractTransportOdeForCoupledPdeSystem::RetrieveReservedCellConcentrationByName(std::string name)
{
    unsigned index = mpTransportReactionSystem -> GetSystemChemistry() -> GetChemicalIndexByName(name);

    return mReservedCellConcentration[index];
}

double AbstractTransportOdeForCoupledPdeSystem::RetrieveChangeCellConcentrationByIndex(unsigned index)
{
    return mChangeCellConc[index];
}

double AbstractTransportOdeForCoupledPdeSystem::RetrieveChangeCellConcentrationByName(std::string name)
{
    unsigned index = mpTransportReactionSystem -> GetSystemChemistry() -> GetChemicalIndexByName(name);

    return mChangeCellConc[index];
}


template<>
void ChemicalOdeSystemInformation<AbstractTransportOdeForCoupledPdeSystem>::Initialise()
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