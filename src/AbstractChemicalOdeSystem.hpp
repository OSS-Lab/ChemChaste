#ifndef ABSTRACTCHEMICALODESYSTEM_HPP_
#define ABSTRACTCHEMICALODESYSTEM_HPP_

#include "AbstractOdeSystem.hpp"
#include "AbstractReactionSystem.hpp"
#include "ChemicalOdeSystemInformation.hpp"

// Ode system for performing chemical reactions of templated chemical kinetic rate laws

class AbstractChemicalOdeSystem: public AbstractOdeSystem
{
private:

    // the chemical reaction system to be modelled 
    AbstractReactionSystem* mpReactionSystem;

    // number of chemical species within the chemical system
    unsigned mNumberOfSpecies;

    // number of chemical reactions within the reaction system to be solved
    unsigned mNumberOfReactions;

    // error threshold 
    // defaullt 1e-6
    double mDelta_error;

    // wether to check the chemical concentrations (ode states) for negative concentration values
    // default true
    bool mIsCheckConcentration;


public:

    AbstractChemicalOdeSystem(AbstractReactionSystem* pReactionSystem = new AbstractReactionSystem());

    virtual ~AbstractChemicalOdeSystem()
    {
    };
    
    AbstractChemicalOdeSystem(const AbstractChemicalOdeSystem&);

    // virtual methods
    virtual void EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY);

    // concrete methods
    void CheckConcentration(const std::vector<double>&);

    // set methods

    void SetNumberOfSpecies(unsigned);

    void SetNumberOfReactions(unsigned);

    void SetReactionSystem(AbstractReactionSystem* );

    void SetDeltaError(double);

    void SetIsCheckConcentration(bool);

    // get methods

    unsigned GetNumberOfSpecies();

    unsigned GetNumberOfReactions();

    AbstractReactionSystem* GetReactionSystem();

    double GetDeltaError();

    bool GetIsCheckConcentration();
};

AbstractChemicalOdeSystem::AbstractChemicalOdeSystem(AbstractReactionSystem* pReactionSystem)
        :   AbstractOdeSystem(pReactionSystem -> GetSystemChemistry() -> GetNumberChemicals())
            
{
    // set up the chemical reaciton sytem
    SetReactionSystem(pReactionSystem);
    SetNumberOfSpecies(pReactionSystem -> GetSystemChemistry() -> GetNumberChemicals());
    SetNumberOfReactions(pReactionSystem -> GetNumberOfReactions());

    // fulfil the Chaste ode information setup by initialising
    mpSystemInfo.reset(new ChemicalOdeSystemInformation<AbstractChemicalOdeSystem>(pReactionSystem));

    SetStateVariables(GetInitialConditions());

    // set the chemical concentration error checking threhsold
    SetIsCheckConcentration(true);
    SetDeltaError(1e-6);
}

AbstractChemicalOdeSystem::AbstractChemicalOdeSystem(const AbstractChemicalOdeSystem& existingSystem)
    : AbstractOdeSystem(existingSystem.mNumberOfSpecies)
{
    mpReactionSystem = existingSystem.mpReactionSystem;

    mNumberOfSpecies = existingSystem.mNumberOfSpecies;

    mNumberOfReactions = existingSystem.mNumberOfReactions;

    mDelta_error = existingSystem.mDelta_error;

    mIsCheckConcentration = existingSystem.mIsCheckConcentration;
}

// virtual methods

void AbstractChemicalOdeSystem::EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
{
    if(mIsCheckConcentration)
    {
        CheckConcentration(rY);
    }
    // reset rDY
    for(unsigned i=0; i<mNumberOfSpecies; i++)
    {
        rDY[i] = 0.0;
    }
    // perform the reaction, updating rDY with the resultant chemical change due to the reaction
    mpReactionSystem -> ReactSystem(rY, rDY);
}  

// concrete methods

void AbstractChemicalOdeSystem::CheckConcentration(const std::vector<double>& rY)
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

// set methods

void AbstractChemicalOdeSystem::SetNumberOfSpecies(unsigned numberOfSpecies)
{
    mNumberOfSpecies = numberOfSpecies;
}

void AbstractChemicalOdeSystem::SetNumberOfReactions(unsigned numberOfReactions)
{
    mNumberOfReactions = numberOfReactions;
}

void AbstractChemicalOdeSystem::SetReactionSystem(AbstractReactionSystem* p_reactionSystem)
{
    mpReactionSystem = p_reactionSystem;
}

void AbstractChemicalOdeSystem::SetDeltaError(double delta_error)
{
    mDelta_error = delta_error;
}

void AbstractChemicalOdeSystem::SetIsCheckConcentration(bool IsCheckConcentration)
{
    mIsCheckConcentration = IsCheckConcentration;
}

// get methods

unsigned AbstractChemicalOdeSystem::GetNumberOfSpecies()
{
    return mNumberOfSpecies;
}

unsigned AbstractChemicalOdeSystem::GetNumberOfReactions()
{
    return mNumberOfReactions;
}

AbstractReactionSystem* AbstractChemicalOdeSystem::GetReactionSystem()
{
    return mpReactionSystem;
}

bool AbstractChemicalOdeSystem::GetIsCheckConcentration()
{
    return mIsCheckConcentration;
}

double AbstractChemicalOdeSystem::GetDeltaError()
{
    return mDelta_error;
}


// system information template
template<>
void ChemicalOdeSystemInformation<AbstractChemicalOdeSystem>::Initialise()
{

    for( unsigned i=0; i<mp_reaction_system -> GetSystemChemistry() -> GetNumberChemicals(); i++)
    {
        this->mVariableNames.push_back(mp_reaction_system -> GetSystemChemistry() -> GetChemicalNamesByIndex(i));
        this->mVariableUnits.push_back(mp_reaction_system -> GetSystemChemistry() -> GetChemicalDimensionsByIndex(i));
        this->mInitialConditions.push_back(1.0); // will be overridden
    }
    this->mInitialised = true;
}

#endif 