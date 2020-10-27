#ifndef AbstractChemicalOdeForCoupledPdeSystem_HPP_
#define AbstractChemicalOdeForCoupledPdeSystem_HPP_

#include "AbstractOdeSystemForCoupledPdeSystem.hpp"
#include "AbstractReactionSystem.hpp"
#include "ChemicalOdeSystemInformation.hpp"

// chemcial ode system to be coupled to a pde system, coantins the emthods to evaluate a chemical reaction system 
// also checks the output of the chemical reaction system to ensure non-negative concentrations (states) are removed

// not used in the reaction-diffusion system as the chemcial ode is not required to be solved

class AbstractChemicalOdeForCoupledPdeSystem: public AbstractOdeSystemForCoupledPdeSystem
{
private:

    // reaction system object which contains the data structures and functions regarding the reactions
    AbstractReactionSystem* mpReactionSystem;

    // number of chemical species (state varibales) in the pde system 
    unsigned mNumberOfSpecies;

    // numebr of the reacitons in the reaction system which need to be iterated through to fully react the system
    unsigned mNumberOfReactions;
    
    // error tolerence for checking whether the concentrations are >0, i.e =0 if <tolerence
    // default 1e-6
    double mDelta_error; 
    
    // whether to check for negative concentrations
    // default true
    bool mIsCheckConcentration;


public:

    AbstractChemicalOdeForCoupledPdeSystem(AbstractReactionSystem* pReactionSystem = new AbstractReactionSystem());

    virtual ~AbstractChemicalOdeForCoupledPdeSystem()
    {
    };
    
    AbstractChemicalOdeForCoupledPdeSystem(const AbstractChemicalOdeForCoupledPdeSystem&);

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

AbstractChemicalOdeForCoupledPdeSystem::AbstractChemicalOdeForCoupledPdeSystem(AbstractReactionSystem* pReactionSystem)
        :   AbstractOdeSystemForCoupledPdeSystem(pReactionSystem -> GetSystemChemistry() -> GetNumberChemicals(),pReactionSystem -> GetSystemChemistry() -> GetNumberChemicals())            
{
    // set the reaction insformation
    SetReactionSystem(pReactionSystem);
    SetNumberOfSpecies(pReactionSystem -> GetSystemChemistry() -> GetNumberChemicals());
    SetNumberOfReactions(pReactionSystem -> GetNumberOfReactions());

    // fulfil the Chaste Ode system requirement of initialising the ode system information
    mpSystemInfo.reset(new ChemicalOdeSystemInformation<AbstractChemicalOdeForCoupledPdeSystem>(pReactionSystem));

    SetStateVariables(GetInitialConditions());

    // setup the data needed for checking the concentration checking procedures
    SetIsCheckConcentration(true);
    SetDeltaError(1e-6);
}

AbstractChemicalOdeForCoupledPdeSystem::AbstractChemicalOdeForCoupledPdeSystem(const AbstractChemicalOdeForCoupledPdeSystem& existingSystem)
    : AbstractOdeSystemForCoupledPdeSystem(existingSystem.mNumberOfSpecies,existingSystem.mNumberOfSpecies)
{
    mpReactionSystem = existingSystem.mpReactionSystem;

    mNumberOfSpecies = existingSystem.mNumberOfSpecies;

    mNumberOfReactions = existingSystem.mNumberOfReactions;
}

// virtual methods

void AbstractChemicalOdeForCoupledPdeSystem::EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
{
    std::vector<double>& sol_pde  = this->rGetPdeSolution();

    for(unsigned i=0; i<mNumberOfSpecies; i++)
    {
        rDY[i] = 0.0;
    }
    
    // perform the reactions in the system updating rDY with the resultant change in concentrations
    mpReactionSystem -> ReactSystem(rY, rDY);

    for(unsigned i=0; i<mNumberOfSpecies; i++)
    {
        // due to the discrete nature occasionally rY can evaluate to <0
        // ensure rY >= 0
        if(rDY[i]<mDelta_error && rDY[i]>mDelta_error)
        {
            rDY[i] = 0;
        }
    }

    for(unsigned i=0; i<mNumberOfSpecies; i++)
    {
        const_cast<double&>(rY[i]) = rDY[i]; // this can't be correct?
    }

    for(unsigned i=0; i<mNumberOfSpecies; i++)
    {
        rDY[i] = 0.0;
    }

}  

// concrete methods

void AbstractChemicalOdeForCoupledPdeSystem::CheckConcentration(const std::vector<double>& rY)
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

void AbstractChemicalOdeForCoupledPdeSystem::SetNumberOfSpecies(unsigned numberOfSpecies)
{
    mNumberOfSpecies = numberOfSpecies;
}

void AbstractChemicalOdeForCoupledPdeSystem::SetNumberOfReactions(unsigned numberOfReactions)
{
    mNumberOfReactions = numberOfReactions;
}

void AbstractChemicalOdeForCoupledPdeSystem::SetReactionSystem(AbstractReactionSystem* p_reactionSystem)
{
    mpReactionSystem = p_reactionSystem;
}

void AbstractChemicalOdeForCoupledPdeSystem::SetDeltaError(double delta_error)
{
    mDelta_error = delta_error;
}

void AbstractChemicalOdeForCoupledPdeSystem::SetIsCheckConcentration(bool IsCheckConcentration)
{
    mIsCheckConcentration = IsCheckConcentration;
}

// get methods

unsigned AbstractChemicalOdeForCoupledPdeSystem::GetNumberOfSpecies()
{
    return mNumberOfSpecies;
}

unsigned AbstractChemicalOdeForCoupledPdeSystem::GetNumberOfReactions()
{
    return mNumberOfReactions;
}

AbstractReactionSystem* AbstractChemicalOdeForCoupledPdeSystem::GetReactionSystem()
{
    return mpReactionSystem;
}

bool AbstractChemicalOdeForCoupledPdeSystem::GetIsCheckConcentration()
{
    return mIsCheckConcentration;
}

double AbstractChemicalOdeForCoupledPdeSystem::GetDeltaError()
{
    return mDelta_error;
}

template<>
void ChemicalOdeSystemInformation<AbstractChemicalOdeForCoupledPdeSystem>::Initialise()
{
    // initilaise the ode information required by Chaste for all ode systems
    // derive the states from the chemical vectors of the system chemistry
    for( unsigned i=0; i<mp_reaction_system -> GetSystemChemistry() -> GetNumberChemicals(); i++)
    {
        this->mVariableNames.push_back(mp_reaction_system -> GetSystemChemistry() -> GetChemicalNamesByIndex(i));
        this->mVariableUnits.push_back(mp_reaction_system -> GetSystemChemistry() -> GetChemicalDimensionsByIndex(i));
        this->mInitialConditions.push_back(1.0); // will be overridden
    }
    this->mInitialised = true;
}

#endif 