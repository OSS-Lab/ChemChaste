#ifndef AbstractInhomogenousChemicalOdeSystemForCoupledPdeSystem_HPP_
#define AbstractInhomogenousChemicalOdeSystemForCoupledPdeSystem_HPP_

#include "AbstractInhomogenousOdeSystemForCoupledPdeSystem.hpp"
#include "AbstractReactionSystem.hpp"
#include "ChemicalOdeSystemInformation.hpp"
#include "StateVariableRegister.hpp"


class AbstractInhomogenousChemicalOdeSystemForCoupledPdeSystem: public AbstractInhomogenousOdeSystemForCoupledPdeSystem
{
private:

    // pointer to the reaction system which determines the chemical dynamics of the cell
    AbstractReactionSystem* mpReactionSystem;

    unsigned mNumberOfSpecies;

    unsigned mNumberOfReactions;

    double mDelta_error;

    bool mIsCheckConcentration;


public:

    AbstractInhomogenousChemicalOdeSystemForCoupledPdeSystem(AbstractReactionSystem* pReactionSystem = new AbstractReactionSystem());

    virtual ~AbstractInhomogenousChemicalOdeSystemForCoupledPdeSystem()
    {
    };
    
    AbstractInhomogenousChemicalOdeSystemForCoupledPdeSystem(const AbstractInhomogenousChemicalOdeSystemForCoupledPdeSystem&);

    virtual void EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY);

    unsigned GetNumberOfSpecies();

    unsigned GetNumberOfReactions();

    AbstractReactionSystem* GetReactionSystem();

    void SetNumberOfSpecies(unsigned);

    void SetNumberOfReactions(unsigned);

    void SetReactionSystem(AbstractReactionSystem* );

    void SetDeltaError(double);

    double GetDeltaError();

    void CheckConcentration(const std::vector<double>&);

    void SetIsCheckConcentration(bool);

    bool GetIsCheckConcentration();

};

AbstractInhomogenousChemicalOdeSystemForCoupledPdeSystem::AbstractInhomogenousChemicalOdeSystemForCoupledPdeSystem(AbstractReactionSystem* pReactionSystem)
    : AbstractInhomogenousOdeSystemForCoupledPdeSystem(pReactionSystem -> GetSystemChemistry() -> GetNumberChemicals(),pReactionSystem -> GetSystemChemistry() -> GetNumberChemicals())
{

    //mpSystemInfo -> ChemicalOdeSystemInformation<AbstractChemicalOdeSystem>::InstanceWrapper(pReactionSystem);
    SetReactionSystem(pReactionSystem);
    SetNumberOfSpecies(pReactionSystem -> GetSystemChemistry() -> GetNumberChemicals());
    SetNumberOfReactions(pReactionSystem -> GetNumberOfReactions());



    mpSystemInfo.reset(new ChemicalOdeSystemInformation<AbstractInhomogenousChemicalOdeSystemForCoupledPdeSystem>(pReactionSystem));
    //mpSystemInfo -> SetReactionSystem(pReactionSystem);
    SetStateVariables(GetInitialConditions());

    SetIsCheckConcentration(true);
    SetDeltaError(1e-6);
}

AbstractInhomogenousChemicalOdeSystemForCoupledPdeSystem::AbstractInhomogenousChemicalOdeSystemForCoupledPdeSystem(const AbstractInhomogenousChemicalOdeSystemForCoupledPdeSystem& existingSystem)
    : AbstractInhomogenousOdeSystemForCoupledPdeSystem(existingSystem.mNumberOfSpecies,existingSystem.mNumberOfSpecies)
{
    mpReactionSystem = existingSystem.mpReactionSystem;

    mNumberOfSpecies = existingSystem.mNumberOfSpecies;

    mNumberOfReactions = existingSystem.mNumberOfReactions;
}


void AbstractInhomogenousChemicalOdeSystemForCoupledPdeSystem::EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
{

    if(mIsCheckConcentration)
    {
        CheckConcentration(rY);
    }

    for(unsigned i=0; i<mNumberOfSpecies; i++)
    {
        rDY[i] = 0.0;
    }
    

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
        const_cast<double&>(rY[i]) = rDY[i];
    }

    for(unsigned i=0; i<mNumberOfSpecies; i++)
    {
        rDY[i] = 0.0;
    }
}  

unsigned AbstractInhomogenousChemicalOdeSystemForCoupledPdeSystem::GetNumberOfSpecies()
{
    return mNumberOfSpecies;
}

unsigned AbstractInhomogenousChemicalOdeSystemForCoupledPdeSystem::GetNumberOfReactions()
{
    return mNumberOfReactions;
}

AbstractReactionSystem* AbstractInhomogenousChemicalOdeSystemForCoupledPdeSystem::GetReactionSystem()
{
    return mpReactionSystem;
}

void AbstractInhomogenousChemicalOdeSystemForCoupledPdeSystem::SetNumberOfSpecies(unsigned numberOfSpecies)
{
    mNumberOfSpecies = numberOfSpecies;
}

void AbstractInhomogenousChemicalOdeSystemForCoupledPdeSystem::SetNumberOfReactions(unsigned numberOfReactions)
{
    mNumberOfReactions = numberOfReactions;
}

void AbstractInhomogenousChemicalOdeSystemForCoupledPdeSystem::SetReactionSystem(AbstractReactionSystem* p_reactionSystem)
{
    mpReactionSystem = p_reactionSystem;
}

void AbstractInhomogenousChemicalOdeSystemForCoupledPdeSystem::SetDeltaError(double delta_error)
{
    mDelta_error = delta_error;
}

double AbstractInhomogenousChemicalOdeSystemForCoupledPdeSystem::GetDeltaError()
{
    return mDelta_error;
}

void AbstractInhomogenousChemicalOdeSystemForCoupledPdeSystem::CheckConcentration(const std::vector<double>& rY)
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

void AbstractInhomogenousChemicalOdeSystemForCoupledPdeSystem::SetIsCheckConcentration(bool IsCheckConcentration)
{
    mIsCheckConcentration = IsCheckConcentration;
}

bool AbstractInhomogenousChemicalOdeSystemForCoupledPdeSystem::GetIsCheckConcentration()
{
    return mIsCheckConcentration;
}


template<>
void ChemicalOdeSystemInformation<AbstractInhomogenousChemicalOdeSystemForCoupledPdeSystem>::Initialise()
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