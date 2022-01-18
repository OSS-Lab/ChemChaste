#ifndef ABSTRACTREACTIONSYSTEM_HPP_
#define ABSTRACTREACTIONSYSTEM_HPP_

// general includes
#include <string>
#include <tuple>
#include <vector>
#include <iostream>

// reaction includes
#include "AbstractReaction.hpp"
#include "AbstractChemistry.hpp"
#include "Cell.hpp"

// class to hold the information and data structures defining a system of chemical reactions. 

class AbstractReactionSystem
{
protected:

    AbstractChemistry* mpSystemChemistry;

    std::vector<AbstractReaction*> mpReactionVector;

    unsigned mNumberOfReactions;

    CellPtr mpCell;

public:
    AbstractReactionSystem( AbstractChemistry* systemChemistry = new AbstractChemistry(), 
                            std::vector<AbstractReaction*> reactionVector = std::vector<AbstractReaction*>());

    virtual ~AbstractReactionSystem()
    {
    };

    AbstractReactionSystem(const AbstractReactionSystem&);

    virtual void ReactSystem(const std::vector<double>&, std::vector<double>&);

    virtual void UpdateReactionSystem(const std::vector<double>&);

    void DistributeCellPtr();

    std::vector<AbstractReaction*> GetReactionVector();

    AbstractReaction* GetReactionByIndex(unsigned);

    void SetReactionVector(std::vector<AbstractReaction*>);

    void SetReactionByIndex(AbstractReaction*, unsigned);

    AbstractChemistry* GetSystemChemistry();

    void SetSystemChemistry(AbstractChemistry*);

    unsigned GetNumberOfReactions();

    void SetNumberOfReactions(unsigned);

    CellPtr GetCell();

    void SetCell(CellPtr);


};

AbstractReactionSystem::AbstractReactionSystem(AbstractChemistry* systemChemistry, 
                            std::vector<AbstractReaction*> reactionVector)
    :   mpSystemChemistry(systemChemistry),
        mpReactionVector(reactionVector),
        mNumberOfReactions(reactionVector.size())
{
}

AbstractReactionSystem::AbstractReactionSystem(const AbstractReactionSystem& existingReactionSystem)
{
    mpReactionVector = existingReactionSystem.mpReactionVector;
    mNumberOfReactions = existingReactionSystem.mNumberOfReactions;
    mpSystemChemistry = existingReactionSystem.mpSystemChemistry;
}

void AbstractReactionSystem::ReactSystem(const std::vector<double>& currentChemistryConc, std::vector<double>& changeChemistryConc)
{
    unsigned number_of_species = currentChemistryConc.size();

    // update the reaction system if any variables depend on the current system concentrations
    UpdateReactionSystem(currentChemistryConc);
    
    std::vector<double> deltaConcentration(number_of_species, 0.0);
    for(std::vector<AbstractReaction*>::iterator reaction_iter = mpReactionVector.begin();
            reaction_iter != mpReactionVector.end();
            ++reaction_iter)
    {
        // iterate through the reactions in the system, modify a temporary change of concentration
        // update the system change in concentration
        
        AbstractReaction *p_system_reaction = dynamic_cast<AbstractReaction*>(*reaction_iter);
        
        p_system_reaction -> GiveCell(mpCell); // some reactions may need the cell to extract properties 
        p_system_reaction -> React(mpSystemChemistry, currentChemistryConc, deltaConcentration);
       
        // update the change in concentration
        
        for(unsigned i=0; i<number_of_species; i++)
        {
            changeChemistryConc[i] += deltaConcentration[i];

            deltaConcentration[i] =0;
        }
    }
}

void AbstractReactionSystem::UpdateReactionSystem(const std::vector<double>& currentChemistryConc)
{
    // virtual
}

void AbstractReactionSystem::DistributeCellPtr()
{
    assert(mpCell != nullptr);
    // distribute the cell ptr to the individual reactions
    for(std::vector<AbstractReaction*>::iterator reaction_iter = mpReactionVector.begin();
            reaction_iter != mpReactionVector.end();
            ++reaction_iter)
    {
        AbstractReaction *p_system_reaction = dynamic_cast<AbstractReaction*>(*reaction_iter);
        //std::cout<<"AbstractReactionSystem::DistributeCellPtr()"<<std::endl;
        p_system_reaction -> GiveCell(mpCell); // some reactions may need the cell to extract properties 

    }

}

std::vector<AbstractReaction*> AbstractReactionSystem::GetReactionVector()
{
    return mpReactionVector;
}

AbstractReaction* AbstractReactionSystem::GetReactionByIndex(unsigned index)
{
    if(index < mNumberOfReactions)
    {
        return mpReactionVector[index];
    }else{
        return new AbstractReaction();
    }
}

void AbstractReactionSystem::SetReactionVector(std::vector<AbstractReaction*> reactionVector)
{
    mpReactionVector = reactionVector;
}

void AbstractReactionSystem::SetReactionByIndex(AbstractReaction* reaction, unsigned index)
{
    if(index < mNumberOfReactions)
    {
        mpReactionVector[index] = reaction;
    }
}

AbstractChemistry* AbstractReactionSystem::GetSystemChemistry()
{
    return mpSystemChemistry;
}

void AbstractReactionSystem::SetSystemChemistry(AbstractChemistry* systemChemistry)
{
    mpSystemChemistry = systemChemistry;
}

unsigned AbstractReactionSystem::GetNumberOfReactions()
{
    return mNumberOfReactions;
}

void AbstractReactionSystem::SetNumberOfReactions(unsigned numberOfReactions)
{
    mNumberOfReactions = numberOfReactions;
}

void AbstractReactionSystem::SetCell(CellPtr pCell)
{
    mpCell = pCell;
}

CellPtr AbstractReactionSystem::GetCell()
{
    assert(mpCell != nullptr);
    return mpCell;
}


#endif