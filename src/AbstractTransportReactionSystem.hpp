#ifndef ABSTRACTTRANSPORTREACTIONSYSTEM_HPP_
#define ABSTRACTTRANSPORTREACTIONSYSTEM_HPP_

// general includes
#include <string>
#include <tuple>
#include <vector>
#include <iostream>

// reaction includes
#include "AbstractTransportReaction.hpp"
#include "AbstractChemistry.hpp"
#include "Cell.hpp"

// class to contain and handle a system of transport reactions. Transport reactions are defined across 
// a membrane where either side have their own concentration vectors

class AbstractTransportReactionSystem
{
protected:

    AbstractChemistry* mpBulkChemistry;

    AbstractChemistry* mpCellChemistry;

    std::vector<AbstractTransportReaction*> mpReactionVector;

    unsigned mNumberOfReactions;

    unsigned mNumberOfBulkStates=0;

    unsigned mNumberOfCellStates=0;

    CellPtr mpCell;

public:
    AbstractTransportReactionSystem(AbstractChemistry* bulkChemistry = new AbstractChemistry(), 
                                    AbstractChemistry* cellChemistry = new AbstractChemistry(), 
                                    std::vector<AbstractTransportReaction*> reactionVector = std::vector<AbstractTransportReaction*>());

    virtual ~AbstractTransportReactionSystem()
    {
    };

    AbstractTransportReactionSystem(const AbstractTransportReactionSystem&);

    // virtual methods

    virtual void ReactSystem(const std::vector<double>&, const std::vector<double>&, std::vector<double>&, std::vector<double>&);

    virtual void UpdateReactionSystem(const std::vector<double>&, const std::vector<double>&);

    void DistributeCellPtr();

    // set methods

    void SetReactionVector(std::vector<AbstractTransportReaction*>);

    void SetReactionByIndex(AbstractTransportReaction*, unsigned);

    void SetBulkChemistry(AbstractChemistry*);

    void SetCellChemistry(AbstractChemistry*);

    void SetNumberOfReactions(unsigned);


    // get methods 

    std::vector<AbstractTransportReaction*> GetReactionVector();

    AbstractTransportReaction* GetReactionByIndex(unsigned);

    AbstractChemistry* GetBulkChemistry();

    AbstractChemistry* GetCellChemistry();

    unsigned GetNumberOfReactions();

    unsigned GetNumberOfBulkStates();

    unsigned GetNumberOfCellStates();

    CellPtr GetCell();

    void SetCell(CellPtr);

};

AbstractTransportReactionSystem::AbstractTransportReactionSystem(   AbstractChemistry* bulkChemistry,
                                                                    AbstractChemistry* cellChemistry, 
                                                                    std::vector<AbstractTransportReaction*> reactionVector)
    :   mpBulkChemistry(bulkChemistry),
        mpCellChemistry(cellChemistry),
        mpReactionVector(reactionVector),
        mNumberOfReactions(reactionVector.size()),
        mNumberOfBulkStates(bulkChemistry -> GetNumberChemicals()),
        mNumberOfCellStates(cellChemistry -> GetNumberChemicals())
{
}

AbstractTransportReactionSystem::AbstractTransportReactionSystem(const AbstractTransportReactionSystem& existingReactionSystem)
{
    mpBulkChemistry = existingReactionSystem.mpBulkChemistry;
    mpCellChemistry = existingReactionSystem.mpCellChemistry;
    mpReactionVector = existingReactionSystem.mpReactionVector;
    mNumberOfReactions = existingReactionSystem.mNumberOfReactions;
    mNumberOfCellStates = existingReactionSystem.mNumberOfCellStates;
    mNumberOfBulkStates = existingReactionSystem.mNumberOfBulkStates;
}


// virtual methods

void AbstractTransportReactionSystem::ReactSystem(const std::vector<double>& currentBulkConcentration, const std::vector<double>& currentCellConcentration, std::vector<double>& changeBulkConc,std::vector<double>& changeCellConc)
{
    // update the reaction system if any variables depend on the current system concentrations
    UpdateReactionSystem(currentBulkConcentration,currentCellConcentration);

    std::vector<double> deltaBulkConcentration(mNumberOfBulkStates, 0.0);
    std::vector<double> deltaCellConcentration(mNumberOfCellStates, 0.0);

    for(std::vector<AbstractTransportReaction*>::iterator reaction_iter = mpReactionVector.begin();
            reaction_iter != mpReactionVector.end();
            ++reaction_iter)
    {
        // iterate through the reactions in the system, modify a temporary change of concentration
        // update the system change in concentration
        deltaBulkConcentration.resize(mNumberOfBulkStates, 0.0);
        deltaCellConcentration.resize(mNumberOfCellStates, 0.0);

        AbstractTransportReaction *p_system_reaction = static_cast<AbstractTransportReaction*>(*reaction_iter);

        p_system_reaction -> React(mpBulkChemistry, mpCellChemistry, currentBulkConcentration, currentCellConcentration, deltaBulkConcentration, deltaCellConcentration);

        // update the change in concentrations
        for(unsigned i=0; i<mNumberOfBulkStates; i++)
        {
            changeBulkConc.at(i) += deltaBulkConcentration.at(i);
        }

        for(unsigned i=0; i<mNumberOfCellStates; i++)
        {
            changeCellConc.at(i) += deltaCellConcentration.at(i);
        }

    }

}

void AbstractTransportReactionSystem::UpdateReactionSystem(const std::vector<double>& currentBulkConc, const std::vector<double>& currentCellConc)
{
    // virtual
    return;
}

void AbstractTransportReactionSystem::DistributeCellPtr()
{
    assert(mpCell != nullptr);
    // distribute the cell ptr to the individual reactions
    for(std::vector<AbstractTransportReaction*>::iterator reaction_iter = mpReactionVector.begin();
            reaction_iter != mpReactionVector.end();
            ++reaction_iter)
    {
        AbstractTransportReaction *p_system_reaction = static_cast<AbstractTransportReaction*>(*reaction_iter);
        
        p_system_reaction -> GiveCell(mpCell); // some reactions may need the cell to extract properties 

    }

}

// set methods

void AbstractTransportReactionSystem::SetReactionVector(std::vector<AbstractTransportReaction*> reactionVector)
{
    mpReactionVector = reactionVector;
}

void AbstractTransportReactionSystem::SetReactionByIndex(AbstractTransportReaction* reaction, unsigned index)
{
    if(index < mNumberOfReactions)
    {
        mpReactionVector[index] = reaction;
    }
}

void AbstractTransportReactionSystem::SetBulkChemistry(AbstractChemistry* bulkChemistry)
{
    mpBulkChemistry = bulkChemistry;
    mNumberOfBulkStates = bulkChemistry -> GetNumberChemicals();
}

void AbstractTransportReactionSystem::SetCellChemistry(AbstractChemistry* cellChemistry)
{
    mpCellChemistry = cellChemistry;
    mNumberOfCellStates = cellChemistry -> GetNumberChemicals();
}

void AbstractTransportReactionSystem::SetNumberOfReactions(unsigned numberOfReactions)
{
    mNumberOfReactions = numberOfReactions;
}


// get methods

std::vector<AbstractTransportReaction*> AbstractTransportReactionSystem::GetReactionVector()
{
    return mpReactionVector;
}

AbstractTransportReaction* AbstractTransportReactionSystem::GetReactionByIndex(unsigned index)
{
    if(index < mNumberOfReactions)
    {
        return mpReactionVector[index];
    }else{
        return new AbstractTransportReaction();
    }
}

AbstractChemistry* AbstractTransportReactionSystem::GetBulkChemistry()
{
    return mpBulkChemistry;
}

AbstractChemistry* AbstractTransportReactionSystem::GetCellChemistry()
{
    return mpCellChemistry;
}

unsigned AbstractTransportReactionSystem::GetNumberOfReactions()
{
    return mNumberOfReactions;
}

unsigned AbstractTransportReactionSystem::GetNumberOfBulkStates()
{
    return mNumberOfBulkStates;
}

unsigned AbstractTransportReactionSystem::GetNumberOfCellStates()
{
    return mNumberOfCellStates;
}

void AbstractTransportReactionSystem::SetCell(CellPtr pCell)
{
    mpCell = pCell;
}

CellPtr AbstractTransportReactionSystem::GetCell()
{
    assert(mpCell != nullptr);
    return mpCell;
}

#endif