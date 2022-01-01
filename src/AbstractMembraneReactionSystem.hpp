#ifndef ABSTRACTMEMBRANEREACTIONSYSTEM_HPP_
#define ABSTRACTMEMBRANEREACTIONSYSTEM_HPP_

// general includes
#include <string>
#include <tuple>
#include <vector>
#include <iostream>

// reaction includes
#include "AbstractMembraneReaction.hpp"
#include "AbstractChemistry.hpp"
#include "Cell.hpp"

// class to hold a set of membrane reactions and handle the parsing and reforming of the different chemical concentration (state variable)
// vectors. These vectors being the concentrations external and internal to the infinitesimal thickness membrane. 

class AbstractMembraneReactionSystem
{
protected:

    AbstractChemistry* mpBulkChemistry;

    AbstractChemistry* mpCellChemistry;

    std::vector<AbstractMembraneReaction*> mpReactionVector;

    unsigned mNumberOfReactions;

    unsigned mNumberOfBulkStates=0;

    unsigned mNumberOfCellStates=0;

    CellPtr mpCell;

public:
    AbstractMembraneReactionSystem( AbstractChemistry* bulkChemistry = new AbstractChemistry(), 
                                    AbstractChemistry* cellChemistry = new AbstractChemistry(), 
                                    std::vector<AbstractMembraneReaction*> reactionVector = std::vector<AbstractMembraneReaction*>());

    virtual ~AbstractMembraneReactionSystem()
    {
    };

    AbstractMembraneReactionSystem(const AbstractMembraneReactionSystem&);

    // virtual methods

    virtual void ReactSystem(const std::vector<double>&, const std::vector<double>&, std::vector<double>&, std::vector<double>&);

    virtual void UpdateReactionSystem(const std::vector<double>&, const std::vector<double>&);

    void DistributeCellPtr();

    // set methods

    void SetReactionVector(std::vector<AbstractMembraneReaction*>);

    void SetReactionByIndex(AbstractMembraneReaction*, unsigned);

    void SetBulkChemistry(AbstractChemistry*);

    void SetCellChemistry(AbstractChemistry*);

    void SetNumberOfReactions(unsigned);


    // get methods 

    std::vector<AbstractMembraneReaction*> GetReactionVector();

    AbstractMembraneReaction* GetReactionByIndex(unsigned);

    AbstractChemistry* GetBulkChemistry();

    AbstractChemistry* GetCellChemistry();

    unsigned GetNumberOfReactions();

    unsigned GetNumberOfBulkStates();

    unsigned GetNumberOfCellStates();

    CellPtr GetCell();

    void SetCell(CellPtr);

};

AbstractMembraneReactionSystem::AbstractMembraneReactionSystem( AbstractChemistry* bulkChemistry,
                                                                AbstractChemistry* cellChemistry, 
                                                                std::vector<AbstractMembraneReaction*> reactionVector)
    :   mpBulkChemistry(bulkChemistry),
        mpCellChemistry(cellChemistry),
        mpReactionVector(reactionVector),
        mNumberOfReactions(reactionVector.size()),
        mNumberOfBulkStates(bulkChemistry -> GetNumberChemicals()),
        mNumberOfCellStates(cellChemistry -> GetNumberChemicals())
{
}

AbstractMembraneReactionSystem::AbstractMembraneReactionSystem(const AbstractMembraneReactionSystem& existingReactionSystem)
{
    mpBulkChemistry = existingReactionSystem.mpBulkChemistry;
    mpCellChemistry = existingReactionSystem.mpCellChemistry;
    mpReactionVector = existingReactionSystem.mpReactionVector;
    mNumberOfReactions = existingReactionSystem.mNumberOfReactions;
    mNumberOfCellStates = existingReactionSystem.mNumberOfCellStates;
    mNumberOfBulkStates = existingReactionSystem.mNumberOfBulkStates;
}


// virtual methods

void AbstractMembraneReactionSystem::ReactSystem(const std::vector<double>& currentBulkConcentration, const std::vector<double>& currentCellConcentration, std::vector<double>& changeBulkConc,std::vector<double>& changeCellConc)
{
    // update the reaction system if any variables depend on the current system concentrations
    UpdateReactionSystem(currentBulkConcentration,currentCellConcentration);

    std::vector<double> deltaBulkConcentration(mNumberOfBulkStates, 0.0);
    std::vector<double> deltaCellConcentration(mNumberOfCellStates, 0.0);


    for(std::vector<AbstractMembraneReaction*>::iterator reaction_iter = mpReactionVector.begin();
            reaction_iter != mpReactionVector.end();
            ++reaction_iter)
    {
        // iterate through the reactions in the system, modify a temporary change of concentration
        // update the system change in concentration
        deltaBulkConcentration.resize(mNumberOfBulkStates, 0.0);
        deltaCellConcentration.resize(mNumberOfCellStates, 0.0);

        AbstractMembraneReaction *p_system_reaction = static_cast<AbstractMembraneReaction*>(*reaction_iter);

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

void AbstractMembraneReactionSystem::UpdateReactionSystem(const std::vector<double>& currentBulkConc, const std::vector<double>& currentCellConc)
{
    // virtual
    return;
}

void AbstractMembraneReactionSystem::DistributeCellPtr()
{
    assert(mpCell != nullptr);
    // distribute the cell ptr to the individual reactions
    for(std::vector<AbstractMembraneReaction*>::iterator reaction_iter = mpReactionVector.begin();
            reaction_iter != mpReactionVector.end();
            ++reaction_iter)
    {
        AbstractMembraneReaction *p_system_reaction = static_cast<AbstractMembraneReaction*>(*reaction_iter);
        
        p_system_reaction -> GiveCell(mpCell); // some reactions may need the cell to extract properties 

    }

}



// set methods

void AbstractMembraneReactionSystem::SetReactionVector(std::vector<AbstractMembraneReaction*> reactionVector)
{
    mpReactionVector = reactionVector;
}

void AbstractMembraneReactionSystem::SetReactionByIndex(AbstractMembraneReaction* reaction, unsigned index)
{
    if(index < mNumberOfReactions)
    {
        mpReactionVector[index] = reaction;
    }
}

void AbstractMembraneReactionSystem::SetBulkChemistry(AbstractChemistry* bulkChemistry)
{
    mpBulkChemistry = bulkChemistry;
    mNumberOfBulkStates = bulkChemistry -> GetNumberChemicals();
}

void AbstractMembraneReactionSystem::SetCellChemistry(AbstractChemistry* cellChemistry)
{
    mpCellChemistry = cellChemistry;
    mNumberOfCellStates = cellChemistry -> GetNumberChemicals();
}

void AbstractMembraneReactionSystem::SetNumberOfReactions(unsigned numberOfReactions)
{
    mNumberOfReactions = numberOfReactions;
}


// get methods

std::vector<AbstractMembraneReaction*> AbstractMembraneReactionSystem::GetReactionVector()
{
    return mpReactionVector;
}

AbstractMembraneReaction* AbstractMembraneReactionSystem::GetReactionByIndex(unsigned index)
{
    if(index < mNumberOfReactions)
    {
        return mpReactionVector[index];
    }else{
        return new AbstractMembraneReaction();
    }
}

AbstractChemistry* AbstractMembraneReactionSystem::GetBulkChemistry()
{
    return mpBulkChemistry;
}

AbstractChemistry* AbstractMembraneReactionSystem::GetCellChemistry()
{
    return mpCellChemistry;
}

unsigned AbstractMembraneReactionSystem::GetNumberOfReactions()
{
    return mNumberOfReactions;
}

unsigned AbstractMembraneReactionSystem::GetNumberOfBulkStates()
{
    return mNumberOfBulkStates;
}

unsigned AbstractMembraneReactionSystem::GetNumberOfCellStates()
{
    return mNumberOfCellStates;
}

void AbstractMembraneReactionSystem::SetCell(CellPtr pCell)
{
    mpCell = pCell;
}

CellPtr AbstractMembraneReactionSystem::GetCell()
{
    assert(mpCell != nullptr);
    return mpCell;
}

#endif