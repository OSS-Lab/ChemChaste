#ifndef ABSTRACTREVERSIBLEMEMBRANEREACTION_HPP_
#define ABSTRACTREVERSIBLEMEMBRANEREACTION_HPP_

// general includes
#include <string>
#include <stdlib.h> 
#include <vector>
#include <iostream>
// custom includes
#include "AbstractChemistry.hpp"
#include "AbstractChemical.hpp"

// reaction includes
#include "AbstractMembraneReaction.hpp"

// abstract property to contain information about the interactions of chemical species at a cell boundary
// define two coupled reactions; first the reaction in the bulk, second the reaction on the cell side
// bulk -> bulk | cell -> cell


class AbstractReversibleMembraneReaction : public AbstractMembraneReaction
{
protected:
    double mReverseReactionRate;
    // forward reaction rate is stored in the upcast as double mReactionRate;

    std::string mReversibleDelimiter = "<->";

    std::string mReversibleName = "kr =";

private:
    using AbstractMembraneReaction::React;
    using AbstractMembraneReaction::UpdateReaction;
    using AbstractMembraneReaction::UpdateReactionRate;
    using AbstractMembraneReaction::GetReactionType;
    using AbstractMembraneReaction::ParseReactionInformation;

public:

    // constructor
    AbstractReversibleMembraneReaction(   
                                std::vector<AbstractChemical*> bulkSubstrates = std::vector<AbstractChemical*>(),
                                std::vector<AbstractChemical*> bulkProducts = std::vector<AbstractChemical*>(),
                                std::vector<AbstractChemical*> cellSubstrates = std::vector<AbstractChemical*>(),
                                std::vector<AbstractChemical*> cellProducts = std::vector<AbstractChemical*>(),
                                std::vector<unsigned> stoichBulkSubstrates = std::vector<unsigned>(),
                                std::vector<unsigned> stoichBulkProducts = std::vector<unsigned>(),
                                std::vector<unsigned> stoichCellSubstrates = std::vector<unsigned>(),
                                std::vector<unsigned> stoichCellProducts = std::vector<unsigned>(),
                                double reactionRate = 1.0,
                                double reverseReactionRate =0.0
    );

    // copy constructor
    AbstractReversibleMembraneReaction(const AbstractReversibleMembraneReaction&);

    // destructor
    virtual ~AbstractReversibleMembraneReaction()
    {
    };

    virtual double GetForwardReactionRate();

    virtual double GetReverseReactionRate();

    virtual void SetForwardReactionRate(double);

    virtual void SetReverseReactionRate(double);

    virtual void React(AbstractChemistry*, AbstractChemistry*, const std::vector<double>&, const std::vector<double>&, std::vector<double>&, std::vector<double>&);

    virtual void UpdateReaction();

    virtual void UpdateReactionRate(AbstractChemistry*, AbstractChemistry*, const std::vector<double>&, const std::vector<double>&);

    virtual std::string GetReactionType();

    virtual void ParseReactionInformation(std::string, bool);

    void SetReversibleDelimiter(std::string);
    
    std::string GetReversibleDelimiter();

    void SetReversibleRateName(std::string);

    std::string GetReversibleRateName();

};

// constructor
AbstractReversibleMembraneReaction::AbstractReversibleMembraneReaction(   
                            std::vector<AbstractChemical*> bulkSubstrates,
                            std::vector<AbstractChemical*> bulkProducts,
                            std::vector<AbstractChemical*> cellSubstrates,
                            std::vector<AbstractChemical*> cellProducts,
                            std::vector<unsigned> stoichBulkSubstrates,
                            std::vector<unsigned> stoichBulkProducts,
                            std::vector<unsigned> stoichCellSubstrates,
                            std::vector<unsigned> stoichCellProducts,
                            double reactionRate,
                            double reverseReactionRate)
        :   AbstractMembraneReaction(bulkSubstrates,bulkProducts,cellSubstrates,cellProducts,stoichBulkSubstrates,stoichBulkProducts,stoichCellSubstrates,stoichCellProducts,reactionRate),
            mReverseReactionRate(reverseReactionRate)
{
}


// copy constructor
AbstractReversibleMembraneReaction::AbstractReversibleMembraneReaction(const AbstractReversibleMembraneReaction& existingReaction)
    : AbstractMembraneReaction(existingReaction)
{
    mReverseReactionRate = existingReaction.mReverseReactionRate;
}

double AbstractReversibleMembraneReaction::GetForwardReactionRate()
{
    return mReactionRate;
}

double AbstractReversibleMembraneReaction::GetReverseReactionRate()
{
    return mReverseReactionRate;
}

void AbstractReversibleMembraneReaction::SetForwardReactionRate(double reactionRate)
{
    SetReactionRate(reactionRate);
}

void AbstractReversibleMembraneReaction::SetReverseReactionRate(double reactionRate)
{
    if(mIsRateCheck)
    {
        reactionRate = CheckRate(reactionRate);
    }
    mReverseReactionRate = reactionRate;
}

void AbstractReversibleMembraneReaction::React(AbstractChemistry* bulkChemistry, AbstractChemistry* cellChemistry, const std::vector<double>& currentBulkConcentration, const std::vector<double>& currentCellConcentration, std::vector<double>& changeBulkConc, std::vector<double>& changeCellConc)
{
    std::vector<AbstractChemical*> p_bulk_chemical_vector = bulkChemistry -> rGetChemicalVector();
    
    std::vector<AbstractChemical*> p_cell_chemical_vector = cellChemistry -> rGetChemicalVector();
    
    UpdateReactionRate(bulkChemistry, cellChemistry, currentBulkConcentration, currentCellConcentration);
    // perform the reaction
    // run through the bulk species
    unsigned index=0;
    for(std::vector<AbstractChemical*>::iterator chem_iter = p_bulk_chemical_vector.begin();
            chem_iter != p_bulk_chemical_vector.end();
            ++chem_iter, ++index)
    {
        
        AbstractChemical *p_system_chemical = dynamic_cast<AbstractChemical*>(*chem_iter);

        // for each bulk chemical, parse whether it is involved in this reaction.
        for(unsigned j=0; j<mNumberOfBulkSubstrates; j++)
        {
            if(mpBulkSubstrates[j] -> GetChemicalName()==p_system_chemical -> GetChemicalName())
            {
                changeBulkConc[index] -= mStoichBulkSubstrates[j]*GetForwardReactionRate();
                changeBulkConc[index] += mStoichBulkSubstrates[j]*GetReverseReactionRate();
                break;
            }
        }
        // a reactant may be present on both sides of the reaction and may convert at different functional rates
        for(unsigned j=0; j<mNumberOfBulkProducts; j++)
        {
            if(mpBulkProducts[j] -> GetChemicalName()==p_system_chemical -> GetChemicalName())
            {
                changeBulkConc[index] += mStoichBulkProducts[j]*GetForwardReactionRate();
                changeBulkConc[index] -= mStoichBulkProducts[j]*GetReverseReactionRate();
                break;
            }
        }
    }
    // run through the cell species
    index=0;
    for(std::vector<AbstractChemical*>::iterator chem_iter = p_cell_chemical_vector.begin();
            chem_iter != p_cell_chemical_vector.end();
            ++chem_iter, ++index)
    {
        AbstractChemical *p_system_chemical = dynamic_cast<AbstractChemical*>(*chem_iter);

        // for each bulk chemical, parse whether it is involved in this reaction.
        for(unsigned j=0; j<mNumberOfCellSubstrates; j++)
        {
            if(mpCellSubstrates[j] -> GetChemicalName()==p_system_chemical -> GetChemicalName())
            {
                changeCellConc[index] -= mStoichCellSubstrates[j]*GetForwardReactionRate();
                changeCellConc[index] += mStoichCellSubstrates[j]*GetReverseReactionRate();
                break;
            }
        }
        // a reactant may be present on both sides of the reaction and may convert at different functional rates
        for(unsigned j=0; j<mNumberOfCellProducts; j++)
        {
            if(mpCellProducts[j] -> GetChemicalName()==p_system_chemical -> GetChemicalName())
            {
                changeCellConc[index] += mStoichCellProducts[j]*GetForwardReactionRate();
                changeCellConc[index] -= mStoichCellProducts[j]*GetReverseReactionRate();
                break;
            }
        }
    }
};


void AbstractReversibleMembraneReaction::UpdateReaction()
{
    return;
}

void AbstractReversibleMembraneReaction::UpdateReactionRate(AbstractChemistry* bulkChemistry, AbstractChemistry* cellChemistry, const std::vector<double>& currentBulkConc, const std::vector<double>& currentCellConc)
{
    // default to zeroth order reaction
    SetForwardReactionRate(mReactionRate);
    SetReverseReactionRate(mReverseReactionRate);
}

std::string AbstractReversibleMembraneReaction::GetReactionType()
{
    // virtual function to be overriden in derived classes, used to idnetify reaction inheritance 
    return "ZerothOrderReversibleMembrane";
}

void AbstractReversibleMembraneReaction::ParseReactionInformation(std::string reaction_information, bool IsReversible = true)
{
     

    if(!IsReversible)
    {
        AbstractMembraneReaction::ParseReactionInformation(reaction_information, false);
        SetReverseReactionRate(0.0);
    }
    else
    {
        
        bool IsForward = reaction_information.find(mIrreversibleRateName);
        bool IsReverse = reaction_information.find(mReversibleName);
        size_t posForward;
        size_t posReverse;
        unsigned length_string_forward =0;
        unsigned length_string_reverse =0;

        if(IsForward)
        {
            posForward = reaction_information.find(mIrreversibleRateName);

            length_string_forward = reaction_information.substr(posForward,std::string::npos).length();
        }
        if(IsReverse)
        {
            posReverse = reaction_information.find(mReversibleName);

            length_string_reverse = reaction_information.substr(posReverse,std::string::npos).length();
        }

        if( length_string_forward<length_string_reverse )
        {
            SetForwardReactionRate(atof(reaction_information.substr(posForward+mIrreversibleRateName.size()+1,std::string::npos).c_str()));
            reaction_information.erase(posForward,std::string::npos);
            SetReverseReactionRate(atof(reaction_information.substr(posReverse+mReversibleName.size()+1,std::string::npos).c_str()));
 
        }
        else
        {
            SetReverseReactionRate(atof(reaction_information.substr(posReverse+mReversibleName.size()+1,std::string::npos).c_str()));
            reaction_information.erase(posReverse,std::string::npos);
            SetForwardReactionRate(atof(reaction_information.substr(posForward+mIrreversibleRateName.size()+1,std::string::npos).c_str()));
        }
    }
}

// file read functions
void AbstractReversibleMembraneReaction::SetReversibleDelimiter(std::string delim)
{
    mReversibleDelimiter = delim;
}

std::string AbstractReversibleMembraneReaction::GetReversibleDelimiter()
{
    return mReversibleDelimiter;
}

void AbstractReversibleMembraneReaction::SetReversibleRateName(std::string rateName)
{
    mReversibleDelimiter = rateName;
}

std::string AbstractReversibleMembraneReaction::GetReversibleRateName()
{
    return mReversibleDelimiter;
}

#endif