#ifndef ABSTRACTREVERSIBLEREACTION_HPP_
#define ABSTRACTREVERSIBLEREACTION_HPP_

// general includes
#include <string>
#include <tuple>
#include <vector>

// reaction includes
#include "AbstractReaction.hpp"

// reaction class formed from two individual AbstractReactions, handles running both reactions
// in their coupled sense.

class AbstractReversibleReaction : public AbstractReaction
{
protected:
    double mReverseReactionRate;
    // forward reaction rate is stored in the upcast as double mReactionRate;

    std::string mReversibleDelimiter = "<->";

    std::string mReversibleName = "kr =";

private:
    using AbstractReaction::React;
    using AbstractReaction::UpdateReaction;
    using AbstractReaction::UpdateReactionRate;
    using AbstractReaction::GetReactionType;
    using AbstractReaction::ParseReactionInformation;

public:

    // constructor
    AbstractReversibleReaction( std::vector<AbstractChemical*> substrates = std::vector<AbstractChemical*>(),
                                std::vector<AbstractChemical*> products = std::vector<AbstractChemical*>(),
                                std::vector<unsigned> stoichSubstrates = std::vector<unsigned>(),
                                std::vector<unsigned> stoichProducts = std::vector<unsigned>(),
                                double reactionRate = 1.0,
                                double reverseReactionRate = 0.0);

    // copy constructor
    AbstractReversibleReaction(const AbstractReversibleReaction&);

    virtual ~AbstractReversibleReaction()
    {
    };

    virtual double GetForwardReactionRate();

    virtual double GetReverseReactionRate();

    virtual void SetForwardReactionRate(double);

    virtual void SetReverseReactionRate(double);

    virtual void React(AbstractChemistry*, const std::vector<double>&, std::vector<double>&);

    virtual void UpdateReaction();

    virtual void UpdateReactionRate(AbstractChemistry*, const std::vector<double>&);

    virtual std::string GetReactionType();

    virtual void ParseReactionInformation(std::string, bool);

    void SetReversibleDelimiter(std::string);
    
    std::string GetReversibleDelimiter();

    void SetReversibleRateName(std::string);

    std::string GetReversibleRateName();

};


AbstractReversibleReaction::AbstractReversibleReaction(std::vector<AbstractChemical*> substrates,
                                    std::vector<AbstractChemical*> products,
                                    std::vector<unsigned> stoichSubstrates,
                                    std::vector<unsigned> stoichProducts,
                                    double reactionRate,
                                    double reverseReactionRate)
        :   AbstractReaction(substrates,products,stoichSubstrates,stoichProducts, reactionRate),
            mReverseReactionRate(reverseReactionRate)
{
}

// copy constructor
AbstractReversibleReaction::AbstractReversibleReaction(const AbstractReversibleReaction& existingReaction)
    : AbstractReaction(existingReaction) // should reduce the explicitly copied elements?
{
    mpSubstrates = existingReaction.mpSubstrates;
    mpProducts = existingReaction.mpProducts;
    mStoichSubstrates = existingReaction.mStoichSubstrates;
    mStoichProducts = existingReaction.mStoichProducts;
    mReactionRate = existingReaction.mReactionRate;
    mReverseReactionRate = existingReaction.mReverseReactionRate;
    mNumberOfProducts = existingReaction.mNumberOfProducts;
    mNumberOfSubstrates = existingReaction.mNumberOfSubstrates;
}


double AbstractReversibleReaction::GetForwardReactionRate()
{
    return mReactionRate;
}

double AbstractReversibleReaction::GetReverseReactionRate()
{
    return mReverseReactionRate;
}

void AbstractReversibleReaction::SetForwardReactionRate(double reactionRate)
{
    SetReactionRate(reactionRate);
}

void AbstractReversibleReaction::SetReverseReactionRate(double reactionRate)
{
    if(mIsRateCheck)
    {
        reactionRate = CheckRate(reactionRate);
    }
    mReverseReactionRate = reactionRate;
}

void AbstractReversibleReaction::React(AbstractChemistry* systemChemistry, const std::vector<double>& currentChemistryConc, std::vector<double>& changeChemistryConc)
{
    std::vector<AbstractChemical*> p_chemical_vector = systemChemistry -> rGetChemicalVector();
    
    UpdateReactionRate(systemChemistry, currentChemistryConc);
    
    // perform the reaction
    unsigned index = 0;
    for(std::vector<AbstractChemical*>::iterator chem_iter = p_chemical_vector.begin();
            chem_iter != p_chemical_vector.end();
            ++chem_iter, ++index)
    {
        AbstractChemical *p_system_chemical = dynamic_cast<AbstractChemical*>(*chem_iter);

        // for each system chemical, parse whether it is involved in this reaction.
        for(unsigned j=0; j<mNumberOfSubstrates; j++)
        {
            if(mpSubstrates[j] -> GetChemicalName()==p_system_chemical -> GetChemicalName())
            {
                changeChemistryConc[index] -= mStoichSubstrates[j]*GetForwardReactionRate();
                changeChemistryConc[index] += mStoichSubstrates[j]*GetReverseReactionRate();
                break;
            }
        }
        // a reactant may be present on both sides of the reaction and may convert at different functional rates
        for(unsigned j=0; j<mNumberOfProducts; j++)
        {
            if(mpProducts[j] -> GetChemicalName()==p_system_chemical -> GetChemicalName())
            {
                changeChemistryConc[index] += mStoichProducts[j]*GetForwardReactionRate();
                changeChemistryConc[index] -= mStoichProducts[j]*GetReverseReactionRate();
                break;
            }
        }
    }
};


void AbstractReversibleReaction::UpdateReaction()
{
    return;
}

void AbstractReversibleReaction::UpdateReactionRate(AbstractChemistry* systemChemistry, const std::vector<double>& currentChemistryConc)
{
    // default to zeroth order reaction

    SetForwardReactionRate(mReactionRate);
    SetReverseReactionRate(mReverseReactionRate);
}

std::string AbstractReversibleReaction::GetReactionType()
{
    // virtual function to be overriden in derived classes, used to idnetify reaction inheritance 
    return "ZerothOrderReversibleReaction";
}

void AbstractReversibleReaction::ParseReactionInformation(std::string reaction_information, bool IsReversible = true)
{
     

    if(!IsReversible)
    {
        AbstractReaction::ParseReactionInformation(reaction_information, false);
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
void AbstractReversibleReaction::SetReversibleDelimiter(std::string delim)
{
    mReversibleDelimiter = delim;
}

std::string AbstractReversibleReaction::GetReversibleDelimiter()
{
    return mReversibleDelimiter;
}

void AbstractReversibleReaction::SetReversibleRateName(std::string rateName)
{
    mReversibleDelimiter = rateName;
}

std::string AbstractReversibleReaction::GetReversibleRateName()
{
    return mReversibleDelimiter;
}

#endif