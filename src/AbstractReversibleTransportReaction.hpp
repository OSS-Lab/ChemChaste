#ifndef ABSTRACTREVERSIBLETRANSPORTREACTION_HPP_
#define ABSTRACTREVERSIBLETRANSPORTREACTION_HPP_

// general includes
#include <string>
#include <tuple>
#include <vector>

// reaction includes
#include "AbstractTransportReaction.hpp"

// reaction class formed from two individual AbstractTransportReaction, handles running both reactions
// in their coupled sense.

// abstract property to contain information about the interactions of chemical species at a cell boundary
// bulk <-> cell
// reaction substrates denote the species in the bulk while products denote species in the cell

// ReversibleTransportReaction

class AbstractReversibleTransportReaction : public AbstractTransportReaction
{
protected:
    double mReverseReactionRate;
    // forward reaction rate is stored in the upcast as double mReactionRate;

    std::string mReversibleDelimiter = "<->";

    std::string mReversibleName = "kr =";

private:
    using AbstractTransportReaction::React;
    using AbstractTransportReaction::UpdateReaction;
    using AbstractTransportReaction::UpdateReactionRate;
    using AbstractTransportReaction::GetReactionType;
    using AbstractTransportReaction::ParseReactionInformation;

public:

    // constructor
    AbstractReversibleTransportReaction(
                        std::vector<AbstractChemical*> bulkReactionSpecies = std::vector<AbstractChemical*>(),
                        std::vector<AbstractChemical*> cellReactionSpecies = std::vector<AbstractChemical*>(),
                        std::vector<unsigned> stoichBulk = std::vector<unsigned>(),
                        std::vector<unsigned> stoichCell = std::vector<unsigned>(),
                        double reactionRate = 1.0,
                        double reverseReactionRate = 0.0);

    // copy constructor
    AbstractReversibleTransportReaction(const AbstractReversibleTransportReaction&);

    virtual ~AbstractReversibleTransportReaction()
    {
    };

    virtual double GetForwardReactionRate();

    virtual double GetReverseReactionRate();

    virtual void SetForwardReactionRate(double);

    virtual void SetReverseReactionRate(double);

    virtual void React(AbstractChemistry*, AbstractChemistry*, const std::vector<double>&, const std::vector<double>&, std::vector<double>&, std::vector<double>&);

    virtual void UpdateReaction();

    virtual void UpdateReactionRate(AbstractChemistry*, AbstractChemistry*, const std::vector<double>&, const std::vector<double>&);

    virtual std::string GetReactionType()
    {
        // virtual function to be overriden in derived classes, used to idnetify reaction inheritance 
        return "ZerothOrderReversibleTransport";
    }


    virtual void ParseReactionInformation(std::string, bool);

    void SetReversibleDelimiter(std::string);
    
    std::string GetReversibleDelimiter();

    void SetReversibleRateName(std::string);

    std::string GetReversibleRateName();

};


AbstractReversibleTransportReaction::AbstractReversibleTransportReaction(
                                    std::vector<AbstractChemical*> bulkReactionSpecies,
                                    std::vector<AbstractChemical*> cellReactionSpecies,
                                    std::vector<unsigned> stoichBulk,
                                    std::vector<unsigned> stoichCell,
                                    double reactionRate,
                                    double reverseReactionRate)
        :   AbstractTransportReaction(bulkReactionSpecies,cellReactionSpecies,stoichBulk,stoichCell, reactionRate),
            mReverseReactionRate(reverseReactionRate)
{
}

// copy constructor
AbstractReversibleTransportReaction::AbstractReversibleTransportReaction(const AbstractReversibleTransportReaction& existingReaction)
    : AbstractTransportReaction(existingReaction)
{
    mReverseReactionRate = existingReaction.mReverseReactionRate;
}


double AbstractReversibleTransportReaction::GetForwardReactionRate()
{
    return mReactionRate;
}

double AbstractReversibleTransportReaction::GetReverseReactionRate()
{
    return mReverseReactionRate;
}

void AbstractReversibleTransportReaction::SetForwardReactionRate(double reactionRate)
{
    SetReactionRate(reactionRate);
}

void AbstractReversibleTransportReaction::SetReverseReactionRate(double reactionRate)
{
    if(mIsRateCheck)
    {
        reactionRate = CheckRate(reactionRate);
    }
    mReverseReactionRate = reactionRate;
}

void AbstractReversibleTransportReaction::React(AbstractChemistry* bulkChemistry, AbstractChemistry* cellChemistry, const std::vector<double>& currentBulkConcentration, const std::vector<double>& currentCellConcentration, std::vector<double>& changeBulkConc, std::vector<double>& changeCellConc)
{   

    std::vector<AbstractChemical*> p_bulk_chemical_vector = bulkChemistry -> rGetChemicalVector();
    std::vector<AbstractChemical*> p_cell_chemical_vector = cellChemistry -> rGetChemicalVector();

    
    UpdateReactionRate(bulkChemistry, cellChemistry, currentBulkConcentration, currentCellConcentration);
    
    // perform the reaction
    
    // bulk species
    unsigned index = 0;
    for(std::vector<AbstractChemical*>::iterator chem_iter = p_bulk_chemical_vector.begin();
            chem_iter != p_bulk_chemical_vector.end();
            ++chem_iter, ++index)
    {
        AbstractChemical *p_bulk_chemical = dynamic_cast<AbstractChemical*>(*chem_iter);

        // for each system chemical, parse whether it is involved in this reaction.
        for(unsigned j=0; j<mNumberOfBulkSpecies; j++)
        {
            if(mpBulkReactionSpecies[j] -> GetChemicalName()==p_bulk_chemical -> GetChemicalName())
            {
                changeBulkConc[index] -= mStoichBulk[j]*GetForwardReactionRate();
                changeBulkConc[index] += mStoichBulk[j]*GetReverseReactionRate();
                break;
            }
        } 
    }

    // cell species 
    index = 0;
    for(std::vector<AbstractChemical*>::iterator chem_iter = p_cell_chemical_vector.begin();
            chem_iter != p_cell_chemical_vector.end();
            ++chem_iter, ++index)
    {
        AbstractChemical *p_cell_chemical = dynamic_cast<AbstractChemical*>(*chem_iter);

        for(unsigned j=0; j<mNumberOfCellSpecies; j++)
        {
            if(mpCellReactionSpecies[j] -> GetChemicalName()==p_cell_chemical -> GetChemicalName())
            {
                changeCellConc[index] += mStoichCell[j]*GetForwardReactionRate();
                changeCellConc[index] -= mStoichCell[j]*GetReverseReactionRate();
                break;
            }
        }
    }
    
};


void AbstractReversibleTransportReaction::UpdateReaction()
{
    return;
}

void AbstractReversibleTransportReaction::UpdateReactionRate(AbstractChemistry* bulkChemistry, AbstractChemistry* cellChemistry, const std::vector<double>& currentBulkConc, const std::vector<double>& currentCellConc)
{
    // default to zeroth order reaction
    SetForwardReactionRate(mReactionRate);
    SetReverseReactionRate(mReverseReactionRate);
}

void AbstractReversibleTransportReaction::ParseReactionInformation(std::string reaction_information, bool IsReversible=true)
{
    //std::cout<<"AbstractReversibleTransportReaction::ParseReactionInformation - start"<<std::endl;

    if(!IsReversible)
    {
        AbstractTransportReaction::ParseReactionInformation(reaction_information, false);
        SetReverseReactionRate(0.0);
    }
    else
    {
        // check which rate, kf or kr, comes first. first one will have bool value 0
        bool IsForward = reaction_information.find(mIrreversibleRateName);
        bool IsReverse = reaction_information.find(mReversibleName);
     
        size_t posForward =0;
        size_t posReverse =0;
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
            // kf defined first

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
    //std::cout<<"AbstractReversibleTransportReaction::ParseReactionInformation - end"<<std::endl;
}

// file read functions
void AbstractReversibleTransportReaction::SetReversibleDelimiter(std::string delim)
{
    mReversibleDelimiter = delim;
}

std::string AbstractReversibleTransportReaction::GetReversibleDelimiter()
{
    return mReversibleDelimiter;
}

void AbstractReversibleTransportReaction::SetReversibleRateName(std::string rateName)
{
    mReversibleDelimiter = rateName;
}

std::string AbstractReversibleTransportReaction::GetReversibleRateName()
{
    return mReversibleDelimiter;
}

#endif