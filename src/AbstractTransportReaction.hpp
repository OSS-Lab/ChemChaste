#ifndef ABSTRACTTRANSPORTREACTION_HPP_
#define ABSTRACTTRANSPORTREACTION_HPP_

// general includes
#include <string>
#include <stdlib.h> 
#include <vector>
#include <iostream>
// custom includes
#include "AbstractChemistry.hpp"
#include "AbstractChemical.hpp"

// abstract property to contain information about the interactions of chemical species at a cell boundary
// bulk -> cell
// reaction substrates denote the species in the bulk while products denote species in the cell

// ZerothOrderTransportIntoCell

class AbstractTransportReaction 
{
protected:

    // vector with AbstractChemical of bulk species, in the forward reaction sense.  To be used to track system concentrations
    std::vector<AbstractChemical*> mpBulkReactionSpecies;

    // vector with the AbstractChemical cell species, in the forward reaction sense.  To be used to track system concentrations
    std::vector<AbstractChemical*> mpCellReactionSpecies;

    unsigned mNumberOfBulkSpecies;

    unsigned mNumberOfCellSpecies;


    // vector containing the stoichmetry of the substrates, bulk species.  Will be "-ve" in reaction
    std::vector<unsigned> mStoichBulk;

    // vector containing the stoichmetry of the products, cell species.  Will be "+ve" in reaction
    std::vector<unsigned> mStoichCell;

    // container for the reaction rate constant.  '+ve' value. May be overriden in derived classes to account for functional rates
    double mReactionRate;

    bool mIsRateCheck;

    double mDelta_rate_min;

    double mDelta_rate_max;

    std::string mIrreversibleDelimiter = "->";
  
    std::string mIrreversibleRateName = "kf =";

public:

    // constructor
    AbstractTransportReaction(  std::vector<AbstractChemical*> bulkReactionSpecies = std::vector<AbstractChemical*>(),
                                std::vector<AbstractChemical*> cellReactionSpecies = std::vector<AbstractChemical*>(),
                                std::vector<unsigned> stoichBulk = std::vector<unsigned>(),
                                std::vector<unsigned> stoichCell = std::vector<unsigned>(),
                                double reactionRate = 1.0
    );
    

    // copy constructor
    AbstractTransportReaction(const AbstractTransportReaction&);

    // destructor
    virtual ~AbstractTransportReaction()
    {
    };


    // virtual methods

    // function to take in pointer to current concentration state vector of the state vector for change in cocnentration.  Basic update via multiplying by constant reaction rate
    virtual void React(AbstractChemistry*, AbstractChemistry*, const std::vector<double>&, const std::vector<double>&, std::vector<double>&, std::vector<double>&);

    virtual void UpdateReaction();

    virtual void UpdateReactionRate(AbstractChemistry*, AbstractChemistry*, const std::vector<double>&, const std::vector<double>&);

    virtual void SetReactionRate(double);

    virtual double GetReactionRate();

    virtual std::string GetReactionType();

    virtual void ParseReactionInformation(std::string, bool);


    // Chemical handeling functions
    std::vector<AbstractChemical*> GetBulkSpecies();

    AbstractChemical* GetBulkSpeciesByIndex(unsigned);

    std::vector<AbstractChemical*> GetCellSpecies();

    AbstractChemical* GetCellSpeciesByIndex(unsigned);

    void SetBulkSpecies(std::vector<AbstractChemical*>);

    void SetCellSpecies(std::vector<AbstractChemical*>);

    std::vector<unsigned> GetStoichBulk();

    unsigned GetStoichBulkByIndex(unsigned);

    void SetStoichBulk(std::vector<unsigned>);

    std::vector<unsigned> GetStoichCell();

    unsigned GetStoichCellByIndex(unsigned);

    void SetStoichCell(std::vector<unsigned>);

    void SetNumberOfBulkSpecies(unsigned);

    unsigned GetNumberOfBulkSpecies();

    void SetNumberOfCellSpecies(unsigned);

    unsigned GetNumberOfCellSpecies();


    // reaction concentration checking functions
    void SetIsRateCheck(bool);

    bool GetIsRateCheck();

    void SetDeltaErrorRateMin(double);

    double GetDeltaErrorRateMin();

    void SetDeltaErrorRateMax(double);

    double GetDeltaErrorRateMax();

    double CheckRate(double);



    // file read functions
    void SetIrreversibleDelimiter(std::string);

    std::string GetIrreversibleDelimiter();

    void SetIrreversibleRateName(std::string);

    std::string GetIrreversibleRateName();

};

// constructor
AbstractTransportReaction::AbstractTransportReaction(   std::vector<AbstractChemical*> bulkReactionSpecies,
                                                        std::vector<AbstractChemical*> cellReactionSpecies,
                                                        std::vector<unsigned> stoichBulk,
                                                        std::vector<unsigned> stoichCell,
                                                        double reactionRate)
        :   mpBulkReactionSpecies(bulkReactionSpecies),
            mpCellReactionSpecies(cellReactionSpecies),
            mStoichBulk(stoichBulk),
            mStoichCell(stoichCell),
            mReactionRate(reactionRate)
{
    mNumberOfCellSpecies = cellReactionSpecies.size();
    mNumberOfBulkSpecies = bulkReactionSpecies.size();

    SetIsRateCheck(true);
    SetDeltaErrorRateMin(1e-8);
    SetDeltaErrorRateMax(1e8);
}


// copy constructor
AbstractTransportReaction::AbstractTransportReaction(const AbstractTransportReaction& existingReaction)
{
    
    mpBulkReactionSpecies = existingReaction.mpBulkReactionSpecies;
    mpCellReactionSpecies = existingReaction.mpCellReactionSpecies;
    mStoichBulk = existingReaction.mStoichBulk;
    mStoichCell = existingReaction.mStoichCell;
    mReactionRate = existingReaction.mReactionRate;
    mNumberOfCellSpecies = existingReaction.mNumberOfCellSpecies;
    mNumberOfBulkSpecies = existingReaction.mNumberOfBulkSpecies;
    mIsRateCheck = existingReaction.mIsRateCheck;
    mDelta_rate_min = existingReaction.mDelta_rate_min;
    mDelta_rate_max = existingReaction.mDelta_rate_max;
}

void AbstractTransportReaction::React(AbstractChemistry* bulkChemistry, AbstractChemistry* cellChemistry, const std::vector<double>& currentBulkConcentration, const std::vector<double>& currentCellConcentration, std::vector<double>& changeBulkConc, std::vector<double>& changeCellConc)
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
        AbstractChemical *p_bulk_chemical = static_cast<AbstractChemical*>(*chem_iter);

        // for each bulk chemical, parse whether it is involved in this reaction.
        for(unsigned j=0; j<mNumberOfBulkSpecies; j++)
        {
            if(mpBulkReactionSpecies.at(j) -> GetChemicalName()==p_bulk_chemical -> GetChemicalName())
            {
                changeBulkConc.at(index)  -= mStoichBulk.at(j)*GetReactionRate();
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
        AbstractChemical *p_cell_chemical = static_cast<AbstractChemical*>(*chem_iter);
   
        for(unsigned j=0; j<mNumberOfCellSpecies; j++)
        {
            if(mpCellReactionSpecies.at(j) -> GetChemicalName()==p_cell_chemical -> GetChemicalName())
            {
                changeCellConc.at(index) += mStoichCell.at(j) *GetReactionRate();
                break;
            }
        }

    }

}

// vitual function to update member variables and conditions of derived classes
void AbstractTransportReaction::UpdateReaction()
{
    return;
}

void AbstractTransportReaction::UpdateReactionRate(AbstractChemistry* bulkChemistry, AbstractChemistry* cellChemistry, const std::vector<double>& currentBulkConc, const std::vector<double>& currentCellConc)
{
    // default to zeroth order reaction
    SetReactionRate(mReactionRate);
}

void AbstractTransportReaction::SetReactionRate(double reactionRate)
{
    if(mIsRateCheck)
    {
        reactionRate = CheckRate(reactionRate);
    }
    
    mReactionRate = reactionRate;
}

double AbstractTransportReaction::GetReactionRate()
{
    return mReactionRate;
}

std::string AbstractTransportReaction::GetReactionType()
{
    // virtual function to be overriden in derived classes, used to idnetify reaction inheritance 
    return "ZerothOrderTransportIntoCell";
}


void AbstractTransportReaction::ParseReactionInformation(std::string reaction_information, bool IsReversible=false)
{
   
    if(reaction_information.find(mIrreversibleRateName) != std::string::npos)
    {
        size_t pos= reaction_information.find(mIrreversibleRateName);

        SetReactionRate(atof(reaction_information.substr(pos+mIrreversibleRateName.size()+1,std::string::npos).c_str()));
    }

}


// member functions
std::vector<AbstractChemical*> AbstractTransportReaction::GetBulkSpecies()
{
    return mpBulkReactionSpecies;
}

AbstractChemical* AbstractTransportReaction::GetBulkSpeciesByIndex(unsigned index)
{
    if(index < mNumberOfBulkSpecies)
    {
        return mpBulkReactionSpecies[index];
    }
    else
    {
        std::cout<<"Error: AbstractTransportReaction::GetBulkSpeciesByIndex(unsigned index), index out of bounds"<<std::endl;
        return new AbstractChemical();
    } 
}

std::vector<AbstractChemical*> AbstractTransportReaction::GetCellSpecies()
{
    return mpCellReactionSpecies;
}


AbstractChemical* AbstractTransportReaction::GetCellSpeciesByIndex(unsigned index)
{
    if(index < mNumberOfCellSpecies)
    {
        return mpCellReactionSpecies[index];
    }
    else
    {
        std::cout<<"Error: AbstractTransportReaction::GetCellSpeciesByIndex(unsigned index), index out of bounds"<<std::endl;
        return new AbstractChemical();
    } 
}

void AbstractTransportReaction::SetBulkSpecies(std::vector<AbstractChemical*> bulkReactionSpecies)
{
    mpBulkReactionSpecies = bulkReactionSpecies;
}


void AbstractTransportReaction::SetCellSpecies(std::vector<AbstractChemical*> cellReactionSpecies)
{
    mpCellReactionSpecies = cellReactionSpecies;
}


std::vector<unsigned> AbstractTransportReaction::GetStoichBulk()
{
    return mStoichBulk;
}


unsigned AbstractTransportReaction::GetStoichBulkByIndex(unsigned index)
{
    if(index < mNumberOfBulkSpecies)
    {
        return mStoichBulk[index];
    }
    else
    {
        std::cout<<"Error: AbstractTransportReaction::GetStoichBulkByIndex(unsigned index), index out of bounds"<<std::endl;
        return 0;
    } 
}


void AbstractTransportReaction::SetStoichBulk(std::vector<unsigned> stoichBulk)
{
    mStoichBulk = stoichBulk;
}


std::vector<unsigned> AbstractTransportReaction::GetStoichCell()
{
    return mStoichCell;
}


unsigned AbstractTransportReaction::GetStoichCellByIndex(unsigned index)
{
    if(index < mNumberOfCellSpecies)
    {
        return mStoichCell[index];
    }
    else
    {
        std::cout<<"Error: AbstractTransportReaction::GetStoichCellByIndex(unsigned index), index out of bounds"<<std::endl;
        return 0;
    } 
}

void AbstractTransportReaction::SetStoichCell(std::vector<unsigned> StoichCell)
{
    mStoichCell = StoichCell;
}

void AbstractTransportReaction::SetNumberOfBulkSpecies(unsigned numberOfBulkSpecies)
{
    mNumberOfBulkSpecies = numberOfBulkSpecies;
}


unsigned AbstractTransportReaction::GetNumberOfBulkSpecies()
{
    return mNumberOfBulkSpecies;
}


void AbstractTransportReaction::SetNumberOfCellSpecies(unsigned numberOfCellSpecies)
{
    mNumberOfCellSpecies = numberOfCellSpecies;
}


unsigned AbstractTransportReaction::GetNumberOfCellSpecies()
{
    return mNumberOfCellSpecies;
}

void AbstractTransportReaction::SetIsRateCheck(bool IsRateCheck)
{
    mIsRateCheck = IsRateCheck;
}

bool AbstractTransportReaction::GetIsRateCheck()
{
    return mIsRateCheck;
}

void AbstractTransportReaction::SetDeltaErrorRateMin(double delta_rate_min)
{
    mDelta_rate_min = delta_rate_min;
}

double AbstractTransportReaction::GetDeltaErrorRateMin()
{
    return mDelta_rate_min;
}

void AbstractTransportReaction::SetDeltaErrorRateMax(double delta_rate_max)
{
    mDelta_rate_max = delta_rate_max;
}

double AbstractTransportReaction::GetDeltaErrorRateMax()
{
    return mDelta_rate_max;
}

double AbstractTransportReaction::CheckRate(double rate)
{
    // if reaction rate gets too low or high then undefined behaviour can occur

    if(abs(rate) < mDelta_rate_min)
    {
        rate = 0.0;
    }else if(abs(rate) > mDelta_rate_max)
    {
        rate = mDelta_rate_max;
    }
    return rate;
}


// file read functions
void AbstractTransportReaction::SetIrreversibleDelimiter(std::string delim)
{
    mIrreversibleDelimiter = delim;
}

std::string AbstractTransportReaction::GetIrreversibleDelimiter()
{
    return mIrreversibleDelimiter;
}

void AbstractTransportReaction::SetIrreversibleRateName(std::string rateName)
{
    mIrreversibleRateName = rateName;
}

std::string AbstractTransportReaction::GetIrreversibleRateName()
{
    return mIrreversibleRateName;
}

#endif