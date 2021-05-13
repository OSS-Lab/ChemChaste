#ifndef ABSTRACTMEMBRANEREACTION_HPP_
#define ABSTRACTMEMBRANEREACTION_HPP_

// general includes
#include <string>
#include <stdlib.h> 
#include <vector>
#include <iostream>
// custom includes
#include "AbstractChemistry.hpp"
#include "AbstractChemical.hpp"

// abstract property to contain information about the interactions of chemical species at a cell boundary
// define two coupled reactions; first the reaction in the bulk, second the reaction on the cell side
// bulk -> bulk | cell -> cell


class AbstractMembraneReaction 
{
protected:
    
    std::vector<AbstractChemical*> mpBulkSubstrates;

    std::vector<AbstractChemical*> mpBulkProducts;

    std::vector<AbstractChemical*> mpCellSubstrates;

    std::vector<AbstractChemical*> mpCellProducts;

    unsigned mNumberOfBulkProducts;

    unsigned mNumberOfBulkSubstrates;

    unsigned mNumberOfCellProducts;

    unsigned mNumberOfCellSubstrates;

    // vector containing the stoichmetry of the bulk substrates.  Will be "-ve" in reaction
    std::vector<unsigned> mStoichBulkSubstrates;

    // vector containing the stoichmetry of the bulk products.  Will be "+ve" in reaction
    std::vector<unsigned> mStoichBulkProducts;

    // vector containing the stoichmetry of the cell substrates.  Will be "-ve" in reaction
    std::vector<unsigned> mStoichCellSubstrates;

    // vector containing the stoichmetry of the cell products.  Will be "+ve" in reaction
    std::vector<unsigned> mStoichCellProducts;


    // container for the reaction rate constant.  '+ve' value. May be overriden in derived classes to account for functional rates
    double mReactionRate;

    bool mIsRateCheck;

    double mDelta_rate_min;

    double mDelta_rate_max;

    std::string mIrreversibleDelimiter = "->";
  
    std::string mIrreversibleRateName = "kf =";

public:

    // constructor
    AbstractMembraneReaction(   std::vector<AbstractChemical*> bulkSubstrates = std::vector<AbstractChemical*>(),
                                std::vector<AbstractChemical*> bulkProducts = std::vector<AbstractChemical*>(),
                                std::vector<AbstractChemical*> cellSubstrates = std::vector<AbstractChemical*>(),
                                std::vector<AbstractChemical*> cellProducts = std::vector<AbstractChemical*>(),
                                std::vector<unsigned> stoichBulkSubstrates = std::vector<unsigned>(),
                                std::vector<unsigned> stoichBulkProducts = std::vector<unsigned>(),
                                std::vector<unsigned> stoichCellSubstrates = std::vector<unsigned>(),
                                std::vector<unsigned> stoichCellProducts = std::vector<unsigned>(),
                                double reactionRate = 1.0
    );

    // copy constructor
    AbstractMembraneReaction(const AbstractMembraneReaction&);

    // destructor
    virtual ~AbstractMembraneReaction()
    {
    };

    // function to take in pointer to current concentration state vector of the state vector for change in concentration.  Basic update via multiplying by constant reaction rate
    virtual void React(AbstractChemistry*, AbstractChemistry*, const std::vector<double>&, const std::vector<double>&, std::vector<double>&, std::vector<double>&);

    // pure vitual function to update member variables and conditions of derived classes
    virtual void UpdateReaction();

    virtual void UpdateReactionRate(AbstractChemistry*, AbstractChemistry*, const std::vector<double>&, const std::vector<double>&);

    virtual void SetReactionRate(double);

    virtual double GetReactionRate();

    virtual std::string GetReactionType();

    virtual void ParseReactionInformation(std::string, bool);

    // Chemical handeling functions
    std::vector<AbstractChemical*> GetBulkSubstrates();

    AbstractChemical* GetBulkSubstratesByIndex(unsigned);

    std::vector<AbstractChemical*> GetCellSubstrates();

    AbstractChemical* GetCellSubstratesByIndex(unsigned);

    std::vector<AbstractChemical*> GetBulkProducts();

    AbstractChemical* GetBulkProductsByIndex(unsigned);

    std::vector<AbstractChemical*> GetCellProducts();

    AbstractChemical* GetCellProductsByIndex(unsigned);

    void SetBulkSubstrates(std::vector<AbstractChemical*>);

    void SetCellSubstrates(std::vector<AbstractChemical*>);

    void SetBulkProducts(std::vector<AbstractChemical*>);

    void SetCellProducts(std::vector<AbstractChemical*>);

    std::vector<unsigned> GetBulkStoichSubstrates();

    unsigned GetBulkStoichSubstratesByIndex(unsigned);

    void SetBulkStoichSubstrates(std::vector<unsigned>);

    std::vector<unsigned> GetCellStoichSubstrates();

    unsigned GetCellStoichSubstratesByIndex(unsigned);

    void SetCellStoichSubstrates(std::vector<unsigned>);

    std::vector<unsigned> GetBulkStoichProducts();

    unsigned GetBulkStoichProductsByIndex(unsigned);

    void SetBulkStoichProducts(std::vector<unsigned>);

    std::vector<unsigned> GetCellStoichProducts();

    unsigned GetCellStoichProductsByIndex(unsigned);

    void SetCellStoichProducts(std::vector<unsigned>);

    void SetNumberOfBulkSubstrates(unsigned);

    unsigned GetNumberOfBulkSubstrates();

    void SetNumberOfCellSubstrates(unsigned);

    unsigned GetNumberOfCellSubstrates();

    void SetNumberOfBulkProducts(unsigned);

    unsigned GetNumberOfBulkProducts();

    void SetNumberOfCellProducts(unsigned);

    unsigned GetNumberOfCellProducts();

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
AbstractMembraneReaction::AbstractMembraneReaction(   
                            std::vector<AbstractChemical*> bulkSubstrates,
                            std::vector<AbstractChemical*> bulkProducts,
                            std::vector<AbstractChemical*> cellSubstrates,
                            std::vector<AbstractChemical*> cellProducts,
                            std::vector<unsigned> stoichBulkSubstrates,
                            std::vector<unsigned> stoichBulkProducts,
                            std::vector<unsigned> stoichCellSubstrates,
                            std::vector<unsigned> stoichCellProducts,
                            double reactionRate)
        :   mpBulkSubstrates(bulkSubstrates),
            mpBulkProducts(bulkProducts),
            mpCellSubstrates(cellSubstrates),
            mpCellProducts(cellProducts),
            mStoichBulkSubstrates(stoichBulkSubstrates),
            mStoichBulkProducts(stoichBulkProducts),
            mStoichCellSubstrates(stoichCellSubstrates),
            mStoichCellProducts(stoichCellProducts),
            mReactionRate(reactionRate)
{
    mNumberOfBulkProducts = bulkProducts.size();
    mNumberOfBulkSubstrates = bulkSubstrates.size();
    mNumberOfCellProducts = cellProducts.size();
    mNumberOfCellSubstrates = cellSubstrates.size();

    SetIsRateCheck(true);
    SetDeltaErrorRateMin(1e-8);
    SetDeltaErrorRateMax(1e8);
}

// copy constructor
AbstractMembraneReaction::AbstractMembraneReaction(const AbstractMembraneReaction& existingReaction)
{
    mpBulkSubstrates = existingReaction.mpBulkSubstrates;
    mpBulkProducts = existingReaction.mpBulkProducts;
    mpCellSubstrates = existingReaction.mpCellSubstrates;
    mpCellProducts = existingReaction.mpCellProducts;
    mNumberOfBulkProducts = existingReaction.mNumberOfBulkProducts;
    mNumberOfBulkSubstrates = existingReaction.mNumberOfBulkSubstrates;
    mNumberOfCellProducts = existingReaction.mNumberOfCellProducts;
    mNumberOfCellSubstrates = existingReaction.mNumberOfCellSubstrates;
    mStoichBulkSubstrates = existingReaction.mStoichBulkSubstrates;
    mStoichBulkProducts = existingReaction.mStoichBulkProducts;
    mStoichCellSubstrates = existingReaction.mStoichCellSubstrates;
    mStoichCellProducts = existingReaction.mStoichCellProducts;
    mReactionRate = existingReaction.mReactionRate;
    mIsRateCheck = existingReaction.mIsRateCheck;
    mDelta_rate_min = existingReaction.mDelta_rate_min;
    mDelta_rate_max = existingReaction.mDelta_rate_max;
}


// function to take in pointer to current concentration state vector of the state vector for change in concentration.  Basic update via multiplying by constant reaction rate
void AbstractMembraneReaction::React(AbstractChemistry* bulkChemistry, AbstractChemistry* cellChemistry, const std::vector<double>& currentBulkConcentration, const std::vector<double>& currentCellConcentration, std::vector<double>& changeBulkConc, std::vector<double>& changeCellConc)
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
                changeBulkConc[index] -= mStoichBulkSubstrates[j]*GetReactionRate();
                break;
            }
        }
        // a reactant may be present on both sides of the reaction and may convert at different functional rates
        for(unsigned j=0; j<mNumberOfBulkProducts; j++)
        {
            if(mpBulkProducts[j] -> GetChemicalName()==p_system_chemical -> GetChemicalName())
            {
                changeBulkConc[index] += mStoichBulkProducts[j]*GetReactionRate();
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
                changeCellConc[index] -= mStoichCellSubstrates[j]*GetReactionRate();
                break;
            }
        }
        // a reactant may be present on both sides of the reaction and may convert at different functional rates
        for(unsigned j=0; j<mNumberOfCellProducts; j++)
        {
            if(mpCellProducts[j] -> GetChemicalName()==p_system_chemical -> GetChemicalName())
            {
                changeCellConc[index] += mStoichCellProducts[j]*GetReactionRate();
                break;
            }
        }
    }

}

// pure vitual function to update member variables and conditions of derived classes
void AbstractMembraneReaction::UpdateReaction()
{
    return;
}

void AbstractMembraneReaction::UpdateReactionRate(AbstractChemistry* bulkChemistry, AbstractChemistry* cellChemistry, const std::vector<double>& currentBulkConc, const std::vector<double>& currentCellConc)
{
    // default to zeroth order reaction
    SetReactionRate(mReactionRate);
}

void AbstractMembraneReaction::SetReactionRate(double reactionRate)
{
    if(mIsRateCheck)
    {
        reactionRate = CheckRate(reactionRate);
    }
    
    mReactionRate = reactionRate;
}

double AbstractMembraneReaction::GetReactionRate()
{
    return mReactionRate;
}

std::string AbstractMembraneReaction::GetReactionType()
{
    // virtual function to be overriden in derived classes, used to idnetify reaction inheritance 
    return "ZerothOrderCoupledMembrane";
}

void AbstractMembraneReaction::ParseReactionInformation(std::string reaction_information, bool IsReversible=false)
{
    if(reaction_information.find(mIrreversibleRateName) != std::string::npos)
    {
        if(reaction_information.find(";") != std::string::npos)
        {
            size_t pos= reaction_information.find(";");
            reaction_information = reaction_information.substr(0,pos);
        }
        size_t pos= reaction_information.find(mIrreversibleRateName);

        SetReactionRate(atof(reaction_information.substr(pos+mIrreversibleRateName.size()+1,std::string::npos).c_str()));
    }
}

// Chemical handeling functions
std::vector<AbstractChemical*> AbstractMembraneReaction::GetBulkSubstrates()
{
    return mpBulkSubstrates;
}

AbstractChemical* AbstractMembraneReaction::GetBulkSubstratesByIndex(unsigned index)
{
    if(index < mNumberOfBulkSubstrates)
    {
        return mpBulkSubstrates[index];
    }
    else
    {
        std::cout<<"Error: AbstractMembraneReaction::GetBulkSubstratesByIndex(unsigned index), index out of bounds"<<std::endl;
        return new AbstractChemical();
    } 
}

std::vector<AbstractChemical*> AbstractMembraneReaction::GetCellSubstrates()
{
    return mpCellSubstrates;
}

AbstractChemical* AbstractMembraneReaction::GetCellSubstratesByIndex(unsigned index)
{
    if(index < mNumberOfCellSubstrates)
    {
        return mpCellSubstrates[index];
    }
    else
    {
        std::cout<<"Error: AbstractMembraneReaction::GetCellSubstratesByIndex(unsigned index), index out of bounds"<<std::endl;
        return new AbstractChemical();
    } 
}

std::vector<AbstractChemical*> AbstractMembraneReaction::GetBulkProducts()
{
    return mpBulkProducts;
}

AbstractChemical* AbstractMembraneReaction::GetBulkProductsByIndex(unsigned index)
{
    if(index < mNumberOfBulkProducts)
    {
        return mpBulkProducts[index];
    }
    else
    {
        std::cout<<"Error: AbstractMembraneReaction::GetBulkProductsByIndex(unsigned index), index out of bounds"<<std::endl;
        return new AbstractChemical();
    } 
}


std::vector<AbstractChemical*> AbstractMembraneReaction::GetCellProducts()
{
    return mpCellProducts;
}

AbstractChemical* AbstractMembraneReaction::GetCellProductsByIndex(unsigned index)
{
    if(index < mNumberOfCellProducts)
    {
        return mpCellProducts[index];
    }
    else
    {
        std::cout<<"Error: AbstractMembraneReaction::GetCellProductsByIndex(unsigned index), index out of bounds"<<std::endl;
        return new AbstractChemical();
    } 
}

void AbstractMembraneReaction::SetBulkSubstrates(std::vector<AbstractChemical*> substrates)
{
    mpBulkSubstrates = substrates;
}

void AbstractMembraneReaction::SetCellSubstrates(std::vector<AbstractChemical*> substrates)
{
    mpCellSubstrates = substrates;
}

void AbstractMembraneReaction::SetBulkProducts(std::vector<AbstractChemical*> products)
{
    mpBulkProducts = products;
}

void AbstractMembraneReaction::SetCellProducts(std::vector<AbstractChemical*> products)
{
    mpCellProducts = products;
}

std::vector<unsigned> AbstractMembraneReaction::GetBulkStoichSubstrates()
{
    return mStoichBulkSubstrates;
}

unsigned AbstractMembraneReaction::GetBulkStoichSubstratesByIndex(unsigned index)
{
    if(index < mNumberOfBulkSubstrates)
    {
        return mStoichBulkSubstrates[index];
    }
    else
    {
        std::cout<<"Error: AbstractMembraneReaction::GetBulkStoichSubstratesByIndex(unsigned index), index out of bounds"<<std::endl;
        return 0;
    } 
}

void AbstractMembraneReaction::SetBulkStoichSubstrates(std::vector<unsigned> stoich)
{
    mStoichBulkSubstrates = stoich;
}

std::vector<unsigned> AbstractMembraneReaction::GetCellStoichSubstrates()
{
    return mStoichCellSubstrates;
}

unsigned AbstractMembraneReaction::GetCellStoichSubstratesByIndex(unsigned index)
{
    if(index < mNumberOfCellSubstrates)
    {
        return mStoichCellSubstrates[index];
    }
    else
    {
        std::cout<<"Error: AbstractMembraneReaction::GetCellStoichSubstratesByIndex(unsigned index), index out of bounds"<<std::endl;
        return 0;
    } 
}

void AbstractMembraneReaction::SetCellStoichSubstrates(std::vector<unsigned> stoich)
{
    mStoichCellSubstrates = stoich;
}

std::vector<unsigned> AbstractMembraneReaction::GetBulkStoichProducts()
{
    return mStoichBulkProducts;
}

unsigned AbstractMembraneReaction::GetBulkStoichProductsByIndex(unsigned index)
{
    if(index < mNumberOfBulkProducts)
    {
        return mStoichBulkProducts[index];
    }
    else
    {
        std::cout<<"Error: AbstractMembraneReaction::GetBulkStoichProductsByIndex(unsigned index), index out of bounds"<<std::endl;
        return 0;
    } 
}

void AbstractMembraneReaction::SetBulkStoichProducts(std::vector<unsigned> stoich)
{
    mStoichBulkProducts = stoich;
}

std::vector<unsigned> AbstractMembraneReaction::GetCellStoichProducts()
{
    return mStoichCellProducts;
}

unsigned AbstractMembraneReaction::GetCellStoichProductsByIndex(unsigned index)
{
    if(index < mNumberOfCellProducts)
    {
        return mStoichCellProducts[index];
    }
    else
    {
        std::cout<<"Error: AbstractMembraneReaction::GetCellStoichProductsByIndex(unsigned index), index out of bounds"<<std::endl;
        return 0;
    } 
}

void AbstractMembraneReaction::SetCellStoichProducts(std::vector<unsigned> stoich)
{
    mStoichCellProducts = stoich;
}

void AbstractMembraneReaction::SetNumberOfBulkSubstrates(unsigned numberSubstrates)
{
    mNumberOfBulkSubstrates = numberSubstrates;
}

unsigned AbstractMembraneReaction::GetNumberOfBulkSubstrates()
{
    return mNumberOfBulkSubstrates;
}

void AbstractMembraneReaction::SetNumberOfCellSubstrates(unsigned numberSubstrates)
{
    mNumberOfCellSubstrates = numberSubstrates;
}

unsigned AbstractMembraneReaction::GetNumberOfCellSubstrates()
{
    return mNumberOfCellSubstrates;
}

void AbstractMembraneReaction::SetNumberOfBulkProducts(unsigned numberProducts)
{
    mNumberOfBulkProducts = numberProducts;
}

unsigned AbstractMembraneReaction::GetNumberOfBulkProducts()
{
    return mNumberOfBulkProducts;
}

void AbstractMembraneReaction::SetNumberOfCellProducts(unsigned numberProducts)
{
    mNumberOfCellProducts = numberProducts;
}

unsigned AbstractMembraneReaction::GetNumberOfCellProducts()
{
    return mNumberOfCellProducts;
}


// reaction concentration checking functions
void AbstractMembraneReaction::SetIsRateCheck(bool IsRateCheck)
{
    mIsRateCheck = IsRateCheck;
}

bool AbstractMembraneReaction::GetIsRateCheck()
{
    return mIsRateCheck;
}

void AbstractMembraneReaction::SetDeltaErrorRateMin(double delta_rate_min)
{
    mDelta_rate_min = delta_rate_min;
}

double AbstractMembraneReaction::GetDeltaErrorRateMin()
{
    return mDelta_rate_min;
}

void AbstractMembraneReaction::SetDeltaErrorRateMax(double delta_rate_max)
{
    mDelta_rate_max = delta_rate_max;
}

double AbstractMembraneReaction::GetDeltaErrorRateMax()
{
    return mDelta_rate_max;
}

double AbstractMembraneReaction::CheckRate(double rate)
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
void AbstractMembraneReaction::SetIrreversibleDelimiter(std::string delim)
{
    mIrreversibleDelimiter = delim;
}

std::string AbstractMembraneReaction::GetIrreversibleDelimiter()
{
    return mIrreversibleDelimiter;
}

void AbstractMembraneReaction::SetIrreversibleRateName(std::string rateName)
{
    mIrreversibleRateName = rateName;
}

std::string AbstractMembraneReaction::GetIrreversibleRateName()
{
    return mIrreversibleRateName;
}

#endif