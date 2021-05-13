#ifndef ABSTRACTTRANSPORTOUTREACTION_HPP_
#define ABSTRACTTRANSPORTOUTREACTION_HPP_

// general includes
#include <string>
#include <stdlib.h> 
#include <vector>
#include <iostream>
// custom includes
#include "AbstractChemistry.hpp"
#include "AbstractChemical.hpp"
#include "AbstractTransportReaction.hpp"

// abstract property to contain information about the interactions of chemical species at a cell boundary
// bulk <- cell
// reaction substrates denote the species in the bulk while products denote species in the cell

// ZerothOrderTransportOutOfCell

class AbstractTransportOutReaction : public AbstractTransportReaction
{
private:

    using AbstractTransportReaction::React;
    using AbstractTransportReaction::GetReactionType;
    using AbstractTransportReaction::ParseReactionInformation;

public:

    // constructor
    AbstractTransportOutReaction(   std::vector<AbstractChemical*> bulkReactionSpecies = std::vector<AbstractChemical*>(),
                                    std::vector<AbstractChemical*> cellReactionSpecies = std::vector<AbstractChemical*>(),
                                    std::vector<unsigned> stoichBulk = std::vector<unsigned>(),
                                    std::vector<unsigned> stoichCell = std::vector<unsigned>(),
                                    double reactionRate = 1.0
    );
    

    // destructor
    virtual ~AbstractTransportOutReaction()
    {
    };


    // function to take in pointer to current concentration state vector of the state vector for change in cocnentration.  Basic update via multiplying by constant reaction rate
    virtual void React(AbstractChemistry*, AbstractChemistry*,const std::vector<double>&, const std::vector<double>&, std::vector<double>&, std::vector<double>&);

    virtual std::string GetReactionType();

    virtual void ParseReactionInformation(std::string, bool);

    // file read functions
    void SetIrreversibleDelimiter(std::string );


    std::string GetIrreversibleDelimiter();

    void SetIrreversibleRateName(std::string );

    std::string GetIrreversibleRateName();

};

// constructor
AbstractTransportOutReaction::AbstractTransportOutReaction(     std::vector<AbstractChemical*> bulkReactionSpecies,
                                                                std::vector<AbstractChemical*> cellReactionSpecies,
                                                                std::vector<unsigned> stoichBulk,
                                                                std::vector<unsigned> stoichCell,
                                                                double reactionRate)
        :   AbstractTransportReaction(
                            bulkReactionSpecies,
                            cellReactionSpecies,
                            stoichBulk,
                            stoichCell,
                            reactionRate
                            )
{
    mIrreversibleDelimiter = "<-";
    mIrreversibleRateName = "kr =";
}


void AbstractTransportOutReaction::React(AbstractChemistry* bulkChemistry, AbstractChemistry* cellChemistry, const std::vector<double>& currentBulkConcentration, const std::vector<double>& currentCellConcentration, std::vector<double>& changeBulkConc, std::vector<double>& changeCellConc)
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
        AbstractChemical *p_bulk_chemical = dynamic_cast<AbstractChemical*>(*chem_iter);

        // for each system chemical, parse whether it is involved in this reaction.
        for(unsigned j=0; j<mNumberOfBulkSpecies; j++)
        {
            if(mpBulkReactionSpecies[j] -> GetChemicalName()==p_bulk_chemical -> GetChemicalName())
            {
                changeBulkConc[index] += mStoichBulk[j]*GetReactionRate();
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
        AbstractChemical *p_cell_chemical = dynamic_cast<AbstractChemical*>(*chem_iter);

        // for each cell chemical, parse whether it is involved in this reaction.
        for(unsigned j=0; j<mNumberOfCellSpecies; j++)
        {
            if(mpCellReactionSpecies[j] -> GetChemicalName()==p_cell_chemical -> GetChemicalName())
            {
                changeCellConc[index] -= mStoichCell[j]*GetReactionRate();
                break;
            }
        }
    }
    
}


std::string AbstractTransportOutReaction::GetReactionType()
{
    // virtual function to be overriden in derived classes, used to idnetify reaction inheritance 
    return "ZerothOrderTransportOutOfCell";
}


void AbstractTransportOutReaction::ParseReactionInformation(std::string reaction_information, bool IsReversible=false)
{
    if(reaction_information.find(mIrreversibleRateName) != std::string::npos)
    {

        size_t pos= reaction_information.find(mIrreversibleRateName);

        SetReactionRate(atof(reaction_information.substr(pos+mIrreversibleRateName.size()+1,std::string::npos).c_str()));
    }
}


// file read functions
void AbstractTransportOutReaction::SetIrreversibleDelimiter(std::string delim)
{
    mIrreversibleDelimiter = delim;
}

std::string AbstractTransportOutReaction::GetIrreversibleDelimiter()
{
    return mIrreversibleDelimiter;
}

void AbstractTransportOutReaction::SetIrreversibleRateName(std::string rateName)
{
    mIrreversibleRateName = rateName;
}

std::string AbstractTransportOutReaction::GetIrreversibleRateName()
{
    return mIrreversibleRateName;
}

#endif