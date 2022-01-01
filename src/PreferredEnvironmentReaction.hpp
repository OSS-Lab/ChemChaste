#ifndef PREFERREDENVIRONMENTREACTION_HPP_
#define PREFERREDENVIRONMENTREACTION_HPP_

// general includes
#include <string>
#include <tuple>
#include <vector>
#include <cmath>

// custom includes
#include "AbstractChemical.hpp"
#include "AbstractChemistry.hpp"
#include "AbstractReaction.hpp"
#include "EnvironmentDependentReaction.hpp"
#include "Cell.hpp"


// reaction class which takes in chemical data for the whole system and uses
// species not involved directly in the reaction to change the reaction rate 

class PreferredEnvironmentReaction : public EnvironmentDependentReaction
{
protected:

    CellPtr mpCell; // ptr to the cell to check environment cell property

    std::vector<std::string> mEnvironmentChemicalNames;

    unsigned mNumberOfEnvironmentChemicals;

    std::string mEnvironmentDelimiter = "Environment =";

    double mReactionRateConstant;


private:
    using AbstractReaction::UpdateReactionRate;
    using AbstractReaction::GetReactionType;
    using AbstractReaction::ParseReactionInformation;

public:

    // constructor
    PreferredEnvironmentReaction( std::vector<AbstractChemical*> substrates = std::vector<AbstractChemical*>(),
                                std::vector<AbstractChemical*> products = std::vector<AbstractChemical*>(),
                                std::vector<unsigned> stoichSubstrates = std::vector<unsigned>(),
                                std::vector<unsigned> stoichProducts = std::vector<unsigned>(),
                                AbstractChemistry* p_systemChemistry = new AbstractChemistry(),
                                double reactionRate = 1.0
                                );

    // destructor
    virtual ~PreferredEnvironmentReaction()
    {
    };


    virtual void UpdateReactionRate(AbstractChemistry*, const std::vector<double>&);

    virtual std::string GetReactionType();

    virtual void ParseReactionInformation(std::string, bool);

    virtual void GiveCell(CellPtr);

    // class specific methods

    void SetEnvironmentChemicalNames(std::vector<std::string>);

    std::vector<std::string> GetEnvironmentChemicalNames();

    void SetReactionRateConstant(double);

    double GetReactionRateConstant();

    void SetNumberOfEnvironmentChemicals(unsigned);

    unsigned GetNumberOfEnvironmentChemicals();

    void SetEnvironmentDelimiter(std::string);
    
    std::string GetEnvironmentDelimiter();

};

void PreferredEnvironmentReaction::UpdateReactionRate(AbstractChemistry* systemChemistry, const std::vector<double>& currentSystemConc)
{
    // multiply the reaction rate constant by the product of the species set concentrations

    double environmentConcentration =1.0;
    
    assert(mpCell->rGetCellPropertyCollection().HasProperty<EnvironmentCellProperty>());

    boost::shared_ptr<EnvironmentCellProperty> environmentProp = boost::static_pointer_cast<EnvironmentCellProperty>(mpCell->rGetCellPropertyCollection().GetPropertiesType<TransportCellProperty>().GetProperty());
    
    StateVariableRegister* mpEnvironmentStateVariableRegister = environmentProp -> GetEnvironmentStateVariableRegister();

    std::vector<double> environmentConcentrationVector = environmentProp -> GetEnvironmentVector();

    std::vector<double> preferredEnvironmentConcentrationVector = environmentProp -> GetPreferredEnvironmentVector();

    unsigned index = 0;
    for(unsigned chem=0; chem<mNumberOfEnvironmentChemicals; chem++)
    {
        index = mpEnvironmentStateVariableRegister->RetrieveStateVariableIndex(mEnvironmentChemicalNames[chem]);
        
        // calculate the product of environment concentrations
        environmentConcentration *=  abs(environmentConcentrationVector[index] - preferredEnvironmentConcentrationVector[index]);
        
    };
    SetReactionRate(mReactionRateConstant*exp(-environmentConcentration));
}


std::string PreferredEnvironmentReaction::GetReactionType()
{
    return "PreferredEnvironmentReaction";
}

void PreferredEnvironmentReaction::ParseReactionInformation(std::string, bool)
{
    // read in from the reaction information string the necessary information for the basic multiplicative constant
    // and the species, may be multiple species but only one reaction rate, they can be in any order

    // assume that there must be a species
    size_t positionSpecies = reaction_information.find(mEnvironmentDelimiter);
    // initialse a temporary position, used if there are more than one species
    size_t tempPositionSpecies = reaction_information.find(mEnvironmentDelimiter);
    if(positionSpecies==std::string::npos)
    {
        // species not found
        std::cout<<"Error PreferredEnvironmentReaction::ParseReactionInformation: environment species not found"<<std::endl;
    }

    // assume that there must be reaction rate constant
    size_t positionRate = reaction_information.find(mIrreversibleRateName);

    if(positionRate==std::string::npos)
    {
        // rate not found
        std::cout<<"Error PreferredEnvironmentReaction::ParseReactionInformation: reaction rate not found"<<std::endl;
    }

    double reactionRate;
    std::string species;
    std::vector<std::string> environmentNames;
    std::vector<unsigned> speciesIndices;
    std::string tempString;
    

    while(positionSpecies != std::string::npos)
    {
        if(positionRate<positionSpecies)
        {
            // rate value is first in the string
            reactionRate = atof(reaction_information.substr(positionRate+mIrreversibleRateName.size()+1,positionSpecies).c_str());
            reaction_information = reaction_information.substr(positionRate+mIrreversibleRateName.size()+1,std::string::npos);

        }
        else
        {
            // next data element is a environment species, test if the next value is an environment species or rate
            tempString = reaction_information.substr(positionSpecies+mEnvironmentDelimiter.size()+1,std::string::npos);

            // find the new positions of the delimiters
            tempPositionSpecies = tempString.find(mEnvironmentDelimiter); // if no more species, should be npos
            positionRate = tempString.find(mIrreversibleRateName); // if rate already taken, should be npos

            if(positionRate>tempPositionSpecies)
            {
                // next value is species 
                species = reaction_information.substr(positionSpecies+mEnvironmentDelimiter.size()+1,tempPositionSpecies);

            }
            else if(positionRate<tempPositionSpecies)
            {
                // next value is rate
                species = reaction_information.substr(positionSpecies+mEnvironmentDelimiter.size()+1,positionRate-1);
            }
            else
            {
                // both are npos so no other value
                species = reaction_information.substr(positionSpecies+mEnvironmentDelimiter.size()+1,std::string::npos);
            }
            
            environmentNames.push_back(species);
            
            reaction_information = reaction_information.substr(positionSpecies+mEnvironmentDelimiter.size()+1,std::string::npos);

        }

        // update the delimiter positions from the new substrings
        positionSpecies = reaction_information.find(mEnvironmentDelimiter);
        positionRate = reaction_information.find(mIrreversibleRateName);
        
        if(positionSpecies == std::string::npos && positionRate !=std::string::npos)
        {
            // there are no more species but the rate is still unknown and in the string
            reactionRate = atof(reaction_information.substr(positionRate+mIrreversibleRateName.size()+1,std::string::npos).c_str());

        }
    }

    SetReactionRateConstant(reactionRate);
    SetEnvironmentChemicalNames(environmentNames);
    SetNumberOfEnvironmentSpecies(environmentNames.size());
}

void PreferredEnvironmentReaction::GiveCell(CellPtr p_cell)
{
    mpCell = p_cell;
}

void PreferredEnvironmentReaction::SetEnvironmentChemicalNames(std::vector<std::string> names)
{
    mEnvironmentChemicalNames = names;
}


std::vector<std::string> PreferredEnvironmentReaction::GetEnvironmentChemicalNames()
{
    return mEnvironmentChemicalNames;
}


void PreferredEnvironmentReaction::SetReactionRateConstant(double reactionRate)
{
    mReactionRateConstant = reactionRate;
}

double PreferredEnvironmentReaction::GetReactionRateConstant()
{
    return mReactionRateConstant;
}

void PreferredEnvironmentReaction::SetNumberOfEnvironmentChemicals(unsigned numberOfChemicals)
{
    mNumberOfEnvironmentChemicals = numberOfChemicals;
}

unsigned PreferredEnvironmentReaction::GetNumberOfEnvironmentChemicals()
{
    return mNumberOfEnvironmentChemicals;
}

void PreferredEnvironmentReaction::SetEnvironmentDelimiter(std::string delim)
{
    mEnvironmentDelimiter = delim;
}
    
std::string PreferredEnvironmentReaction::GetEnvironmentDelimiter()
{
    return mEnvironmentDelimiter;
}




#endif