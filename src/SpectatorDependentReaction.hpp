#ifndef SPECTATORDEPENDENTREACTION_HPP_
#define SPECTATORDEPENDENTREACTION_HPP_

// general includes
#include <string>
#include <tuple>
#include <vector>
#include <cmath>

// custom includes
#include "AbstractChemical.hpp"
#include "AbstractChemistry.hpp"
#include "AbstractReaction.hpp"

// reaction class which takes in chemical data for the whole system and uses
// species not involved directly in the reaction to change the reaction rate 

class SpectatorDependentReaction : public AbstractReaction
{
protected:

    std::vector<std::string> mSpectatorNames;

    unsigned mNumberOfSpectators;

    std::string mSpectatorDelimiter = "Spectator =";

    AbstractChemistry* mpSystemChemistry;

    double mReactionRateConstant;
private:
    using AbstractReaction::UpdateReactionRate;
    using AbstractReaction::GetReactionType;
    using AbstractReaction::ParseReactionInformation;

public:

    // constructor
    SpectatorDependentReaction( std::vector<AbstractChemical*> substrates = std::vector<AbstractChemical*>(),
                                std::vector<AbstractChemical*> products = std::vector<AbstractChemical*>(),
                                std::vector<unsigned> stoichSubstrates = std::vector<unsigned>(),
                                std::vector<unsigned> stoichProducts = std::vector<unsigned>(),
                                AbstractChemistry* p_systemChemistry = new AbstractChemistry(),
                                double reactionRate = 1.0
                                );

    // destructor
    virtual ~SpectatorDependentReaction()
    {
    };

    virtual void UpdateReactionRate(AbstractChemistry*, const std::vector<double>&);

    virtual std::string GetReactionType();

    virtual void ParseReactionInformation(std::string, bool);

    // class specific methods

    void UpdateSystemChemistry(std::vector<std::string>);

    void SetSpectatorNames(std::vector<std::string>);

    std::vector<std::string> GetSpectatorNames();

    void SetReactionRateConstant(double);

    double GetReactionRateConstant();

    void SetNumberOfSpectators(unsigned);

    unsigned GetNumberOfSpectators();
};

SpectatorDependentReaction::SpectatorDependentReaction( 
                                std::vector<AbstractChemical*> substrates,
                                std::vector<AbstractChemical*> products,
                                std::vector<unsigned> stoichSubstrates,
                                std::vector<unsigned> stoichProducts,
                                AbstractChemistry* p_systemChemistry,
                                double reactionRate)
                            : AbstractReaction(
                                substrates,
                                products,
                                stoichSubstrates,
                                stoichProducts,
                                reactionRate
                            ),
                            mpSystemChemistry(p_systemChemistry)
                            {
                            }


void SpectatorDependentReaction::UpdateReactionRate(AbstractChemistry* systemChemistry, const std::vector<double>& currentSystemConc)
{
    // multiply the reaction rate constant by the product of the spectator set concentrations

    double spectatorConcentration =1.0;

    unsigned spectatorCount=0;

    std::vector<AbstractChemical*> p_chemical_vector = systemChemistry -> rGetChemicalVector();
    unsigned index = 0;
    for(std::vector<AbstractChemical*>::iterator chem_iter = p_chemical_vector.begin();
            chem_iter != p_chemical_vector.end();
            ++chem_iter, ++index)
    {
        if(spectatorCount == mNumberOfSpectators)
        {
            // don't need to continue processing the species
            break;
        }
        else
        {
            AbstractChemical *p_system_chemical = dynamic_cast<AbstractChemical*>(*chem_iter);

            for(unsigned i=0; i<mNumberOfSpectators; i++)
            {
                if(mSpectatorNames[i]==p_system_chemical -> GetChemicalName())
                {
                    // calculate the product of spectator concentrations
                    spectatorConcentration *=  currentSystemConc[index];
                    // increase the count of spectators processed
                    spectatorCount = spectatorCount+1;
                    break;
                }
            }
        }
    };
    if(spectatorCount==0)
    {
        // if no spectators are in the system, default to zero rate
        spectatorConcentration=0.0;
    }
    SetReactionRate(mReactionRateConstant*spectatorConcentration);
}

std::string SpectatorDependentReaction::GetReactionType()
{
    return "SpectatorDependentReaction";
}

void SpectatorDependentReaction::ParseReactionInformation(std::string reaction_information, bool IsReversible=false)
{
    // read in from the reaction information string the necessary information for the basic multiplicative constant
    // and the spectator species, may be multiple spectator species but only one reaction rate, they can be in any order

    // assume that there must be a spectator
    size_t positionSpectator = reaction_information.find(mSpectatorDelimiter);
    // initialse a temporary position, used if there are more than one spectator species
    size_t tempPositionSpectator = reaction_information.find(mSpectatorDelimiter);
    if(positionSpectator==std::string::npos)
    {
        // spectator not found
        std::cout<<"Error SpectatorDependentReaction::ParseReactionInformation: spectator species not found"<<std::endl;
    }

    // assume that there must be reaction rate constant
    size_t positionRate = reaction_information.find(mIrreversibleRateName);

    if(positionRate==std::string::npos)
    {
        // spectator not found
        std::cout<<"Error SpectatorDependentReaction::ParseReactionInformation: reaction rate not found"<<std::endl;
    }

    double reactionRate;
    std::string spectator;
    std::vector<std::string> spectatorNames;
    std::vector<unsigned> spectatorIndices;
    std::string tempString;
    

    while(positionSpectator != std::string::npos)
    {
        if(positionRate<positionSpectator)
        {
            // rate value is first in the string
            reactionRate = atof(reaction_information.substr(positionRate+mIrreversibleRateName.size()+1,positionSpectator).c_str());
            reaction_information = reaction_information.substr(positionRate+mIrreversibleRateName.size()+1,std::string::npos);

        }
        else
        {
            // next data element is a spectator species, test if the next value is a spectator or rate
            tempString = reaction_information.substr(positionSpectator+mSpectatorDelimiter.size()+1,std::string::npos);

            // find the new positions of the delimiters
            tempPositionSpectator = tempString.find(mSpectatorDelimiter); // if no more spectators, should be npos
            positionRate = tempString.find(mIrreversibleRateName); // if rate already taken, should be npos

            if(positionRate>tempPositionSpectator)
            {
                // next value is spectator 
                spectator = reaction_information.substr(positionSpectator+mSpectatorDelimiter.size()+1,tempPositionSpectator);

            }
            else if(positionRate<tempPositionSpectator)
            {
                // next value is rate
                spectator = reaction_information.substr(positionSpectator+mSpectatorDelimiter.size()+1,positionRate-1);
            }
            else
            {
                // both are npos so no other value
                spectator = reaction_information.substr(positionSpectator+mSpectatorDelimiter.size()+1,std::string::npos);
            }
            
            spectatorNames.push_back(spectator);
            
            reaction_information = reaction_information.substr(positionSpectator+mSpectatorDelimiter.size()+1,std::string::npos);

        }

        // update the delimiter positions from the new substrings
        positionSpectator = reaction_information.find(mSpectatorDelimiter);
        positionRate = reaction_information.find(mIrreversibleRateName);
        
        if(positionSpectator == std::string::npos && positionRate !=std::string::npos)
        {
            // there are no more species but the rate is still unknown and in the string
            reactionRate = atof(reaction_information.substr(positionRate+mIrreversibleRateName.size()+1,std::string::npos).c_str());

        }
    }

    SetReactionRateConstant(reactionRate);
    SetSpectatorNames(spectatorNames);
    SetNumberOfSpectators(spectatorNames.size());
    UpdateSystemChemistry(spectatorNames);
}

// class specific methods

void SpectatorDependentReaction::UpdateSystemChemistry(std::vector<std::string> spectatorNames)
{
    // if spectator species are not involved within the reactions of the reaction system then they will not
    // be present in the system chemistry.  As the chemical ODE system derives the state variable vector,
    // which controls the concentrations data flow, from the reaction system chemistry the spectators need 
    // to be in the chemistry for the concentrations to be avaliable

    for(unsigned i=0; i<mNumberOfSpectators; i++)
    {
        mpSystemChemistry -> AddChemical(new AbstractChemical(spectatorNames[i]));
    }
}

void SpectatorDependentReaction::SetSpectatorNames(std::vector<std::string> names)
{
    mSpectatorNames = names;
}


std::vector<std::string> SpectatorDependentReaction::GetSpectatorNames()
{
    return mSpectatorNames;
}


void SpectatorDependentReaction::SetReactionRateConstant(double reactionRate)
{
    mReactionRateConstant = reactionRate;
}

double SpectatorDependentReaction::GetReactionRateConstant()
{
    return mReactionRateConstant;
}

void SpectatorDependentReaction::SetNumberOfSpectators(unsigned numberOfSpectators)
{
    mNumberOfSpectators = numberOfSpectators;
}

unsigned SpectatorDependentReaction::GetNumberOfSpectators()
{
    return mNumberOfSpectators;
}

#endif