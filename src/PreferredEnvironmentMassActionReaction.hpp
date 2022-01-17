#ifndef PREFERREDENVIRONMENTMASSACTIONREACTION_HPP_
#define PREFERREDENVIRONMENTMASSACTIONREACTION_HPP_

class ComplexCell;
class AbstractChemical;
class AbstractChemistry;
// general includes
#include <string>
#include <tuple>
#include <vector>
#include <cmath>

// custom includes
//#include "ComplexCell.hpp"
//#include "AbstractChemical.hpp"
//#include "AbstractChemistry.hpp"
//#include "AbstractReaction.hpp"
//#include "EnvironmentDependentReaction.hpp"

//#include "EnvironmentCellProperty.hpp"

// reaction class which takes in chemical data for the whole system and uses
// species not involved directly in the reaction to change the reaction rate 


class PreferredEnvironmentMassActionReaction : public EnvironmentDependentReaction
{
private:
    using EnvironmentDependentReaction::UpdateReactionRate;
    using EnvironmentDependentReaction::GetReactionType;
    using EnvironmentDependentReaction::ParseReactionInformation;

public:

    // constructor
    PreferredEnvironmentMassActionReaction( 
                                std::vector<AbstractChemical*> substrates = std::vector<AbstractChemical*>(),
                                std::vector<AbstractChemical*> products = std::vector<AbstractChemical*>(),
                                std::vector<unsigned> stoichSubstrates = std::vector<unsigned>(),
                                std::vector<unsigned> stoichProducts = std::vector<unsigned>(),
                                double reactionRate = 1.0
                                );

    // destructor
    virtual ~PreferredEnvironmentMassActionReaction()
    {
    };


    virtual void UpdateReactionRate(AbstractChemistry*, const std::vector<double>&);

    virtual std::string GetReactionType();

    virtual void ParseReactionInformation(std::string, bool);

    //virtual void GiveCell(CellPtr);

    // class specific methods


};


PreferredEnvironmentMassActionReaction::PreferredEnvironmentMassActionReaction( 
                                std::vector<AbstractChemical*> substrates,
                                std::vector<AbstractChemical*> products,
                                std::vector<unsigned> stoichSubstrates,
                                std::vector<unsigned> stoichProducts,
                                double reactionRate)
                            : EnvironmentDependentReaction(
                                substrates,
                                products,
                                stoichSubstrates,
                                stoichProducts,
                                reactionRate
                            )
                            {
                            }


void PreferredEnvironmentMassActionReaction::UpdateReactionRate(AbstractChemistry* systemChemistry, const std::vector<double>& currentSystemConc)
{
    // multiply the reaction rate constant by the product of the species set concentrations

    double environmentConcentration =1.0;

    double forwardFlux=1.0;

    boost::shared_ptr<ComplexCell> cellPtr = boost::dynamic_pointer_cast<ComplexCell>(mpCell);
    assert(mpCell->rGetCellPropertyCollection().HasProperty<EnvironmentCellProperty>());

    CellPropertyCollection cPCollection = cellPtr->rGetCellPropertyCollection();
    

    boost::shared_ptr<EnvironmentCellProperty> environmentProp = boost::static_pointer_cast<EnvironmentCellProperty>(cPCollection.GetPropertiesType<EnvironmentCellProperty>().GetProperty());

    StateVariableRegister* pEnvironmentStateVariableRegister = environmentProp -> GetEnvironmentStateVariableRegister();

    std::vector<double> environmentConcentrationVector = environmentProp -> GetEnvironmentVector();

    std::vector<double> preferredEnvironmentConcentrationVector = environmentProp -> GetPreferredEnvironmentVector();
    unsigned index = 0;
    for(unsigned chem=0; chem<this->mNumberOfEnvironmentChemicals; chem++)
    {
        index = pEnvironmentStateVariableRegister->RetrieveStateVariableIndex(this->mEnvironmentChemicalNames[chem]);
        // calculate the product of environment concentrations
        environmentConcentration *=  abs(environmentConcentrationVector[index] - preferredEnvironmentConcentrationVector[index]);
    };

    std::vector<AbstractChemical*> p_chemical_vector = systemChemistry -> rGetChemicalVector();
    index = 0;
    for(std::vector<AbstractChemical*>::iterator chem_iter = p_chemical_vector.begin();
            chem_iter != p_chemical_vector.end();
            ++chem_iter, ++index)
    {
        AbstractChemical *p_system_chemical = dynamic_cast<AbstractChemical*>(*chem_iter);

        for(unsigned j=0; j<this->mNumberOfSubstrates; j++)
        {
            if(mpSubstrates[j] -> GetChemicalName()==p_system_chemical -> GetChemicalName())
            {
                forwardFlux *=  std::pow(currentSystemConc[index],mStoichSubstrates[j]);
                break;
            }
        }
    }

    SetReactionRate(mReactionRateConstant*forwardFlux*exp(-environmentConcentration));

}


std::string PreferredEnvironmentMassActionReaction::GetReactionType()
{
    return "PreferredEnvironmentMassActionReaction";
}

void PreferredEnvironmentMassActionReaction::ParseReactionInformation(std::string reaction_information, bool IsReversible=false)
{
    // read in from the reaction information string the necessary information for the basic multiplicative constant
    // and the species, may be multiple species but only one reaction rate, they can be in any order

    // assume that there must be a species
    size_t positionSpecies = reaction_information.find(this->mEnvironmentDelimiter);
    // initialse a temporary position, used if there are more than one species
    size_t tempPositionSpecies = reaction_information.find(this->mEnvironmentDelimiter);
    if(positionSpecies==std::string::npos)
    {
        // species not found
        std::cout<<"Error PreferredEnvironmentMassActionReaction::ParseReactionInformation: environment species not found"<<std::endl;
    }

    // assume that there must be reaction rate constant
    size_t positionRate = reaction_information.find(this->mIrreversibleRateName);

    if(positionRate==std::string::npos)
    {
        // rate not found
        std::cout<<"Error PreferredEnvironmentMassActionReaction::ParseReactionInformation: reaction rate not found"<<std::endl;
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
            reactionRate = atof(reaction_information.substr(positionRate+this->mIrreversibleRateName.size()+1,positionSpecies).c_str());
            reaction_information = reaction_information.substr(positionRate+this->mIrreversibleRateName.size()+1,std::string::npos);

        }
        else
        {
            // next data element is a environment species, test if the next value is an environment species or rate
            tempString = reaction_information.substr(positionSpecies+this->mEnvironmentDelimiter.size()+1,std::string::npos);

            // find the new positions of the delimiters
            tempPositionSpecies = tempString.find(this->mEnvironmentDelimiter); // if no more species, should be npos
            positionRate = tempString.find(this->mIrreversibleRateName); // if rate already taken, should be npos

            if(positionRate>tempPositionSpecies)
            {
                // next value is species 
                species = reaction_information.substr(positionSpecies+this->mEnvironmentDelimiter.size()+1,tempPositionSpecies);

            }
            else if(positionRate<tempPositionSpecies)
            {
                // next value is rate
                species = reaction_information.substr(positionSpecies+this->mEnvironmentDelimiter.size()+1,positionRate-1);
            }
            else
            {
                // both are npos so no other value
                species = reaction_information.substr(positionSpecies+this->mEnvironmentDelimiter.size()+1,std::string::npos);
            }
            
            environmentNames.push_back(species);
            
            reaction_information = reaction_information.substr(positionSpecies+this->mEnvironmentDelimiter.size()+1,std::string::npos);

        }

        // update the delimiter positions from the new substrings
        positionSpecies = reaction_information.find(this->mEnvironmentDelimiter);
        positionRate = reaction_information.find(this->mIrreversibleRateName);
        
        if(positionSpecies == std::string::npos && positionRate !=std::string::npos)
        {
            // there are no more species but the rate is still unknown and in the string
            reactionRate = atof(reaction_information.substr(positionRate+this->mIrreversibleRateName.size()+1,std::string::npos).c_str());

        }
    }

    this->SetReactionRateConstant(reactionRate);
    this->SetEnvironmentChemicalNames(environmentNames);
    this->SetNumberOfEnvironmentChemicals(environmentNames.size());
}

//void PreferredEnvironmentMassActionReaction::GiveCell(CellPtr p_cell)
//{
//    std::cout<<"[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[ mpCell update"<<std::endl;
//    mpCell = p_cell;
//}


#endif