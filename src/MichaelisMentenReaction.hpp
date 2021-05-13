#ifndef MICHAELISMENTENREACTION_HPP_
#define MICHAELISMENTENREACTION_HPP_

// general includes
#include <string>
#include <tuple>
#include <vector>
#include <cmath>

// custom includes
#include "AbstractChemical.hpp"
#include "AbstractChemistry.hpp"
#include "AbstractReaction.hpp"
#include "SpectatorDependentReaction.hpp"

// reaction class for performing a Michaelis-Menten kinetic reation
// E + S <==> ES -> E + x
// rate(x) = dx/dt = Vmax*(x/(Km + x)) = kcat*[Enzyme]*(x/(Km + x))
// There may be more than one enzyme, should there be a mechanical need,
// and herein enzymes are referred to as "Spectators" due to class hierarchy 
// referring to the enzymes not being consumed/produced over the full mechanism.

class MichaelisMentenReaction : public SpectatorDependentReaction
{
protected:

    std::string mMichaelisConstantDelimiter = "Km =";

    double mKm = 0;

    std::string mCatalyticConstantDelimiter = "kcat =";
    // kcat is set as the reaction rate

    std::string mEnzymeDelimiter = "Enzyme =";

private:
    using SpectatorDependentReaction::UpdateReactionRate;
    using SpectatorDependentReaction::GetReactionType;
    using SpectatorDependentReaction::ParseReactionInformation;

public:

    // constructor
    MichaelisMentenReaction(    std::vector<AbstractChemical*> substrates = std::vector<AbstractChemical*>(),
                                std::vector<AbstractChemical*> products = std::vector<AbstractChemical*>(),
                                std::vector<unsigned> stoichSubstrates = std::vector<unsigned>(),
                                std::vector<unsigned> stoichProducts = std::vector<unsigned>(),
                                AbstractChemistry* p_systemChemistry = new AbstractChemistry(),
                                double reactionRate = 1.0
                                );

    // destructor
    virtual ~MichaelisMentenReaction()
    {
    };

    virtual void UpdateReactionRate(AbstractChemistry*, const std::vector<double>&);

    virtual std::string GetReactionType();

    virtual void ParseReactionInformation(std::string, bool);

    // class specific methods

    void SetKm(double); 

    double GetKm();
};

MichaelisMentenReaction::MichaelisMentenReaction( 
                                std::vector<AbstractChemical*> substrates,
                                std::vector<AbstractChemical*> products,
                                std::vector<unsigned> stoichSubstrates,
                                std::vector<unsigned> stoichProducts,
                                AbstractChemistry* p_systemChemistry,
                                double reactionRate)
                            : SpectatorDependentReaction(
                                substrates,
                                products,
                                stoichSubstrates,
                                stoichProducts,
                                p_systemChemistry,
                                reactionRate
                            )
                            {
                                // recast upstream class text delimiter to version more in line with Michaeli-Menten 
                                // terminology
                                SetSpectatorDelimiter(mEnzymeDelimiter);
                                SetIrreversibleDelimiter(mCatalyticConstantDelimiter);
                            }


void MichaelisMentenReaction::UpdateReactionRate(AbstractChemistry* systemChemistry, const std::vector<double>& currentSystemConc)
{
    // multiply the reaction rate constant by the product of the spectator set concentrations
    SpectatorDependentReaction::UpdateReactionRate(systemChemistry,currentSystemConc);
    double reactionRate = GetReactionRate();
    double substrateConcentration=0.0;

    // apply non-linear term (x/(Km + x))
    std::vector<AbstractChemical*> p_chemical_vector = systemChemistry -> rGetChemicalVector();
    unsigned index = 0;
    for(std::vector<AbstractChemical*>::iterator chem_iter = p_chemical_vector.begin();
            chem_iter != p_chemical_vector.end();
            ++chem_iter, ++index)
    {
        AbstractChemical *p_system_chemical = dynamic_cast<AbstractChemical*>(*chem_iter);

        for(unsigned j=0; j<mNumberOfSubstrates; j++)
        {
            if(mpSubstrates[j] -> GetChemicalName()==p_system_chemical -> GetChemicalName())
            {
                substrateConcentration *=  std::pow(currentSystemConc[index],mStoichSubstrates[j]);
                break;
            }
        }         
    }  

    reactionRate *= substrateConcentration/(mKm + substrateConcentration);

    SetReactionRate(reactionRate);
}

std::string MichaelisMentenReaction::GetReactionType()
{
    return "MichaelisMentenReaction";
}

void MichaelisMentenReaction::ParseReactionInformation(std::string reaction_information, bool IsReversible=false)
{
    // read in from the reaction information string the necessary information for the Michaelis-Menten kinetics;
    // a basic multiplicative rate constant, Michaelis-Menten constant and enzyme (catalytic) identifiers (names)
    // There may be multiple spectator species but only one catalytic reaction rate or 
    // Michaelis-Menten constant, they can be in any order
    //std::cout<<"MichaelisMentenReaction::ParseReactionInformation( - start"<<std::endl;
    // form the vector of delimiters
    unsigned numberOfDelimiters = 3;
    std::vector<std::string> delimiterVector(numberOfDelimiters,"");
    delimiterVector[0] = mSpectatorDelimiter;
    delimiterVector[1] = mCatalyticConstantDelimiter;
    delimiterVector[2] = mMichaelisConstantDelimiter;

    // assume that there must be a spectator
    size_t positionSpectator = reaction_information.find(mSpectatorDelimiter);
    // assume that there must be reaction rate constant
    size_t positionRate = reaction_information.find(mCatalyticConstantDelimiter);
    // assume that there must be a michaelis menten constant
    size_t positionMichaelisMenten = reaction_information.find(mMichaelisConstantDelimiter);

    // check that each of the delimiters exists
    if(positionSpectator==std::string::npos)
    {
        // spectator not found
        std::cout<<"Error MichaelisMentenReaction::ParseReactionInformation: Enzyme not found"<<std::endl;
    }

    if(positionRate==std::string::npos)
    {
        // spectator not found
        std::cout<<"Error MichaelisMentenReaction::ParseReactionInformation: reaction rate not found: "<< mIrreversibleRateName<<std::endl;
    }

    if(positionMichaelisMenten==std::string::npos)
    {
        // spectator not found
        std::cout<<"Error MichaelisMentenReaction::ParseReactionInformation: Michaelis-Menten constant not found"<<std::endl;
    }


    double reactionRate;
    double Km;
    std::string spectator;
    std::vector<std::string> spectatorNames;
    size_t delimiterPosition = 0;
    size_t nextDelimiterPosition = 0;
    std::string tempString = reaction_information;
    std::string valueString = "";

    unsigned delimiterIndex=0;

    while(delimiterIndex != numberOfDelimiters)
    {
        // find the position of the first delimiter
        delimiterIndex = FindIndexOfLastDelimiterPosition(delimiterVector, tempString);
        if(delimiterIndex == numberOfDelimiters)
        {
            break;
        }
        delimiterPosition = tempString.rfind(delimiterVector[delimiterIndex]);
    
        valueString = reaction_information.substr(delimiterPosition+delimiterVector[delimiterIndex].size()+1,std::string::npos);
        tempString = tempString.erase(delimiterPosition,std::string::npos);

        switch(delimiterIndex)
        {
            // part specific to the paresing the particular reaction string

            //delimiterVector[0] = mSpectatorDelimiter;
            //delimiterVector[1] = mIrreversibleRateName;
            //delimiterVector[2] = mMichaelisConstantDelimiter;
            
            case 0:
                // spectator
                spectatorNames.push_back(valueString);
                break;
            
            case 1:
                // reaction rate
                reactionRate = atof(valueString.c_str());
                break;

            case 2:
                // Michaelis-Menten constant
                Km = atof(valueString.c_str());
                break;

            default:
                break;
        }
        
        // replace the string with the substring
        reaction_information = tempString;
        //delimiterIndex++;
    }
    

    SetReactionRateConstant(reactionRate);
    SetKm(Km);
    SetSpectatorNames(spectatorNames);
    SetNumberOfSpectators(spectatorNames.size());
    UpdateSystemChemistry(spectatorNames);
    //std::cout<<"MichaelisMentenReaction::ParseReactionInformation( - end"<<std::endl;
}

// class specific methods

void MichaelisMentenReaction::SetKm(double Km)
{
    mKm = Km;
}

double MichaelisMentenReaction::GetKm()
{
    return mKm;
}



#endif