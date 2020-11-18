#ifndef MASSACTIONREACTION_HPP_
#define MASSACTIONREACTION_HPP_

// general includes
#include <string>
#include <tuple>
#include <vector>
#include <cmath>

// custom includes
#include "AbstractChemical.hpp"
#include "AbstractReversibleReaction.hpp"

// mass action reactions are generally reversible

class MassActionReaction : public AbstractReversibleReaction
{
private:

    using AbstractReversibleReaction::UpdateReaction;
    using AbstractReversibleReaction::UpdateReactionRate;
    using AbstractReversibleReaction::GetReactionType;
    using AbstractReversibleReaction::ParseReactionInformation;

    double mReactionGibbs;
    
    double mRkj = 8.3144598e-3;
    double mTemp = 300; // kelvin

    bool mIsGibbs;

    bool mIsReversible;

    double mGibbsFreeEnergy = 0.0;

    double mForwardReactionRateConstant;

    double mReverseReactionRateConstant;

    std::string mGibbsDelimiter = "deltaG =";

public:
    MassActionReaction( std::vector<AbstractChemical*> substrates = std::vector<AbstractChemical*>(),
                        std::vector<AbstractChemical*> products = std::vector<AbstractChemical*>(),
                        std::vector<unsigned> stoichSubstrates = std::vector<unsigned>(),
                        std::vector<unsigned> stoichProducts = std::vector<unsigned>(),
                        bool IsGibbs = false,
                        bool IsReversible = true,
                        double ForwardReactionRateConstant = 1.0,
                        double ReverseReactionRateConstant = 1.0);

    virtual ~MassActionReaction()
    {
    };

    virtual void UpdateReaction();

    virtual void UpdateReactionRate(AbstractChemistry*, const std::vector<double>&);

    virtual std::string GetReactionType();

    virtual void ParseReactionInformation(std::string, bool);

    double GetGibbsFreeEnergy();

    void SetGibbsFreeEnergy(double);

    double CalculateReactionQuotient(AbstractChemistry* , const std::vector<double>& );

    double DeltaGtoKf(double deltaG, double Kr = 1.0);

    double KftodeltaG(double Kf, double Kr = 1.0);

    double CalculateGibbsFromQuotient(double G0, double Q);

    void SetReactionTemperature(double);

    double GetReactionTemperature();


    double GetForwardReactionRateConstant();

    void SetForwardReactionRateConstant(double);
    
    double GetReverseReactionRateConstant();

    void SetReverseReactionRateConstant(double);

    void SetGibbsDelmiter(std::string);

    std::string GetGibbsDelimiter();

};

// constructor
MassActionReaction::MassActionReaction(std::vector<AbstractChemical*> substrates,
                        std::vector<AbstractChemical*> products,
                        std::vector<unsigned> stoichSubstrates,
                        std::vector<unsigned> stoichProducts,
                        bool IsGibbs,
                        bool IsReversible,
                        double ForwardReactionRateConstant, // takes the form of the gibbs energy
                        double ReverseReactionRateConstant)
    :   AbstractReversibleReaction(substrates,
                                   products,
                                   stoichSubstrates,
                                   stoichProducts,
                                   ForwardReactionRateConstant,
                                   ReverseReactionRateConstant),
        mIsGibbs(IsGibbs),
        mIsReversible(IsReversible),
        mForwardReactionRateConstant(ForwardReactionRateConstant),
        mReverseReactionRateConstant(ReverseReactionRateConstant)
{
    //std::cout<<"MasActionConstructor: "<<ForwardReactionRateConstant<<" "<<ReverseReactionRateConstant<<std::endl;
    if(mIsGibbs)
    {
        mGibbsFreeEnergy = ForwardReactionRateConstant;
    }
    if(!mIsReversible)
    {
        mReverseReactionRateConstant = 0.0;
    }
}

void MassActionReaction::UpdateReaction()
{
    return;
}

void MassActionReaction::UpdateReactionRate(AbstractChemistry* systemChemistry, const std::vector<double>& currentSystemConc)
{

    double kf = mForwardReactionRateConstant;
    double kr = mReverseReactionRateConstant;

    if(mIsGibbs)
    {
        kf = DeltaGtoKf(mGibbsFreeEnergy, kr);
    }

    //double reaction_quotient = CalculateReactionQuotient(systemChemistry, currentSystemConc);

    double forwardFlux=1.0;
    double reverseFlux=1.0;

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
                forwardFlux *=  std::pow(currentSystemConc[index],mStoichSubstrates[j]);
                break;
            }
        }
        for(unsigned j=0; j<mNumberOfProducts; j++)
        {
            if(mpProducts[j] -> GetChemicalName()==p_system_chemical -> GetChemicalName())
            {
                reverseFlux *=  std::pow(currentSystemConc[index],mStoichProducts[j]);
                break;
            }
        }
    }  

    SetForwardReactionRate(kf*forwardFlux);
    SetReverseReactionRate(kr*reverseFlux);
}

std::string MassActionReaction::GetReactionType()
{
    return "MassActionReaction";
}


void MassActionReaction::ParseReactionInformation(std::string reaction_information, bool IsReversible=false)
{
    
    mIsReversible = IsReversible;

    if(!mIsReversible)
    {
       
        if(reaction_information.find(mIrreversibleRateName) != std::string::npos)
        {

            size_t pos= reaction_information.find(mIrreversibleRateName);
            mForwardReactionRateConstant = atof(reaction_information.substr(pos+mIrreversibleRateName.size()+1,std::string::npos).c_str());
            mReverseReactionRateConstant = 0.0;

        }
    }
    else
    {  
        if(reaction_information.find(mGibbsDelimiter) != std::string::npos)
        {
            size_t pos= reaction_information.find(mGibbsDelimiter);
            std::cout<<"Gibbs raw: "<<reaction_information.substr(pos+mGibbsDelimiter.size()+1,std::string::npos).c_str()<<std::endl;
            mGibbsFreeEnergy = atof(reaction_information.substr(pos+mGibbsDelimiter.size()+1,std::string::npos).c_str());
            mIsGibbs = true;
            std::cout<<"Gibbs translated: "<<mGibbsFreeEnergy<<std::endl;
        }
        else
        {
            
            //abstractReversibleReaction set the reaction rates not constants so update
            AbstractReversibleReaction::ParseReactionInformation(reaction_information,mIsReversible);
            SetForwardReactionRateConstant(GetForwardReactionRate());
            SetReverseReactionRateConstant(GetReverseReactionRate());
        }
        
    }
    //std::cout<<"MasActionCorrected: "<<mForwardReactionRateConstant<<" "<<mReverseReactionRateConstant<<std::endl;
}



// member function sspecific to this reaction class

double MassActionReaction::CalculateReactionQuotient(AbstractChemistry* systemChemistry, const std::vector<double>& currentSystemConc)
{
    double quotient = 1.0;
    
    if(mNumberOfSubstrates ==0 || mNumberOfProducts==0)
    {
        quotient = 0.0;
    }else{
        double products_concentrations = 1.0;
        double substrates_concentrations = 1.0;
        // need to check against the concentration of each chemical in the system
        
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
                    substrates_concentrations *=  std::pow(currentSystemConc[index],mStoichSubstrates[j]);
                    break;
                }
            }
            for(unsigned j=0; j<mNumberOfProducts; j++)
            {
                if(mpProducts[j] -> GetChemicalName()==p_system_chemical -> GetChemicalName())
                {
                    products_concentrations *=  std::pow(currentSystemConc[index],mStoichProducts[j]);
                    break;
                }
            }
        }    

        quotient=products_concentrations/substrates_concentrations;
    }

    return quotient;
}

double MassActionReaction::GetGibbsFreeEnergy()
{
    return mGibbsFreeEnergy;
}

void MassActionReaction::SetGibbsFreeEnergy(double gibbsEnergy)
{
    mGibbsFreeEnergy = gibbsEnergy;
}

double MassActionReaction::DeltaGtoKf(double deltaG, double Kr)
{
    return Kr*exp(-deltaG/(mRkj*mTemp));
}

double MassActionReaction::KftodeltaG(double Kf, double Kr)
{
    return -mRkj*mTemp*log(Kf/Kr);
}

double MassActionReaction::CalculateGibbsFromQuotient(double G0, double Q)
{
    // Q reaction quotient
    return G0 - 8.3144598e-3*mTemp*log(Q);
}

void MassActionReaction::SetReactionTemperature(double temp)
{
    mTemp = temp;
}

double MassActionReaction::GetReactionTemperature()
{
    return mTemp;
}

double MassActionReaction::GetForwardReactionRateConstant()
{
    return mForwardReactionRateConstant;
}

void MassActionReaction::SetForwardReactionRateConstant(double forwardReactionRateConstant)
{
    mForwardReactionRateConstant = forwardReactionRateConstant;
}

double MassActionReaction::GetReverseReactionRateConstant()
{
    return mReverseReactionRateConstant;
}

void MassActionReaction::SetReverseReactionRateConstant(double reverseReactionRateConstant)
{
    mReverseReactionRateConstant = reverseReactionRateConstant;
}

void MassActionReaction::SetGibbsDelmiter(std::string gibbsDelimiter)
{
    mGibbsDelimiter = gibbsDelimiter;
}

std::string MassActionReaction::GetGibbsDelimiter()
{
    return mGibbsDelimiter;
}


#endif