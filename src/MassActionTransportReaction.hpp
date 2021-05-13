#ifndef MASSACTIONTRANSPORTREACTION_HPP_
#define MASSACTIONTRANSPORTREACTION_HPP_

// general includes
#include <string>
#include <tuple>
#include <vector>
#include <cmath>

// custom includes
#include "AbstractChemical.hpp"
#include "AbstractReversibleTransportReaction.hpp"

// mass action reaction where the concentrations are derived from either side of a membrane
// splitting the compartments. The mass action kinetics defined as being proportional
// to the product of the species concentrations either side of the membrane.  

class MassActionTransportReaction : public AbstractReversibleTransportReaction
{
private:

    using AbstractReversibleTransportReaction::UpdateReaction;
    using AbstractReversibleTransportReaction::UpdateReactionRate;
    using AbstractReversibleTransportReaction::GetReactionType;
    using AbstractReversibleTransportReaction::ParseReactionInformation;

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
    MassActionTransportReaction( 
                        std::vector<AbstractChemical*> bulkReactionSpecies = std::vector<AbstractChemical*>(),
                        std::vector<AbstractChemical*> cellReactionSpecies = std::vector<AbstractChemical*>(),
                        std::vector<unsigned> stoichBulk = std::vector<unsigned>(),
                        std::vector<unsigned> stoichCell = std::vector<unsigned>(),
                        bool IsGibbs = false,
                        bool IsReversible = true,
                        double ForwardReactionRateConstant = 1.0,
                        double ReverseReactionRateConstant = 1.0);

    virtual ~MassActionTransportReaction()
    {
    };

    virtual void UpdateReaction();

    virtual void UpdateReactionRate(AbstractChemistry*, AbstractChemistry*, const std::vector<double>&, const std::vector<double>&);

    virtual std::string GetReactionType()
    {
        return "MassActionTransportReaction";
    }

    virtual void ParseReactionInformation(std::string, bool);
    
    double GetGibbsFreeEnergy();

    void SetGibbsFreeEnergy(double);

    double CalculateReactionQuotient(AbstractChemistry*, AbstractChemistry*, const std::vector<double>&, const std::vector<double>&);

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
MassActionTransportReaction::MassActionTransportReaction(
                        std::vector<AbstractChemical*> bulkReactionSpecies,
                        std::vector<AbstractChemical*> cellReactionSpecies,
                        std::vector<unsigned> stoichBulk,
                        std::vector<unsigned> stoichCell,
                        bool IsGibbs,
                        bool IsReversible,
                        double ForwardReactionRateConstant, // takes the form of the gibbs energy
                        double ReverseReactionRateConstant)
    :   AbstractReversibleTransportReaction(
                                    bulkReactionSpecies,
                                    cellReactionSpecies,
                                    stoichBulk,
                                    stoichCell,
                                    ForwardReactionRateConstant,
                                    ReverseReactionRateConstant),
        mIsGibbs(IsGibbs),
        mIsReversible(IsReversible),
        mForwardReactionRateConstant(ForwardReactionRateConstant),
        mReverseReactionRateConstant(ReverseReactionRateConstant)
{

    if(mIsGibbs)
    {
        mGibbsFreeEnergy = ForwardReactionRateConstant;
    }
    if(!mIsReversible)
    {
        mReverseReactionRateConstant = 0.0;
    }
}

void MassActionTransportReaction::UpdateReaction()
{
    return;
}

void MassActionTransportReaction::UpdateReactionRate(AbstractChemistry* bulkChemistry, AbstractChemistry* cellChemistry, const std::vector<double>& currentBulkConc, const std::vector<double>& currentCellConc)
{
    // custom inheritence based on reversibility?

    double kf = mForwardReactionRateConstant;
    double kr = mReverseReactionRateConstant;
 
    if(mIsGibbs)
    {
        kf = DeltaGtoKf(mGibbsFreeEnergy, kr);
    }

    //double reaction_quotient = CalculateReactionQuotient(systemChemistry, currentSystemConc);

    double forwardFlux=1.0;
    double reverseFlux=1.0;

    std::vector<AbstractChemical*> p_chemical_vector = bulkChemistry -> rGetChemicalVector();
    unsigned index = 0;
    for(std::vector<AbstractChemical*>::iterator chem_iter = p_chemical_vector.begin();
            chem_iter != p_chemical_vector.end();
            ++chem_iter, ++index)
    {
        AbstractChemical *p_system_chemical = dynamic_cast<AbstractChemical*>(*chem_iter);

        for(unsigned j=0; j<mNumberOfBulkSpecies; j++)
        {
            if(mpBulkReactionSpecies[j] -> GetChemicalName()==p_system_chemical -> GetChemicalName())
            {
                forwardFlux *=  std::pow(currentBulkConc[index],mStoichBulk[j]);
                break;
            }
        }
    }

    p_chemical_vector = cellChemistry -> rGetChemicalVector();
    index = 0;
    for(std::vector<AbstractChemical*>::iterator chem_iter = p_chemical_vector.begin();
            chem_iter != p_chemical_vector.end();
            ++chem_iter, ++index)
    {
        AbstractChemical *p_system_chemical = dynamic_cast<AbstractChemical*>(*chem_iter);

        for(unsigned j=0; j<mNumberOfCellSpecies; j++)
        {
            if(mpCellReactionSpecies[j] -> GetChemicalName()==p_system_chemical -> GetChemicalName())
            {
                reverseFlux *=  std::pow(currentCellConc[index],mStoichCell[j]);
                break;
            }
        }
    }  

    SetForwardReactionRate(kf*forwardFlux);
    SetReverseReactionRate(kr*reverseFlux);
}

void MassActionTransportReaction::ParseReactionInformation(std::string reaction_information, bool IsReversible=false)
{
    //std::cout<<"MassActionTransportReaction::ParseReactionInformation - start"<<std::endl;
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
            AbstractReversibleTransportReaction::ParseReactionInformation(reaction_information,IsReversible);
            SetForwardReactionRateConstant(GetForwardReactionRate());
            SetReverseReactionRateConstant(GetReverseReactionRate());
        }
        
    }
    //std::cout<<"MassActionTransportReaction::ParseReactionInformation - end"<<std::endl;
}



// member function sspecific to this reaction class
double MassActionTransportReaction::CalculateReactionQuotient(AbstractChemistry* bulkChemistry, AbstractChemistry* cellChemistry, const std::vector<double>& currentBulkConc, const std::vector<double>& currentCellConc)
{
    double quotient = 1.0;
    
    if(mNumberOfBulkSpecies ==0 || mNumberOfCellSpecies==0)
    {
        quotient = 0.0;
    }else{
        double products_concentrations = 1.0;
        double substrates_concentrations = 1.0;
        // need to check against the concentration of each chemical in the system
        
        std::vector<AbstractChemical*> p_chemical_vector = bulkChemistry -> rGetChemicalVector();
        unsigned index = 0;
        for(std::vector<AbstractChemical*>::iterator chem_iter = p_chemical_vector.begin();
                chem_iter != p_chemical_vector.end();
                ++chem_iter, ++index)
        {
            AbstractChemical *p_system_chemical = dynamic_cast<AbstractChemical*>(*chem_iter);

            for(unsigned j=0; j<mNumberOfBulkSpecies; j++)
            {
                if(mpBulkReactionSpecies[j] -> GetChemicalName()==p_system_chemical -> GetChemicalName())
                {
                    substrates_concentrations *=  std::pow(currentBulkConc[index],mStoichBulk[j]);
                    break;
                }
            }
        }    
        p_chemical_vector = cellChemistry -> rGetChemicalVector();
        index = 0;
        for(std::vector<AbstractChemical*>::iterator chem_iter = p_chemical_vector.begin();
                chem_iter != p_chemical_vector.end();
                ++chem_iter, ++index)
        {
            AbstractChemical *p_system_chemical = dynamic_cast<AbstractChemical*>(*chem_iter);
            for(unsigned j=0; j<mNumberOfCellSpecies; j++)
            {
                if(mpCellReactionSpecies[j] -> GetChemicalName()==p_system_chemical -> GetChemicalName())
                {
                    products_concentrations *=  std::pow(currentCellConc[index],mStoichCell[j]);
                    break;
                }
            }
        }    

        quotient=products_concentrations/substrates_concentrations;
    }

    return quotient;
}

double MassActionTransportReaction::GetGibbsFreeEnergy()
{
    return mGibbsFreeEnergy;
}

void MassActionTransportReaction::SetGibbsFreeEnergy(double gibbsEnergy)
{
    mGibbsFreeEnergy = gibbsEnergy;
}

double MassActionTransportReaction::DeltaGtoKf(double deltaG, double Kr)
{
    return Kr*exp(-deltaG/(mRkj*mTemp));
}

double MassActionTransportReaction::KftodeltaG(double Kf, double Kr)
{
    return -mRkj*mTemp*log(Kf/Kr);
}

double MassActionTransportReaction::CalculateGibbsFromQuotient(double G0, double Q)
{
    // Q reaction quotient
    return G0 - 8.3144598e-3*mTemp*log(Q);
}

void MassActionTransportReaction::SetReactionTemperature(double temp)
{
    mTemp = temp;
}

double MassActionTransportReaction::GetReactionTemperature()
{
    return mTemp;
}

double MassActionTransportReaction::GetForwardReactionRateConstant()
{
    return mForwardReactionRateConstant;
}

void MassActionTransportReaction::SetForwardReactionRateConstant(double forwardReactionRateConstant)
{
    mForwardReactionRateConstant = forwardReactionRateConstant;
}

double MassActionTransportReaction::GetReverseReactionRateConstant()
{
    return mReverseReactionRateConstant;
}

void MassActionTransportReaction::SetReverseReactionRateConstant(double reverseReactionRateConstant)
{
    mReverseReactionRateConstant = reverseReactionRateConstant;
}

void MassActionTransportReaction::SetGibbsDelmiter(std::string gibbsDelimiter)
{
    mGibbsDelimiter = gibbsDelimiter;
}

std::string MassActionTransportReaction::GetGibbsDelimiter()
{
    return mGibbsDelimiter;
}

#endif