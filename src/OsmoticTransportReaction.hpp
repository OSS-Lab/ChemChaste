#ifndef OSMOTICTRANSPORTREACTION_HPP_
#define OSMOTICTRANSPORTREACTION_HPP_

// general includes
#include <string>
#include <tuple>
#include <vector>
#include <cmath>

// custom includes
#include "AbstractChemical.hpp"
#include "AbstractReversibleTransportReaction.hpp"

// mass action reactions are genreally reversible

class OsmoticTransportReaction : public AbstractReversibleTransportReaction
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
    OsmoticTransportReaction( 
                        std::vector<AbstractChemical*> bulkReactionSpecies = std::vector<AbstractChemical*>(),
                        std::vector<AbstractChemical*> cellReactionSpecies = std::vector<AbstractChemical*>(),
                        std::vector<unsigned> stoichBulk = std::vector<unsigned>(),
                        std::vector<unsigned> stoichCell = std::vector<unsigned>(),
                        bool IsGibbs = false,
                        bool IsReversible = true,
                        double ForwardReactionRateConstant = 1.0,
                        double ReverseReactionRateConstant = 1.0);

    virtual ~OsmoticTransportReaction()
    {
    };

    virtual void UpdateReaction();

    virtual void UpdateReactionRate(AbstractChemistry*, const std::vector<double>&, const std::vector<double>&);

    virtual std::string GetReactionType()
    {
        return "OsmoticTransportReaction";
    }

    virtual void ParseReactionInformation(std::string, bool);

    double GetGibbsFreeEnergy();

    void SetGibbsFreeEnergy(double);

    double CalculateReactionQuotient(AbstractChemistry*, const std::vector<double>&, const std::vector<double>&);

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
OsmoticTransportReaction::OsmoticTransportReaction(
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

void OsmoticTransportReaction::UpdateReaction()
{
    return;
}

void OsmoticTransportReaction::UpdateReactionRate(AbstractChemistry* systemChemistry, const std::vector<double>& currentBulkConc, const std::vector<double>& currentCellConc)
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

    std::vector<AbstractChemical*> p_chemical_vector = systemChemistry -> rGetChemicalVector();
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

void OsmoticTransportReaction::ParseReactionInformation(std::string reaction_information, bool IsReversible=false)
{
    //std::cout<<"Parse mass action"<<std::endl;
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

}



// member function sspecific to this reaction class

double OsmoticTransportReaction::CalculateReactionQuotient(AbstractChemistry* systemChemistry, const std::vector<double>& currentBulkConc, const std::vector<double>& currentCellConc)
{
    double quotient = 1.0;
    
    if(mNumberOfBulkSpecies ==0 || mNumberOfCellSpecies==0)
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

            for(unsigned j=0; j<mNumberOfBulkSpecies; j++)
            {
                if(mpBulkReactionSpecies[j] -> GetChemicalName()==p_system_chemical -> GetChemicalName())
                {
                    substrates_concentrations *=  std::pow(currentBulkConc[index],mStoichBulk[j]);
                    break;
                }
            }
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

double OsmoticTransportReaction::GetGibbsFreeEnergy()
{
    return mGibbsFreeEnergy;
}

void OsmoticTransportReaction::SetGibbsFreeEnergy(double gibbsEnergy)
{
    mGibbsFreeEnergy = gibbsEnergy;
}

double OsmoticTransportReaction::DeltaGtoKf(double deltaG, double Kr)
{
    return Kr*exp(-deltaG/(mRkj*mTemp));
}

double OsmoticTransportReaction::KftodeltaG(double Kf, double Kr)
{
    return -mRkj*mTemp*log(Kf/Kr);
}

double OsmoticTransportReaction::CalculateGibbsFromQuotient(double G0, double Q)
{
    // Q reaction quotient
    return G0 - 8.3144598e-3*mTemp*log(Q);
}

void OsmoticTransportReaction::SetReactionTemperature(double temp)
{
    mTemp = temp;
}

double OsmoticTransportReaction::GetReactionTemperature()
{
    return mTemp;
}

double OsmoticTransportReaction::GetForwardReactionRateConstant()
{
    return mForwardReactionRateConstant;
}

void OsmoticTransportReaction::SetForwardReactionRateConstant(double forwardReactionRateConstant)
{
    mForwardReactionRateConstant = forwardReactionRateConstant;
}

double OsmoticTransportReaction::GetReverseReactionRateConstant()
{
    return mReverseReactionRateConstant;
}

void OsmoticTransportReaction::SetReverseReactionRateConstant(double reverseReactionRateConstant)
{
    mReverseReactionRateConstant = reverseReactionRateConstant;
}

void OsmoticTransportReaction::SetGibbsDelmiter(std::string gibbsDelimiter)
{
    mGibbsDelimiter = gibbsDelimiter;
}

std::string OsmoticTransportReaction::GetGibbsDelimiter()
{
    return mGibbsDelimiter;
}



#endif