#ifndef MASSACTIONCOUPLEDMEMBRANEREACTION_HPP_
#define MASSACTIONCOUPLEDMEMBRANEREACTION_HPP_

// general includes
#include <string>
#include <tuple>
#include <vector>
#include <cmath>

// custom includes
#include "AbstractChemical.hpp"
#include "AbstractReversibleMembraneReaction.hpp"

// membrane reaction class wherein the reaction rate is determined as proportional to the product of 
// all species concentrations both inside and outside of the membrane.

class MassActionCoupledMembraneReaction : public AbstractReversibleMembraneReaction
{
private:

    using AbstractReversibleMembraneReaction::UpdateReaction;
    using AbstractReversibleMembraneReaction::UpdateReactionRate;
    using AbstractReversibleMembraneReaction::GetReactionType;
    using AbstractReversibleMembraneReaction::ParseReactionInformation;

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

    MassActionCoupledMembraneReaction( 
                        std::vector<AbstractChemical*> bulkSubstrates = std::vector<AbstractChemical*>(),
                        std::vector<AbstractChemical*> bulkProducts = std::vector<AbstractChemical*>(),
                        std::vector<AbstractChemical*> cellSubstrates = std::vector<AbstractChemical*>(),
                        std::vector<AbstractChemical*> cellProducts = std::vector<AbstractChemical*>(),
                        std::vector<unsigned> stoichBulkSubstrates = std::vector<unsigned>(),
                        std::vector<unsigned> stoichBulkProducts = std::vector<unsigned>(),
                        std::vector<unsigned> stoichCellSubstrates = std::vector<unsigned>(),
                        std::vector<unsigned> stoichCellProducts = std::vector<unsigned>(),
                        bool IsGibbs = false,
                        bool IsReversible = true,
                        double ForwardReactionRateConstant = 1.0,
                        double ReverseReactionRateConstant = 1.0);

    virtual ~MassActionCoupledMembraneReaction()
    {
    };

    virtual void UpdateReaction();

    virtual void UpdateReactionRate(AbstractChemistry*, AbstractChemistry*, const std::vector<double>&, const std::vector<double>&);

    virtual std::string GetReactionType();

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
MassActionCoupledMembraneReaction::MassActionCoupledMembraneReaction(
                        std::vector<AbstractChemical*> bulkSubstrates,
                        std::vector<AbstractChemical*> bulkProducts,
                        std::vector<AbstractChemical*> cellSubstrates,
                        std::vector<AbstractChemical*> cellProducts,
                        std::vector<unsigned> stoichBulkSubstrates,
                        std::vector<unsigned> stoichBulkProducts,
                        std::vector<unsigned> stoichCellSubstrates,
                        std::vector<unsigned> stoichCellProducts,
                        bool IsGibbs,
                        bool IsReversible,
                        double ForwardReactionRateConstant, // takes the form of the gibbs energy
                        double ReverseReactionRateConstant)
    :   AbstractReversibleMembraneReaction(
                                bulkSubstrates,
                                bulkProducts,
                                cellSubstrates,
                                cellProducts,
                                stoichBulkSubstrates,
                                stoichBulkProducts,
                                stoichCellSubstrates,
                                stoichCellProducts,
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

void MassActionCoupledMembraneReaction::UpdateReaction()
{
    return;
}

void MassActionCoupledMembraneReaction::UpdateReactionRate(AbstractChemistry* bulkChemistry, AbstractChemistry* cellChemistry, const std::vector<double>& currentBulkConc, const std::vector<double>& currentCellConc)
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
        
        for(unsigned j=0; j<mNumberOfBulkSubstrates; j++)
        {
            if(mpBulkSubstrates[j] -> GetChemicalName()==p_system_chemical -> GetChemicalName())
            {
                forwardFlux *=  std::pow(currentBulkConc[index],mStoichBulkSubstrates[j]);
                break;
            }
        }
        for(unsigned j=0; j<mNumberOfBulkProducts; j++)
        {
            if(mpBulkProducts[j] -> GetChemicalName()==p_system_chemical -> GetChemicalName())
            {
                reverseFlux *=  std::pow(currentBulkConc[index],mStoichBulkProducts[j]);
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

        for(unsigned j=0; j<mNumberOfCellSubstrates; j++)
        {
            if(mpCellSubstrates[j] -> GetChemicalName()==p_system_chemical -> GetChemicalName())
            {
                forwardFlux *=  std::pow(currentCellConc[index],mStoichCellSubstrates[j]);
                break;
            }
        }
        for(unsigned j=0; j<mNumberOfCellProducts; j++)
        {
            if(mpCellProducts[j] -> GetChemicalName()==p_system_chemical -> GetChemicalName())
            {
                reverseFlux *=  std::pow(currentCellConc[index],mStoichCellProducts[j]);
                break;
            }
        }

    }  

    SetForwardReactionRate(kf*forwardFlux);
    SetReverseReactionRate(kr*reverseFlux);
}

std::string MassActionCoupledMembraneReaction::GetReactionType()
{
    return "MassActionCoupledMembrane";
}


void MassActionCoupledMembraneReaction::ParseReactionInformation(std::string reaction_information, bool IsReversible=false)
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
            AbstractReversibleMembraneReaction::ParseReactionInformation(reaction_information,IsReversible);
            SetForwardReactionRateConstant(GetForwardReactionRate());
            SetReverseReactionRateConstant(GetReverseReactionRate());
        }
        
    }

}

// member function sspecific to this reaction class

double MassActionCoupledMembraneReaction::CalculateReactionQuotient(AbstractChemistry* bulkChemistry, AbstractChemistry* cellChemistry, const std::vector<double>& currentBulkConc, const std::vector<double>& currentCellConc)
{
    double quotient = 1.0;
    
    if(mNumberOfBulkSubstrates ==0 || mNumberOfBulkProducts==0 || mNumberOfCellSubstrates==0 || mNumberOfCellProducts==0)
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

            for(unsigned j=0; j<mNumberOfBulkSubstrates; j++)
            {
                if(mpBulkSubstrates[j] -> GetChemicalName()==p_system_chemical -> GetChemicalName())
                {
                    substrates_concentrations *=  std::pow(currentBulkConc[index],mStoichBulkSubstrates[j]);
                    break;
                }
            }
            for(unsigned j=0; j<mNumberOfBulkProducts; j++)
            {
                if(mpBulkProducts[j] -> GetChemicalName()==p_system_chemical -> GetChemicalName())
                {
                    products_concentrations *=  std::pow(currentBulkConc[index],mStoichBulkProducts[j]);
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

            for(unsigned j=0; j<mNumberOfCellSubstrates; j++)
            {
                if(mpCellSubstrates[j] -> GetChemicalName()==p_system_chemical -> GetChemicalName())
                {
                    substrates_concentrations *=  std::pow(currentCellConc[index],mStoichCellSubstrates[j]);
                    break;
                }
            }
            for(unsigned j=0; j<mNumberOfCellProducts; j++)
            {
                if(mpCellProducts[j] -> GetChemicalName()==p_system_chemical -> GetChemicalName())
                {
                    products_concentrations *=  std::pow(currentCellConc[index],mStoichCellProducts[j]);
                    break;
                }
            }
        }

        quotient=products_concentrations/substrates_concentrations;
    }

    return quotient;
}

double MassActionCoupledMembraneReaction::GetGibbsFreeEnergy()
{
    return mGibbsFreeEnergy;
}

void MassActionCoupledMembraneReaction::SetGibbsFreeEnergy(double gibbsEnergy)
{
    mGibbsFreeEnergy = gibbsEnergy;
}

double MassActionCoupledMembraneReaction::DeltaGtoKf(double deltaG, double Kr)
{
    return Kr*exp(-deltaG/(mRkj*mTemp));
}

double MassActionCoupledMembraneReaction::KftodeltaG(double Kf, double Kr)
{
    return -mRkj*mTemp*log(Kf/Kr);
}

double MassActionCoupledMembraneReaction::CalculateGibbsFromQuotient(double G0, double Q)
{
    // Q reaction quotient
    return G0 - 8.3144598e-3*mTemp*log(Q);
}

void MassActionCoupledMembraneReaction::SetReactionTemperature(double temp)
{
    mTemp = temp;
}

double MassActionCoupledMembraneReaction::GetReactionTemperature()
{
    return mTemp;
}

double MassActionCoupledMembraneReaction::GetForwardReactionRateConstant()
{
    return mForwardReactionRateConstant;
}

void MassActionCoupledMembraneReaction::SetForwardReactionRateConstant(double forwardReactionRateConstant)
{
    mForwardReactionRateConstant = forwardReactionRateConstant;
}

double MassActionCoupledMembraneReaction::GetReverseReactionRateConstant()
{
    return mReverseReactionRateConstant;
}

void MassActionCoupledMembraneReaction::SetReverseReactionRateConstant(double reverseReactionRateConstant)
{
    mReverseReactionRateConstant = reverseReactionRateConstant;
}

void MassActionCoupledMembraneReaction::SetGibbsDelmiter(std::string gibbsDelimiter)
{
    mGibbsDelimiter = gibbsDelimiter;
}

std::string MassActionCoupledMembraneReaction::GetGibbsDelimiter()
{
    return mGibbsDelimiter;
}

#endif