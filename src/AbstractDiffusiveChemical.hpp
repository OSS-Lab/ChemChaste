#ifndef ABSTRACTDIFFUSIVECHEMICAL_HPP_
#define ABSTRACTDIFFUSIVECHEMICAL_HPP_

// general includes
#include <string>
#include <tuple>
#include <vector>

// custom includes
#include "AbstractChemical.hpp"

// abstract property to contain information about diffusive chemical species

class AbstractDiffusiveChemical : public AbstractChemical
{
private:

    // vector containing the labelled sub-domains in which the chemical species diffuses
    std::vector<std::string> mDiffusiveDomains = std::vector<std::string>();
    // the respective valeus for the diffusive subdomains
    std::vector<double> mDiffusivities = std::vector<double>();

    unsigned mNumberOfDiffusiveDomains =0;

    bool mIsDiffusionSetup = false; 

public:

    AbstractDiffusiveChemical(  std::string chemicalName = "", 
                                double size = 0.0,
                                double mass = 0.0,
                                int valence = 0,
                                std::string diffusiveDomain = "",
                                double diffusivity = 0.0);

    virtual ~AbstractDiffusiveChemical()
    {
    };

    // update methods

    void SetupDiffusiveProperties();

    void AddDiffusiveDomain(std::string newDomain = "Bulk", double newDiffusivity =0.0);

    void AddDiffusiveDomainVector(std::vector<std::string>, std::vector<double>);

    // set methods

    void SetDiffusiveDomainVector(std::vector<std::string>);

    void SetChemicalDiffusivityVector(std::vector<double>);

    void SetNumberOfDiffusiveDomains(unsigned); 

    void SetIsDiffusionSetUp(bool);

    // get methods

    std::vector<std::string> GetDiffusiveDomainVector();

    std::vector<double> GetChemicalDiffusivityVector();

    unsigned GetNumberOfDiffusiveDomains();

    std::string GetDiffusiveDomainByIndex(unsigned);

    double GetChemicalDiffusivityByIndex(unsigned);

    // tracker method

    virtual std::string GetChemicalType()
    {
        return "AbstractDiffusiveChemical";
    };

};

AbstractDiffusiveChemical::AbstractDiffusiveChemical(  
                                std::string chemicalName, 
                                double size,
                                double mass,
                                int valence,
                                std::string diffusiveDomain,
                                double diffusivity)
    : AbstractChemical( chemicalName,
                        size,
                        mass,
                        valence)
{
    // set up the member variables
    std::vector<std::string> mDiffusiveDomain =std::vector<std::string>();
    std::vector<double> mDiffusivity = std::vector<double>();
    mIsDiffusionSetup =true;
    if(diffusiveDomain != "")
    {
        // add the first input domain as an elements of the data vectors
        mDiffusiveDomain.push_back(diffusiveDomain);
        mDiffusivity.push_back(diffusivity);
        mNumberOfDiffusiveDomains=mNumberOfDiffusiveDomains+1;
    }
    
}

void AbstractDiffusiveChemical::SetupDiffusiveProperties()
{
    // if constructor is not called, method to set up the member variables
    std::vector<std::string> diffusiveDomain =std::vector<std::string>();
    std::vector<double> diffusivity = std::vector<double>();

    // store these new memeber variables
    SetDiffusiveDomainVector(diffusiveDomain);
    SetChemicalDiffusivityVector(diffusivity);
    SetIsDiffusionSetUp(true);
    
}

void AbstractDiffusiveChemical::AddDiffusiveDomain(std::string newDomain, double newDiffusivity)
{
    // add a new domain and diffusivity to the chemical's diffusion properties, check for duplicate domains
    unsigned numberOfDiffusiveDomains = GetNumberOfDiffusiveDomains();

    std::vector<std::string> diffusiveDomains = GetDiffusiveDomainVector();
    std::vector<double> diffusivities = GetChemicalDiffusivityVector();

    // if the domains are empty, add the domain as a new property
    if(numberOfDiffusiveDomains==0)
    {

        diffusiveDomains.push_back(newDomain);
        diffusivities.push_back(newDiffusivity);
        numberOfDiffusiveDomains = numberOfDiffusiveDomains + 1;
        SetDiffusiveDomainVector(diffusiveDomains);

        SetChemicalDiffusivityVector(diffusivities);

        SetNumberOfDiffusiveDomains(numberOfDiffusiveDomains); 
    }
    else
    {
        // there are already domains in which the chemical diffuses, check that the new domain to add is not a duplicate
        bool isNewDomain = true;
        //check if domain is new
        for(unsigned domain_index=0; domain_index < numberOfDiffusiveDomains; domain_index++)
        {
            if(diffusiveDomains[domain_index] == newDomain)
            {
                isNewDomain = false;
            }
        }
        if(isNewDomain)
        {
            // if the candidate domain is new the update diffusion properties
            diffusiveDomains.push_back(newDomain);
            diffusivities.push_back(newDiffusivity);
            numberOfDiffusiveDomains = mNumberOfDiffusiveDomains + 1;
            SetDiffusiveDomainVector(diffusiveDomains);

            SetChemicalDiffusivityVector(diffusivities);

            SetNumberOfDiffusiveDomains(numberOfDiffusiveDomains); 
        }
    }
}

void AbstractDiffusiveChemical::AddDiffusiveDomainVector(std::vector<std::string> domainLabelVector, std::vector<double> diffusivityVector)
{
    if(domainLabelVector.size() == diffusivityVector.size())
    {
        for(unsigned domain_index = 0; domain_index<domainLabelVector.size(); domain_index++)
        {
            AddDiffusiveDomain(domainLabelVector[domain_index],diffusivityVector[domain_index]);
        }

    }
    else
    {
        std::cout<<"Error: AbstractDiffusiveChemical::AddDiffusiveDomainVector:   Domain vector and diffusivity vector not of equla size"<<std::endl;
    }
    
}


std::vector<std::string> AbstractDiffusiveChemical::GetDiffusiveDomainVector()
{
    return mDiffusiveDomains;
}

std::vector<double> AbstractDiffusiveChemical::GetChemicalDiffusivityVector()
{
    return mDiffusivities;
}


void AbstractDiffusiveChemical::SetDiffusiveDomainVector(std::vector<std::string> diffusiveDomains)
{
    mDiffusiveDomains = diffusiveDomains;
}

void AbstractDiffusiveChemical::SetChemicalDiffusivityVector(std::vector<double> diffusivities)
{
    mDiffusivities = diffusivities;
}

void AbstractDiffusiveChemical::SetNumberOfDiffusiveDomains(unsigned numberOfDomains)
{
    mNumberOfDiffusiveDomains = numberOfDomains;
}

void AbstractDiffusiveChemical::SetIsDiffusionSetUp(bool isSetup)
{
    mIsDiffusionSetup = isSetup;
}

unsigned AbstractDiffusiveChemical::GetNumberOfDiffusiveDomains()
{
    return mNumberOfDiffusiveDomains;
}

std::string AbstractDiffusiveChemical::GetDiffusiveDomainByIndex(unsigned index)
{
    if(index<mNumberOfDiffusiveDomains)
    {
        return mDiffusiveDomains[index];
    }
    else
    {    
        std::cout<<"Error: AbstractDiffusiveChemical::GetDiffusiveDomainByIndex:   index out of bounds"<<std::endl;
        return "Error";
    }
    
}

double AbstractDiffusiveChemical::GetChemicalDiffusivityByIndex(unsigned index)
{
    if(index<mNumberOfDiffusiveDomains)
    {
        return mDiffusivities[index];
    }
    else
    {    
        std::cout<<"Error: AbstractDiffusiveChemical::GetChemicalDiffusivityByIndex:   index out of bounds"<<std::endl;
        return 0.0;
    }
}

#endif