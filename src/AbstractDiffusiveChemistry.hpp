#ifndef ABSTRACTDIFFUSIVECHEMISTRY_HPP_
#define ABSTRACTDIFFUSIVECHEMISTRY_HPP_

#include "AbstractChemistry.hpp"
#include "AbstractDiffusiveChemical.hpp"

#include <string>
#include <tuple>
#include <vector>

// abstract property to contain information about the collection of chemical species in the system

class AbstractDiffusiveChemistry : public AbstractChemistry
{
private:
    unsigned mNumberDiffusiveChemicals;
    
    using AbstractChemistry::CheckChemical;
    using AbstractChemistry::AddChemical;
    using AbstractChemistry::UpdateChemicalVectors;
    using AbstractChemistry::AddChemistry;

protected:
    std::vector<std::string> mDiffusionChemicalNames;

    // maybe these aren't necessary and the data structure of the AbstractDiffusiveChemicals is all that is necessary?
    std::vector<std::vector<double>> mDiffusivityMatrix;
    std::vector<std::vector<std::string>> mDiffusionDomainMatrix;

public:

    AbstractDiffusiveChemistry();

    virtual ~AbstractDiffusiveChemistry()
    {
    };

    virtual void AddChemical(AbstractDiffusiveChemical*);

    virtual bool CheckChemical(AbstractDiffusiveChemical*); 

    virtual void UpdateChemicalVectors(AbstractDiffusiveChemical*);

    virtual void UpdateDomainVector(AbstractDiffusiveChemical*);

    virtual void AddChemistry(AbstractDiffusiveChemistry*);

    virtual std::string GetChemistryType();

    void SetDiffusivityMatrix(std::vector<std::vector<double>> );

    std::vector<std::vector<double>> GetDiffusivityMatrix();

    std::vector<double> GetDiffusivityVectorByIndex(unsigned);

    double GetDiffusivityValueByChemicalAndDomainName(std::string, std::string);

    double GetDiffusivityValueByChemicalName(std::string);

    void SetDiffusionDomainsMatrix(std::vector<std::vector<std::string>>);

    std::vector<std::vector<std::string>> GetDiffusionDomainsMatrix();

    std::vector<std::string> GetDiffusionDomainsVectorByIndex(unsigned);

    unsigned GetNumberDiffusiveChemicals();

    std::string GetDiffusiveChemicalNamesByIndex(unsigned);

    unsigned GetDiffusiveChemicalIndexByName(std::string);


};

AbstractDiffusiveChemistry::AbstractDiffusiveChemistry()
{
    // chemical species (domain vector)
    mDiffusivityMatrix = std::vector<std::vector<double>>();
    mDiffusionDomainMatrix = std::vector<std::vector<std::string>>();
    mNumberDiffusiveChemicals = 0;
}

void  AbstractDiffusiveChemistry::AddChemical(AbstractDiffusiveChemical* chemical)
{
    // add a chemical to the diffusive chemistry, check whether the diffusive domain of the chemical is new 
    // if new then add to the chemical diffusion properties
    bool isNewChem = false;

    if(mpChemicalVector.empty())
    {
        isNewChem = true;
    }
    else
    {
        isNewChem = CheckChemical(chemical);
    }
    if(isNewChem)
    {
        // if candidate is a new chemical then naturally the domain would be new for the chemical
        UpdateChemicalVectors(chemical);
    }
    else
    {
        // check whether the existing chemical is in a new domain
        UpdateDomainVector(chemical);
        
    }
}

// same structure as the AbstractChemical, remove?
bool AbstractDiffusiveChemistry::CheckChemical(AbstractDiffusiveChemical* chemical)
{
    // function to return true if chemical is already in class chemical vector

    bool isNewChem  = true;

    std::vector<AbstractChemical*> :: iterator chemical_iter;
    for(chemical_iter = mpChemicalVector.begin();
        chemical_iter != mpChemicalVector.end();
        ++chemical_iter )
    {
        AbstractDiffusiveChemical* query_chemical = dynamic_cast<AbstractDiffusiveChemical*>(*chemical_iter);
        if(chemical -> GetChemicalName() ==  query_chemical -> GetChemicalName())
        {   
            isNewChem = false;
            break;
        }
    }
    return isNewChem;
}

void AbstractDiffusiveChemistry::UpdateDomainVector(AbstractDiffusiveChemical* chemical)
{
    // function to return true if chemical is already in class chemical vector
    std::vector<AbstractChemical*> :: iterator chemical_iter;
    for(chemical_iter = mpChemicalVector.begin();
        chemical_iter != mpChemicalVector.end();
        ++chemical_iter )
    {
        // assume the new chemical is known in the system
        AbstractDiffusiveChemical* query_chemical = dynamic_cast<AbstractDiffusiveChemical*>(*chemical_iter);
        if(chemical -> GetChemicalName() ==  query_chemical -> GetChemicalName())
        {   
            unsigned number_candidate_domains = chemical -> GetNumberOfDiffusiveDomains();

            // chemical whose domain is to be tested has been found in the chemical data structure
            for(unsigned candidate_domain_index=0; candidate_domain_index<number_candidate_domains; candidate_domain_index++)
            {
                // test each of the potential new domains
                // assume to be new
                bool isNewDomain  = true;
                unsigned indexRecord=0;
                
                for(unsigned domain_index=0; domain_index<query_chemical ->GetNumberOfDiffusiveDomains(); domain_index++)
                {
                    // for each of the current existing domains
                    std::string existing_query_domain = query_chemical->GetDiffusiveDomainByIndex(domain_index);
                
                    if(chemical -> GetDiffusiveDomainByIndex(candidate_domain_index) ==  existing_query_domain)
                    {
                        // new domain has already be recorded
                        isNewDomain = false;
                        indexRecord = candidate_domain_index;
                        break;
                    }
                }
                if(isNewDomain)
                {
                    // then domain is new and add to DiffusiveChemical 
                    query_chemical ->AddDiffusiveDomain(chemical -> GetDiffusiveDomainByIndex(indexRecord),chemical -> GetChemicalDiffusivityByIndex(indexRecord));
                }
            }
            
            break;
        }
    }

}

void AbstractDiffusiveChemistry::UpdateChemicalVectors(AbstractDiffusiveChemical* chemical)
{
    // virtual function to overide with further updates
    // additional diffusion properties
    AbstractChemistry::UpdateChemicalVectors(chemical); // implicit upcasting to AbstractChemcial

    UpdateDomainVector(chemical);

    mDiffusivityMatrix.push_back(chemical -> GetChemicalDiffusivityVector());
    mDiffusionDomainMatrix.push_back(chemical -> GetDiffusiveDomainVector());
    mDiffusionChemicalNames.push_back(chemical -> GetChemicalName());
    mNumberDiffusiveChemicals +=1;
}

void AbstractDiffusiveChemistry::AddChemistry(AbstractDiffusiveChemistry* newChemistry)
{
    // function to combine the chemistries in the effort to form an overall union of  chemistries, for use as domain
    // domain chemistry determines the solver state varaible ordering.  chemistry classes need not be of the same type, 
    // hence derived classes will need to dertermine the "highest class" define the highest class to be the parent chemistry

    // want to add the chemical vectors while preventing duplicates

    std::vector<AbstractChemical*> :: iterator chemical_iter;
    std::vector<AbstractChemical*> chemicalVector = newChemistry->rGetChemicalVector();

    for(chemical_iter = chemicalVector.begin(); chemical_iter != chemicalVector.end(); ++chemical_iter )
    {
        // add each of the chemicals from the additional chemistry in turn.  This checks for duplicates implicitly
        AddChemical(dynamic_cast<AbstractDiffusiveChemical*>(*chemical_iter));
    }
}

std::string AbstractDiffusiveChemistry::GetChemistryType()
{
    // virtual function to be overriden in derived classes, used to identify properties to compy when adding chemistries
    return "AbstractDiffusiveChemistry";
}

void AbstractDiffusiveChemistry::SetDiffusivityMatrix(std::vector<std::vector<double>> diffusivityMatrix)
{
    mDiffusivityMatrix=diffusivityMatrix;
}

std::vector<std::vector<double>> AbstractDiffusiveChemistry::GetDiffusivityMatrix()
{
    return mDiffusivityMatrix;
}

std::vector<double> AbstractDiffusiveChemistry::GetDiffusivityVectorByIndex(unsigned index)
{
    if(index < mNumberDiffusiveChemicals)
    {
        return mDiffusivityMatrix[index];
    }
    else
    {
        std::cout<<index<<std::endl;
        std::cout<<"Error: AbstractDiffusiveChemistry::GetDiffusivityByIndex(unsigned index), index out of bounds"<<std::endl;
        std::vector<double> returnVec(1,0.0);
        return returnVec;
    } 
}

double AbstractDiffusiveChemistry::GetDiffusivityValueByChemicalAndDomainName(std::string chemical_name, std::string domain_name)
{
    // return the diffusivity value for a given species with a given domain name
    unsigned index = GetDiffusiveChemicalIndexByName(chemical_name);

    std::vector<std::string>    domainsVector = GetDiffusionDomainsVectorByIndex(index);
    std::vector<double>     diffusivityVector = GetDiffusivityVectorByIndex(index);

    bool IsDomainFound = false;
    for(unsigned domain_index=0; domain_index<domainsVector.size();domain_index++)
    {
        if(domain_name == domainsVector[domain_index])
        {
            IsDomainFound = true;
            return diffusivityVector[domain_index];
        }
    }
    if(!IsDomainFound)
    {
        // the domain is unaccounted for in the AbstractDiffusiveChemical data structure
        return 0.0; // i.e not diffusive
    }
    return 0.0;
}

double AbstractDiffusiveChemistry::GetDiffusivityValueByChemicalName(std::string chemical_name)
{
    unsigned index = GetDiffusiveChemicalIndexByName(chemical_name);
    return GetDiffusivityVectorByIndex(index)[0];
}


void AbstractDiffusiveChemistry::SetDiffusionDomainsMatrix(std::vector<std::vector<std::string>> diffusionDomains)
{
    mDiffusionDomainMatrix=diffusionDomains;
}

std::vector<std::vector<std::string>> AbstractDiffusiveChemistry::GetDiffusionDomainsMatrix()
{
    return mDiffusionDomainMatrix;
}

std::vector<std::string> AbstractDiffusiveChemistry::GetDiffusionDomainsVectorByIndex(unsigned index)
{

    if(index < mNumberDiffusiveChemicals)
    {
        
        return mDiffusionDomainMatrix[index];
    }
    else
    {
        std::cout<<"Error: AbstractDiffusiveChemistry::GetDiffusionDomainsByIndex(unsigned index), index out of bounds"<<std::endl;
        std::vector<std::string> returnVec(1,"Null");
        return returnVec;
    } 
}

unsigned AbstractDiffusiveChemistry::GetNumberDiffusiveChemicals()
{
    return mNumberDiffusiveChemicals;
}

std::string AbstractDiffusiveChemistry::GetDiffusiveChemicalNamesByIndex(unsigned index)
{
    if(index<mNumberDiffusiveChemicals)
    {
        return mDiffusionChemicalNames[index];
    }
    else
    {
        std::cout<<"Error: AbstractDiffusiveChemistry::GetDiffusiveChemicalNamesByIndex(unsigned index), index out of bounds"<<std::endl;
        return "Error";
    }
    
    
}

unsigned AbstractDiffusiveChemistry::GetDiffusiveChemicalIndexByName(std::string chemicalName)
{
    bool IsFound=false;
    for(unsigned index=0; index<mNumberDiffusiveChemicals; index++)
    {
        if(chemicalName == mDiffusionChemicalNames[index])
        {
            IsFound = true;
            return index;
        }

    }
    if(!IsFound)
    {
        return mNumberDiffusiveChemicals; // should act out of bounds for future methods, yielding non-diffusive chemical
    }
    return mNumberDiffusiveChemicals;
}

#endif