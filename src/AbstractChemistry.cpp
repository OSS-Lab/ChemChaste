#include "AbstractChemistry.hpp"

AbstractChemistry::AbstractChemistry()
{
    mpChemicalVector = std::vector<AbstractChemical*>();
    mChemicalNames = std::vector<std::string>();
    mChemicalDimensions  =std::vector<std::string>();
    mNumberChemicals =0;
}

void AbstractChemistry::AddChemistry(AbstractChemistry *newChemistry)
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
        AddChemical(dynamic_cast<AbstractChemical*>(*chemical_iter));
    }
}

void AbstractChemistry::AddChemical(AbstractChemical *chemical)
{
    // function to check whether a chemical is present, if not add to mChemicalVector

    bool newChem = false;
    if(mpChemicalVector.empty())
    {
        newChem = true;
    }
    else
    {
    
        newChem = CheckChemical(chemical);
    }

    if(newChem)
    {
    
        UpdateChemicalVectors(chemical);
    }

}

bool AbstractChemistry::CheckChemical(AbstractChemical* chemical)
{
    // function to return true if chemical is already in class chemical vector

    bool newChem  = true;

    for(std::vector<AbstractChemical*> :: iterator chemical_iter = mpChemicalVector.begin(); 
        chemical_iter != mpChemicalVector.end();
         ++chemical_iter)
    {
    
        if(chemical -> GetChemicalName() ==  dynamic_cast<AbstractChemical*>(*chemical_iter) -> GetChemicalName())
        {

            newChem = false;
            break;
        }
    }
    return newChem;
}


void AbstractChemistry::UpdateChemicalVectors(AbstractChemical* chemical)
{
    mpChemicalVector.push_back(chemical);
    // properties needed for general function
    mChemicalNames.push_back(chemical -> GetChemicalName()); // may remove this, will lead to duplicates
    mChemicalDimensions.push_back(chemical -> GetChemicalDimensions()); // may remove this, will lead to duplicates
    // further properties, to be overridded in derived classes
    mNumberChemicals +=1;
    return;
}

void AbstractChemistry::SetChemicalVector(std::vector<AbstractChemical*> chemicalVector)
{
    mpChemicalVector=chemicalVector;
}

void AbstractChemistry::SetChemicalNames(std::vector<std::string> chemicalNames)
{
    mChemicalNames=chemicalNames;
}

void AbstractChemistry::SetChemicalDimensions(std::vector<std::string> chemicalDimensions)
{
    mChemicalDimensions=chemicalDimensions;
}

std::vector<AbstractChemical*> AbstractChemistry::rGetChemicalVector()
{
    return mpChemicalVector;
}

std::vector<std::string> AbstractChemistry::GetChemicalNames()
{
    return mChemicalNames;
}

std::string AbstractChemistry::GetChemicalNamesByIndex(unsigned index)
{
    if(index < mNumberChemicals)
    {
        return mChemicalNames[index];
    }
    else
    {
        std::cout<<"Error: AbstractChemistry::GetChemicalNamesByIndex(unsigned index), index out of bounds"<<std::endl;
        return "Null";
    } 
}

unsigned AbstractChemistry::GetChemicalIndexByName(std::string name)
{
    unsigned index=0;
    for(std::vector<AbstractChemical*> :: iterator chemical_iter = mpChemicalVector.begin(); 
        chemical_iter != mpChemicalVector.end();
         ++chemical_iter)
    {
    
        if(name ==  dynamic_cast<AbstractChemical*>(*chemical_iter) -> GetChemicalName())
        {
            break;
        }
        index +=1;
    }

    // if name is not one of the chemicals then index== mNumberChemicals; would throw error later on
    return index;
}

std::vector<std::string> AbstractChemistry::GetChemicalDimensions()
{
    return mChemicalDimensions;
}

std::string AbstractChemistry::GetChemicalDimensionsByIndex(unsigned index)
{
    if(index < mNumberChemicals)
    {
        return mChemicalDimensions[index];
    }
    else
    {
        std::cout<<"Error: AbstractChemistry::GetChemicalDimensionsByIndex(unsigned index), index out of bounds"<<std::endl;
       return "Null";
    } 
    
}

std::string AbstractChemistry::GetChemistryType()
{
    // virtual function to be overriden in derived classes, used to identify properties to compy when adding chemistries
    return "AbstractChemistry";
}

unsigned AbstractChemistry::GetNumberChemicals()
{
    return mNumberChemicals;
}
