#ifndef ABSTRACTCHEMISTRY_HPP_
#define ABSTRACTCHEMISTRY_HPP_

#include "AbstractChemical.hpp"

#include <string>
#include <tuple>
#include <vector>
#include <iostream>

// abstract property to contain information about the collection of chemical species in the system
// container operator for adding two AbstractReactionSystems to form union for bulk 
// copy contructor for mutating reaction network

// inherit? : public AbstractChemical

class AbstractChemistry 
{
private:
    unsigned mNumberChemicals;
protected:
    std::vector<AbstractChemical*> mpChemicalVector;

    std::vector<std::string> mChemicalNames;

    std::vector<std::string> mChemicalDimensions;
    
public:
    AbstractChemistry();

    virtual ~AbstractChemistry()
    {
    };

    virtual void AddChemical(AbstractChemical*);

    virtual bool CheckChemical(AbstractChemical*);

    virtual void UpdateChemicalVectors(AbstractChemical*);

    void SetChemicalVector(std::vector<AbstractChemical*>);

    void SetChemicalNames(std::vector<std::string>);

    void SetChemicalDimensions(std::vector<std::string>);

    std::vector<AbstractChemical*> rGetChemicalVector();

    std::vector<std::string> GetChemicalNames();

    std::string GetChemicalNamesByIndex(unsigned);

    unsigned GetChemicalIndexByName(std::string);

    std::vector<std::string> GetChemicalDimensions();

    std::string GetChemicalDimensionsByIndex(unsigned);

    virtual std::string GetChemistryType();

    virtual void AddChemistry(AbstractChemistry*);

    unsigned GetNumberChemicals();

};

#endif