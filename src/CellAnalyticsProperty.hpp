#ifndef CELLANALYTICSPROPERTY_HPP
#define CELLANALYTICSPROPERTY_HPP

//general includes
#include <vector>
#include <string>
#include <boost/shared_ptr.hpp>

// chaste includes
#include "Cell.hpp"
#include "AbstractCellProperty.hpp"



class CellAnalyticsProperty : public AbstractCellProperty
{

protected:

    // CellPtr to a given cell for access to cell data and for properties
    CellPtr mThis_cellPtr;

    unsigned mCellID;

    std::string mCellTypeName = "";

    std::vector<unsigned> mNeighbourTypes;

    std::vector<std::string> mPopulationCellTypeNames;

public:

    CellAnalyticsProperty();


    virtual ~CellAnalyticsProperty();

    CellAnalyticsProperty(const CellAnalyticsProperty&);

    // virtual methods

    virtual void SetUp(CellPtr,unsigned, std::string);

    virtual void PreparePostDivisionParent(double);
    
    virtual void PreparePostDivisionDaughter(const CellAnalyticsProperty&, double);

    void SetCellPtr(CellPtr);

    void SetCellID(unsigned);

    void SetCellTypeName(std::string);

    void SetNeighbourTypes(std::vector<unsigned>);

    void SetPopulationCellTypeNames(std::vector<std::string>);

    CellPtr GetCellPtr();

    unsigned GetCellID();

    std::string GetCellTypeName();

    std::vector<unsigned> GetNeighbourTypes();

    std::vector<std::string> GetPopulationCellTypeNames();
};




#endif