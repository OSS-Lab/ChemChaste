#ifndef CELLANALYTICSPROPERTY_HPP
#define CELLANALYTICSPROPERTY_HPP

//general includes
#include <vector>
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


public:

    CellAnalyticsProperty();


    virtual ~CellAnalyticsProperty();

    CellAnalyticsProperty(const CellAnalyticsProperty&);

    // virtual methods

    virtual void SetUp(CellPtr,unsigned);

    virtual void PreparePostDivisionParent(double);
    
    virtual void PreparePostDivisionDaughter(const CellAnalyticsProperty&, double);

    void SetCellPtr(CellPtr);

    void SetCellID(unsigned);


    CellPtr GetCellPtr();

    unsigned GetCellID();

};




#endif