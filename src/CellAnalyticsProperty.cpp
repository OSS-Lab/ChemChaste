#include "CellAnalyticsProperty.hpp"



CellAnalyticsProperty::CellAnalyticsProperty()
    : AbstractCellProperty()
{
}


CellAnalyticsProperty::~CellAnalyticsProperty()
{
}

CellAnalyticsProperty::CellAnalyticsProperty(const CellAnalyticsProperty& existingProperty)
{
    mCellID = existingProperty.mCellID;
}

// virtual methods

void CellAnalyticsProperty::SetUp(CellPtr this_cellPtr,unsigned cellID)
{
    SetCellID(cellID);
    SetCellPtr(this_cellPtr);

}

void CellAnalyticsProperty::PreparePostDivisionParent(double splitRatio)
{
    // split any properties that are shared
}
    
void  CellAnalyticsProperty::PreparePostDivisionDaughter(const CellAnalyticsProperty& parentProperty, double splitRatio)
{
    // split any properties that are shared
}


void CellAnalyticsProperty::SetCellPtr(CellPtr this_cellPtr)
{
    mThis_cellPtr = this_cellPtr;
}

void CellAnalyticsProperty::SetCellID(unsigned cellID)
{
    mCellID = cellID;
}


CellPtr CellAnalyticsProperty::GetCellPtr()
{
    return mThis_cellPtr;
}

unsigned CellAnalyticsProperty::GetCellID()
{
    return mCellID;
}
