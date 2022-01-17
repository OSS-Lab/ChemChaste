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
    mCellTypeName = existingProperty.mCellTypeName;
}

// virtual methods

void CellAnalyticsProperty::SetUp(CellPtr this_cellPtr,unsigned cellID, std::string cellTypeName)
{
    SetCellID(cellID);
    SetCellPtr(this_cellPtr);
    SetCellTypeName(cellTypeName);

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

void CellAnalyticsProperty::SetCellTypeName(std::string cellTypeName)
{
    mCellTypeName = cellTypeName;
}

void CellAnalyticsProperty::SetNeighbourTypes(std::vector<unsigned> neighbourTypes)
{
    mNeighbourTypes = neighbourTypes;
}

void CellAnalyticsProperty::SetNeighbourTypeProbability(std::vector<double> neighbourProbability)
{
    mNeighbourTypeProbability = neighbourProbability;
}

void CellAnalyticsProperty::SetPopulationCellTypeNames(std::vector<std::string> cellTypeNames)
{
    mPopulationCellTypeNames = cellTypeNames;
}

void CellAnalyticsProperty::SetCellDiversity(double cellDiversity)
{
    mCellDiversity = cellDiversity;
}

void CellAnalyticsProperty::SetMeanCellDiversity(double cellDiversity)
{
    mMeanCellDiversity = cellDiversity;
}

CellPtr CellAnalyticsProperty::GetCellPtr()
{
    return mThis_cellPtr;
}

unsigned CellAnalyticsProperty::GetCellID()
{
    return mCellID;
}

std::string CellAnalyticsProperty::GetCellTypeName()
{
    return mCellTypeName;
}

std::vector<unsigned> CellAnalyticsProperty::GetNeighbourTypes()
{
    return mNeighbourTypes;
}

std::vector<double> CellAnalyticsProperty::GetNeighbourTypeProbability()
{
    return mNeighbourTypeProbability;
}

std::vector<std::string> CellAnalyticsProperty::GetPopulationCellTypeNames()
{
    return mPopulationCellTypeNames;
}

double CellAnalyticsProperty::GetCellDiversity()
{
    return mCellDiversity;
}

double CellAnalyticsProperty::GetMeanCellDiversity()
{
    return mMeanCellDiversity;
}
