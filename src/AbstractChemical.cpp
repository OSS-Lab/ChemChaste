#include "AbstractChemical.hpp"

AbstractChemical::AbstractChemical(std::string chemicalName, 
                    double size,
                    double mass,
                    int valence,
                    std::string chemicalDimensions)
    :   mChemicalName(chemicalName),
        mSize(size),
        mMass(mass),
        mValence(valence),
        mChemicalDimensions(chemicalDimensions)
{
    // default to not knowning the formation energy
    SetChemicalFormationGibbs(0.0);
    SetChemicalFormationKnown(false);
}

// set methods

void AbstractChemical::SetChemicalName(std::string ChemicalName)
{
    mChemicalName = ChemicalName;
}

void AbstractChemical::SetChemicalSize(double Size)
{
    mSize = Size;
}

void AbstractChemical::SetChemicalMass(double Mass)
{
    mMass  =Mass;
}

void AbstractChemical::SetChemicalValence(int Valence)
{
    mValence = Valence;
}

void AbstractChemical::SetChemicalFormationGibbs(double formationGibbs)
{
    mFormationGibbs = formationGibbs;
}

void AbstractChemical::SetChemicalFormationKnown(bool formationKnown)
{
    mFormationKnown = formationKnown;
}

void AbstractChemical::SetChemicalDimensions(std::string dimensions)
{
    mChemicalDimensions = dimensions;
}

// get methods

std::string AbstractChemical::GetChemicalName()
{
    return mChemicalName;
}

double AbstractChemical::GetChemicalSize()
{
    return mSize;
}

double AbstractChemical::GetChemicalMass()
{
    return mMass;
}

int AbstractChemical::GetChemicalValence()
{
    return mValence;
}

double AbstractChemical::GetChemicalFormationGibbs()
{
    return mFormationKnown;
}

std::string AbstractChemical::GetChemicalDimensions()
{
    return mChemicalDimensions;
}

bool AbstractChemical::IsChemicalFormationKnown()
{
    return mFormationKnown;
}