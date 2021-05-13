#ifndef ABSTRACTCHEMICAL_HPP_
#define ABSTRACTCHEMICAL_HPP_

// general includes
#include <string>
#include <tuple>
#include <vector>

// abstract property to contain information about a chemical species used in ChemChaste
// conatains the name of the chemical and basic chemcial properties

class AbstractChemical
{
private:

    std::string mChemicalName;

    // chemical properties list for accessing by inherited properties

    // default value 1.0
    double mSize;

    // default value 1.0
    double mMass;

    // molecular charge for case of electrical reactions
    // default value 0
    int mValence;

    // numerical dimensions for the chemical traits, used in the ode information to define the states
    // default "non-dim"
    std::string mChemicalDimensions;

    // energetic information for gibbs formation energies and reaction rate calculation
    // no default set
    bool mFormationKnown;

    // no default set
    double mFormationGibbs;


public:

    AbstractChemical(   std::string chemicalName = "", 
                        double size = 1.0,
                        double mass = 1.0,
                        int valence = 0,
                        std::string chemicalDimensions = "non-dim"
                        );

    virtual ~AbstractChemical()
    {
    };

    // set methods

    void SetChemicalName(std::string);

    void SetChemicalSize(double);

    void SetChemicalMass(double);

    void SetChemicalValence(int);

    void SetChemicalFormationGibbs(double);

    void SetChemicalFormationKnown(bool);

    void SetChemicalDimensions(std::string);

    // get methods

    std::string GetChemicalName();

    double GetChemicalSize();

    double GetChemicalMass();

    int GetChemicalValence();

    double GetChemicalFormationGibbs();

    std::string GetChemicalDimensions();

    bool IsChemicalFormationKnown();

    virtual std::string GetChemicalType()
    {
        return "AbstractChemical";
    };

};

#endif