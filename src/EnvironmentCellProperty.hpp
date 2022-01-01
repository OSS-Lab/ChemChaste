#ifndef ENVIRONMENTCELLPROPERTY_HPP
#define ENVIRONMENTCELLPROPERTY_HPP

//general includes
#include <vector>
#include <boost/shared_ptr.hpp>

// chaste includes
#include "Cell.hpp"
#include "AbstractCellProperty.hpp"
#include "StateVariableRegister.hpp"
#include "AbstractChemistry.hpp"
#include "ChastePoint.hpp"


class EnvironmentCellProperty : public AbstractCellProperty
{
protected:

    // CellPtr to a given cell for access to cell data and properties
    CellPtr mThis_cellPtr;

    // register for variables in the environment for recording
    StateVariableRegister* mpEnvironmentStateVariableRegister;

    // the values for the environment parameters
    std::vector<double> mEnvironmentVector;

    // the preferred values for the environment parameters
    std::vector<double> mPreferredEnvironmentVector;



public:

    EnvironmentCellProperty();

    virtual ~EnvironmentCellProperty();

    EnvironmentCellProperty(const EnvironmentCellProperty&);

    // virtual methods

    virtual void SetUp(CellPtr, std::vector<std::string>);

    virtual void UpdateEnvironmentVector(std::vector<double>);
    
    virtual void UpdatePreferredEnvironmentVector(std::vector<double>);

    virtual void PreparePostDivisionParent(double);
    
    virtual void PreparePostDivisionDaughter(const EnvironmentCellProperty&, double);

    // concrete methods

    void SetUpEnvironmentVectors(unsigned);

    // set methods

    void SetEnvironmentStateVariableRegister(StateVariableRegister*);

    void SetEnvironmentVector(std::vector<double>);

    void SetEnvironmentValueByIndex(unsigned,double);

    void SetEnvironmentValueByName(std::string,double);

    void SetPreferredEnvironmentVector(std::vector<double>);

    void SetCellPtr(CellPtr);

    // get methods

    StateVariableRegister* GetEnvironmentStateVariableRegister();

    std::vector<double> GetEnvironmentVector();

    std::vector<double> GetPreferredEnvironmentVector();

    double GetEnvironmentValueByIndex(unsigned);

    double GetPreferredEnvironmentValueByIndex(unsigned);

    double GetEnvironmentValueByName(std::string);

    double GetPreferredEnvironmentValueByName(std::string);

    CellPtr GetCellPtr();
};


#endif