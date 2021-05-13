#ifndef CHEMICALCELLPROPERTY_HPP
#define CHEMICALCELLPROPERTY_HPP

//general includes
#include <vector>
#include <string>
#include <boost/shared_ptr.hpp>

// chaste includes
#include "AbstractCellProperty.hpp"
#include "StateVariableRegister.hpp"

// cell porperty to handle the basic chemical property required for cell simulations

class ChemicalCellProperty : public AbstractCellProperty
{
protected:

    // register for the state variables present in the cell

    StateVariableRegister* mpStateVariableRegister;

    // concentration vector of the cell, would be derived from the CellData property
    std::vector<double> mConcentrationVector;


public:

    ChemicalCellProperty();

    virtual ~ChemicalCellProperty();

    virtual void InitialiseCell(std::vector<std::string>, std::vector<double>);

    virtual void InitialiseCell(StateVariableRegister*, std::vector<double>);

    virtual void UpdateCellConcentrationVector(std::vector<double>&);

    void SetStateVariableRegister(StateVariableRegister*);

    StateVariableRegister* GetStateVariableRegister();

    std::vector<double> GetCellConcentrationVector();

    double GetCellConcentrationByIndex(unsigned);

    double GetCellConcentrationByName(std::string);

};

#endif