#ifndef TRANSPORTCELLPROPERTY_HPP
#define TRANSPORTCELLPROPERTY_HPP

//general includes
#include <vector>
#include <boost/shared_ptr.hpp>

// chaste includes
#include "Cell.hpp"
#include "AbstractCellProperty.hpp"
#include "StateVariableRegister.hpp"
#include "AbstractTransportReactionSystem.hpp"
#include "AbstractTransportOdeSystem.hpp"
#include "AbstractIvpOdeSolver.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "ChastePoint.hpp"


class TransportCellProperty : public AbstractCellProperty
{
protected:

    // CellPtr to a given cell for access to cell data and for properties that relate ot the transport ode
    CellPtr mThis_cellPtr;

    // register for variables that cross the membrane; either side
    StateVariableRegister* mpBulkStateVariableRegister;
    StateVariableRegister* mpCellStateVariableRegister;

    // transport reaction system for the passage of state variables across the boundary
    AbstractTransportReactionSystem* mpTransportReactionSystem;

    // transport reaction system Ode corresponding to the mpTransportReactionSystem
    AbstractTransportOdeSystem* mpTransportOdeSystem;

    // transport Ode solver
    boost::shared_ptr<AbstractIvpOdeSolver> mpSolver;

    // concentration of state variables at the internal boundary of the cell if this is to remain constant
    bool mIsConstantCellConcentration = false; // for use in averaged source for example?

    std::vector<double> mBulkBoundaryConcentrationVector;

    std::vector<double> mChangeBulkBoundaryConcentrationVector;

    std::vector<double> mCellBoundaryConcentrationVector;

    std::vector<double> mChangeCellBoundaryConcentrationVector;

    bool mIncludeOdeInterpolationOnBoundary=false;


    //15/10/2020
    std::vector<double> mInitBulkBoundaryConcentrationVector;

    unsigned mNumberOCalls_this_reaction_step =0;

public:

    TransportCellProperty();

    virtual ~TransportCellProperty();

    TransportCellProperty(const TransportCellProperty&);

    // virtual methods

    virtual void SetUp(AbstractTransportReactionSystem*, CellPtr);

    virtual void UpdateCellConcentrationVector(std::vector<double>);

    virtual void UpdateBulkConcentrationVector(std::vector<double>);

    virtual void PerformTransportSystem(const std::vector<double>&, const std::vector<double>&, std::vector<double>&, std::vector<double>&);

    virtual void UpdateTransportReactionSystem(AbstractTransportReactionSystem*);

    virtual void UpdateTransportOdeSystem(AbstractTransportOdeSystem*);

    virtual void PreparePostDivisionParent(double);
    
    virtual void PreparePostDivisionDaughter(const TransportCellProperty&, double);

    // concrete methods
    //double RetrieveBoundarySourceByStateName(std::string,ChastePoint<SPACE_DIM> ); // moved to extended cell property

    void SetUpCellConcentrationVector(unsigned);

    void SetUpBulkConcentrationVector(unsigned);

    void SetUpChangeCellConcentrationVector(unsigned);

    void SetUpChangeBulkConcentrationVector(unsigned);

    double RetrieveBoundarySourceByStateName(std::string); // overload for non-extended cell property case

    double RetrieveChangeBoundarySourceByStateName(std::string); // overload for non-extended cell property case


    void AppendInternalCellBoundaryConcentrations(std::vector<double>&);

    void ReplaceBoundaryStateVariables(std::vector<double>&);

    void ReplaceChangeBoundaryStateVariables(std::vector<double>&);

    void ResetReactionCalls();

    unsigned GetReactionCalls();

    // set methods

    void SetBulkStateVariableRegister(StateVariableRegister*);

    void SetCellStateVariableRegister(StateVariableRegister*);

    void SetConstantCellConcentrationBool(bool);

    void SetCellBoundaryConcentrationVector(std::vector<double>);

    void SetBulkBoundaryConcentrationVector(std::vector<double>);

    void SetInitBulkBoundaryConcentrationVector(std::vector<double>);

    void SetTransportOdeSolver(boost::shared_ptr<AbstractIvpOdeSolver>);

    void SetIncludeOdeInterpolationOnBoundary(bool);

    void SetCellPtr(CellPtr);

    // get methods

    StateVariableRegister* GetBulkStateVariableRegister();

    StateVariableRegister* GetCellStateVariableRegister();

    bool GetConstantCellConcentrationBool();

    std::vector<double> GetInternalCellBoundaryConcentrationVector();

    double GetInternalCellBoundaryConcentrationByIndex(unsigned);

    double GetInternalCellBoundaryConcentrationByName(std::string);

    double GetChangeInternalCellBoundaryConcentrationByName(std::string);

    std::vector<double> GetExternalCellBoundaryConcentrationVector();

    double GetExternalCellBoundaryConcentrationByIndex(unsigned);

    double GetChangeExternalCellBoundaryConcentrationByIndex(unsigned);

    double GetExternalCellBoundaryConcentrationByName(std::string);

    AbstractTransportReactionSystem* GetTransportReactionSystem();

    AbstractTransportOdeSystem* GetTransportOdeSystem();

    boost::shared_ptr<AbstractIvpOdeSolver> GetTransportOdeSolver();

    bool GetIncludeOdeInterpolationOnBoundary();

    CellPtr GetCellPtr();
};


#endif