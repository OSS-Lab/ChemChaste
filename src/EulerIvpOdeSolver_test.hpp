#ifndef _EULERIVPODESOLVERTEST_HPP_
#define _EULERIVPODESOLVERTEST_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <iostream>
#include "AbstractOneStepIvpOdeSolver.hpp"

/**
 * A concrete one step ODE solver class that employs the forward Euler
 * method. This numerical method is explicit.
 */
class EulerIvpOdeSolver_test : public AbstractOneStepIvpOdeSolver
{
private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the abstract IVP Solver, never used directly - boost uses this.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // This calls serialize on the base class.
        archive & boost::serialization::base_object<AbstractOneStepIvpOdeSolver>(*this);
    }

protected:

    /**
     * Calculate the solution to the ODE system at the next timestep.
     *
     * A usage example:
     *     EulerIvpOdeSolver_test mySolver;
     *     OdeSolution solution = mySolver.Solve(pMyOdeSystem, yInit, StartTime, EndTime, TimeStep, SamplingTime);
     *
     * @param pAbstractOdeSystem  the ODE system to solve
     * @param timeStep  dt
     * @param time  the current time
     * @param rCurrentYValues  the current (initial) state
     * @param rNextYValues  the state at the next timestep
     */
    void CalculateNextYValue(AbstractOdeSystem* pAbstractOdeSystem,
                             double timeStep,
                             double time,
                             std::vector<double>& rCurrentYValues,
                             std::vector<double>& rNextYValues);

public:

    /**
     * Constructor.
     */
    EulerIvpOdeSolver_test()
    {}

    /**
     * Destructor.
     */
    virtual ~EulerIvpOdeSolver_test()
    {}
};

#endif //_EulerIvpOdeSolver_test_HPP_
