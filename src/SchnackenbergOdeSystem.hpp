
#ifndef SCHNACKENBERGODESYSTEM_HPP_
#define SCHNACKENBERGODESYSTEM_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include <cmath>
#include <iostream>

#include "AbstractOdeSystem.hpp"

class SchnackenbergOdeSystem : public AbstractOdeSystem
{

public:

    /**
     * Default constructor.
     *
     * @param stateVariables optional initial conditions for state variables (only used in archiving)
     */
    SchnackenbergOdeSystem(std::vector<double> stateVariables=std::vector<double>());

    /**
     * Destructor.
     */
    ~SchnackenbergOdeSystem();

    /**
     * Compute the RHS of the  Collier et al. system of ODEs.
     *
     * Returns a vector representing the RHS of the ODEs at each time step, y' = [y1' ... yn'].
     * An ODE solver will call this function repeatedly to solve for y = [y1 ... yn].
     *
     * @param time used to evaluate the RHS.
     * @param rY value of the solution vector used to evaluate the RHS.
     * @param rDY filled in with the resulting derivatives (using  Collier et al. system of equations).
     */
    void EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY);
};


#endif 
