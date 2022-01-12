#ifndef FIXEDCELLBOUNDARYCONDITION_HPP
#define FIXEDCELLBOUNDARYCONDITION_HPP

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/shared_ptr.hpp>
#include <vector>
#include <string>
#include "AbstractCellPopulationBoundaryCondition.hpp"



template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class FixedCellBoundaryCondition : public AbstractCellPopulationBoundaryCondition<ELEMENT_DIM, SPACE_DIM>
{
private:

    // maximum values for bounadries xmax, ymax, zmax
    std::vector<double> mBoundaryMax;

    // minimum values for bounadries xmin, ymin, zmin
    std::vector<double> mBoundaryMin;

public:

    FixedCellBoundaryCondition(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation,
                           std::vector<double> boundaryMax,
                           std::vector<double> boundaryMin);

    /**
     * Overridden ImposeBoundaryCondition() method.
     *
     * Apply the cell population boundary conditions.
     *
     * @param rOldLocations the node locations before any boundary conditions are applied
     */
    void ImposeBoundaryCondition(const std::map<Node<SPACE_DIM>*, c_vector<double, SPACE_DIM> >& rOldLocations);

    /**
     * Overridden VerifyBoundaryCondition() method.
     * Verify the boundary conditions have been applied.
     * This is called after ImposeBoundaryCondition() to ensure the condition is still satisfied.
     *
     * @return whether the boundary conditions are satisfied.
     */
    bool VerifyBoundaryCondition();

    /**
     * Overridden OutputCellPopulationBoundaryConditionParameters() method.
     * Output cell population boundary condition parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile);


    void SetBoundaryMax(std::vector<double> boundaryMax);

    void SetBoundaryMin(std::vector<double> boundaryMin);

    std::vector<double> GetBoundaryMax();

    std::vector<double> GetBoundaryMin();

    
};



#endif