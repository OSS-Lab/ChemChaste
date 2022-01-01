#ifndef INHOMOGENOUSCOUPLEDPDEODECOUPLEDCELLSOLVER_HPP_
#define INHOMOGENOUSCOUPLEDPDEODECOUPLEDCELLSOLVER_HPP_

#include "AbstractAssemblerSolverHybrid.hpp"
#include "AbstractDynamicLinearPdeSolver.hpp"
#include "AbstractLinearParabolicPdeSystemForCoupledOdeSystem.hpp"
#include "TetrahedralMesh.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "AbstractInhomogenousOdeSystemForCoupledPdeSystem.hpp"
#include "AbstractOdeSystemForCoupledPdeSystem.hpp"
#include "CvodeAdaptor.hpp"
#include "BackwardEulerIvpOdeSolver.hpp"
#include "Warnings.hpp"
#include "VtkMeshWriter.hpp"
#include "StateVariableRegister.hpp"
#include "InhomogenousCoupledPdeOdeSolver_templated.hpp"
#include "TransportCellProperty.hpp"

#include <boost/shared_ptr.hpp>

// Same as other class but the numberOfStateVariables to the ode terms are not all equal.
// The ode's change on a nodal basis

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM=ELEMENT_DIM, unsigned PROBLEM_DIM=2>
class InhomogenousCoupledPdeOdeCoupledCellSolver
    : public AbstractAssemblerSolverHybrid<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM, NORMAL>,
      public AbstractDynamicLinearPdeSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>
{
private:

    /** Pointer to the mesh. */
    AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>* mpFeMesh;

    /** The PDE system to be solved. */ 
    //CellSourceInhomogenousParabolicPdeOdeSystem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>* mpPdeSystem;
    InhomogenousParabolicPdeForCoupledOdeSystemTemplated<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>* mpPdeSystem;

    /** Vector of pointers to ODE systems, defined at nodes. */
    std::vector<AbstractInhomogenousOdeSystemForCoupledPdeSystem*> mOdeSystemsAtNodes;

    /** The values of the ODE system state variables, interpolated at a quadrature point. */
    std::vector<double> mInterpolatedOdeStateVariables;

    /** The ODE solvers. */
    std::vector<boost::shared_ptr<AbstractIvpOdeSolver>> mpOdeSolvers;

    bool mCellTransportOdeSystemsPresent;

    bool mCellMembraneOdeSystemsPresent;

    bool mConditionsInterpolated=false;

    AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& mrCellPopulation;

    ChastePoint<SPACE_DIM> mX;

    c_vector<double,PROBLEM_DIM> mU;

    unsigned mInterpolationCount=0;

    /**
     * A sampling timestep for writing results to file. Set to
     * PdeSimulationTime::GetPdeTimeStep() in the constructor;
     * may be overwritten using the SetSamplingTimeStep() method.
     */
    double mSamplingTimeStep;

    /** Whether ODE systems are present (if not, then the system comprises coupled PDEs only). */
    bool mOdeSystemsPresent;

    /** Meta results file for VTK. */
    out_stream mpVtkMetaFile;

    /**
     * Whether the output directory should be cleared before solve or not. False by default.
     * Can be changed when setting the output directory
     */
    bool mClearOutputDirectory;

    /**
     * Write the current results to mpVtkMetaFile.
     */
    void WriteVtkResultsToFile();

    /**
     * @return the term to be added to the element stiffness matrix.
     *
     * @param rPhi The basis functions, rPhi(i) = phi_i, i=1..numBases
     * @param rGradPhi Basis gradients, rGradPhi(i,j) = d(phi_j)/d(X_i)
     * @param rX The point in space
     * @param rU The unknown as a vector, u(i) = u_i
     * @param rGradU The gradient of the unknown as a matrix, rGradU(i,j) = d(u_i)/d(X_j)
     * @param pElement Pointer to the element
     */
    c_matrix<double, PROBLEM_DIM*(ELEMENT_DIM+1), PROBLEM_DIM*(ELEMENT_DIM+1)> ComputeMatrixTerm(
        c_vector<double, ELEMENT_DIM+1>& rPhi,
        c_matrix<double, SPACE_DIM, ELEMENT_DIM+1>& rGradPhi,
        ChastePoint<SPACE_DIM>& rX,
        c_vector<double,PROBLEM_DIM>& rU,
        c_matrix<double, PROBLEM_DIM, SPACE_DIM>& rGradU,
        Element<ELEMENT_DIM, SPACE_DIM>* pElement);

    /**
     * @return the term to be added to the element stiffness vector.
     *
     * @param rPhi The basis functions, rPhi(i) = phi_i, i=1..numBases
     * @param rGradPhi Basis gradients, rGradPhi(i,j) = d(phi_j)/d(X_i)
     * @param rX The point in space
     * @param rU The unknown as a vector, u(i) = u_i
     * @param rGradU The gradient of the unknown as a matrix, rGradU(i,j) = d(u_i)/d(X_j)
     * @param pElement Pointer to the element
     */
    c_vector<double, PROBLEM_DIM*(ELEMENT_DIM+1)> ComputeVectorTerm(
        c_vector<double, ELEMENT_DIM+1>& rPhi,
        c_matrix<double, SPACE_DIM, ELEMENT_DIM+1>& rGradPhi,
        ChastePoint<SPACE_DIM>& rX,
        c_vector<double,PROBLEM_DIM>& rU,
        c_matrix<double,PROBLEM_DIM,SPACE_DIM>& rGradU,
        Element<ELEMENT_DIM, SPACE_DIM>* pElement);

    /**
     * Reset the member variable mInterpolatedOdeStateVariables.
     */
    void ResetInterpolatedQuantities();

    /**
     * Update the member variable mInterpolatedOdeStateVariables by computing the
     * interpolated value of each ODE state variable at each Gauss point.
     *
     * @param phiI
     * @param pNode pointer to a Node
     */
    void IncrementInterpolatedQuantities(double phiI, const Node<SPACE_DIM>* pNode);

    /**
     * Initialise method: sets up the linear system (using the mesh to
     * determine the number of unknowns per row to preallocate) if it is not
     * already set up. Can use an initial solution as PETSc template,
     * or base it on the mesh size.
     *
     * @param initialSolution Initial solution (defaults to NULL) for PETSc to use as a template.
     */
    void InitialiseForSolve(Vec initialSolution=NULL);

    /**
     * Completely set up the linear system that has to be solved each timestep.
     *
     * @param currentSolution The current solution which can be used in setting up
     *  the linear system if needed (NULL if there isn't a current solution)
     * @param computeMatrix Whether to compute the LHS matrix of the linear system
     *   (mainly for dynamic solves).
     */
    void SetupLinearSystem(Vec currentSolution, bool computeMatrix);

public:

    /**
     * Constructor.
     *
     * @param pMesh pointer to the mesh
     * @param pPdeSystem pointer to the PDE system
     * @param pBoundaryConditions pointer to the boundary conditions.
     * @param odeSystemsAtNodes optional vector of pointers to ODE systems, defined at nodes
     * @param pOdeSolver optional pointer to an ODE solver (defaults to NULL)
     */
    InhomogenousCoupledPdeOdeCoupledCellSolver(TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>* pMesh,
                                    InhomogenousParabolicPdeForCoupledOdeSystemTemplated<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>* pPdeSystem,
                                    BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>* pBoundaryConditions,
                                    AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation,
                                    std::vector<AbstractInhomogenousOdeSystemForCoupledPdeSystem*> odeSystemsAtNodes=std::vector<AbstractInhomogenousOdeSystemForCoupledPdeSystem*>(),
                                    std::vector<boost::shared_ptr<AbstractIvpOdeSolver>> pOdeSolvers=std::vector<boost::shared_ptr<AbstractIvpOdeSolver>>(),
                                    bool conditionsInterpolated =false,
                                    ChastePoint<SPACE_DIM> point = ChastePoint<SPACE_DIM>(),
                                    c_vector<double,PROBLEM_DIM> stateVector = zero_vector<double>(PROBLEM_DIM));

    /**
     * Destructor.
     * If an ODE system is present, the pointers to the ODE system objects are deleted here.
     */
    ~InhomogenousCoupledPdeOdeCoupledCellSolver();

    /**
     * Overridden PrepareForSetupLinearSystem() method.
     * Pass the current solution to the PDE system to the ODE system and solve it over the next timestep.
     *
     * @param currentPdeSolution the solution to the PDE system at the current time
     */
    void PrepareForSetupLinearSystem(Vec currentPdeSolution);

    /**
     * Set mOutputDirectory.
     *
     * @param outputDirectory the output directory to use
     * @param clearDirectory whether to clear outputDirectory or not. Note that the actual clearing happens when you call SolveAndWriteResultsToFile().
     *                       False by default.
     */
    void SetOutputDirectory(std::string outputDirectory, bool clearDirectory=false);

    /**
     * Set mSamplingTimeStep.
     *
     * @param samplingTimeStep the sampling timestep to use
     */
    void SetSamplingTimeStep(double samplingTimeStep);

    /**
     * Solve the coupled PDE/ODE system over the pre-specified time interval,
     * and record results using mSamplingTimeStep.
     */
    void SolveAndWriteResultsToFile();

    /**
     * Write the solution to VTK. Called by SolveAndWriteResultsToFile().
     *
     * @param solution the solution of the coupled PDE/ODE system
     * @param numTimeStepsElapsed the number of timesteps that have elapsed
     */
    void WriteVtkResultsToFile(Vec solution, unsigned numTimeStepsElapsed);

    /**
     * Get a pointer to the ODE system defined at a given node.
     *
     * @param index the global index of a node in the mpFeMesh
     * @return mOdeSystemsAtNodes[index]
     */
    AbstractOdeSystemForCoupledPdeSystem* GetOdeSystemAtNode(unsigned index);

    void SetInterpolationPoint(ChastePoint<SPACE_DIM> );

    void SetCurrentStateVector(c_vector<double,PROBLEM_DIM>);

    bool CheckChastePointsForEquality(ChastePoint<SPACE_DIM>,ChastePoint<SPACE_DIM>);
};

///////////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////////

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
c_matrix<double, PROBLEM_DIM*(ELEMENT_DIM+1), PROBLEM_DIM*(ELEMENT_DIM+1)> InhomogenousCoupledPdeOdeCoupledCellSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::ComputeMatrixTerm(
    c_vector<double, ELEMENT_DIM+1>& rPhi,
    c_matrix<double, SPACE_DIM, ELEMENT_DIM+1>& rGradPhi,
    ChastePoint<SPACE_DIM>& rX,
    c_vector<double,PROBLEM_DIM>& rU,
    c_matrix<double, PROBLEM_DIM, SPACE_DIM>& rGradU,
    Element<ELEMENT_DIM, SPACE_DIM>* pElement)
{
 //   std::cout<<"InhomogenousCoupledPdeOdeCoupledCellSolver - ComputeMatrixTerm"<<std::endl;
    double timestep_inverse = PdeSimulationTime::GetPdeTimeStepInverse();
    c_matrix<double, PROBLEM_DIM*(ELEMENT_DIM+1), PROBLEM_DIM*(ELEMENT_DIM+1)> matrix_term = zero_matrix<double>(PROBLEM_DIM*(ELEMENT_DIM+1), PROBLEM_DIM*(ELEMENT_DIM+1));

    // Loop over PDEs and populate matrix_term
    for (unsigned pde_index=0; pde_index<PROBLEM_DIM; pde_index++)
    {
        double this_dudt_coefficient = mpPdeSystem->ComputeDuDtCoefficientFunction(rX, pde_index);

        // in general this should be looking to the domain field to determine diffusion properties

        c_matrix<double, SPACE_DIM, SPACE_DIM> this_pde_diffusion_term = mpPdeSystem->ComputeDiffusionTerm(rX, pde_index, pElement);
        
        c_matrix<double, 1*(ELEMENT_DIM+1), 1*(ELEMENT_DIM+1)> this_stiffness_matrix =
            prod(trans(rGradPhi), c_matrix<double, SPACE_DIM, ELEMENT_DIM+1>(prod(this_pde_diffusion_term, rGradPhi)) )
                + timestep_inverse * this_dudt_coefficient * outer_prod(rPhi, rPhi);

        for (unsigned i=0; i<ELEMENT_DIM+1; i++)
        {
            for (unsigned j=0; j<ELEMENT_DIM+1; j++)
            {
                matrix_term(i*PROBLEM_DIM + pde_index, j*PROBLEM_DIM + pde_index) = this_stiffness_matrix(i,j);
            }
        }
    }
  //  std::cout<<"InhomogenousCoupledPdeOdeCoupledCellSolver - ComputeMatrixTerm - end"<<std::endl;
    return matrix_term;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
c_vector<double, PROBLEM_DIM*(ELEMENT_DIM+1)> InhomogenousCoupledPdeOdeCoupledCellSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::ComputeVectorTerm(
    c_vector<double, ELEMENT_DIM+1>& rPhi,
    c_matrix<double, SPACE_DIM, ELEMENT_DIM+1>& rGradPhi,
    ChastePoint<SPACE_DIM>& rX,
    c_vector<double,PROBLEM_DIM>& rU,
    c_matrix<double,PROBLEM_DIM,SPACE_DIM>& rGradU,
    Element<ELEMENT_DIM, SPACE_DIM>* pElement)
{ //  std::cout<<"InhomogenousCoupledPdeOdeCoupledCellSolver - ComputeVectorTerm - start"<<std::endl;
    double timestep_inverse = PdeSimulationTime::GetPdeTimeStepInverse();
    c_vector<double, PROBLEM_DIM*(ELEMENT_DIM+1)> vector_term;
    vector_term = zero_vector<double>(PROBLEM_DIM*(ELEMENT_DIM+1));

    // Loop over PDEs and populate vector_term
    for (unsigned pde_index=0; pde_index<PROBLEM_DIM; pde_index++)
    {
   
        // property of the pde not the ode systems
        double this_dudt_coefficient = mpPdeSystem->ComputeDuDtCoefficientFunction(rX, pde_index);

        // property of the ode systems; returns indexed state of the interpolated state variables, odes already solved
        double this_source_term = mpPdeSystem->ComputeSourceTerm(rX, rU, mInterpolatedOdeStateVariables, pde_index);

        double this_constant_source_term =0.0;

        // check whether adding cell contribution
        if(mCellTransportOdeSystemsPresent || mCellMembraneOdeSystemsPresent)
        {  
          
            // at least one of the cells has the transport property and can run an transport ode system which needs to be considered
            // check whether point rX is associated with a cell
            bool isContributed=false;
            unsigned cell_count=0;
            for (AbstractCellPopulation<2>::Iterator cell_iter = mrCellPopulation.Begin();
             cell_iter != mrCellPopulation.End();
             ++cell_iter)
            {   

                unsigned cell_location_index = mrCellPopulation.GetLocationIndexUsingCell(*cell_iter);
                const ChastePoint<SPACE_DIM>& cellCentrePoint = mrCellPopulation.GetLocationOfCellCentre(*cell_iter);
                
                if (cell_iter->rGetCellPropertyCollection().HasProperty<TransportCellProperty>())
                {
                   
                    boost::shared_ptr<TransportCellProperty> transport_cell_property = boost::static_pointer_cast<TransportCellProperty>(cell_iter->rGetCellPropertyCollection().GetPropertiesType<TransportCellProperty>().GetProperty());
 
                    if(false)//cell_iter->HasCellProperty<ExtendedCellProperty<SPACE_DIM>>())
                    {
                        boost::shared_ptr<ExtendedCellProperty<SPACE_DIM>> extended_cell_property = boost::static_pointer_cast<ExtendedCellProperty<SPACE_DIM>>(cell_iter->rGetCellPropertyCollection().GetPropertiesType<ExtendedCellProperty<SPACE_DIM>>().GetProperty());

                        if(extended_cell_property -> IsPointInCell(cellCentrePoint, rX))
                        {
                            // point is in cell, test if it is a cell boundary point
                            if(extended_cell_property -> CheckCellBoundary(cellCentrePoint, rX))
                            {
                                // if so then add the transport contribution
                                std::string name_of_pde_index = mpPdeSystem->GetStateVariableRegister()->RetrieveStateVariableName(pde_index);
                                

                                // if cell is on boundary add the result of the transport Ode system
                                this_constant_source_term = extended_cell_property ->RetrieveBoundarySourceByStateName(name_of_pde_index,mX);

                                if(!transport_cell_property ->GetIncludeOdeInterpolationOnBoundary())
                                {
                                    // override the ode source term
                                    //this_source_term =0.0;
                                }
                            }
                            else
                            {
                                // point rX is in cell
                                std::string name_of_pde_index = mpPdeSystem->GetStateVariableRegister()->RetrieveStateVariableName(pde_index);
              //                  this_constant_source_term += (extended_cell_property->RetrieveInternalCellSourceByStateName(name_of_pde_index)/(double)(ELEMENT_DIM+1));
                                if(!extended_cell_property -> GetIncludeOdeInterpolationInCell())
                                {
                                    // source is null
                                    //this_source_term = 0.0;
                                }
                                // else source term is the ode term is present; see above
                            }
                            // point rX can only be associated with one cell at a time, so here it has been found so break
                            break;
                        }
                    
                    }
                    else
                    {
                        // for case that cell has no extended property, can check that rX is close to cell_location
                        // then use the concentrations, bulk and cell, with the transport property
                        if(CheckChastePointsForEquality(cellCentrePoint,rX))
                        {
                            // if so then add the transport contribution
                            std::string name_of_pde_index = mpPdeSystem->GetStateVariableRegister()->RetrieveStateVariableName(pde_index);
         
                            // if cell is on boundary add the result of the transport Ode system
                            // scale the constant source term to remove the additive affect of the multiple element nodes 
                            //this_constant_source_term += (transport_cell_property ->RetrieveBoundarySourceByStateName(name_of_pde_index)/(double)(ELEMENT_DIM+1));
                            
                            this_constant_source_term += (transport_cell_property ->RetrieveChangeBoundarySourceByStateName(name_of_pde_index)/(double)(ELEMENT_DIM+1));
                            
                            if(!transport_cell_property ->GetIncludeOdeInterpolationOnBoundary())
                            {
                                // override the ode source term
                                //this_source_term =0.0;
                            }
                            isContributed =true; // as already found the only possible cell which may act as a source
                        }
                        // else not this cell
                    }
                } 
                
                // check for membrane system
                if (cell_iter->rGetCellPropertyCollection().HasProperty<MembraneCellProperty>())
                {
                   
                    boost::shared_ptr<MembraneCellProperty> membrane_cell_property = boost::static_pointer_cast<MembraneCellProperty>(cell_iter->rGetCellPropertyCollection().GetPropertiesType<MembraneCellProperty>().GetProperty());
                    // check that rX is close to cell_location
                    // then use the concentrations, bulk and cell, with the membrane property
                    if(CheckChastePointsForEquality(cellCentrePoint,rX))
                    {
                        // if so then add the membrane contribution
                        std::string name_of_pde_index = mpPdeSystem->GetStateVariableRegister()->RetrieveStateVariableName(pde_index);
                  
                        // if cell is on boundary add the result of the membrane Ode system
                        // scale the constant source term to remove the additive affect of the multiple element nodes 
                        this_constant_source_term += (membrane_cell_property ->RetrieveChangeBoundarySourceByStateName(name_of_pde_index)/(double)(ELEMENT_DIM+1));

                        if(!membrane_cell_property ->GetIncludeMembraneOdeInterpolationOnBoundary())
                        {
                            // override the ode source term
                            //this_source_term =0.0;
                        }
                        isContributed =true; // as already found the only possible cell which may act as a source
                    }
                    // else not this cell
                }



                
                
                if(isContributed)
                {
                    // only one cell may contribute at a single point
                    //std::cout<<this_constant_source_term<<std::endl;
                    break;
                }
                // else ignore cell; the cell has no transport property
                cell_count++;
            }
        }


        // include the cell contribution alongside the interpolated ode contribution
        c_vector<double, ELEMENT_DIM+1> this_vector_term;
        
        this_vector_term = (this_source_term + this_constant_source_term + timestep_inverse*this_dudt_coefficient*rU(pde_index))* rPhi;

        for (unsigned i=0; i<ELEMENT_DIM+1; i++)
        {
            if(this_vector_term(i)<0)
            {
                this_vector_term(i)=0.0;
            }
            vector_term(i*PROBLEM_DIM + pde_index) = this_vector_term(i);
        }
    }
 //   std::cout<<"InhomogenousCoupledPdeOdeCoupledCellSolver - ComputeVectorTerm - end"<<std::endl;
    return vector_term;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void InhomogenousCoupledPdeOdeCoupledCellSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::ResetInterpolatedQuantities()                
{  // std::cout<<"InhomogenousCoupledPdeOdeCoupledCellSolver - ResetInterpolatedQuantities"<<std::endl;
    mInterpolatedOdeStateVariables.clear();

    if (mOdeSystemsPresent)
    {
        unsigned num_state_variables = mpPdeSystem->GetStateVariableRegister()->GetNumberOfStateVariables();
        
        mInterpolatedOdeStateVariables.resize(num_state_variables, 0.0);
    }

    if(mCellTransportOdeSystemsPresent)
    {
        // reset the point and current state vector
        ChastePoint<SPACE_DIM> x(0,0,0);

        SetInterpolationPoint(x);

        c_vector<double,PROBLEM_DIM> u = zero_vector<double>(PROBLEM_DIM);

        SetCurrentStateVector(u);

        mInterpolationCount=0;
        for (AbstractCellPopulation<2>::Iterator cell_iter = mrCellPopulation.Begin();
        cell_iter != mrCellPopulation.End();
        ++cell_iter)
        {
            if (false)//cell_iter->HasCellProperty<ExtendedCellProperty<SPACE_DIM>>())
            {
                boost::shared_ptr<ExtendedCellProperty<SPACE_DIM>> extended_cell_property = boost::static_pointer_cast<ExtendedCellProperty<SPACE_DIM>>(cell_iter->rGetCellPropertyCollection().GetPropertiesType<ExtendedCellProperty<SPACE_DIM>>().GetProperty());

                extended_cell_property -> ResetVectorOfBoundaryStateVariables();
                extended_cell_property -> ResetVectorOfInternalBoundaryStateVariables();
                extended_cell_property -> ResetVectorOfBoundaryLocations();

            }
            // if cell has only the transport property then there is nothing needed to be reset
        }

    }
  //  std::cout<<"InhomogenousCoupledPdeOdeCoupledCellSolver - ResetInterpolatedQuantities - end"<<std::endl;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void InhomogenousCoupledPdeOdeCoupledCellSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::IncrementInterpolatedQuantities(double phiI, const Node<SPACE_DIM>* pNode)
{   //std::cout<<"InhomogenousCoupledPdeOdeCoupledCellSolver - IncrementInterpolatedQuantities"<<std::endl;
    // interploates a quantity from a node location to point x through the basis function phi associated with the node 
    if (mOdeSystemsPresent)
    {
        unsigned matchedIndex;
        unsigned num_state_variables = mpPdeSystem->GetStateVariableRegister()->GetNumberOfStateVariables();
        std::vector<std::string> pde_variable_register = mpPdeSystem->GetStateVariableRegister() -> GetStateVariableRegisterVector();

        for (unsigned i=0; i<num_state_variables; i++)
        {
            // each node may not have the total number of states in the ode system
            // select using the state variable registers at the nodes, if state variable isn't present the value is 0

            if(mOdeSystemsAtNodes[pNode->GetIndex()] -> GetStateVariableRegister() -> IsStateVariablePresent(pde_variable_register[i]))
            {
                matchedIndex = mOdeSystemsAtNodes[pNode->GetIndex()] -> GetStateVariableRegister() -> RetrieveStateVariableIndex(pde_variable_register[i]);

                //rGetStateVariables() returns the current values of the state variables, index
                mInterpolatedOdeStateVariables[i] += phiI * mOdeSystemsAtNodes[pNode->GetIndex()]->rGetStateVariables()[matchedIndex];
                // rGetStateVariables if the value of the ode states variables at the nodes, size < num_state_vars in the pde system
            }
            else
            {
                // if the state variable isn't a part of the ODE then interpolated the pde solution at the node
                mInterpolatedOdeStateVariables[i] += phiI * 0.0;

            }
                                
        }

    }

    // check whether adding cell contribution
    if(mCellTransportOdeSystemsPresent || mCellMembraneOdeSystemsPresent)
    {

        // interpolate mX and mU
        mX.rGetLocation() += phiI*pNode->rGetLocation();
        for(unsigned prob_dim=0; prob_dim<PROBLEM_DIM; prob_dim++)
        {
            mU(prob_dim) += phiI*this->GetCurrentSolutionOrGuessValue(pNode->GetIndex(), prob_dim);
        }
        
        // store mX and mU over the interpolation; assume ELEMENT_DIM + 1 nodes for interpolation
     //   std::cout<<"mInterpolationCount = "<<mInterpolationCount<<std::endl;
        if(mInterpolationCount==2)
        {
            // check for -ve concentrations
            for(unsigned prob_dim=0; prob_dim<PROBLEM_DIM; prob_dim++)
            {
                if(mU(prob_dim)<0.0)
                {
                    mU(prob_dim) = 0.0;
                }
            }


            // the mX and mU are fully interpolated
     //       std::cout<<"hit interpolation"<<std::endl;
            // check whether point mX is associated with a cell
            bool this_point_found = false;
            for (AbstractCellPopulation<2>::Iterator cell_iter = mrCellPopulation.Begin();
            cell_iter != mrCellPopulation.End();
            ++cell_iter)
            {
               
                if(cell_iter->rGetCellPropertyCollection().HasProperty<TransportCellProperty>())
                {
          //          std::cout<<"transport"<<std::endl;
                    const ChastePoint<SPACE_DIM>& cellCentrePoint = mrCellPopulation.GetLocationOfCellCentre(*cell_iter);

                    //ChastePoint<SPACE_DIM> cellCentrePoint(cell_location);

                    if (false)//cell_iter->HasCellProperty<ExtendedCellProperty<SPACE_DIM>>())
                    {

                        boost::shared_ptr<ExtendedCellProperty<SPACE_DIM>> extended_cell_property = boost::static_pointer_cast<ExtendedCellProperty<SPACE_DIM>>(cell_iter->rGetCellPropertyCollection().GetPropertiesType<ExtendedCellProperty<SPACE_DIM>>().GetProperty());

                        if(extended_cell_property -> IsPointInCell(cellCentrePoint, mX))
                        {
                            // ignore point if it is in the cell and not on boundary for time being. The internal cell points are handled differently.
                            if(extended_cell_property -> IsPointOnCellBoundary(cellCentrePoint, mX))
                            {
                                // need mU as std::vector<double>
                                std::vector<double> mUstd(PROBLEM_DIM,0.0);
                                for(unsigned i=0; i<PROBLEM_DIM;i++)
                                {
                                    mUstd[i] = mU[i];
                                }
                                extended_cell_property -> RecordLocationAndStateVariable(mX,mUstd);
                                mUstd.clear();
                            }
                            this_point_found = true; // found where the point rX is associated
                        }
                    }
                    else
                    {
                        boost::shared_ptr<TransportCellProperty> transport_cell_property = boost::static_pointer_cast<TransportCellProperty>(cell_iter->rGetCellPropertyCollection().GetPropertiesType<TransportCellProperty>().GetProperty());

                        // check whether the interpolated point is on the cell centre
             //           std::cout<<"Check for equality"<<std::endl;
                        if(CheckChastePointsForEquality(cellCentrePoint, mX))
                        {
              //              std::cout<<"IncrementInterpolatedQuantities - CheckChastePointsForEquality "<<std::endl;
                            // set the value for rU as the mBulkBoundaryConcentrationVector
                            // need mU as std::vector<double>
                            std::vector<double> mUstd(PROBLEM_DIM,0.0);
                            for(unsigned i=0; i<PROBLEM_DIM;i++)
                            {
                                mUstd[i] = mU[i];
                               
                            }
                            transport_cell_property -> SetBulkBoundaryConcentrationVector(mUstd);
                            transport_cell_property -> SetInitBulkBoundaryConcentrationVector(mUstd); //15/10/2020
                            transport_cell_property -> ResetReactionCalls();
                            mUstd.clear();

                            this_point_found = true; // as have found the cell associated with the point interpolated
                        }

                    }
                    
                }
          
                if(cell_iter->rGetCellPropertyCollection().HasProperty<MembraneCellProperty>())
                {

                    const ChastePoint<SPACE_DIM>& cellCentrePoint = mrCellPopulation.GetLocationOfCellCentre(*cell_iter);

                    boost::shared_ptr<MembraneCellProperty> membrane_cell_property = boost::static_pointer_cast<MembraneCellProperty>(cell_iter->rGetCellPropertyCollection().GetPropertiesType<MembraneCellProperty>().GetProperty());

                    // check whether the interpolated point is on the cell centre

                    if(CheckChastePointsForEquality(cellCentrePoint, mX))
                    {

                        // set the value for rU as the mBulkBoundaryConcentrationVector
                        // need mU as std::vector<double>
                        std::vector<double> mUstd(PROBLEM_DIM,0.0);
                        for(unsigned i=0; i<PROBLEM_DIM;i++)
                        {
                            mUstd[i] = mU[i];
                            
                        }
                        membrane_cell_property -> SetBulkBoundaryConcentrationVector(mUstd);
                        membrane_cell_property -> SetInitBulkBoundaryConcentrationVector(mUstd); //15/10/2020
                        membrane_cell_property -> ResetReactionCalls();
                        mUstd.clear();

                        this_point_found = true; // as have found the cell associated with the point interpolated
                    } 
                }

                // update the cell data based on any environment properties
                if (cell_iter-> rGetCellPropertyCollection().HasProperty<EnvironmentCellProperty>())
                {
                    boost::shared_ptr<EnvironmentCellProperty> environment_cell_property = boost::static_pointer_cast<EnvironmentCellProperty>(cell_iter-> rGetCellPropertyCollection(). GetPropertiesType<EnvironmentCellProperty>().GetProperty());

                    StateVariableRegister* pEnvironmentRegister = environment_cell_property -> GetEnvironmentStateVariableRegister();

                    const ChastePoint<SPACE_DIM>& cellCentrePoint = mrCellPopulation.GetLocationOfCellCentre(*cell_iter);


                    if(CheckChastePointsForEquality(cellCentrePoint,mX))
                    {
                        std::string name_of_pde_index ="";
                        for(unsigned pde_index=0; pde_index<PROBLEM_DIM;pde_index++)
                        {
                            name_of_pde_index = mpPdeSystem->GetStateVariableRegister()->RetrieveStateVariableName(pde_index);
                            if(pEnvironmentRegister->IsStateVariablePresent(name_of_pde_index))
                            {
                                std::cout<<"Environment: "<<name_of_pde_index<<" : "<<mU[pde_index]<<std::endl;
                                std::cout<<"Preferred: "<<name_of_pde_index<<" : "<<environment_cell_property->GetPreferredEnvironmentValueByName(name_of_pde_index)<<std::endl;
                                environment_cell_property->SetEnvironmentValueByIndex(pEnvironmentRegister->RetrieveStateVariableIndex(name_of_pde_index),mU[pde_index]);
                            }
                            
                        }
                        this_point_found = true; // as have found the cell associated with the point interpolated
                    }


                }
                if(this_point_found)
                {
                    break;
                }
            }
            // reset the interpolation count
            mInterpolationCount=0;
        }
        else
        {
            // not fully interpolated the point and state vector yet
            mInterpolationCount +=1;
        }
    }
    //std::cout<<"InhomogenousCoupledPdeOdeCoupledCellSolver - IncrementInterpolatedQuantities - end"<<std::endl;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void InhomogenousCoupledPdeOdeCoupledCellSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::InitialiseForSolve(Vec initialSolution)
{  // std::cout<<"InhomogenousCoupledPdeOdeCoupledCellSolver - InitialiseForSolve"<<std::endl;
    if (this->mpLinearSystem == NULL)
    {
        unsigned preallocation = mpFeMesh->CalculateMaximumContainingElementsPerProcess() + ELEMENT_DIM;
        if (ELEMENT_DIM > 1)
        {
            // Highest connectivity is closed
            preallocation--;
        }
        preallocation *= PROBLEM_DIM;

        /*
         * Use the current solution (ie the initial solution) as the
         * template in the alternative constructor of LinearSystem.
         * This is to avoid problems with VecScatter.
         */
        this->mpLinearSystem = new LinearSystem(initialSolution, preallocation);
    }

    assert(this->mpLinearSystem);
    this->mpLinearSystem->SetMatrixIsSymmetric(true);
    this->mpLinearSystem->SetKspType("cg");

    // set up each of the cells
    if(mCellTransportOdeSystemsPresent)
    {
       
        // reset the point and current state vector
        ChastePoint<SPACE_DIM> x(0,0,0);

        SetInterpolationPoint(x);

        c_vector<double,PROBLEM_DIM> u = zero_vector<double>(PROBLEM_DIM);

        SetCurrentStateVector(u);
     
        for (AbstractCellPopulation<2>::Iterator cell_iter = mrCellPopulation.Begin();
        cell_iter != mrCellPopulation.End();
        ++cell_iter)
        {   
            if (false)//cell_iter -> rGetCellPropertyCollection().HasProperty<ExtendedCellProperty>())
            {
                boost::shared_ptr<ExtendedCellProperty<SPACE_DIM>> extended_cell_property = boost::static_pointer_cast<ExtendedCellProperty<SPACE_DIM>>(cell_iter->rGetCellPropertyCollection().GetPropertiesType<ExtendedCellProperty<SPACE_DIM>>().GetProperty());

                unsigned num_voxels = extended_cell_property->GetTotalNumberMeshVoxels();
                if(num_voxels == 0)
                {
                    num_voxels=1;
                }
                std::vector<double> initialCellConcentration = extended_cell_property->GetCellConcentrationVector();

                for(unsigned state=0; state<initialCellConcentration.size(); state++)
                {
                    // scale by the number of voxels
                    initialCellConcentration[state] /= num_voxels;
                }

                extended_cell_property -> ResetVectorOfBoundaryStateVariables();
                extended_cell_property -> ResetVectorOfInternalBoundaryStateVariables();
                extended_cell_property -> ResetVectorOfBoundaryLocations();
                extended_cell_property -> ResetNextTimestepConcentrationVector(PROBLEM_DIM);

                for(unsigned boundary_index=0; boundary_index<extended_cell_property->GetNumberMeshVoxelsOnBoundary(); boundary_index++)
                {
                    extended_cell_property -> RecordLocationAndStateVariable(x,initialCellConcentration);

                }
            }
        }

    }

  //  std::cout<<"InhomogenousCoupledPdeOdeCoupledCellSolver - InitialiseForSolve -end"<<std::endl;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void InhomogenousCoupledPdeOdeCoupledCellSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetupLinearSystem(Vec currentSolution, bool computeMatrix)
{  // std::cout<<"InhomogenousCoupledPdeOdeCoupledCellSolver - SetupLinearSystem"<<std::endl;
    this->SetupGivenLinearSystem(currentSolution, computeMatrix, this->mpLinearSystem);
  //  std::cout<<"InhomogenousCoupledPdeOdeCoupledCellSolver - SetupLinearSystem -end"<<std::endl;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
InhomogenousCoupledPdeOdeCoupledCellSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::InhomogenousCoupledPdeOdeCoupledCellSolver(
        TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>* pMesh,
        InhomogenousParabolicPdeForCoupledOdeSystemTemplated<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>* pPdeSystem,
        BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>* pBoundaryConditions,
        AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation,
        std::vector<AbstractInhomogenousOdeSystemForCoupledPdeSystem*> odeSystemsAtNodes,
        std::vector<boost::shared_ptr<AbstractIvpOdeSolver>> pOdeSolvers,
        bool conditionsInterpolated,
        ChastePoint<SPACE_DIM> point,
        c_vector<double,PROBLEM_DIM> stateVector
        )
    : AbstractAssemblerSolverHybrid<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM, NORMAL>(pMesh, pBoundaryConditions),
      AbstractDynamicLinearPdeSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>(pMesh),
      mpFeMesh(pMesh),
      mpPdeSystem(pPdeSystem),
      mrCellPopulation(rCellPopulation),
      mOdeSystemsAtNodes(odeSystemsAtNodes),
      mpOdeSolvers(pOdeSolvers),
      mConditionsInterpolated(conditionsInterpolated),
      mX(point),
      mU(stateVector),
      mSamplingTimeStep(DOUBLE_UNSET),
      mOdeSystemsPresent(false),
      mCellTransportOdeSystemsPresent(false),
      mCellMembraneOdeSystemsPresent(false),
      mClearOutputDirectory(false)
{  // std::cout<<"InhomogenousCoupledPdeOdeCoupledCellSolver - InhomogenousCoupledPdeOdeCoupledCellSolver - start"<<std::endl;
    this->mpBoundaryConditions = pBoundaryConditions;

    /*
     * If any ODE systems are passed in to the constructor, then we aren't just
     * solving a coupled PDE system, in which case the number of ODE system objects
     * must match the number of nodes in the finite element mesh.
     */
    
    if (!mOdeSystemsAtNodes.empty())
    {
      
        mOdeSystemsPresent = true;
        assert(mOdeSystemsAtNodes.size() == mpFeMesh->GetNumNodes());

        /*
         * In this case, if an ODE solver is not explicitly passed into the 
         * constructor, then we create a default solver.
         */
        if (!mpOdeSolvers[0])
        {
#ifdef CHASTE_CVODE
            for(unsigned i=0; i<mpOdeSolvers.size(); i++)
            {
                mpOdeSolvers[i].reset(new CvodeAdaptor);
            }
            
#else
            for(unsigned i=0; i<mOdeSystemsAtNodes.size(); i++)
            {
                mpOdeSolvers.push_back(new BackwardEulerIvpOdeSolver(mOdeSystemsAtNodes[i]->GetNumberOfStateVariables()));
            }
            
#endif //CHASTE_CVODE
        }
    }

//std::cout<<"coupledPdeCellSolver set membrane property"<<std::endl;
    // run through the cell popualtion for the occurance of a cell with the membrane property defined
    for (AbstractCellPopulation<2>::Iterator cell_iter = rCellPopulation.Begin();
        cell_iter != rCellPopulation.End();
        ++cell_iter)
    {
   //     std::cout<<"test membrane "<<cell_iter->HasCellProperty<MembraneCellProperty>()<<std::endl;

        //if (cell_iter-> template HasCellProperty<MembraneCellProperty>())
        
        if (cell_iter-> rGetCellPropertyCollection().HasProperty<MembraneCellProperty>())
        {
            // if there exists at least one cell with the membrane property then solve transport odes
            mCellMembraneOdeSystemsPresent=true;
            break;
        }
    }

  //  std::cout<<"coupledPdeCellSolver set transport property"<<std::endl;
    // run through the cell popualtion for the occurance of a cell with the transport property defined
    for (AbstractCellPopulation<2>::Iterator cell_iter = rCellPopulation.Begin();
        cell_iter != rCellPopulation.End();
        ++cell_iter)
    {
  
        if (cell_iter-> rGetCellPropertyCollection().HasProperty<TransportCellProperty>())
        {
  //          std::cout<<"has transport"<<std::endl;
            // if there exists at least one cell with the transport property then solve transport odes
            mCellTransportOdeSystemsPresent=true;
            break;
        }
    }

    //std::cout<<"InhomogenousCoupledPdeOdeCoupledCellSolver - InhomogenousCoupledPdeOdeCoupledCellSolver - end"<<std::endl;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
InhomogenousCoupledPdeOdeCoupledCellSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::~InhomogenousCoupledPdeOdeCoupledCellSolver()
{
    //std::cout<<"InhomogenousCoupledPdeOdeCoupledCellSolver - ~InhomogenousCoupledPdeOdeCoupledCellSolver"<<std::endl;

    /*  don't delete as solver will be called multiple times on same odes, assume that the ode FeMesh remains unchanged
    if (mOdeSystemsPresent)
    {
        for (unsigned i=0; i<mOdeSystemsAtNodes.size(); i++)
        {
            delete mOdeSystemsAtNodes[i];
        }
    }
    */
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void InhomogenousCoupledPdeOdeCoupledCellSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::PrepareForSetupLinearSystem(Vec currentPdeSolution)
{   
// std::cout<<"InhomogenousCoupledPdeOdeCoupledCellSolver - PrepareForSetupLinearSystem - strat"<<std::endl;
    if (mOdeSystemsPresent)
    {
        double time = PdeSimulationTime::GetTime();
        double next_time = PdeSimulationTime::GetNextTime();
        double dt = PdeSimulationTime::GetPdeTimeStep();

        ReplicatableVector soln_repl(currentPdeSolution);
        std::vector<double> current_soln_this_node(PROBLEM_DIM);

        // Loop over nodes
        for (unsigned node_index=0; node_index<mpFeMesh->GetNumNodes(); node_index++)
        { 
            // Store the current solution to the PDE system at this node
            for (unsigned pde_index=0; pde_index<PROBLEM_DIM; pde_index++)
            {   
                double current_soln_this_pde_this_node = soln_repl[PROBLEM_DIM*node_index + pde_index];

                current_soln_this_node[pde_index] = current_soln_this_pde_this_node;

            }

            // Pass it into the ODE system at this node, of full state space dimensions
            mOdeSystemsAtNodes[node_index]->SetPdeSolution(current_soln_this_node);

            // Solve ODE system at this node
            mpOdeSolvers[node_index]->SolveAndUpdateStateVariable(mOdeSystemsAtNodes[node_index], time, next_time, dt);

        }
    }

    // run cell transport chemical ode

    if (mConditionsInterpolated && (mCellTransportOdeSystemsPresent||mCellMembraneOdeSystemsPresent))
    {
        //std::cout<<"run cell transport"<<std::endl;
        double time = PdeSimulationTime::GetTime();
        double next_time = PdeSimulationTime::GetNextTime();
        double dt = PdeSimulationTime::GetPdeTimeStep();
        //std::cout<<"Here - 4 -time: "<<time<<std::endl;
        for (AbstractCellPopulation<2>::Iterator cell_iter = mrCellPopulation.Begin();
             cell_iter != mrCellPopulation.End();
             ++cell_iter)
        {
             //std::cout<<"for cell:"<<std::endl;
            if (cell_iter->rGetCellPropertyCollection().HasProperty<TransportCellProperty>())
            {
                
                boost::shared_ptr<TransportCellProperty> transport_property = boost::static_pointer_cast<TransportCellProperty>(cell_iter->rGetCellPropertyCollection().GetPropertiesType<TransportCellProperty>().GetProperty());
                if(cell_iter->rGetCellPropertyCollection().HasProperty<ExtendedCellProperty<SPACE_DIM>>())
                //if(false)//cell_iter->HasCellProperty<ExtendedCellProperty<SPACE_DIM>>())
                {
                    boost::shared_ptr<ExtendedCellProperty<SPACE_DIM>> extended_cell_property = boost::static_pointer_cast<ExtendedCellProperty<SPACE_DIM>>(cell_iter->rGetCellPropertyCollection().GetPropertiesType<ExtendedCellProperty<SPACE_DIM>>().GetProperty());
    
                    for(unsigned boundary_index=0; boundary_index<extended_cell_property->GetNumberMeshVoxelsOnBoundary(); boundary_index++)
                    {
                        // rY for bulk state variables at the location specified by the boundary_index
                        std::vector<double> rY = extended_cell_property -> GetVectorOfBoundaryStateVariablesByLocationIndex(boundary_index);
                        
                        // append the cell state variables
                        extended_cell_property -> AppendInternalCellBoundaryConcentrations(rY,boundary_index);
                        
                        // rY needs to be both bulk and cell state variables
                        transport_property->GetTransportOdeSolver()->Solve(transport_property->GetTransportOdeSystem(), rY, time, next_time, dt);

                        extended_cell_property -> ReplaceBoundaryStateVariables(boundary_index, rY);
                    }
                }
                else
                {   
                    // cell is not exetended and as such the cell state variables are not stored in the extended cell property 
                    // no voxels to iterate over
                    //std::cout<<"transport get external"<<std::endl;
                    // rY for bulk state variables at the location specified by the boundary_index
                    std::vector<double> rY = transport_property -> GetExternalCellBoundaryConcentrationVector();
                    //std::cout<<"transport append internal"<<std::endl;
                    // append the cell state variables
                    transport_property -> AppendInternalCellBoundaryConcentrations(rY);

                    // need to combine the outer nad inner values of concentrations either side of the cell into one rY vector
                    std::vector<double> rDY(rY.size(),0.0);
                    // whether we need to solve the ode or just perform the reaction system
                    //transport_property->GetTransportOdeSolver()->Solve(transport_property->GetTransportOdeSystem(), rY, time, next_time, dt);
                    transport_property->GetTransportOdeSystem()->SetPdeStateRegister(mpPdeSystem->GetStateVariableRegister());
                    //std::cout<<"transport evalualte"<<std::endl;
                    transport_property->GetTransportOdeSystem()->EvaluateYDerivatives(time, rY, rDY);
                    for(unsigned i=0; i<rDY.size(); i++)
                    {
                        rDY[i] = rDY[i]*dt;
                    }
                    
                    transport_property -> ReplaceBoundaryStateVariables(rY);
                    transport_property -> ReplaceChangeBoundaryStateVariables(rDY);
                }
            }
   
            if (cell_iter->rGetCellPropertyCollection().HasProperty<MembraneCellProperty>())
            {
                    
                boost::shared_ptr<MembraneCellProperty> membrane_property = boost::static_pointer_cast<MembraneCellProperty>(cell_iter->rGetCellPropertyCollection().GetPropertiesType<MembraneCellProperty>().GetProperty());
                // rY for bulk state variables at the location specified by the boundary_index
                std::vector<double> rY = membrane_property -> GetExternalCellBoundaryConcentrationVector();
 
                membrane_property -> AppendInternalCellBoundaryConcentrations(rY);
                std::vector<double> rDY(rY.size(),0.0);
                membrane_property->GetMembraneOdeSystem()->SetPdeStateRegister(mpPdeSystem->GetStateVariableRegister());
      
                membrane_property->GetMembraneOdeSystem()->EvaluateYDerivatives(time, rY, rDY);
         
                for(unsigned i=0; i<rDY.size(); i++)
                {
                    rDY[i] = rDY[i]*dt;
                }

                membrane_property -> ReplaceBoundaryStateVariables(rY);

                membrane_property -> ReplaceChangeBoundaryStateVariables(rDY);
            }
        }
    }
//std::cout<<"InhomogenousCoupledPdeOdeCoupledCellSolver - PrepareForSetupLinearSystem - end"<<std::endl;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void InhomogenousCoupledPdeOdeCoupledCellSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetOutputDirectory(std::string outputDirectory, bool clearDirectory)
{
    mClearOutputDirectory = clearDirectory;
    this->mOutputDirectory = outputDirectory;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void InhomogenousCoupledPdeOdeCoupledCellSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetSamplingTimeStep(double samplingTimeStep)
{
    assert(samplingTimeStep >= this->mIdealTimeStep);
    mSamplingTimeStep = samplingTimeStep;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void InhomogenousCoupledPdeOdeCoupledCellSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SolveAndWriteResultsToFile()
{
  //  std::cout<<"InhomogenousCoupledPdeOdeCoupledCellSolver - SolveAndWriteResultsToFile - start"<<std::endl;
   
    // A number of methods must have been called prior to this method
    if (this->mOutputDirectory == "")
    {
        EXCEPTION("SetOutputDirectory() must be called prior to SolveAndWriteResultsToFile()");
    }
    if (this->mTimesSet == false)
    {
        EXCEPTION("SetTimes() must be called prior to SolveAndWriteResultsToFile()");
    }
    if (this->mIdealTimeStep <= 0.0)
    {
        EXCEPTION("SetTimeStep() must be called prior to SolveAndWriteResultsToFile()");
    }
    if (mSamplingTimeStep == DOUBLE_UNSET)
    {
        EXCEPTION("SetSamplingTimeStep() must be called prior to SolveAndWriteResultsToFile()");
    }
    if (!this->mInitialCondition)
    {
        EXCEPTION("SetInitialCondition() must be called prior to SolveAndWriteResultsToFile()");
    }

#ifdef CHASTE_VTK
    // Create a .pvd output file

    OutputFileHandler output_file_handler(this->mOutputDirectory, mClearOutputDirectory);
    mpVtkMetaFile = output_file_handler.OpenOutputFile("results.pvd");
    *mpVtkMetaFile << "<?xml version=\"1.0\"?>\n";
    *mpVtkMetaFile << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">\n";
    *mpVtkMetaFile << "    <Collection>\n";

    

    // Write initial condition to VTK
    Vec initial_condition = this->mInitialCondition;
    WriteVtkResultsToFile(initial_condition, 0);

    // The helper class TimeStepper deals with issues such as small final timesteps so we don't have to
    TimeStepper stepper(this->mTstart, this->mTend, mSamplingTimeStep);

    // Main time loop
    while (!stepper.IsTimeAtEnd())
    {
        // Reset start and end times
        this->SetTimes(stepper.GetTime(), stepper.GetNextTime());

        // Solve the system up to the new end time
        Vec soln = this->Solve();

        // Reset the initial condition for the next timestep
        if (this->mInitialCondition != initial_condition)
        {
            PetscTools::Destroy(this->mInitialCondition);
        }
        this->mInitialCondition = soln;

        // Move forward in time
        stepper.AdvanceOneTimeStep();

        // Write solution to VTK
        WriteVtkResultsToFile(soln, stepper.GetTotalTimeStepsTaken());
    }

    // Restore saved initial condition to avoid user confusion!
    if (this->mInitialCondition != initial_condition)
    {
        PetscTools::Destroy(this->mInitialCondition);
    }
    this->mInitialCondition = initial_condition;

    // Close .pvd output file
    *mpVtkMetaFile << "    </Collection>\n";
    *mpVtkMetaFile << "</VTKFile>\n";
    mpVtkMetaFile->close();
#else //CHASTE_VTK
// LCOV_EXCL_START // We only test this in weekly builds
    WARNING("VTK is not installed and is required for this functionality");
// LCOV_EXCL_STOP
#endif //CHASTE_VTK

//    std::cout<<"InhomogenousCoupledPdeOdeCoupledCellSolver - SolveAndWriteResultsToFile - end"<<std::endl;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void InhomogenousCoupledPdeOdeCoupledCellSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::WriteVtkResultsToFile(Vec solution, unsigned numTimeStepsElapsed)
{
 //   std::cout<<"InhomogenousCoupledPdeOdeCoupledCellSolver - WriteVtkResultsToFile - start"<<std::endl;
#ifdef CHASTE_VTK


    // Create a new VTK file for this time step
    std::stringstream time;
    time << numTimeStepsElapsed;
    VtkMeshWriter<ELEMENT_DIM, SPACE_DIM> mesh_writer(this->mOutputDirectory, "results_"+time.str(), false);
    // need to ensure StateVariableRegister is defined
    std::vector<std::string> p_pde_stateVariableNames = mpPdeSystem -> GetStateVariableRegister() ->GetStateVariableRegisterVector();
    
    /*
     * We first loop over PDEs. For each PDE we store the solution
     * at each node in a vector, then pass this vector to the mesh
     * writer.
     */
    ReplicatableVector solution_repl(solution);
    unsigned num_nodes = mpFeMesh->GetNumNodes();
    for (unsigned pde_index=0; pde_index<PROBLEM_DIM; pde_index++)
    {
        // Store the solution of this PDE at each node
        std::vector<double> pde_index_data;
        pde_index_data.resize(num_nodes, 0.0);
        for (unsigned node_index=0; node_index<num_nodes; node_index++)
        {
            pde_index_data[node_index] = solution_repl[PROBLEM_DIM*node_index + pde_index];
        }

        // Add this data to the mesh writer
        std::stringstream data_name;
        data_name << "PDE variable " << p_pde_stateVariableNames[pde_index];
        mesh_writer.AddPointData(data_name.str(), pde_index_data);
    }
    
    if (mOdeSystemsPresent)
    {
        /*
         * We cannot loop over ODEs like PDEs, since the solutions are not
         * stored in one place. Therefore we build up a large 'vector of
         * vectors', then pass each component of this vector to the mesh
         * writer.
         */


        std::vector<std::vector<double> > ode_data;
        unsigned num_state_vars = mpPdeSystem -> GetStateVariableRegister() ->GetNumberOfStateVariables();
        ode_data.resize(num_state_vars);
        for (unsigned state_var_index=0; state_var_index<num_state_vars; state_var_index++)
        {
            ode_data[state_var_index].resize(num_nodes, 0.0);
        }

        for (unsigned node_index=0; node_index<num_nodes; node_index++)
        {
            std::vector<double> all_odes_this_node = mOdeSystemsAtNodes[node_index]->rGetStateVariables();
            // this could be of variable size, is of only the states that the ode modifies


            std::vector<std::string> ode_var_names = mOdeSystemsAtNodes[node_index]->GetStateVariableRegister() -> GetStateVariableRegisterVector();
            for (unsigned ode_index=0; ode_index<ode_var_names.size(); ode_index++)
            {
                // for each state variable in the ode system, find and update the corresponding variable in the pde system 
                ode_data[mpPdeSystem -> GetStateVariableRegister() -> RetrieveStateVariableIndex(ode_var_names[ode_index])][node_index] = all_odes_this_node[ode_index];
                
            }


            //---------------------------------------------------------------------------------
            // do we need to use the rPDEsolution for the other nodes or is 0 the correct (current) choice?
            //---------------------------------------------------------------------------------



        }

        for (unsigned ode_index=0; ode_index<num_state_vars; ode_index++)
        {
            // ode_index is the state variable
            std::vector<double> ode_index_data = ode_data[ode_index];

            // Add this data to the mesh writer
            std::stringstream data_name;
            data_name << "ODE variable " << p_pde_stateVariableNames[ode_index];
            mesh_writer.AddPointData(data_name.str(), ode_index_data);
        }
    }

    mesh_writer.WriteFilesUsingMesh(*mpFeMesh);
    *mpVtkMetaFile << "        <DataSet timestep=\"";
    *mpVtkMetaFile << numTimeStepsElapsed;
    *mpVtkMetaFile << "\" group=\"\" part=\"0\" file=\"results_";
    *mpVtkMetaFile << numTimeStepsElapsed;
    *mpVtkMetaFile << ".vtu\"/>\n";
#endif // CHASTE_VTK
  //  std::cout<<"InhomogenousCoupledPdeOdeCoupledCellSolver - SolveAndWriteResultsToFile - end"<<std::endl;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
AbstractOdeSystemForCoupledPdeSystem* InhomogenousCoupledPdeOdeCoupledCellSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::GetOdeSystemAtNode(unsigned index)
{
    return mOdeSystemsAtNodes[index];
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void InhomogenousCoupledPdeOdeCoupledCellSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetInterpolationPoint(ChastePoint<SPACE_DIM> point)
{
    mX = point;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void InhomogenousCoupledPdeOdeCoupledCellSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetCurrentStateVector(c_vector<double,PROBLEM_DIM> stateVector)
{
    mU = stateVector;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
bool InhomogenousCoupledPdeOdeCoupledCellSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::CheckChastePointsForEquality(ChastePoint<SPACE_DIM> rX1,ChastePoint<SPACE_DIM> rX2)
{   //std::cout<<"InhomogenousCoupledPdeOdeCoupledCellSolver - CheckChastePointsForEquality"<<std::endl;
    // function to check whether two chaste points are congruent
    // use half the scale dimension as an error

    // half of the voxel dimensions of the tetrahedral mesh, a value to which a point should be congruent
    std::vector<double> mesh_single_voxel_dimensions(SPACE_DIM,0.4); 

    bool equality=true;

    for(unsigned i=0; i<SPACE_DIM; i++)
    {
        if(rX1[i] < rX2[i] + mesh_single_voxel_dimensions[i] && rX1[i] >= rX2[i] - mesh_single_voxel_dimensions[i])
        {
            // both points are congruent within 0.5 of mesh scale in the ith dimension
            equality=true;
        }
        else
        {
            // as soon as a dimension where the points are disaparate is found end the function return false
            equality=false;
            break;
        }
        
    }
    return equality;
}


#endif