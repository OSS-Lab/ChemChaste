#ifndef INHOMOGENOUSCOUPLEDPDEODESOLVER_TEMPLATED_HPP_
#define INHOMOGENOUSCOUPLEDPDEODESOLVER_TEMPLATED_HPP_

#include "ChemChasteFeAssemblerCommon.hpp"
#include "ChemChasteVolumeAssembler.hpp"
#include "ChemChasteSurfaceAssembler.hpp"
#include "AbstractAssemblerSolverHybrid.hpp"
#include "AbstractDynamicLinearPdeSolver.hpp"
#include "AbstractLinearParabolicPdeSystemForCoupledOdeSystem.hpp"
#include "TetrahedralMesh.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "AbstractInhomogenousOdeSystemForCoupledPdeSystem.hpp"
#include "InhomogenousParabolicPdeForCoupledOdeSystem_templated.hpp"
#include "AbstractOdeSystemForCoupledPdeSystem.hpp"
#include "CvodeAdaptor.hpp"
#include "BackwardEulerIvpOdeSolver.hpp"
#include "Warnings.hpp"
#include "VtkMeshWriter.hpp"
#include "StateVariableRegister.hpp"


#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <boost/shared_ptr.hpp>

// Same as other class but the numberOfStateVariables to the ode terms are not all equal.
// The ode's change on a nodal basis

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM=ELEMENT_DIM, unsigned PROBLEM_DIM=1>
class InhomogenousCoupledPdeOdeSolverTemplated
    : public AbstractAssemblerSolverHybrid<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM, NORMAL>,
      public AbstractDynamicLinearPdeSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>
{
private:

    /** Pointer to the mesh. */
    AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>* mpMesh;

    /** The PDE system to be solved. */
    InhomogenousParabolicPdeForCoupledOdeSystemTemplated<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>* mpPdeSystem;

    /** Vector of pointers to ODE systems, defined at nodes. */
    std::vector<AbstractInhomogenousOdeSystemForCoupledPdeSystem*> mOdeSystemsAtNodes;

    /** The values of the ODE system state variables, interpolated at a quadrature point. */
    std::vector<double> mInterpolatedOdeStateVariables;

    /** The ODE solvers. */
    std::vector<boost::shared_ptr<AbstractIvpOdeSolver>> mpOdeSolvers;

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


    // for calculating the total mass, interpolate the nodal solution and add to sum
    bool mCalculateTotalMass = false;
    
    std::vector<double> mTotalMass = {0.0};

    std::vector<double> mOldTotalMass = {0.0};


    bool mInterpolateStateVariable = false;

 //   std::vector<double> mValueAtX(PROBLEM_DIM,0.0);

   

    // for calculating the trace over a line segment
    bool mCalculateSliceMass = false;

    unsigned mInterpolationNodeCount =0;

    unsigned mNodesInElement =3;

    bool mInterpolatePosition=false;

  //  ChastePoint<SPACE_DIM> mX(0,0,0);

    std::vector<std::vector<double>> mSliceMass;

    std::vector<std::vector<double>> mSlicePositions;

    std::vector<double> mDomainPositionsX={0.0};

    std::vector<double> mDomainPositionsY={0.0};

    std::vector<std::vector<double>> mPositionMass = std::vector<std::vector<double>> ();


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
    InhomogenousCoupledPdeOdeSolverTemplated(TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>* pMesh,
                                    InhomogenousParabolicPdeForCoupledOdeSystemTemplated<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>* pPdeSystem,
                                    BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>* pBoundaryConditions,
                                    std::vector<AbstractInhomogenousOdeSystemForCoupledPdeSystem*> odeSystemsAtNodes=std::vector<AbstractInhomogenousOdeSystemForCoupledPdeSystem*>(),
                                    std::vector<boost::shared_ptr<AbstractIvpOdeSolver>> pOdeSolvers=std::vector<boost::shared_ptr<AbstractIvpOdeSolver>>());

    /**
     * Destructor.
     * If an ODE system is present, the pointers to the ODE system objects are deleted here.
     */
    ~InhomogenousCoupledPdeOdeSolverTemplated();

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
     * @param index the global index of a node in the mpMesh
     * @return mOdeSystemsAtNodes[index]
     */
    AbstractOdeSystemForCoupledPdeSystem* GetOdeSystemAtNode(unsigned index);

    void ResetInterpolatedQuantitiesOnNewTimeStep();

    void OutputInterpolatedQuantitiesOnTimestep();
};

///////////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////////

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
c_matrix<double, PROBLEM_DIM*(ELEMENT_DIM+1), PROBLEM_DIM*(ELEMENT_DIM+1)> InhomogenousCoupledPdeOdeSolverTemplated<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::ComputeMatrixTerm(
    c_vector<double, ELEMENT_DIM+1>& rPhi,
    c_matrix<double, SPACE_DIM, ELEMENT_DIM+1>& rGradPhi,
    ChastePoint<SPACE_DIM>& rX,
    c_vector<double,PROBLEM_DIM>& rU,
    c_matrix<double, PROBLEM_DIM, SPACE_DIM>& rGradU,
    Element<ELEMENT_DIM, SPACE_DIM>* pElement)
{
//std::cout<<"InhomogenousCoupledPdeOdeSolverTemplated<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::ComputeMatrixTerm( - start"<<std::endl;
    double timestep_inverse = PdeSimulationTime::GetPdeTimeStepInverse();
    c_matrix<double, PROBLEM_DIM*(ELEMENT_DIM+1), PROBLEM_DIM*(ELEMENT_DIM+1)> matrix_term = zero_matrix<double>(PROBLEM_DIM*(ELEMENT_DIM+1), PROBLEM_DIM*(ELEMENT_DIM+1));

    // Loop over PDEs and populate matrix_term
    for (unsigned pde_index=0; pde_index<PROBLEM_DIM; pde_index++)
    {
        double this_dudt_coefficient = mpPdeSystem->ComputeDuDtCoefficientFunction(rX, pde_index);

        // in general this should be looking to the domain field to determine diffusion properties, interpolation required?

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
//std::cout<<"InhomogenousCoupledPdeOdeSolverTemplated<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::ComputeMatrixTerm( - end"<<std::endl;
    return matrix_term;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
c_vector<double, PROBLEM_DIM*(ELEMENT_DIM+1)> InhomogenousCoupledPdeOdeSolverTemplated<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::ComputeVectorTerm(
    c_vector<double, ELEMENT_DIM+1>& rPhi,
    c_matrix<double, SPACE_DIM, ELEMENT_DIM+1>& rGradPhi,
    ChastePoint<SPACE_DIM>& rX,
    c_vector<double,PROBLEM_DIM>& rU,
    c_matrix<double,PROBLEM_DIM,SPACE_DIM>& rGradU,
    Element<ELEMENT_DIM, SPACE_DIM>* pElement)
{
    //std::cout<<"InhomogenousCoupledPdeOdeSolverTemplated<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::ComputeVectorTerm( - start"<<std::endl;
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


        c_vector<double, ELEMENT_DIM+1> this_vector_term;
        this_vector_term = (this_source_term + timestep_inverse*this_dudt_coefficient*rU(pde_index))* rPhi;

        for (unsigned i=0; i<ELEMENT_DIM+1; i++)
        {
            vector_term(i*PROBLEM_DIM + pde_index) = this_vector_term(i);
        }
    }
    //std::cout<<"InhomogenousCoupledPdeOdeSolverTemplated<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::ComputeVectorTerm( - end"<<std::endl;
    return vector_term;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void InhomogenousCoupledPdeOdeSolverTemplated<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::ResetInterpolatedQuantities()
{
    mInterpolatedOdeStateVariables.clear();

    if (mOdeSystemsPresent)
    {
        unsigned num_state_variables = mpPdeSystem->GetStateVariableRegister()->GetNumberOfStateVariables();
        mInterpolatedOdeStateVariables.resize(num_state_variables, 0.0);
    }

    


    // reset for interpolation of slice mass
  //  ChastePoint<SPACE_DIM> tempX(0,0,0)
  //  mX = tempX;
    mInterpolationNodeCount = 0;

  //  std::vector<double> temp(PROBLEM_DIM,0.0);

   // mValueAtX = temp;


}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void InhomogenousCoupledPdeOdeSolverTemplated<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::IncrementInterpolatedQuantities(double phiI, const Node<SPACE_DIM>* pNode)
{
  //  std::cout<<"InhomogenousCoupledPdeOdeSolverTemplated::IncrementInterpolatedQuantities - start"<<std::endl;
    if (mOdeSystemsPresent)
    {   
        unsigned matchedIndex;
        unsigned num_state_variables = mpPdeSystem->GetStateVariableRegister()->GetNumberOfStateVariables();
        std::vector<std::string> pde_variable_register = mpPdeSystem->GetStateVariableRegister() -> GetStateVariableRegisterVector();

        for (unsigned i=0; i<num_state_variables; i++)
        {
            // each node may not have the total number of states in the ode system
            // select using the statevariable registers at the nodes, if state variable isn't present the value is 0

            if(mOdeSystemsAtNodes[pNode->GetIndex()] -> GetStateVariableRegister() -> IsStateVariablePresent(pde_variable_register[i]))
            {
                matchedIndex = mOdeSystemsAtNodes[pNode->GetIndex()] -> GetStateVariableRegister() -> RetrieveStateVariableIndex(pde_variable_register[i]);

                //rGetStateVariables() returns the current values of the state variables, index
                mInterpolatedOdeStateVariables[i] += phiI * mOdeSystemsAtNodes[pNode->GetIndex()]->rGetStateVariables()[matchedIndex];
                // rGetStateVariables if the value of the ode states variables at the nodes, size < num_state_vars in the pde system
            }
                                
        }

    }

    if(mCalculateTotalMass)
    {
        // interpolate the previous timestep solution at nodes
        // sum the interpolated value 

        unsigned num_state_variables = mpPdeSystem->GetStateVariableRegister()->GetNumberOfStateVariables();

        ReplicatableVector current_nodal_values(this->mInitialCondition);

        for (unsigned i=0; i<num_state_variables; i++)
        {
            mOldTotalMass[i] += phiI * current_nodal_values[PROBLEM_DIM*pNode->GetIndex()+i];
        }
    }

    if(mInterpolatePosition)
    {
     //   mX.rGetLocation() += phiI*pNode->rGetLocation();

    }

    if(mInterpolateStateVariable)
    {
    //    ReplicatableVector current_nodal_values(this->mInitialCondition);

    //    for (unsigned i=0; i<num_state_variables; i++)
     //   {
      //      mValueAtX[i] += phiI * current_nodal_values[PROBLEM_DIM*pNode->GetIndex()+i];
       // }

    }

    mInterpolationNodeCount++;
/*
    if(mInterpolationNodeCount == mNodesInElement)
    {
        //std::cout<<"here1"<<std::endl;
        // position and state variable vector are fully known
        ChastePoint<SPACE_DIM> X = AbstractFeVolumeIntegralAssembler<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM, true, true, NORMAL>::GetPosition();

        c_vector<double,PROBLEM_DIM> U = AbstractFeVolumeIntegralAssembler<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM, true, true, NORMAL>::GetStateVariable();

        //std::cout<<"here2"<<std::endl;
        if(mCalculateTotalMass)
        {
            for (unsigned i=0; i<PROBLEM_DIM; i++)
            {
                // should phiI be here?
                mTotalMass[i] += phiI * U[i];
            }
        }

        mDomainPositionsX.push_back(X[0]);

        mDomainPositionsY.push_back(X[1]);

        //std::cout<<"posiiton mass "<<mPositionMass.size()<<" "<<U.size()<<std::endl;
        for(unsigned j=0; j<PROBLEM_DIM; j++)
        {
            mPositionMass[j].push_back(mOldTotalMass[j]);
            //mPositionMass[j].push_back(U[j]);
        }
        //std::cout<<"edn posiiton mass"<<std::endl;
        if(mCalculateSliceMass)
        {
            //std::cout<<"y: "<<X[1]<<std::endl;
            if(2.85>X[1] && X[1]>2.8)
            {
                // store the x position and the state variable U
                std::vector<double> thisSlicePosition;
                std::vector<double> thisSliceMass;
                
                for(unsigned i=0; i<SPACE_DIM; i++)
                {
                    thisSlicePosition.push_back(X[i]);
                    //std::cout<<"x: "<<X[i]<<std::endl;
                }

                for(unsigned j=0; j<PROBLEM_DIM; j++)
                {
                    thisSliceMass.push_back(U[j]);
                }

                mSlicePositions.push_back(thisSlicePosition);
                mSliceMass.push_back(thisSliceMass);

            }
            // else ignore point
        }

    }
    */
//std::cout<<"InhomogenousCoupledPdeOdeSolverTemplated::IncrementInterpolatedQuantities - end"<<std::endl;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void InhomogenousCoupledPdeOdeSolverTemplated<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::InitialiseForSolve(Vec initialSolution)
{
    if (this->mpLinearSystem == NULL)
    {
        unsigned preallocation = mpMesh->CalculateMaximumContainingElementsPerProcess() + ELEMENT_DIM;
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
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void InhomogenousCoupledPdeOdeSolverTemplated<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetupLinearSystem(Vec currentSolution, bool computeMatrix)
{
    this->SetupGivenLinearSystem(currentSolution, computeMatrix, this->mpLinearSystem);
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
InhomogenousCoupledPdeOdeSolverTemplated<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::InhomogenousCoupledPdeOdeSolverTemplated(
        TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>* pMesh,
        InhomogenousParabolicPdeForCoupledOdeSystemTemplated<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>* pPdeSystem,
        BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>* pBoundaryConditions,
        std::vector<AbstractInhomogenousOdeSystemForCoupledPdeSystem*> odeSystemsAtNodes,
        std::vector<boost::shared_ptr<AbstractIvpOdeSolver>> pOdeSolvers)
    : AbstractAssemblerSolverHybrid<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM, NORMAL>(pMesh, pBoundaryConditions),
      AbstractDynamicLinearPdeSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>(pMesh),
      mpMesh(pMesh),
      mpPdeSystem(pPdeSystem),
      mOdeSystemsAtNodes(odeSystemsAtNodes),
      mpOdeSolvers(pOdeSolvers),
      mSamplingTimeStep(DOUBLE_UNSET),
      mOdeSystemsPresent(false),
      mClearOutputDirectory(false)
{
    this->mpBoundaryConditions = pBoundaryConditions;
    /*
     * If any ODE systems are passed in to the constructor, then we aren't just
     * solving a coupled PDE system, in which case the number of ODE system objects
     * must match the number of nodes in the finite element mesh.
     */
    if (!mOdeSystemsAtNodes.empty())
    {
        mOdeSystemsPresent = true;
        assert(mOdeSystemsAtNodes.size() == mpMesh->GetNumNodes());

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
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
InhomogenousCoupledPdeOdeSolverTemplated<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::~InhomogenousCoupledPdeOdeSolverTemplated()
{
    if (mOdeSystemsPresent)
    {
        for (unsigned i=0; i<mOdeSystemsAtNodes.size(); i++)
        {
            delete mOdeSystemsAtNodes[i];
        }
    }
    
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void InhomogenousCoupledPdeOdeSolverTemplated<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::PrepareForSetupLinearSystem(Vec currentPdeSolution)
{   
    //std::cout<<"InhomogenousCoupledPdeOdeSolverTemplated<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::PrepareForSetupLinearSystem( - start"<<std::endl;
    if (mOdeSystemsPresent)
    {
        double time = PdeSimulationTime::GetTime();
        double next_time = PdeSimulationTime::GetNextTime();
        double dt = PdeSimulationTime::GetPdeTimeStep();

        ReplicatableVector soln_repl(currentPdeSolution);
        std::vector<double> current_soln_this_node(PROBLEM_DIM);

        // Loop over nodes
        for (unsigned node_index=0; node_index<mpMesh->GetNumNodes(); node_index++)
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
    //std::cout<<"InhomogenousCoupledPdeOdeSolverTemplated<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::PrepareForSetupLinearSystem( - end"<<std::endl;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void InhomogenousCoupledPdeOdeSolverTemplated<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetOutputDirectory(std::string outputDirectory, bool clearDirectory)
{
    mClearOutputDirectory = clearDirectory;
    this->mOutputDirectory = outputDirectory;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void InhomogenousCoupledPdeOdeSolverTemplated<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetSamplingTimeStep(double samplingTimeStep)
{
    assert(samplingTimeStep >= this->mIdealTimeStep);
    mSamplingTimeStep = samplingTimeStep;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void InhomogenousCoupledPdeOdeSolverTemplated<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SolveAndWriteResultsToFile()
{
//std::cout<<"InhomogenousCoupledPdeOdeSolverTemplated<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SolveAndWriteResultsToFile( - start"<<std::endl;
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

//#ifdef CHASTE_VTK
    ResetInterpolatedQuantitiesOnNewTimeStep();
    // Create a .pvd output file
    //std::cout<<"outputfilehandeler: "<<this->mOutputDirectory<<std::endl;

    OutputFileHandler output_file_handler(this->mOutputDirectory, mClearOutputDirectory);
   
    mpVtkMetaFile = output_file_handler.OpenOutputFile("results.pvd");
    *mpVtkMetaFile << "<?xml version=\"1.0\"?>\n";
    *mpVtkMetaFile << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">\n";
    *mpVtkMetaFile << "    <Collection>\n";

    std::cout<<"Get output loc: "<<output_file_handler.GetChasteTestOutputDirectory()<<std::endl;
    std::cout<<"Get output dir full path: "<<output_file_handler.GetOutputDirectoryFullPath()<<std::endl;
    std::cout<<"Get relative path: "<<output_file_handler.GetRelativePath()<<std::endl;

    // Write initial condition to VTK
    Vec initial_condition = this->mInitialCondition;

    ReplicatableVector result_repl(this->mInitialCondition);

    WriteVtkResultsToFile(initial_condition, 0);

    // The helper class TimeStepper deals with issues such as small final timesteps so we don't have to
    TimeStepper stepper(this->mTstart, this->mTend, mSamplingTimeStep);


    std::vector<std::string> p_pde_stateVariableNames = mpPdeSystem -> GetStateVariableRegister() ->GetStateVariableRegisterVector();
    
   /* 
    // make a vector of file streams to store slice values per pde dim
    std::vector<std::shared_ptr<ofstream> > sliceMassFiles;
    for(unsigned i=0; i<PROBLEM_DIM; i++)
    {
    //    std::ofstream sliceMassFile;
     //   std::string sliceMassFilename = "/home/chaste/testoutput/"+this->mOutputDirectory+"sliceMass_"+p_pde_stateVariableNames[i]+".csv";
     //   std::cout<<"sliceMassFilename: "<<sliceMassFilename<<std::endl;
     //   sliceMassFile.open(sliceMassFilename);
     //   sliceMassFiles.push_back(sliceMassFile);
        std::shared_ptr<std::ofstream> sliceMassFile(new std::ofstream);
        std::string sliceMassFilename = "/home/chaste/testoutput/"+this->mOutputDirectory+"sliceMass_"+p_pde_stateVariableNames[i]+".csv";
        std::cout<<"sliceMassFilename: "<<sliceMassFilename<<std::endl;
        sliceMassFile -> open(sliceMassFilename);
        sliceMassFiles.push_back(sliceMassFile);
    }

    std::ofstream slicePositionFile;
    std::string slicePositionFilename = "/home/chaste/testoutput/"+this->mOutputDirectory+"slicePosition.csv";
    std::cout<<"slicePositionFilename: "<<slicePositionFilename<<std::endl;
    slicePositionFile.open (slicePositionFilename);

    std::ofstream positionFileX;
    std::string positionXFilename = "/home/chaste/testoutput/"+this->mOutputDirectory+"positionX.csv";
    std::cout<<"positionXFilename: "<<positionXFilename<<std::endl;
    positionFileX.open (positionXFilename);

    std::ofstream positionFileY;
    std::string positionYFilename = "/home/chaste/testoutput/"+this->mOutputDirectory+"positionY.csv";
    std::cout<<"positionYFilename: "<<positionYFilename<<std::endl;
    positionFileY.open (positionYFilename);

    std::vector<std::shared_ptr<ofstream> > positionMassFiles;
    for(unsigned i=0; i<PROBLEM_DIM; i++)
    {
        std::shared_ptr<std::ofstream> positionMassFile(new std::ofstream);
        std::string positionMassFilename = "/home/chaste/testoutput/"+this->mOutputDirectory+"positionMass_"+p_pde_stateVariableNames[i]+".csv";
        std::cout<<"positionMassFilename: "<<positionMassFilename<<std::endl;
        positionMassFile -> open(positionMassFilename);
        positionMassFiles.push_back(positionMassFile);
    }


    std::ofstream sumCoeffFile;
    std::string sumCoeffFilename = "/home/chaste/testoutput/"+this->mOutputDirectory+"sumCoeff.csv";
    std::cout<<"sumCoeffFilename: "<<sumCoeffFilename<<std::endl;
    sumCoeffFile.open (sumCoeffFilename);

    for(unsigned pde_index=0; pde_index<PROBLEM_DIM; pde_index++)
    {
        double sumCoeff=0;
        for(unsigned node_index=0;node_index<mpMesh->GetNumNodes();node_index++)
        {
            sumCoeff += result_repl[PROBLEM_DIM*node_index + pde_index];
        }
        sumCoeffFile << sumCoeff<<",";
    }
    sumCoeffFile <<"\n";


    std::ofstream totalMassFile;
    std::string totalMassFilename = "/home/chaste/testoutput/"+this->mOutputDirectory+"totalMass.csv";
    std::cout<<"totalMassFilename: "<<totalMassFilename<<std::endl;
    totalMassFile.open (totalMassFilename);

    std::ofstream totalMassOldFile;
    std::string totalMassOldFilename = "/home/chaste/testoutput/"+this->mOutputDirectory+"totalMassOld.csv";
    std::cout<<"totalMassFilename: "<<totalMassOldFilename<<std::endl;
    totalMassOldFile.open (totalMassOldFilename);

    std::ofstream timeFile;
    std::string timeFilename = "/home/chaste/testoutput/"+this->mOutputDirectory+"time.csv";
    std::cout<<"timeFilename: "<<timeFilename<<std::endl;
    timeFile.open (timeFilename);


    if(mCalculateTotalMass)
    {
        for(unsigned pdeNum=1; pdeNum<PROBLEM_DIM;pdeNum++)
        {
            mTotalMass.push_back(0.0);
            mOldTotalMass.push_back(0.0);
        }
    }
    
    //std::cout<<"here"<<std::endl;
*/
    // Main time loop, interate the stepper
    while (!stepper.IsTimeAtEnd())
    {   
        // Reset start and end times
        this->SetTimes(stepper.GetTime(), stepper.GetNextTime());

       // timeFile << stepper.GetTime()<<"\n";
        
        //std::cout<<"resetInt"<<std::endl;

     //   ResetInterpolatedQuantitiesOnNewTimeStep();

        
        //std::cout<<"get solve"<<std::endl;
        // Solve the system up to the new end time
        Vec soln = this->Solve();
        ReplicatableVector result_repl(soln);
        //std::cout<<"got solve"<<std::endl;
        // Reset the initial condition for the next timestep
        if (this->mInitialCondition != initial_condition)
        {
            PetscTools::Destroy(this->mInitialCondition);
        }
        this->mInitialCondition = soln;

        // Move forward in time
        stepper.AdvanceOneTimeStep();

    /*
        //std::cout<<"sumCoeff"<<std::endl;
        for(unsigned pde_index=0; pde_index<PROBLEM_DIM; pde_index++)
        {
            double sumCoeff=0;
            for(unsigned node_index=0;node_index<mpMesh->GetNumNodes();node_index++)
            {
                sumCoeff += result_repl[PROBLEM_DIM*node_index + pde_index];
            }
            sumCoeffFile << sumCoeff<<",";
        }
        sumCoeffFile << "\n";
       

        //std::cout<<"totalMass"<<std::endl;
        for(unsigned pde_index=0; pde_index<PROBLEM_DIM;pde_index++)
        {
           totalMassFile << mTotalMass[pde_index]<<",";
        }
        totalMassFile <<"\n";

        
        for(unsigned pde_index=0; pde_index<PROBLEM_DIM;pde_index++)
        {
           totalMassOldFile << mOldTotalMass[pde_index]<<",";
        }
        totalMassOldFile <<"\n";

        // for a given timestep output the slice data
        for(unsigned pos=0; pos<mSlicePositions.size(); pos++)
        {
            // take the x value of the position
            slicePositionFile << mSlicePositions[pos][0]<<",";
            for(unsigned pde_index=0; pde_index<PROBLEM_DIM; pde_index++)
            {
                *sliceMassFiles[pde_index] << mSliceMass[pos][pde_index]<<",";
            }
        }  
        // set the new line for the next time trace
        slicePositionFile << "\n";
        for(unsigned i=0; i<PROBLEM_DIM; i++)
        {
            *sliceMassFiles[i] << "\n";
        }



        //std::cout<<"domainPositions"<<std::endl;
        for(unsigned pos=0; pos<mDomainPositionsX.size(); pos++)
        {
            // take the x value of the position
            positionFileX << mDomainPositionsX[pos]<<",";
            positionFileY << mDomainPositionsY[pos]<<",";
            for(unsigned pde_index=0; pde_index<PROBLEM_DIM; pde_index++)
            {
                *positionMassFiles[pde_index] << mPositionMass[pde_index][pos]<<",";
            }
        }  
        // set the new line for the next time trace
        positionFileX << "\n";
        positionFileY << "\n";
        for(unsigned i=0; i<PROBLEM_DIM; i++)
        {
            *positionMassFiles[i] << "\n";
        }
*/


        // Write solution to VTK
        WriteVtkResultsToFile(soln, stepper.GetTotalTimeStepsTaken());
    }

/*
    // close all files
    timeFile.close();
    sumCoeffFile.close();
    totalMassFile.close();
    totalMassOldFile.close();
    slicePositionFile.close();
    positionFileX.close();
    positionFileY.close();
    for(unsigned i=0; i<PROBLEM_DIM; i++)
    {
        sliceMassFiles[i]->close();
        positionMassFiles[i]->close();
    }
*/

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
//#else //CHASTE_VTK
// LCOV_EXCL_START // We only test this in weekly builds
//    WARNING("VTK is not installed and is required for this functionality");
// LCOV_EXCL_STOP
//#endif //CHASTE_VTK
//std::cout<<"InhomogenousCoupledPdeOdeSolverTemplated<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SolveAndWriteResultsToFile( - end"<<std::endl;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void InhomogenousCoupledPdeOdeSolverTemplated<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::WriteVtkResultsToFile(Vec solution, unsigned numTimeStepsElapsed)
{

#ifdef CHASTE_VTK

    // Create a new VTK file for this time step
    std::stringstream time;
    time << numTimeStepsElapsed;
    //std::cout<<this->mOutputDirectory<<std::endl;
    VtkMeshWriter<ELEMENT_DIM, SPACE_DIM> mesh_writer(this->mOutputDirectory, "results_"+time.str(), false);
    // need to ensure StateVariableRegister is defined
  
    std::vector<std::string> p_pde_stateVariableNames = mpPdeSystem -> GetStateVariableRegister() ->GetStateVariableRegisterVector();
    
    /*
     * We first loop over PDEs. For each PDE we store the solution
     * at each node in a vector, then pass this vector to the mesh
     * writer.
     */
    ReplicatableVector solution_repl(solution);

    unsigned num_nodes = mpMesh->GetNumNodes();
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

    mesh_writer.WriteFilesUsingMesh(*mpMesh);
    *mpVtkMetaFile << "        <DataSet timestep=\"";
    *mpVtkMetaFile << numTimeStepsElapsed;
    *mpVtkMetaFile << "\" group=\"\" part=\"0\" file=\"results_";
    *mpVtkMetaFile << numTimeStepsElapsed;
    *mpVtkMetaFile << ".vtu\"/>\n";
#endif // CHASTE_VTK

}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
AbstractOdeSystemForCoupledPdeSystem* InhomogenousCoupledPdeOdeSolverTemplated<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::GetOdeSystemAtNode(unsigned index)
{
    return mOdeSystemsAtNodes[index];
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void InhomogenousCoupledPdeOdeSolverTemplated<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::ResetInterpolatedQuantitiesOnNewTimeStep()
{
/*
    //std::cout<<"reset: "<<mPositionMass.size()<<std::endl;
    if(mCalculateTotalMass)
    {
        // reset previous timstep value of total mass
        std::fill(mTotalMass.begin(), mTotalMass.end(), 0.0);
        
        std::fill(mOldTotalMass.begin(), mOldTotalMass.end(), 0.0);
    }

    if(mCalculateSliceMass)
    {
        mSliceMass = std::vector<std::vector<double>> ();
        mSlicePositions = std::vector<std::vector<double>> ();
    }

    mPositionMass.clear();
    
    for(unsigned i=0;i<PROBLEM_DIM;i++)
    {
        mPositionMass.push_back(std::vector<double>());
    }
    mDomainPositionsX.clear();

    mDomainPositionsY.clear();

    //std::cout<<"reset: "<<mPositionMass.size()<<std::endl;
*/
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void InhomogenousCoupledPdeOdeSolverTemplated<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::OutputInterpolatedQuantitiesOnTimestep()
{

}


#endif
