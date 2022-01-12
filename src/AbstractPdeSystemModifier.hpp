#ifndef ABSTRACTPDESYSTEMMODIFIER_HPP_
#define ABSTRACTPDESYSTEMMODIFIER_HPP_

//#include <cxxtest/GlobalFixture.h>
//#include "PetscSetupAndFinalize.hpp"
#include "ChasteSerialization.hpp"
#include "ClassIsAbstract.hpp"
#include <boost/shared_ptr.hpp>

#include "AbstractCellBasedSimulationModifier.hpp"
#include "TetrahedralMesh.hpp"
#include "AbstractIvpOdeSolver.hpp"
#include "AbstractOdeSystemForCoupledPdeSystem.hpp"
#include "AbstractLinearParabolicPdeSystemForCoupledOdeSystem.hpp"//mod

#include "AbstractBoundaryCondition.hpp"

#include "VtkMeshWriter.hpp" // shouldn't be needed? for #ifdef CHASTE_VTK bit

#include "ReplicatableVector.hpp"
#include "InhomogenousParabolicPdeForCoupledOdeSystem_templated.hpp"
#include "AbstractInhomogenousOdeSystemForCoupledPdeSystem.hpp"
#include "ChemicalDomainFieldForCellCoupling.hpp"

/**
 * An abstract modifier class containing functionality common to AbstractBoxDomainPdeModifier,
 * AbstractGrowingDomainPdeModifier and their subclasses, which solve a linear elliptic or
 * parabolic PDE coupled to a cell-based simulation.
 */
template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
class AbstractPdeSystemModifier : public AbstractCellBasedSimulationModifier<ELEMENT_DIM,SPACE_DIM>
{
protected:

    ChemicalDomainFieldForCellCoupling<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>* mpCoupledDomainField;

   
    boost::shared_ptr<InhomogenousParabolicPdeForCoupledOdeSystemTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>> mpPdeSystem; 


    /**
     * Shared pointer to a boundary condition object.
     */
    boost::shared_ptr<BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM> > mpBoundaryConditionsContainer; 

    std::vector<AbstractInhomogenousOdeSystemForCoupledPdeSystem*> pOdeSystemsAtNodes; // not needed?

    //boost::shared_ptr<AbstractIvpOdeSolver> pOdeSolver; // not needed?
    std::vector<boost::shared_ptr<AbstractIvpOdeSolver>> pOdeSolvers;

    /** The solution to the PDE problem at the current time step. */
    Vec mSolution; // serialised vector for pdeSystem

    /** Pointer to the finite element mesh on which to solve the PDE. */
    TetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* mpFeMesh;

    /** Store the output directory name. */
    std::string mOutputDirectory;

    /** Whether or not to calculate and output the gradient of the solution. */
    bool mOutputGradient;

    /**
     * Whether to output the PDE solution at each node of the FE mesh at output time steps.
     * Defaults to false.
     */
    bool mOutputSolutionAtPdeNodes;

    /** File that the values of the PDE solution are written out to. */
    out_stream mpVizPdeSolutionResultsFile;

    /**
     * Whether to delete the finite element mesh when we are destroyed.
     */
    bool mDeleteFeMesh;

    bool mIsDomainField = true;

public:


    AbstractPdeSystemModifier(ChemicalDomainFieldForCellCoupling<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>* p_domain_field,
                        Vec solution=nullptr);

    /**
     * Destructor.
     */
    virtual ~AbstractPdeSystemModifier();

    /**
     * @return mpPdeSystem
     */
    boost::shared_ptr<InhomogenousParabolicPdeForCoupledOdeSystemTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>> GetPdeSystem();

    boost::shared_ptr<BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM> > GetBoundaryConditionsContainer();


    /**
     * @return whether the PDE has an averaged source
     */
    bool HasAveragedSourcePde();

    /**
     * In the case where the PDE has an averaged source, set the source terms
     * using the information in the given mesh.
     *
     * @param pMesh Pointer to a tetrahedral mesh
     * @param pCellPdeElementMap map between cells and elements
     */
    void SetUpSourceTermsForAveragedSourcePde(TetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh, std::map<CellPtr, unsigned>* pCellPdeElementMap=nullptr);

    /**
     * @return mSolution.
     */
    Vec GetSolution();

    /**
     * @return mSolution (used in archiving)
     */
    Vec GetSolution() const;

    /**
     * @return mpFeMesh.
     */
    TetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* GetFeMesh() const;

    /**
     * Overridden SetupSolve() method.
     *
     * Set mOutputDirectory and, if mOutputSolutionAtPdeNodes is set to true, open mpVizPdeSolutionResultsFile.
     * This method is overridden in subclasses.
     *
     * @param rCellPopulation reference to the cell population
     * @param outputDirectory the output directory, relative to where Chaste output is stored
     */
    virtual void SetupSolve(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation, std::string outputDirectory);

    /**
     * Overridden UpdateAtEndOfTimeStep() method.
     *
     * As this method is pure virtual, it must be overridden in subclasses.
     *
     * @param rCellPopulation reference to the cell population
     */
    virtual void UpdateAtEndOfTimeStep(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation)=0;

    /**
     * Overridden UpdateAtEndOfOutputTimeStep() method,
     * after UpdateAtEndOfTimeStep() has been called.
     *
     * Output the solution to the PDE at each cell to VTK and, if mOutputSolutionAtPdeNodes is set to true,
     * output the solution to the PDE at each node of mpFeMesh to mpVizPdeSolutionResultsFile.
     *
     * @param rCellPopulation reference to the cell population
     */
    virtual void UpdateAtEndOfOutputTimeStep(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation);

    /**
     * Overridden UpdateAtEndOfSolve() method.
     *
     * If mOutputSolutionAtPdeNodes is set to true, close mpVizPdeSolutionResultsFile.
     *
     * @param rCellPopulation reference to the cell population
     */
    virtual void UpdateAtEndOfSolve(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation);

    /**
     * Set whether to calculate and save the gradient of the solution to CellData.
     *
     * @return mOutputGradient
     */
    bool GetOutputGradient();

    /**
     * Set whether to calculate and save the gradient of the solution to CellData.
     *
     * @param outputGradient whether to output the gradient
     */
    void SetOutputGradient(bool outputGradient);

    /**
     * Set mOutputSolutionAtPdeNodes.
     *
     * @param outputSolutionAtPdeNodes whether to output the PDE solution at each node of the FE mesh at output time steps
     */
    void SetOutputSolutionAtPdeNodes(bool outputSolutionAtPdeNodes);

    /**
     * Overridden OutputSimulationModifierParameters() method.
     *
     * Output any simulation modifier parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputSimulationModifierParameters(out_stream& rParamsFile);

    void SetDomainField(AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>*);

    AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>* GetDomainField();
};

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
AbstractPdeSystemModifier<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::AbstractPdeSystemModifier(ChemicalDomainFieldForCellCoupling<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>* p_domain_field,
                                              Vec solution)
      : AbstractCellBasedSimulationModifier<ELEMENT_DIM,SPACE_DIM>(),
        mpCoupledDomainField(p_domain_field),
        mpPdeSystem(p_domain_field->ReturnSharedPtrPdeSystem()), // passed to AbstractCellBasedSimulation, as is a simulation modifier
        mpBoundaryConditionsContainer(p_domain_field->ReturnSharedPtrBoundaryConditionsContainer()), // passed to AbstractCellBasedSimulation
        mSolution(nullptr),
        mOutputDirectory(""),
        mOutputGradient(false),
        mOutputSolutionAtPdeNodes(false),
        mDeleteFeMesh(false)
{
    if (solution)
    {
        mSolution = solution;
    }
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
AbstractPdeSystemModifier<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::~AbstractPdeSystemModifier()
{
    if (mDeleteFeMesh and mpFeMesh!=nullptr)
    {
        delete mpFeMesh;
    }
    if (mSolution)
    {
        PetscTools::Destroy(mSolution);
    }
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
boost::shared_ptr<InhomogenousParabolicPdeForCoupledOdeSystemTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>> AbstractPdeSystemModifier<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::GetPdeSystem()
{
    return mpPdeSystem;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
boost::shared_ptr<BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM> > AbstractPdeSystemModifier<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::GetBoundaryConditionsContainer()
{
    return mpBoundaryConditionsContainer;
}


template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>  // remove? undecided
bool AbstractPdeSystemModifier<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::HasAveragedSourcePde()
{
    // if the pde function iherits from the averaged source pdes
    return false;//((boost::dynamic_pointer_cast<AveragedSourceInhomogenousParabolicPdeOdeSystem<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM> >(mpPdeSystem) != nullptr));    
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM> //remove? or change to linearParabolicPdeSystem?
void AbstractPdeSystemModifier<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetUpSourceTermsForAveragedSourcePde(TetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh, std::map<CellPtr, unsigned>* pCellPdeElementMap)
{
    
    assert(HasAveragedSourcePde());

    /*
    if (boost::dynamic_pointer_cast<AveragedSourceInhomogenousParabolicPdeOdeSystem<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM> >(mpPdeSystem) != nullptr)
    {
        boost::static_pointer_cast<AveragedSourceInhomogenousParabolicPdeOdeSystem<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM> >(mpPdeSystem)->SetupSourceTerms(*pMesh, pCellPdeElementMap);
    }
    
    if (boost::dynamic_pointer_cast<AveragedSourceEllipticPde<SPACE_DIM> >(mpPdeSystem) != nullptr)
    {
        boost::static_pointer_cast<AveragedSourceEllipticPde<SPACE_DIM> >(mpPdeSystem)->SetupSourceTerms(*pMesh, pCellPdeElementMap);
    }
    else if (boost::dynamic_pointer_cast<AveragedSourceParabolicPde<SPACE_DIM> >(mpPdeSystem) != nullptr)
    {
        boost::static_pointer_cast<AveragedSourceParabolicPde<SPACE_DIM> >(mpPdeSystem)->SetupSourceTerms(*pMesh, pCellPdeElementMap);
    }
    */
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
Vec AbstractPdeSystemModifier<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::GetSolution()
{
    return mSolution;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
Vec AbstractPdeSystemModifier<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::GetSolution() const
{
    return mSolution;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
TetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* AbstractPdeSystemModifier<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::GetFeMesh() const
{
    return mpFeMesh;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
void AbstractPdeSystemModifier<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetupSolve(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation, std::string outputDirectory)
{
    //std::cout<<"AbstractPdeSystemModifier<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetupSolve - start"<<std::endl;
    // Cache the output directory
    this->mOutputDirectory = outputDirectory; 

    if (mOutputSolutionAtPdeNodes)
    {
       
        if (PetscTools::AmMaster())
        {
       
            OutputFileHandler output_file_handler(outputDirectory+"/", false);
         
            mpVizPdeSolutionResultsFile = output_file_handler.OpenOutputFile("results.vizpdesolution");
   
        }
    }
    //std::cout<<"AbstractPdeSystemModifier<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetupSolve - end"<<std::endl;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
void AbstractPdeSystemModifier<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::UpdateAtEndOfOutputTimeStep(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation)
{
    //std::cout<<"AbstractPdeSystemModifier - UpdateAtEndOfOutputTimeStep - start"<<std::endl;
    if (mOutputSolutionAtPdeNodes)
    {
        if (PetscTools::AmMaster())
        {
            (*mpVizPdeSolutionResultsFile) << SimulationTime::Instance()->GetTime() << "\t";

            assert(mpFeMesh != nullptr);

            for (unsigned i=0; i<mpFeMesh->GetNumNodes(); i++)
            {
                (*mpVizPdeSolutionResultsFile) << i << " ";
                const c_vector<double,SPACE_DIM>& r_location = mpFeMesh->GetNode(i)->rGetLocation();

                for (unsigned k=0; k<SPACE_DIM; k++)
                {
                    (*mpVizPdeSolutionResultsFile) << r_location[k] << " ";
                }

                assert(mSolution != nullptr); // remove, test against null vector

                ReplicatableVector solution_repl(mSolution);
                for(unsigned pd=0; pd<PROBLEM_DIM; pd++){
                    (*mpVizPdeSolutionResultsFile) << solution_repl[i+pd] << " ";
                }
            }

            (*mpVizPdeSolutionResultsFile) << "\n";
        }
    
    }

#ifdef CHASTE_VTK

    if (SPACE_DIM > 1)
    {
        std::ostringstream time_string;
        time_string << SimulationTime::Instance()->GetTimeStepsElapsed();

        ReplicatableVector solution_repl(mSolution);

        for(unsigned pdeDim=0; pdeDim<PROBLEM_DIM; pdeDim++)
        {
            std::string results_file = "pde_results_" + mpCoupledDomainField->GetDomainStateVariableRegister()->RetrieveStateVariableName(pdeDim) + "_" + time_string.str();
            VtkMeshWriter<ELEMENT_DIM,SPACE_DIM>* p_vtk_mesh_writer = new VtkMeshWriter<ELEMENT_DIM,SPACE_DIM>(mOutputDirectory, results_file, false);
            std::vector<double> pde_solution;
            //std::cout<<results_file<<" : "<<std::endl;
            for (unsigned i=0; i<mpFeMesh->GetNumNodes(); i++)
            {
                pde_solution.push_back(solution_repl[i*PROBLEM_DIM+pdeDim]); 
                //std::cout<<solution_repl[i*PROBLEM_DIM+pdeDim]<<std::endl;
            }
            p_vtk_mesh_writer->AddPointData(mpCoupledDomainField->GetDomainStateVariableRegister()->RetrieveStateVariableName(pdeDim), pde_solution);

            p_vtk_mesh_writer->WriteFilesUsingMesh(*mpFeMesh);

            delete p_vtk_mesh_writer;
        }
        
    }
#endif //CHASTE_VTK
    //std::cout<<"AbstractPdeSystemModifier - UpdateAtEndOfOutputTimeStep - end"<<std::endl;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
void AbstractPdeSystemModifier<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::UpdateAtEndOfSolve(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation)
{
    if (mOutputSolutionAtPdeNodes)
    {
        if (PetscTools::AmMaster())
        {
            mpVizPdeSolutionResultsFile->close();
        }
    }
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
bool AbstractPdeSystemModifier<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::GetOutputGradient()
{
    return mOutputGradient;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
void AbstractPdeSystemModifier<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetOutputGradient(bool outputGradient)
{
    mOutputGradient = outputGradient;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
void AbstractPdeSystemModifier<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetOutputSolutionAtPdeNodes(bool outputSolutionAtPdeNodes)
{
    mOutputSolutionAtPdeNodes = outputSolutionAtPdeNodes;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
void AbstractPdeSystemModifier<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<ELEMENT_DIM,SPACE_DIM>::OutputSimulationModifierParameters(rParamsFile);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractPdeSystemModifier<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetDomainField(AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>* pDomainField)
{
    mpCoupledDomainField = pDomainField;
    mIsDomainField = true;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>* AbstractPdeSystemModifier<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::GetDomainField()
{
    return mpCoupledDomainField;
}

#endif /*ABSTRACTPDEMODIFIER_HPP_*/
