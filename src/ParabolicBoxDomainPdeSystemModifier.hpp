#ifndef PARABOLICBOXDOMAINPDESYSTEMMODIFIER_HPP_
#define PARABOLICBOXDOMAINPDESYSTEMMODIFIER_HPP_

#include "InhomogenousParabolicPdeForCoupledOdeSystem_templated.hpp"
#include "AbstractBoxDomainPdeSystemModifier.hpp"
#include "AbstractPdeSystemModifier.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "AbstractIvpOdeSolver.hpp"
#include "AbstractOdeSystemForCoupledPdeSystem.hpp"
#include "InhomogenousCoupledPdeOdeCoupledCellSolver.hpp"
#include "AbstractOdeSystemForCoupledPdeSystem.hpp"
#include "InhomogenousCoupledPdeOdeSolver_templated.hpp"
#include "ChemicalDomainFieldForCellCoupling.hpp"


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
class ParabolicBoxDomainPdeSystemModifier : public AbstractBoxDomainPdeSystemModifier<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>
{
protected:
    bool mConditionsInterpolated = false;

public:

    ParabolicBoxDomainPdeSystemModifier(ChemicalDomainFieldForCellCoupling<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>* p_domain_field,
                                  boost::shared_ptr<ChasteCuboid<SPACE_DIM> > pMeshCuboid=boost::shared_ptr<ChasteCuboid<SPACE_DIM> >(),
                                  double stepSize=1.0,
                                  Vec solution=nullptr);

    /**
     * Destructor.
     */
    virtual ~ParabolicBoxDomainPdeSystemModifier();

    /**
     * Overridden UpdateAtEndOfTimeStep() method.
     *
     * Specifies what to do in the simulation at the end of each time step.
     *
     * @param rCellPopulation reference to the cell population
     */
    virtual void UpdateAtEndOfTimeStep(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation);

    /**
     * Overridden SetupSolve() method.
     *
     * Specifies what to do in the simulation before the start of the time loop.
     *
     * @param rCellPopulation reference to the cell population
     * @param outputDirectory the output directory, relative to where Chaste output is stored
     */
    virtual void SetupSolve(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation, std::string outputDirectory);

    /**
     * Helper method to construct the boundary conditions container for the PDE.
     *
     * @param rCellPopulation reference to the cell population
     *
     * @return the full boundary conditions container
     */
    virtual boost::shared_ptr<BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM> > ConstructBoundaryConditionsContainer(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation);

    /**
     * Helper method to initialise the PDE solution using the CellData.
     *
     * Here we assume a homogeneous initial consition.
     *
     * @param rCellPopulation reference to the cell population
     */
    void SetupInitialSolutionVector(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation);


    /**
     * Overridden OutputSimulationModifierParameters() method.
     * Output any simulation modifier parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputSimulationModifierParameters(out_stream& rParamsFile);

    void SetPdeDimension(unsigned pdeDim);

    unsigned GetPdeDimension();

    //void SetNodalInitialConditions(std::vector<double> init_nodal_conditions);

    //std::vector<double> GetNodalInitialConditions();

    //void SetCellInitialConditions(std::vector<double> init_cell_conditions);

    //std::vector<double> GetCellInitialConditions();

    
};

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
ParabolicBoxDomainPdeSystemModifier<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::ParabolicBoxDomainPdeSystemModifier(
                                    ChemicalDomainFieldForCellCoupling<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>* p_domain_field,
                                                                  boost::shared_ptr<ChasteCuboid<SPACE_DIM> > pMeshCuboid,
                                                                  double stepSize,
                                                                  Vec solution)
    : AbstractBoxDomainPdeSystemModifier<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>(p_domain_field,
                                        pMeshCuboid,
                                        stepSize,
                                        solution)
{

}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
ParabolicBoxDomainPdeSystemModifier<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::~ParabolicBoxDomainPdeSystemModifier()
{

}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void ParabolicBoxDomainPdeSystemModifier<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation)
{
    //std::cout<<"ParabolicBoxDomainPdeSystemModifier<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::UpdateAtEndOfTimeStep - start"<<std::endl;
    // Set up boundary conditions, comes from the rCellPopulation rather than the constructor
    boost::shared_ptr<BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM> > p_bcc = ConstructBoundaryConditionsContainer(rCellPopulation);


    this->UpdateCellPdeElementMap(rCellPopulation);

    // When using a PDE mesh which doesn't coincide with the cells, we must set up the source terms before solving the PDE.
    // Pass in already updated CellPdeElementMap to speed up finding cells.
    // this line shoudl be fine for the pOdeSystem

    this->SetUpSourceTermsForAveragedSourcePde(this->mpFeMesh, &this->mCellPdeElementMap);

    InhomogenousCoupledPdeOdeCoupledCellSolver<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM> solver(
                                                this->mpFeMesh, 
                                                this->mpCoupledDomainField->ReturnSharedPtrPdeSystem().get(), 
                                                this->mpCoupledDomainField->ReturnSharedPtrBoundaryConditionsContainer().get(),
                                                rCellPopulation,
                                                this->mpCoupledDomainField->GetNodalOdeSystems(),
                                                this->mpCoupledDomainField->GetNodalOdeSolvers(),
                                                mConditionsInterpolated
                                                );

    ///\todo Investigate more than one PDE time step per spatial step

    SimulationTime* p_simulation_time = SimulationTime::Instance();
    double current_time = p_simulation_time->GetTime();
    double dt = p_simulation_time->GetTimeStep();
    // solver is calling the LinearParabolicSystemWithCOupledOdeSystemSolver
   
    solver.SetTimes(current_time,current_time + dt);
    solver.SetTimeStep(dt);
 
    // Use previous solution as the initial condition
    Vec previous_solution = this->mSolution;
    solver.SetInitialCondition(previous_solution);

    // Note that the linear solver creates a vector, so we have to keep a handle on the old one
    // in order to destroy it

    this->mSolution = solver.Solve();

    PetscTools::Destroy(previous_solution);

    this->UpdateCellData(rCellPopulation);
    mConditionsInterpolated = true;
    //std::cout<<"ParabolicBoxDomainPdeSystemModifier<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::UpdateAtEndOfTimeStep - end"<<std::endl;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void ParabolicBoxDomainPdeSystemModifier<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetupSolve(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation, std::string outputDirectory)
{
    //std::cout<<"ParabolicBoxDomainPdeSystemModifier<ELEMENT_DIM,SPACE_DIM>::SetupSolve - start"<<std::endl;
    AbstractBoxDomainPdeSystemModifier<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetupSolve(rCellPopulation,outputDirectory);

    // Copy the cell data to mSolution (this is the initial condition)
    SetupInitialSolutionVector(rCellPopulation);

    // Output the initial conditions on FeMesh
    this->UpdateAtEndOfOutputTimeStep(rCellPopulation);
    //std::cout<<"ParabolicBoxDomainPdeSystemModifier<ELEMENT_DIM,SPACE_DIM>::SetupSolve - end"<<std::endl;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM> // why does this take in a cell population? if() is true, shrinks the box onto the tissue (cellpopulation) here makes no difference
boost::shared_ptr<BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM> > ParabolicBoxDomainPdeSystemModifier<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::ConstructBoundaryConditionsContainer(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation)
{
    boost::shared_ptr<BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM> > p_bcc(new BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>(false)); // false implies not to delete previous conditions but there shouldn't be any

    if (!this->mSetBcsOnBoxBoundary)
    {
        EXCEPTION("Boundary conditions cannot yet be set on the cell population boundary for a ParabolicBoxDomainPdeSystemModifier");
    }
    else // Apply BC at boundary nodes of box domain FE mesh
    {
        p_bcc = this->mpCoupledDomainField -> ReturnSharedPtrBoundaryConditionsContainer();//ChemicalDomainFieldForCellCoupling<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::ReturnSharedPtrBoundaryConditionsContainer();
    }
 
    return p_bcc;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void ParabolicBoxDomainPdeSystemModifier<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetupInitialSolutionVector(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation)
{ 
    
    // set up the initial conditions for the pde mesh

    // set up cell initial conditions
    for (typename AbstractCellPopulation<SPACE_DIM>::Iterator cell_iter = rCellPopulation.Begin();
    cell_iter != rCellPopulation.End();
    ++cell_iter)
    {

        CellPropertyCollection& prop_collection = cell_iter->rGetCellPropertyCollection();

        if (prop_collection.HasProperty<ChemicalCellProperty>())
        {   
            // the cell has it's own concentration vector that may be related to the domain
            boost::shared_ptr<ChemicalCellProperty> property = boost::static_pointer_cast<ChemicalCellProperty>(prop_collection.GetPropertiesType<ChemicalCellProperty>().GetProperty());

            std::vector<std::string> cell_species_names = property -> GetStateVariableRegister() -> GetStateVariableRegisterVector();
  
            for(unsigned name_index=0; name_index<cell_species_names.size();name_index++)
            {   
                cell_iter->GetCellData()->SetItem(cell_species_names[name_index], property -> GetCellConcentrationByIndex(name_index));
            }
        }
        else
        {   
            // assume zero concentration in cell, set up for all species in the PROBLEM_DIM, that is species diffusing through the domain
            std::vector<double> initial_conditions(PROBLEM_DIM,0.0);
            std::vector<std::string> domain_species_names = this->mpCoupledDomainField -> GetDomainStateVariableRegister() -> GetStateVariableRegisterVector();

            for(unsigned name_index=0; name_index<domain_species_names.size();name_index++)
            {
                cell_iter->GetCellData()->SetItem(domain_species_names[name_index], 0.0);
            }

        }
        
    }

    // set pde serialised nodal initial conditions from domain layer

    // Initialise mSolution
    this->mSolution = PetscTools::CreateVec(this->mpCoupledDomainField -> GetInitialNodeConditions());

}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void ParabolicBoxDomainPdeSystemModifier<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractBoxDomainPdeSystemModifier<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::OutputSimulationModifierParameters(rParamsFile);
}


#endif /*PARABOLICBOXDOMAINPDESYSTEMMODIFIER_HPP_*/
