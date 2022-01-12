#ifndef ABSTRACTBOXDOMAINPDESYSTEMMODIFIER_HPP_
#define ABSTRACTBOXDOMAINPDESYSTEMMODIFIER_HPP_

#include "AbstractPdeSystemModifier.hpp"
#include "ReplicatableVector.hpp"
#include "LinearBasisFunction.hpp"


#include "AbstractCellProperty.hpp"

#include "SmartPointers.hpp"
//This test is always run sequentially (never in parallel)
//#include "FakePetscSetup.hpp"


#include "AbstractDomainField.hpp"
#include "CellPropertyCollection.hpp"
#include "MembraneCellProperty.hpp"
#include "TransportCellProperty.hpp"
#include "ExtendedCellProperty.hpp"

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
class AbstractBoxDomainPdeSystemModifier : public AbstractPdeSystemModifier<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>
{

protected:

    /** Map between cells and the elements of the FE mesh containing them. */
    std::map<CellPtr, unsigned> mCellPdeElementMap;

    /**
     * Pointer to a ChasteCuboid storing the outer boundary for the FE mesh.
     * Stored as a member to facilitate archiving.
     */
    boost::shared_ptr<ChasteCuboid<SPACE_DIM> > mpMeshCuboid;

    /**
     * The step size to be used in the FE mesh.
     * Stored as a member to facilitate archiving.
     */
    double mStepSize;

    bool mIsCenterMesh;

    c_vector<double,SPACE_DIM> mOffset;

    /**
     * Whether to set the boundary condition on the edge of the box domain rather than the cell population.
     * Default to true.
     */
    bool mSetBcsOnBoxBoundary;

    
public:


    AbstractBoxDomainPdeSystemModifier(ChemicalDomainFieldForCellCoupling<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>* p_domain_field,
                                 boost::shared_ptr<ChasteCuboid<SPACE_DIM> > pMeshCuboid=boost::shared_ptr<ChasteCuboid<SPACE_DIM> >(),
                                 double stepSize=1.0,
                                 bool isCenterMesh = true,
                                 Vec solution=nullptr);

    virtual ~AbstractBoxDomainPdeSystemModifier();

    /**
     * @return mStepSize.
     */
    double GetStepSize();


    /**
     * Set mSetBcsOnCoarseBoundary.
     *
     * @param setBcsOnBoxBoundary whether to set the boundary condition on the edge of the box domain rather than the cell population
     */
    void SetBcsOnBoxBoundary(bool setBcsOnBoxBoundary);

    /**
     * @return mSetBcsOnCoarseBoundary.
     */
    bool AreBcsSetOnBoxBoundary();

    /**
     * Overridden SetupSolve() method.
     *
     * Specifies what to do in the simulation before the start of the time loop.
     *
     * Here we just initialize the Cell PDE element map
     *
     * @param rCellPopulation reference to the cell population
     * @param outputDirectory the output directory, relative to where Chaste output is stored
     */
    virtual void SetupSolve(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation, std::string outputDirectory);

    /**
     * Helper method to generate the mesh.
     *
     * @param pMeshCuboid the outer boundary for the FE mesh.
     * @param stepSize the step size to be used in the FE mesh.
     */
    void GenerateFeMesh(boost::shared_ptr<ChasteCuboid<SPACE_DIM> > pMeshCuboid, double stepSize, bool isCenterMesh);

    /**
     * Helper method to copy the PDE solution to CellData
     *
     * Here we need to interpolate from the FE mesh onto the cells.
     *
     * @param rCellPopulation reference to the cell population
     */
    void UpdateCellData(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation);

    /**
     * Initialise mCellPdeElementMap.
     *
     * @param rCellPopulation reference to the cell population
     */
    void InitialiseCellPdeElementMap(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation);

    /**
     * Update the mCellPdeElementMap
     *
     * This method should be called before sending the element map to a PDE class
     * to ensure map is up to date.
     *
     * @param rCellPopulation reference to the cell population
     */
    void UpdateCellPdeElementMap(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation);

    /**
     * Overridden OutputSimulationModifierParameters() method.
     * Output any simulation modifier parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputSimulationModifierParameters(out_stream& rParamsFile);

    c_vector<double,SPACE_DIM> GetMeshOffset();

};


template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
AbstractBoxDomainPdeSystemModifier<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::AbstractBoxDomainPdeSystemModifier(ChemicalDomainFieldForCellCoupling<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>* p_domain_field,
                                                                boost::shared_ptr<ChasteCuboid<SPACE_DIM> > pMeshCuboid,
                                                                double stepSize,
                                                                bool isCenterMesh,
                                                                Vec solution) 
    : AbstractPdeSystemModifier<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>(p_domain_field,
                               solution), 
      mpMeshCuboid(pMeshCuboid),
      mStepSize(stepSize),
      mIsCenterMesh(isCenterMesh),
      mSetBcsOnBoxBoundary(true)
      
{
    if (pMeshCuboid)
    {
        // We only need to generate mpFeMesh once, as it does not vary with time
        this->GenerateFeMesh(mpMeshCuboid, mStepSize, mIsCenterMesh);
        this->mDeleteFeMesh = true;
    }
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
AbstractBoxDomainPdeSystemModifier<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::~AbstractBoxDomainPdeSystemModifier()
{
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
double AbstractBoxDomainPdeSystemModifier<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::GetStepSize()
{
     return mStepSize;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
void AbstractBoxDomainPdeSystemModifier<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetBcsOnBoxBoundary(bool setBcsOnBoxBoundary)
{
    mSetBcsOnBoxBoundary = setBcsOnBoxBoundary;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
bool AbstractBoxDomainPdeSystemModifier<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::AreBcsSetOnBoxBoundary()
{
    return mSetBcsOnBoxBoundary;
}


template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
void AbstractBoxDomainPdeSystemModifier<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetupSolve(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation, std::string outputDirectory)
{
//    std::cout<<"AbstractBoxDomainPdeSystemModifier<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetupSolve - start"<<std::endl;
    AbstractPdeSystemModifier<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetupSolve(rCellPopulation, outputDirectory);

    InitialiseCellPdeElementMap(rCellPopulation);
 //   std::cout<<"AbstractBoxDomainPdeSystemModifier<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetupSolve - end"<<std::endl;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
void AbstractBoxDomainPdeSystemModifier<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::GenerateFeMesh(boost::shared_ptr<ChasteCuboid<SPACE_DIM> > pMeshCuboid, double stepSize, bool isCenterMesh)
{
    //std::cout<<"AbstractBoxDomainPdeSystemModifier<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::GenerateFeMesh - start"<<std::endl;
    // generate mesh and PDE system from the Abstract domain and the mesh cuboid
//cell_population.GetCentroidOfCellPopulation()
    // Create a regular coarse tetrahedral mesh

    // Get centroid of meshCuboid
    ChastePoint<SPACE_DIM> upper = pMeshCuboid->rGetUpperCorner();
    ChastePoint<SPACE_DIM> lower = pMeshCuboid->rGetLowerCorner();
    c_vector<double,SPACE_DIM> centre_of_cuboid = 0.5*(upper.rGetLocation() + lower.rGetLocation());

    //std::cout<<"cuboid upper: x= "<<upper[0]<<" y="<<upper[1]<<std::endl;
    //std::cout<<"cuboid lower: x= "<<lower[0]<<" y="<<lower[1]<<std::endl;

    this->mpCoupledDomainField -> GenerateFeMesh();
    this->mpCoupledDomainField -> SetUpDomainFromFiles();

 //   std::cout<<"AbstractBoxDomain::generateFeMesh ScaleFeMesh"<<std::endl;
    this->mpCoupledDomainField -> ScaleFeMesh();
    
    this->mpFeMesh = this->mpCoupledDomainField -> rGetDomainFeMesh();


    // set the lower point of the cuboid as the origin of the mesh
    std::vector<double> origin(SPACE_DIM,0.0);
    c_vector<double,SPACE_DIM> offset = zero_vector<double>(SPACE_DIM);

    if(isCenterMesh)
    {
        //std::cout<<"AbstractBoxDomain::generateFeMesh IsCenterMesh"<<std::endl;
        // calculate the center of both meshes
        // Find the centre of the PDE mesh
        c_vector<double,SPACE_DIM> centre_of_coarse_mesh = zero_vector<double>(SPACE_DIM);
        for (unsigned i=0; i<this->mpFeMesh->GetNumNodes(); i++)
        {
            centre_of_coarse_mesh += this->mpFeMesh->GetNode(i)->rGetLocation();
        }
        centre_of_coarse_mesh /= this->mpFeMesh->GetNumNodes();

        // Find the centre of the cell mesh
        c_vector<double,SPACE_DIM> centre_of_cell_mesh = zero_vector<double>(SPACE_DIM);
        for (unsigned i=0; i<this->mpCoupledDomainField ->rGetCellMesh()->GetNumNodes(); i++)
        {
            centre_of_cell_mesh += this->mpCoupledDomainField ->rGetCellMesh()->GetNode(i)->rGetLocation();
        }
        centre_of_cell_mesh /= this->mpFeMesh->GetNumNodes();

        //std::cout<<"centre_of_cell_mesh: x="<<centre_of_cell_mesh(0)<<" y="<<centre_of_cell_mesh(1)<<std::endl;
        //std::cout<<"centre_of_coarse_mesh: x="<<centre_of_coarse_mesh(0)<<" y="<<centre_of_coarse_mesh(1)<<std::endl;
        //std::cout<<"difference: x="<<-centre_of_cell_mesh(0) + centre_of_coarse_mesh(0)<<" y="<<centre_of_cell_mesh(1) - centre_of_coarse_mesh(1)<<std::endl;
        for(unsigned dim=0; dim<SPACE_DIM; dim++)
        {
            
            offset[dim] = centre_of_coarse_mesh(dim) - centre_of_cell_mesh(dim);
            origin[dim] = centre_of_coarse_mesh(dim) - centre_of_cell_mesh(dim);
            //std::cout<<"dim: "<<dim<<" offset: "<<offset[dim]<<std::endl;
        }

    }
    else
    {
        for(unsigned dim=0; dim<SPACE_DIM; dim++)
        {
            //std::cout<<"else"<<std::endl;
            offset[dim] = lower[dim];
            origin[dim] = lower[dim];
            //std::cout<<"dim: "<<dim<<" origin: "<<origin[dim]<<std::endl;
        }
    }

    mOffset = offset;    

    // Now move the mesh to the correct location
    //this->mpFeMesh->Translate(centre_of_cuboid - offset);
    //this->mpFeMesh->Translate(offset);
    //this->mpCoupledDomainField ->rGetCellMesh()->Translate(800,800);
    this->mpCoupledDomainField -> SetLabelOrigin(origin);


    // c_vector<double,SPACE_DIM> centre_of_coarse_mesh = zero_vector<double>(SPACE_DIM);
    // for (unsigned i=0; i<this->mpFeMesh->GetNumNodes(); i++)
    // {
    //     centre_of_coarse_mesh += this->mpFeMesh->GetNode(i)->rGetLocation();
    // }
    // centre_of_coarse_mesh /= this->mpFeMesh->GetNumNodes();

    // // Find the centre of the cell mesh
    // c_vector<double,SPACE_DIM> centre_of_cell_mesh = zero_vector<double>(SPACE_DIM);
    // for (unsigned i=0; i<this->mpCoupledDomainField ->rGetCellMesh()->GetNumNodes(); i++)
    // {
    //     centre_of_cell_mesh += this->mpCoupledDomainField ->rGetCellMesh()->GetNode(i)->rGetLocation();
    // }
    // centre_of_cell_mesh /= this->mpFeMesh->GetNumNodes();

    // std::cout<<"centre_of_cell_mesh: x="<<centre_of_cell_mesh(0)<<" y="<<centre_of_cell_mesh(1)<<std::endl;
    // std::cout<<"centre_of_coarse_mesh: x="<<centre_of_coarse_mesh(0)<<" y="<<centre_of_coarse_mesh(1)<<std::endl;
    // std::cout<<"difference: x="<<-centre_of_cell_mesh(0) + centre_of_coarse_mesh(0)<<" y="<<centre_of_cell_mesh(1) - centre_of_coarse_mesh(1)<<std::endl;

    //std::cout<<"AbstractBoxDomainPdeSystemModifier<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::GenerateFeMesh - end"<<std::endl;
}


// need to check cell data
template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
void AbstractBoxDomainPdeSystemModifier<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::UpdateCellData(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation)
{
 //   std::cout<<"AbstractBoxDomainPdeSystemModifier - UpdateCellData - start"<<std::endl;
    // Store the PDE solution in an accessible form
    ReplicatableVector solution_repl(this->mSolution); // nodal solution from pdeSolver

    for (typename AbstractCellPopulation<SPACE_DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        // in both cases, with or without the extended cell property, need to refresh the cell concentration vector avaliable for
        // transport from cell data to the appropriate concentration containers

        CellPropertyCollection& prop_collection = cell_iter->rGetCellPropertyCollection();
        
        unsigned elem_index = mCellPdeElementMap[*cell_iter];
        
        if (prop_collection.HasProperty<ExtendedCellProperty<SPACE_DIM>>())//cell_iter->HasCellProperty<ExtendedCellProperty<SPACE_DIM>>())
        {
            /*

            // the cells have a spatial extent where the state varibales are sorted with location keys in the extended cell property
            // therefore no need to interpolate form the nodes but need to recover and set the CellData
            boost::shared_ptr<ExtendedCellProperty<SPACE_DIM>> extended_cell_property = boost::static_pointer_cast<ExtendedCellProperty<SPACE_DIM>>(prop_collection.GetPropertiesType<ExtendedCellProperty<SPACE_DIM>>().GetProperty());
                    
            // for each state variable in the cell data, retrieve the stored cell data value and add the total internal cell concetration
            
            std::vector<double> this_cell_next_total_boundary_concentration = extended_cell_property -> GetNextTimestepConcentrationVector();
            unsigned number_cell_states =this_cell_next_total_boundary_concentration.size();
            std::vector<double> this_cell_total_concentration(number_cell_states,0.0);
            for(unsigned i=0; i<number_cell_states;i++)
            {
                // add the change in cell concentration due to the effects of the internal end of the cell boundary 
                double next_state_value = cell_iter->GetCellData()->GetItem(extended_cell_property -> GetStateVariableRegister()->RetrieveStateVariableName(i)) + this_cell_next_total_concentration_change[i];
                cell_iter->GetCellData()->SetItem(extended_cell_property -> GetStateVariableRegister()->RetrieveStateVariableName(i), next_state_value);
                this_cell_total_concentration[i] = next_state_value;
            }

            // next determine how much of the cell data may be utilised by the next solver call
            // this is a property of the cell so store in extended cell property
            extended_cell_property -> UpdateCellConcentrationVector(this_cell_total_concentration);

            // remove this reserved concentration vector from the cell data
            for(unsigned i=0; i<number_cell_states;i++)
            {
                double cell_value = cell_iter->GetCellData()->GetItem(extended_cell_property -> GetStateVariableRegister()->RetrieveStateVariableName(i))
                cell_iter->GetCellData()->SetItem(extended_cell_property -> GetStateVariableRegister()->RetrieveStateVariableName(i), cell_value - this_cell_total_concentration[i]);
            }
            */
        }
        else
        {
            // the nodes are not extended so have no saved data regarding the simulation concentrations. Interpolate the nodal values
            // to the location of the cells within their associated FeMesh elements

            // for each cell in the population, find the element it belongs to then for each node in the elemnt (SPACE_DIM+1)
            // sum the weighted nodal values for the pde solutions, then save sum to the cell

            // The cells are not nodes of the mesh, so we must interpolate
            std::vector<double> solution_vector_at_cell(PROBLEM_DIM,0.0);
            //std::vector<double> remnants_of_previous_solution_vector_at_cell(PROBLEM_DIM,0.0);

            Element<ELEMENT_DIM,SPACE_DIM>* p_element = this->mpFeMesh->GetElement(elem_index);

            const ChastePoint<SPACE_DIM>& node_location = rCellPopulation.GetLocationOfCellCentre(*cell_iter);
            // find weights of the nodal values based on the location of the cell relative to nodes
            c_vector<double,SPACE_DIM+1> weights = p_element->CalculateInterpolationWeights(node_location);

            for(unsigned pd=0; pd<PROBLEM_DIM; pd++){
         
                // start with whatever was left over from the transport ode call
                //solution_vector_at_cell[pd] = remnants_of_previous_solution_vector_at_cell[pd];
                for (unsigned i=0; i<SPACE_DIM+1; i++)
                {
                    
                    double nodal_value = solution_repl[p_element->GetNodeGlobalIndex(i) + pd]; // serialised [node number * pde_index]
                    solution_vector_at_cell[pd] += nodal_value * weights(i);
                }
            }
       
            StateVariableRegister* p_bulk_register_pde = this->mpCoupledDomainField -> GetDomainStateVariableRegister();

            if (prop_collection.HasProperty<TransportCellProperty>())
            {
      
                boost::shared_ptr<TransportCellProperty> transport_cell_property = boost::static_pointer_cast<TransportCellProperty>(prop_collection.GetPropertiesType<TransportCellProperty>().GetProperty());
     
                //remnants_of_previous_solution_vector_at_cell = transport_cell_property -> GetInternalCellBoundaryConcentrationVector();
                
                // only a subset of the solution vector at the cell are used by the transport property
                StateVariableRegister* p_bulk_register_cell = transport_cell_property -> GetBulkStateVariableRegister();
                std::vector<double> subset_vector_at_cell(p_bulk_register_cell->GetNumberOfStateVariables(),0.0);
                unsigned domain_index=0;
                for(unsigned i=0; i<p_bulk_register_cell->GetNumberOfStateVariables();i++)
                {
                    if(p_bulk_register_pde->IsStateVariablePresent(p_bulk_register_cell->RetrieveStateVariableName(i)))
                    {
                        domain_index =  p_bulk_register_pde->RetrieveStateVariableIndex(p_bulk_register_cell->RetrieveStateVariableName(i));
                        subset_vector_at_cell[i] = solution_vector_at_cell[domain_index];
                    }
                    // if state not found in domain register then automatically has 0.0 value
                }

                transport_cell_property -> UpdateBulkConcentrationVector(subset_vector_at_cell);
            }

            if (prop_collection.HasProperty<MembraneCellProperty>())
            {
                boost::shared_ptr<MembraneCellProperty> membrane_cell_property = boost::static_pointer_cast<MembraneCellProperty>(prop_collection.GetPropertiesType<MembraneCellProperty>().GetProperty());

                //remnants_of_previous_solution_vector_at_cell = transport_cell_property -> GetInternalCellBoundaryConcentrationVector();
                // only a subset of the solution vector at the cell are used by the transport property
                StateVariableRegister* p_bulk_register_cell = membrane_cell_property -> GetBulkStateVariableRegister();
                std::vector<double> subset_vector_at_cell(p_bulk_register_cell->GetNumberOfStateVariables(),0.0);
                unsigned domain_index=0;
                for(unsigned i=0; i<p_bulk_register_cell->GetNumberOfStateVariables();i++)
                {
                    if(p_bulk_register_pde->IsStateVariablePresent(p_bulk_register_cell->RetrieveStateVariableName(i)))
                    {
                        domain_index =  p_bulk_register_pde->RetrieveStateVariableIndex(p_bulk_register_cell->RetrieveStateVariableName(i));
                        subset_vector_at_cell[i] = solution_vector_at_cell[domain_index];
                    }
                    // if state not found in domain register then automatically has 0.0 value
                }

                membrane_cell_property -> UpdateBulkConcentrationVector(subset_vector_at_cell);
            }
            // Find the element in the FE mesh that contains this cell. CellElementMap has been updated so use this.
            
            /*
            // not everything is transported?

            // save the cell values to their variable names
            for(unsigned pd=0; pd<PROBLEM_DIM; pd++){
                // nodal values as in the pde term
                cell_iter->GetCellData()->SetItem(this->mpPdeSystem -> GetStateVariableRegister()->RetrieveStateVariableName(pd), solution_vector_at_cell[pd]);
            }
            */

        }

        if (this->mOutputGradient)
        {
            
            // Now calculate the gradient of the solution and store this in CellVecData
            for(unsigned pd=0; pd<PROBLEM_DIM; pd++){
                c_vector<double, SPACE_DIM> solution_gradient = zero_vector<double>(SPACE_DIM); // change?
                Element<ELEMENT_DIM,SPACE_DIM>* p_element = this->mpFeMesh->GetElement(elem_index);
                // Calculate the basis functions at any point (e.g. zero) in the element
                c_matrix<double, SPACE_DIM, SPACE_DIM> jacobian, inverse_jacobian; // change? lloks like it's for only PROBLEM_DIM=1 but not nodal
                double jacobian_det;
                this->mpFeMesh->GetInverseJacobianForElement(elem_index, jacobian, jacobian_det, inverse_jacobian);
                const ChastePoint<SPACE_DIM> zero_point;
                c_matrix<double, SPACE_DIM, SPACE_DIM+1> grad_phi;
                LinearBasisFunction<SPACE_DIM>::ComputeTransformedBasisFunctionDerivatives(zero_point, inverse_jacobian, grad_phi);

                for (unsigned node_index=0; node_index<SPACE_DIM+1; node_index++)
                {
                    double nodal_value = solution_repl[p_element->GetNodeGlobalIndex(node_index) + pd];

                    for (unsigned j=0; j<SPACE_DIM; j++)
                    {
                        solution_gradient(j) += nodal_value* grad_phi(j, node_index); 
                    }
                }

                switch (SPACE_DIM)
                {
                    case 1:
                        cell_iter->GetCellData()->SetItem("var_node_"+this->mpPdeSystem -> GetStateVariableRegister()->RetrieveStateVariableName(pd)+"_grad_x", solution_gradient(0));
                        break;
                    case 2:
                        cell_iter->GetCellData()->SetItem("var_node_"+this->mpPdeSystem -> GetStateVariableRegister()->RetrieveStateVariableName(pd)+"_grad_x", solution_gradient(0));
                        cell_iter->GetCellData()->SetItem("var_node_"+this->mpPdeSystem -> GetStateVariableRegister()->RetrieveStateVariableName(pd)+"_grad_y", solution_gradient(1));
                        break;
                    case 3:
                        cell_iter->GetCellData()->SetItem("var_node_"+this->mpPdeSystem -> GetStateVariableRegister()->RetrieveStateVariableName(pd)+"_grad_x", solution_gradient(0));
                        cell_iter->GetCellData()->SetItem("var_node_"+this->mpPdeSystem -> GetStateVariableRegister()->RetrieveStateVariableName(pd)+"_grad_y", solution_gradient(1));
                        cell_iter->GetCellData()->SetItem("var_node_"+this->mpPdeSystem -> GetStateVariableRegister()->RetrieveStateVariableName(pd)+"_grad_z", solution_gradient(2));
                        break;
                    default:
                        NEVER_REACHED;
                }
            }
        }

    }
 //   std::cout<<"AbstractBoxDomainPdeSystemModifier - UpdateCellData - end"<<std::endl;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
void AbstractBoxDomainPdeSystemModifier<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::InitialiseCellPdeElementMap(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation)
{
    mCellPdeElementMap.clear();

    // Find the element of mpFeMesh that contains each cell and populate mCellPdeElementMap
    for (typename AbstractCellPopulation<SPACE_DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        const ChastePoint<SPACE_DIM>& r_position_of_cell = rCellPopulation.GetLocationOfCellCentre(*cell_iter);
        
        unsigned elem_index = this->mpFeMesh->GetContainingElementIndex(r_position_of_cell);
        mCellPdeElementMap[*cell_iter] = elem_index;
    }
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
void AbstractBoxDomainPdeSystemModifier<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::UpdateCellPdeElementMap(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation)
{ //std::cout<<"AbstractBoxDomainPdeSystemModifier - UpdateCellPdeElementMap -start"<<std::endl;
    // Find the element of mpCoarsePdeMesh that contains each cell and populate mCellPdeElementMap
    for (typename AbstractCellPopulation<SPACE_DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        const ChastePoint<SPACE_DIM>& r_position_of_cell = rCellPopulation.GetLocationOfCellCentre(*cell_iter);

        unsigned elem_index = this->mpFeMesh->GetContainingElementIndexWithInitialGuess(r_position_of_cell, mCellPdeElementMap[*cell_iter]);

        mCellPdeElementMap[*cell_iter] = elem_index;
    }
  //  std::cout<<"AbstractBoxDomainPdeSystemModifier - UpdateCellPdeElementMap -end"<<std::endl;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
void AbstractBoxDomainPdeSystemModifier<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractPdeSystemModifier<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::OutputSimulationModifierParameters(rParamsFile);
}


template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
c_vector<double,SPACE_DIM>  AbstractBoxDomainPdeSystemModifier<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::GetMeshOffset()
{
    return mOffset;
}

#endif /*ABSTRACTBOXDOMAINPDESYSTEMMODIFIER_HPP_*/
