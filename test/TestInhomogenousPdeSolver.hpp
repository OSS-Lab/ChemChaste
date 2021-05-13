#ifndef TESTINHOMOGENOUSPDESOLVER_HPP_
#define TESTINHOMOGENOUSPDESOLVER_HPP_

// chaste includes
#include <cxxtest/TestSuite.h>
#include "UblasIncludes.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "PetscSetupAndFinalize.hpp"

// general includes
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <tuple>
#include <cmath>

// Extended boundary conditions container
#include "BoundaryConditionsContainer_extended.hpp"

// chaste PdeOde includes
#include "HoneycombMeshGenerator.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "LinearParabolicPdeSystemWithCoupledOdeSystemSolver.hpp"
//#include "BoundaryConditionsContainer.hpp"
#include "ConstBoundaryCondition.hpp"
#include "OutputFileHandler.hpp"
#include "RandomNumberGenerator.hpp"
#include "TrianglesMeshReader.hpp"

// inhomogenous solver
#include "StateVariableRegister.hpp"
#include "InhomogenousParabolicPdeOdeSystem.hpp"
#include "InhomogenousCoupledPdeOdeSolver.hpp"

#include "InhomogenousOdeSchnackenbergCoupledPdeOdeSystem.hpp"
#include "InhomogenousOdeConsumerProducer.hpp"
#include "AbstractInhomogenousChemicalOdeSystemForCoupledPdeSystem.hpp"

// custom pdeOde includes
#include "PdeSchnackenbergCoupledPdeOdeSystem.hpp"
#include "PdeConsumerProducer.hpp"
#include "InhomogenousParabolicPdeOdeSystem.hpp"

// domain field includes
#include "AbstractDomainField.hpp"
#include "ChemicalDomainField.hpp"




// class to solve pde systems in which the nodes contain different ode systems, affecting different vairables
// use StateVariableVector for a system to control the variable names and indices

class TestInhomogenousPdeSolver : public AbstractCellBasedTestSuite
{
public:

    void TestStateVariableRegister()
    {
        std::cout<<"-----------------------"<<std::endl;
        std::cout<<"StateVariableRegister"<<std::endl;
        std::cout<<"-----------------------"<<std::endl;


        std::vector<std::string> system_1 = {"A", "B"};
        std::vector<std::string> system_2 = {"B", "C"};
        

        StateVariableRegister* register_1_2  = new StateVariableRegister(system_1); 

        std::cout<<"Read state variable register:"<<std::endl;
        for(unsigned i=0; i<register_1_2 -> GetNumberOfStateVariables(); i++)
        {
            std::cout<<register_1_2 -> RetrieveStateVariableName(i)<<std::endl;
        }

        std::cout<<"Add new system to state register"<<std::endl;

        register_1_2 -> AddStateVariableVector(system_2);

        std::cout<<"Read state variable register:"<<std::endl;
        for(unsigned i=0; i<register_1_2 -> GetNumberOfStateVariables(); i++)
        {
            std::cout<<register_1_2 -> RetrieveStateVariableName(i)<<std::endl;
        }

        std::cout<<"Add new variable D to state register"<<std::endl;

        register_1_2 -> AddStateVariable("D");

        std::cout<<"Read state variable register:"<<std::endl;
        for(unsigned i=0; i<register_1_2 -> GetNumberOfStateVariables(); i++)
        {
            std::cout<<register_1_2 -> RetrieveStateVariableName(i)<<std::endl;
        }

        std::cout<<"Read index of variable in state register"<<std::endl;

        std::vector<std::string> variable_set = {"A","C"};
        for(unsigned i=0; i<variable_set.size(); i++)
        {
            std::cout<<variable_set[i]<<": "<<register_1_2 -> RetrieveStateVariableIndex(variable_set[i])<<std::endl;
        }


        std::cout<<"Remove variable B from state register"<<std::endl;

        register_1_2 -> RemoveStateVariable("B");

        std::cout<<"Read state variable register:"<<std::endl;
        for(unsigned i=0; i<register_1_2 -> GetNumberOfStateVariables(); i++)
        {
            std::cout<<register_1_2 -> RetrieveStateVariableName(i)<<std::endl;
        }

        

        std::cout<<"Match registers and crossreference indices: "<<std::endl;
        std::vector<std::string> system_3 = {"D", "C"};
        StateVariableRegister* register_3  = new StateVariableRegister(system_3); 

        std::cout<<"Read state variable register_3:"<<std::endl;
        for(unsigned i=0; i<register_3 -> GetNumberOfStateVariables(); i++)
        {
            std::cout<<register_3 -> RetrieveStateVariableName(i)<<std::endl;
        }
        
        
        

        std::cout<<"Common variables: "<<std::endl;
        std::vector<std::string> matchedNames = register_1_2 ->FindCommonNamesInRegisters(register_3);
        for(unsigned i=0; i<matchedNames.size(); i++)
        {
            std::cout<<matchedNames[i]<<std::endl;
        }
        std::cout<<"System_3 variables locations in systems_1_2: "<<std::endl;
        std::vector<unsigned> matchedIndicesThis = register_1_2 ->FindIndicesInThisRegister(register_3);
        for(unsigned i=0; i<matchedIndicesThis.size(); i++)
        {
            std::cout<<matchedIndicesThis[i]<<std::endl;
        }
        std::cout<<"System_1_2 variables locations in systems_3: "<<std::endl;
        std::vector<unsigned> matchedIndicesThat = register_1_2 ->FindIndicesInThatRegister(register_3);
        for(unsigned i=0; i<matchedIndicesThat.size(); i++)
        {
            std::cout<<matchedIndicesThat[i]<<std::endl;
        }

    }

    void TestHomogenousDiffusionInhomogenousOde()
    {

        // direct call to solver for a hard code mesh 


        // make mixed nodal ode vector of consumerproducer and schnackenberg odes
        // consumerproducer acting on state variables A, U

        // scnackenberg acting on state variables U, V


        // system properties
        const unsigned probDim =3;
        const unsigned spaceDim=2;
        const unsigned elementDim=2;
        std::vector<double> initValues = {1.0, 1.0, 1.0};
        std::vector<double> bcValues = {1.0, 1.0, 1.0};
        std::vector<double> diffusionVector = {1e-1, 1e-2, 1e-3};

        // mesh
        HoneycombMeshGenerator generator(15, 18, 0);
        MutableMesh<elementDim,spaceDim>* p_mesh = generator.GetMesh();

        // Process Boundary Conditions
        std::cout<<"Process Boundary Conditions"<<std::endl;
        BoundaryConditionsContainer<elementDim,spaceDim,probDim> bcc;
        std::vector<bool> areNeumannBoundaryConditions(probDim, false);
        std::vector<ConstBoundaryCondition<spaceDim>*> vectorConstBCs;
        
        for(unsigned pdeDim=0; pdeDim<probDim; pdeDim++){
            vectorConstBCs.push_back(new ConstBoundaryCondition<spaceDim>(bcValues[pdeDim]));
        }
        
        for(unsigned pdeDim=0; pdeDim<probDim; pdeDim++)
        {
            if(areNeumannBoundaryConditions[pdeDim]==false)
            {
                for (TetrahedralMesh<elementDim,spaceDim>::BoundaryNodeIterator node_iter = p_mesh->GetBoundaryNodeIteratorBegin();
                 node_iter != p_mesh->GetBoundaryNodeIteratorEnd();
                 ++node_iter)
                {

                    bcc.AddDirichletBoundaryCondition(*node_iter, vectorConstBCs[pdeDim], pdeDim);
                }
            }else{
                for (TetrahedralMesh<elementDim,spaceDim>::BoundaryElementIterator boundary_iter = p_mesh->GetBoundaryElementIteratorBegin();
                boundary_iter != p_mesh->GetBoundaryElementIteratorEnd();
                boundary_iter++)
                {
                    bcc.AddNeumannBoundaryCondition(*boundary_iter, vectorConstBCs[pdeDim], pdeDim);
                }
            }
        }
        std::cout<<"Initial conditions"<<std::endl;
        // initial conditions
        std::vector<double> init_conds(probDim*p_mesh->GetNumNodes());
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {   // set as being a random perturbation about the boundary values
            for(unsigned pdeDim=0; pdeDim<probDim; pdeDim++)
            {   // serialised for nodes
                init_conds[probDim*i + pdeDim] = fabs(initValues[pdeDim] + RandomNumberGenerator::Instance()->ranf());
            }
        }
        // PETSc Vec
        std::cout<<"PETSc Vec"<<std::endl;
        Vec initial_condition = PetscTools::CreateVec(init_conds);



        // pde system
        std::cout<<"Pde"<<std::endl;
        InhomogenousParabolicPdeOdeSystem<elementDim, spaceDim, probDim> pde(diffusionVector);
        pde.SetStateVariableRegister(new StateVariableRegister({"A","U","V"}));
        // coupled ode system
        std::cout<<"Ode loop"<<std::endl;
        std::vector<AbstractInhomogenousOdeSystemForCoupledPdeSystem*> odeSystem;

        //EulerIvpOdeSolver euler_solver;
        std::vector<boost::shared_ptr<AbstractIvpOdeSolver>> odeSolverSystem;

        //std::vector<boost::shared_ptr<EulerIvpOdeSolver>> &odeSolverSystemr = *odeSolverSystem;

        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            // number of ode system objects must match the number of nodes, i.e the individual odes may be multi-dimensional
            if(i==22 || i==77)
            {
                InhomogenousOdeConsumerProducer* p_ConPro_new = new InhomogenousOdeConsumerProducer(1, 0, 0, 0);
                p_ConPro_new -> SetStateVariableRegister(new StateVariableRegister({"A","U"}));
                odeSystem.push_back(p_ConPro_new);

            }else if(i==27 || i==72)
            {
                InhomogenousOdeConsumerProducer* p_ConPro_new = new InhomogenousOdeConsumerProducer(0, 1, 0, 0);
                p_ConPro_new -> SetStateVariableRegister(new StateVariableRegister({"A","U"}));
                odeSystem.push_back(p_ConPro_new); 
            }else{
                InhomogenousOdeSchnackenbergCoupledPdeOdeSystem* p_Sch_new = new InhomogenousOdeSchnackenbergCoupledPdeOdeSystem(1, 1, 1, 1);
                p_Sch_new -> SetStateVariableRegister(new StateVariableRegister({"U","V"}));
                odeSystem.push_back(p_Sch_new);
            }
            boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);
            odeSolverSystem.push_back(p_solver);
        }
        std::cout<<"Number odesAtNodes: "<<p_mesh->GetNumNodes()<<std::endl;

        // used the explicitly defined EulerIvpOdeSolver rather than the default defined BackwardEuler method

        std::cout<<"Solver"<<std::endl;
        // solver
        InhomogenousCoupledPdeOdeSolver<elementDim,spaceDim,probDim> solver(p_mesh, &pde, &bcc, odeSystem, odeSolverSystem);

        // solver properties
        double t_end = 10;
        solver.SetTimes(0, t_end);
        solver.SetTimeStep(1e-2);
        solver.SetSamplingTimeStep(1e-2);
        solver.SetOutputDirectory("TestHomogenousDiffusionInhomogenousOde_test");
        solver.SetInitialCondition(initial_condition);
        // solve
        std::cout<<"Solve"<<std::endl;
        solver.SolveAndWriteResultsToFile();
        std::cout<<"Clean"<<std::endl;

        // clean
        PetscTools::Destroy(initial_condition); //double free or corruption (out)


    }

    void TestDomainField()
    {
        std::cout<<"TestAbstractDomainField"<<std::endl;

        std::string domainFilename = "/home/chaste/projects/ChemicalChaste/src/Data/Domain.csv";
        std::string domainKeyFilename = "/home/chaste/projects/ChemicalChaste/src/Data/DomainKey.csv";
        std::string odeLabelFilename = "/home/chaste/projects/ChemicalChaste/src/Data/NodeSelector.csv";
        std::string odeKeyFilename = "/home/chaste/projects/ChemicalChaste/src/Data/OdeReactionFileKey.csv";
        std::string diffusionFilename = "/home/chaste/projects/ChemicalChaste/src/Data/DiffusionDatabaseFile.csv";
        std::string initialConditionsFilename = "/home/chaste/projects/ChemicalChaste/src/Data/InitialConditionFile.csv";
        std::string boundaryConditionsFilename = "/home/chaste/projects/ChemicalChaste/src/Data/BoundaryConditionFile.csv";


        AbstractDomainField* p_field = new AbstractDomainField(domainFilename, domainKeyFilename, odeLabelFilename, odeKeyFilename, diffusionFilename);
        p_field -> ParseInitialConditionsFromFile(initialConditionsFilename);
        
        p_field -> ParseBoundaryConditionsFromFile(boundaryConditionsFilename);
   
        std::vector<double> initialConditions = p_field -> GetInitialNodeConditions();
        std::vector<double> boundaryConditions = p_field -> GetBoundaryConditionValues();
        std::vector<std::string> boundaryTypes = p_field -> GetBoundaryConditionTypes();
        
        
        for(unsigned i=0; i<initialConditions.size();i++)
        {
            std::cout<<"init: "<<initialConditions[i]<<std::endl;
        }

        for(unsigned i=0; i<boundaryConditions.size();i++)
        {
            std::cout<<"value: "<<boundaryConditions[i]<<std::endl;
            std::cout<<"type : "<<boundaryTypes[i]<<std::endl;
        }
        
        
    }

    void TestChemicalDomainField()
    {
        
        std::cout<<"----------------------------------------------"<<std::endl;
        std::cout<<"TestChemicalDomainField"<<std::endl;
        std::cout<<"----------------------------------------------"<<std::endl;
        // system properties        
        std::string reactionFileRoot = "/home/chaste/projects/ChemicalChaste/src/Data/";
        std::string domainFilename = "/home/chaste/projects/ChemicalChaste/src/Data/Domain.csv";
        std::string domainKeyFilename = "/home/chaste/projects/ChemicalChaste/src/Data/DomainKey.csv";
        std::string odeLabelFilename = "/home/chaste/projects/ChemicalChaste/src/Data/NodeSelector.csv";
        std::string odeKeyFilename = "/home/chaste/projects/ChemicalChaste/src/Data/OdeReactionFileKey.csv";
        std::string diffusionFilename = "/home/chaste/projects/ChemicalChaste/src/Data/DiffusionDatabaseFile.csv";

        ChemicalDomainField* p_field = new ChemicalDomainField(reactionFileRoot,domainFilename, domainKeyFilename, odeLabelFilename, odeKeyFilename, diffusionFilename);

        // print out the diffusion domain
        p_field -> PrintDiffusionDomain();
        // print out the diffusion domain mapped to the chaste cartesian grid
        p_field -> PrintMappedDiffusionDomain();
        // print out the ODE labels and keys
        p_field -> PrintDomainLabelKeys();
        // print out the ode domain
        p_field -> PrintODEDomain();
        // print out the ode domain mapped to the chaste cartesian grid
        p_field -> PrintMappedODEDomain();
        // print out the ODE labels and keys
        p_field -> PrintODELabelKeys();
        // print out the diffusion database
        p_field -> PrintDiffusionDatabase();


        // get domain variable register
        std::cout<<"Read domain state variable register"<<std::endl;
        StateVariableRegister* chemical_domain_register = p_field -> GetDomainStateVariableRegister();
        std::cout<<"number of state variables: "<<chemical_domain_register->GetNumberOfStateVariables()<<std::endl;
        for(unsigned index=0; index<chemical_domain_register->GetNumberOfStateVariables(); index++)
        {
            std::cout<<"State: "<<index<<" : "<<chemical_domain_register -> RetrieveStateVariableName(index)<<std::endl;
        }
    
        // get diffusion value at point
        ChastePoint<2> chastePoint(1.0,1.0);

        double diffusionValueAtPoint = p_field -> GetDiffusionValueBasedOnPoint(chastePoint,1);

        std::cout<<"Diffusion value at 1.0,1.0: "<<diffusionValueAtPoint<<std::endl;

        // run a singlar ode system
        unsigned node_selector = 92;
        AbstractInhomogenousChemicalOdeSystemForCoupledPdeSystem* p_ode_system = p_field -> GetOdeSystem()[node_selector];
        std::cout<<"number of Species: "<<p_ode_system -> GetNumberOfSpecies()<<std::endl;
        std::cout<<"number of Reactions: "<<p_ode_system -> GetNumberOfReactions()<<std::endl;
        std::vector<double> rY(p_ode_system -> GetNumberOfSpecies(),1.0);
        std::vector<double> rDY(p_ode_system -> GetNumberOfSpecies(),0.0);

        p_ode_system -> EvaluateYDerivatives(1.0,rY,rDY);

        for(unsigned index=0; index < p_ode_system -> GetNumberOfSpecies(); index++)
        {
            std::cout<<"Species "<<index<<": "<<"rDY: "<<rDY[index]<<std::endl;
        }

        
    }

    void TestInhomogenousChemicalPde()
    {
        
        std::cout<<"----------------------------------------------"<<std::endl;
        std::cout<<"TestInhomogenousChemicalPde"<<std::endl;
        std::cout<<"----------------------------------------------"<<std::endl;
        // system properties        
        std::string reactionFileRoot = "/home/chaste/projects/ChemicalChaste/src/Data/";
        std::string domainFilename = "/home/chaste/projects/ChemicalChaste/src/Data/Domain.csv";
        std::string domainKeyFilename = "/home/chaste/projects/ChemicalChaste/src/Data/DomainKey.csv";
        std::string odeLabelFilename = "/home/chaste/projects/ChemicalChaste/src/Data/NodeSelector.csv";
        std::string odeKeyFilename = "/home/chaste/projects/ChemicalChaste/src/Data/OdeReactionFileKey.csv";
        std::string diffusionFilename = "/home/chaste/projects/ChemicalChaste/src/Data/DiffusionDatabaseFileInhomogenousChemicalPde.csv";

        ChemicalDomainField* p_field = new ChemicalDomainField(reactionFileRoot,domainFilename, domainKeyFilename, odeLabelFilename, odeKeyFilename, diffusionFilename);

        std::cout<<"File probDim: "<<p_field ->GetProblemDimensions()<<std::endl;
        MutableMesh<2,2>* p_mesh = p_field ->GetMeshGenerator() -> GetMesh();
        std::cout<<"File mesh number nodes: "<<p_mesh->GetNumNodes()<<std::endl;
        
        // System properties
        const unsigned probDim =3; // need to set manually
        const unsigned spaceDim=2;
        const unsigned elementDim=2;


        std::vector<double> initValues(probDim,1.0);
        std::vector<double> bcValues = initValues;

        // Process Boundary Conditions
        std::cout<<"Process Boundary Conditions"<<std::endl;
        BoundaryConditionsContainer<elementDim,spaceDim,probDim> bcc;
        std::vector<bool> areNeumannBoundaryConditions(probDim, false);
        std::vector<ConstBoundaryCondition<spaceDim>*> vectorConstBCs;
        
        for(unsigned pdeDim=0; pdeDim<probDim; pdeDim++){
            vectorConstBCs.push_back(new ConstBoundaryCondition<spaceDim>(bcValues[pdeDim]));
        }

        for(unsigned pdeDim=0; pdeDim<probDim; pdeDim++)
        {
            if(areNeumannBoundaryConditions[pdeDim]==false)
            {
                for (TetrahedralMesh<elementDim,spaceDim>::BoundaryNodeIterator node_iter = p_mesh->GetBoundaryNodeIteratorBegin();
                 node_iter != p_mesh->GetBoundaryNodeIteratorEnd();
                 ++node_iter)
                {

                    bcc.AddDirichletBoundaryCondition(*node_iter, vectorConstBCs[pdeDim], pdeDim);
                }
            }else{
                for (TetrahedralMesh<elementDim,spaceDim>::BoundaryElementIterator boundary_iter = p_mesh->GetBoundaryElementIteratorBegin();
                boundary_iter != p_mesh->GetBoundaryElementIteratorEnd();
                boundary_iter++)
                {
                    bcc.AddNeumannBoundaryCondition(*boundary_iter, vectorConstBCs[pdeDim], pdeDim);
                }
            }
        }
      
        std::cout<<"Initial conditions"<<std::endl;
        // initial conditions
        std::vector<double> init_conds(probDim*p_field->GetDomainMesh()->GetNumNodes());
        for (unsigned i=0; i<p_field->GetDomainMesh()->GetNumNodes(); i++)
        {   // set as being a random perturbation about the boundary values
            for(unsigned pdeDim=0; pdeDim<probDim; pdeDim++)
            {   // serialised for nodes
                init_conds[probDim*i + pdeDim] =fabs(initValues[pdeDim] + RandomNumberGenerator::Instance()->ranf());
            }
        }

        
        // PETSc Vec
        std::cout<<"PETSc Vec"<<std::endl;
        Vec initial_condition = PetscTools::CreateVec(init_conds);

        // pde system
        std::cout<<"Pde"<<std::endl;
        InhomogenousParabolicPdeOdeSystem<elementDim, spaceDim, probDim> pde(p_field);
  


   
        std::cout<<"Number odesAtNodes: "<<p_mesh->GetNumNodes()<<std::endl;

        // used the explicitly defined EulerIvpOdeSolver rather than the default defined BackwardEuler method
        std::vector<boost::shared_ptr<AbstractIvpOdeSolver>> odeSolverSystem;
        std::vector<AbstractInhomogenousOdeSystemForCoupledPdeSystem*> odeSystem;
        for(unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            odeSystem.push_back(p_field->GetOdeSystem()[i]);
            boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);
            odeSolverSystem.push_back(p_solver);
        }

        std::cout<<"Solver"<<std::endl;
        // solver

        InhomogenousCoupledPdeOdeSolver<elementDim,spaceDim,probDim> solver(p_mesh, &pde, &bcc, odeSystem, odeSolverSystem);

        // solver properties
        double t_end = 10;
        solver.SetTimes(0, t_end);
        solver.SetTimeStep(1e-2);
        solver.SetSamplingTimeStep(1e-2);
        solver.SetOutputDirectory("TestInhomogenousParabolicPdeOdeSystemFromFile");
        solver.SetInitialCondition(initial_condition);
        // solve
        std::cout<<"Solve"<<std::endl;

        solver.SolveAndWriteResultsToFile();
        std::cout<<"Clean"<<std::endl;

        // clean
        PetscTools::Destroy(initial_condition);
        
    }
    
    void TestExtendedBoundaryConditionsContainer()
    {
        // set get test for 2D from TestBoundaryConditionsContainer
        std::cout<<"----------------------------------------------"<<std::endl;
        std::cout<<"TestExtendedBoundaryConditionsContainer"<<std::endl;
        std::cout<<"----------------------------------------------"<<std::endl;

        unsigned num_nodes = 10;
        const unsigned probDim =5;
        std::vector<double> bcValues(probDim,1.0);
        BoundaryConditionsContainer<2,2,probDim> bcc;
        std::vector<ConstBoundaryCondition<2>*> vectorConstBCs;

        std::vector<Node<2>*> nodes(num_nodes);

        for(unsigned pdeDim=0; pdeDim<probDim; pdeDim++){
            vectorConstBCs.push_back(new ConstBoundaryCondition<2>(pdeDim*bcValues[pdeDim]));
        }
        
        for (unsigned i=0; i<num_nodes; i++)
        {
            nodes[i] = new Node<2>(i,true,0,0); //Node(unsigned index, bool isBoundaryNode=false, double v1=0, double v2=0, double v3=0);

            for(unsigned pdeDim=0; pdeDim<probDim; pdeDim++)
            {
                bcc.AddDirichletBoundaryCondition(nodes[i], vectorConstBCs[pdeDim], pdeDim);
            }

        }
            

        for (unsigned i=0; i<num_nodes; i++)
        {
            
            for(unsigned j=0; j<probDim;j++)
            {
                double value = bcc.GetDirichletBCValue(nodes[i],j);
                std::cout<<"value node: "<<i<<" "<<value<<std::endl;
            }
            
        }

        for (unsigned i=0; i<num_nodes; i++)
        {
            delete nodes[i];
        }
    }


};


#endif