#ifndef TESTCHEMICALSPATIALSIMULATION_HPP_
#define TESTCHEMICALSPATIALSIMULATION_HPP_

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
#include "BoundaryConditionsContainer.hpp"
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


class TestChemicalSpatialSimulation : public AbstractCellBasedTestSuite
{
public:

    void TestTemplateForChemicalSimulation()
    {
        
        //------------------------------------------------------------------------------//
        //                          Experiment variables                                //
        //------------------------------------------------------------------------------//

        // Variables for the user modify
        std::string dataFileRoot = "/home/chaste/projects/ChemChaste/src/Data/TemplateChemicalSimulation/";
        std::string domainFilename = "Domain.csv";
        std::string domainKeyFilename = "DomainKey.csv";
        std::string odeLabelFilename = "NodeSelector.csv";
        std::string odeKeyFilename = "OdeReactionFileKey.csv";
        std::string diffusionFilename = "DiffusionDatabaseFile.csv";
        std::string initialConditionsFilename = "InitialConditionFile.csv";
        std::string boundaryConditionsFilename = "BoundaryConditionFile.csv";
        // System properties
        const unsigned probDim =6; // need to set manually to the number of diffusive variables for the pde solver to solve
        // solver properties
        double t_end = 10;
        double simulationTimeStep = 1e-2;
        double samplingTimeStep = 1e-2;
        std::string output_filename = "TestTemplateForChemicalSimulation";
        
        //------------------------------------------------------------------------------//
        //                          Simulation construct                                //
        //------------------------------------------------------------------------------//

        // run the domain field set up and parse files
        ChemicalDomainField* p_field = new ChemicalDomainField(dataFileRoot,dataFileRoot+domainFilename, dataFileRoot+domainKeyFilename, dataFileRoot+odeLabelFilename, dataFileRoot+odeKeyFilename, dataFileRoot+diffusionFilename);

        // check that the file input problem dimension is the same as the user defined problem dimension
        std::cout<<"File probDim: "<<p_field ->GetProblemDimensions()<<std::endl;
        std::cout<<"User probDim: "<<probDim<<std::endl;

        MutableMesh<2,2>* p_mesh = p_field ->GetMeshGenerator() -> GetMesh();
        const unsigned spaceDim=2;
        const unsigned elementDim=2;
        
    
        std::cout<<"Domain numNodes: "<<p_field -> GetDomainMesh() -> GetNumNodes()<<std::endl;
        std::cout<<"Class numNodes: "<<p_mesh -> GetNumNodes()<<std::endl;
        // process initial conditions
        std::cout<<"Process initial conditions"<<std::endl;

        p_field -> ParseInitialConditionsFromFile(dataFileRoot+initialConditionsFilename);
        
        // PETSc Vec
        Vec initial_condition = PetscTools::CreateVec(p_field -> GetInitialNodeConditions());

        

        // process boundary conditions
        p_field -> ParseBoundaryConditionsFromFile(dataFileRoot+boundaryConditionsFilename);

        std::vector<std::string> boundaryConditionTypes = p_field -> GetBoundaryConditionTypes();

        std::vector<double> boundaryConditionValues = p_field -> GetBoundaryConditionValues();

        // Process Boundary Conditions
        std::cout<<"Process Boundary Conditions"<<std::endl;
        BoundaryConditionsContainer<elementDim,spaceDim,probDim> bcc;

        std::vector<ConstBoundaryCondition<spaceDim>*> vectorConstBCs;
        
        for(unsigned pdeDim=0; pdeDim<probDim; pdeDim++){
            vectorConstBCs.push_back(new ConstBoundaryCondition<spaceDim>(boundaryConditionValues[pdeDim]));
        }

        for(unsigned pdeDim=0; pdeDim<probDim; pdeDim++)
        {
            if(boundaryConditionTypes[pdeDim]=="Dirichlet"||boundaryConditionTypes[pdeDim]=="dirichlet"||boundaryConditionTypes[pdeDim]=="D"||boundaryConditionTypes[pdeDim]=="d")
            {
                for (TetrahedralMesh<elementDim,spaceDim>::BoundaryNodeIterator node_iter = p_mesh->GetBoundaryNodeIteratorBegin();
                 node_iter != p_mesh->GetBoundaryNodeIteratorEnd();
                 ++node_iter)
                {
                    bcc.AddDirichletBoundaryCondition(*node_iter, vectorConstBCs[pdeDim], pdeDim);
                }
            }else if(boundaryConditionTypes[pdeDim]=="Neumann"||boundaryConditionTypes[pdeDim]=="neumann"||boundaryConditionTypes[pdeDim]=="N"||boundaryConditionTypes[pdeDim]=="n")
            
            {
                for (TetrahedralMesh<elementDim,spaceDim>::BoundaryElementIterator boundary_iter = p_mesh->GetBoundaryElementIteratorBegin();
                boundary_iter != p_mesh->GetBoundaryElementIteratorEnd();
                boundary_iter++)
                {
                    bcc.AddNeumannBoundaryCondition(*boundary_iter, vectorConstBCs[pdeDim], pdeDim);
                }
            }
        }
      
        

        // pde system
        std::cout<<"Pde"<<std::endl;
        InhomogenousParabolicPdeOdeSystem<elementDim, spaceDim, probDim> pde(p_field);
  

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
        solver.SetTimes(0, t_end);
        solver.SetTimeStep(simulationTimeStep);
        solver.SetSamplingTimeStep(samplingTimeStep);
        solver.SetOutputDirectory(output_filename);
        solver.SetInitialCondition(initial_condition);
        // solve
        std::cout<<"Solve"<<std::endl;

        solver.SolveAndWriteResultsToFile();
        std::cout<<"Clean"<<std::endl;

        // clean
        PetscTools::Destroy(initial_condition);
        */
    }


    void TestMembraneModel()
    {
         //------------------------------------------------------------------------------//
        //                          Experiment variables                                //
        //------------------------------------------------------------------------------//

        // Variables for the user modify
        std::string dataFileRoot = "/home/chaste/projects/ChemChaste/src/Data/MembraneModel/";
        std::string domainFilename = "MembraneDomain.csv";
        std::string domainKeyFilename = "MembraneDomainLabelKey.csv";
        std::string odeLabelFilename = "MembraneOdeSelector.csv";
        std::string odeKeyFilename = "MembraneOdeKey.csv";
        std::string diffusionFilename = "Membrane_DiffusionDatabase.csv";
        std::string initialConditionsFilename = "MembraneInitialConditions.csv";
        std::string boundaryConditionsFilename = "Membrane_Boundary_Conditions.csv";
        // System properties
        const unsigned probDim =4; // need to set manually to the number of diffusive variables for the pde solver to solve
        // solver properties
        double t_end = 10;
        double simulationTimeStep = 1e-3;
        double samplingTimeStep = 1e-3;
        std::string output_filename = "TestMembraneModelOutput_SMB_test";
        
        //------------------------------------------------------------------------------//
        //                          Simulation construct                 
                       //
        //------------------------------------------------------------------------------//

        // run the domain field set up and parse files
        ChemicalDomainField* p_field = new ChemicalDomainField(dataFileRoot,dataFileRoot+domainFilename, dataFileRoot+domainKeyFilename, dataFileRoot+odeLabelFilename, dataFileRoot+odeKeyFilename, dataFileRoot+diffusionFilename);

        // check that the file input problem dimension is the same as the user defined problem dimension
        std::cout<<"File probDim: "<<p_field ->GetProblemDimensions()<<std::endl;
        std::cout<<"User probDim: "<<probDim<<std::endl;

        MutableMesh<2,2>* p_mesh = p_field ->GetMeshGenerator() -> GetMesh();
        const unsigned spaceDim=2;
        const unsigned elementDim=2;
        
    
        std::cout<<"Domain numNodes: "<<p_field -> GetDomainMesh() -> GetNumNodes()<<std::endl;
        std::cout<<"Class numNodes: "<<p_mesh -> GetNumNodes()<<std::endl;
        // process initial conditions
        std::cout<<"Process initial conditions"<<std::endl;

        p_field -> ParseInitialConditionsFromFile(dataFileRoot+initialConditionsFilename);
        
        // PETSc Vec
        Vec initial_condition = PetscTools::CreateVec(p_field -> GetInitialNodeConditions());

        

        // process boundary conditions
        p_field -> ParseBoundaryConditionsFromFile(dataFileRoot+boundaryConditionsFilename);

        std::vector<std::string> boundaryConditionTypes = p_field -> GetBoundaryConditionTypes();

        std::vector<double> boundaryConditionValues = p_field -> GetBoundaryConditionValues();

        // Process Boundary Conditions
        std::cout<<"Process Boundary Conditions"<<std::endl;
        BoundaryConditionsContainer<elementDim,spaceDim,probDim> bcc;

        std::vector<ConstBoundaryCondition<spaceDim>*> vectorConstBCs;
        
        for(unsigned pdeDim=0; pdeDim<probDim; pdeDim++){
            vectorConstBCs.push_back(new ConstBoundaryCondition<spaceDim>(boundaryConditionValues[pdeDim]));
        }

        for(unsigned pdeDim=0; pdeDim<probDim; pdeDim++)
        {
            if(boundaryConditionTypes[pdeDim]=="Dirichlet"||boundaryConditionTypes[pdeDim]=="dirichlet"||boundaryConditionTypes[pdeDim]=="D"||boundaryConditionTypes[pdeDim]=="d")
            {
                for (TetrahedralMesh<elementDim,spaceDim>::BoundaryNodeIterator node_iter = p_mesh->GetBoundaryNodeIteratorBegin();
                 node_iter != p_mesh->GetBoundaryNodeIteratorEnd();
                 ++node_iter)
                {
                    bcc.AddDirichletBoundaryCondition(*node_iter, vectorConstBCs[pdeDim], pdeDim);
                }
            }else if(boundaryConditionTypes[pdeDim]=="Neumann"||boundaryConditionTypes[pdeDim]=="neumann"||boundaryConditionTypes[pdeDim]=="N"||boundaryConditionTypes[pdeDim]=="n")
            
            {
                for (TetrahedralMesh<elementDim,spaceDim>::BoundaryElementIterator boundary_iter = p_mesh->GetBoundaryElementIteratorBegin();
                boundary_iter != p_mesh->GetBoundaryElementIteratorEnd();
                boundary_iter++)
                {
                    bcc.AddNeumannBoundaryCondition(*boundary_iter, vectorConstBCs[pdeDim], pdeDim);
                }
            }
        }
      
        

        // pde system
        std::cout<<"Pde"<<std::endl;
        InhomogenousParabolicPdeOdeSystem<elementDim, spaceDim, probDim> pde(p_field);
        

        // used the explicitly defined EulerIvpOdeSolver rather than the default defined BackwardEuler method
        std::vector<boost::shared_ptr<AbstractIvpOdeSolver>> odeSolverSystem;
        std::vector<AbstractInhomogenousOdeSystemForCoupledPdeSystem*> odeSystem;
        for(unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            std::cout<<"node number: "<<i<<std::endl;
            for(unsigned j=0; j<p_field->GetOdeSystem()[i] -> GetReactionSystem() -> GetReactionVector().size();j++)
            {
                std::cout<< p_field->GetOdeSystem()[i] -> GetReactionSystem() -> GetReactionVector()[j] -> GetReactionType()<<std::endl;
            }
            
            odeSystem.push_back(p_field->GetOdeSystem()[i]);
            boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);
            odeSolverSystem.push_back(p_solver);
        }

        std::cout<<"Solver"<<std::endl;
        // solver
        InhomogenousCoupledPdeOdeSolver<elementDim,spaceDim,probDim> solver(p_mesh, &pde, &bcc, odeSystem, odeSolverSystem);

        // solver properties
        solver.SetTimes(0, t_end);
        solver.SetTimeStep(simulationTimeStep);
        solver.SetSamplingTimeStep(samplingTimeStep);
        solver.SetOutputDirectory("TestMembraneModelOutput_SMB_test_oldWay");
        solver.SetInitialCondition(initial_condition);
        // solve
        std::cout<<"Solve"<<std::endl;

        solver.SolveAndWriteResultsToFile();
        std::cout<<"Clean"<<std::endl;

        // clean
        PetscTools::Destroy(initial_condition);
    }

};


#endif