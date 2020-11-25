#ifndef TESTCHEMCHASTE_HPP_
#define TESTCHEMCHASTE_HPP_

#include "Cell_virtual.hpp"
#include "BoundaryConditionsContainer_extended.hpp"
// chaste includes
#include "ChasteHeaders.hpp"

// general includes
#include "GeneralHeaders.hpp"

//ChemChaste includes
#include "ChemChasteHeaders.hpp"
#include "ChemicalCellFromFile.hpp"

// test specific include
#include "ChemicalStructuresForTests.hpp"

#include "InhomogenousFisherPde.hpp"
#include "InhomogenousFisherDiffusiveInhibitionPde.hpp"

struct ControlStruct {
    bool ReactionSystemWithoutCells = false;
    bool Fisher = false;
    bool FisherDiffusiveInhibition = false;
    bool ReactionSystemWithoutCellsUsingInhomogenousSolver = false;
    bool ReactionSystemInhibitedDiffusionWithoutCellsUsingInhomogenousSolver = false;
    bool ReactionSystemWithCells = false;
    bool ReactionSystemWithNodeBasedCells = false;
} control;

class TestChemChaste : public AbstractCellBasedTestSuite
{
public:

    void TestReactionSystemWithoutCells()
    {
        if(control.ReactionSystemWithoutCells)
        {
            std::cout<<"Mass action kinetics of Schnackenberg using the AbstractChemicalOdeSystem Class"<<std::endl;
            std::cout<<"-----------------------------"<<std::endl;

            ChemicalStructuresForTests* chemical_structure = new ChemicalStructuresForTests();
            chemical_structure->SetUpPdeChemicalReactionSystem();
            AbstractReactionSystem* chemicalReactionSystem = chemical_structure -> rGetPtrPdeChemicalReactionSystem();

            // form ode system
            std::cout<<"SchnackenbergCoupledPdeOdeSystem  -As AbstractChemicalODeSystem"<<std::endl;
            // system properties
            const unsigned probDim =2;
            const unsigned spaceDim=2;
            const unsigned elementDim=2;
            std::vector<double> initValues = {2.0, 0.75};
            std::vector<double> bcValues = {0.0, 0.0};

            // mesh
            HoneycombMeshGenerator generator(100, 100, 0);
            MutableMesh<elementDim,spaceDim>* p_mesh = generator.GetMesh();


            std::cout<<"number mesh nodes: "<<p_mesh->GetNumNodes()<<std::endl;

            // Process Boundary Conditions
            BoundaryConditionsContainer<elementDim,spaceDim,probDim> bcc;
            std::vector<bool> areNeumannBoundaryConditions(probDim, true);
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
                    init_conds[probDim*i + pdeDim] = fabs(initValues[pdeDim]+ RandomNumberGenerator::Instance()->ranf());
                }
            }
            // PETSc Vec
            std::cout<<"PETSc Vec"<<std::endl;
            Vec initial_condition = PetscTools::CreateVec(init_conds);

            // coupled ode system
            std::cout<<"Ode loop"<<std::endl;
            std::vector<AbstractOdeSystemForCoupledPdeSystem*> odeSystem;
            for (unsigned i=0; i<p_mesh->GetNumNodes(); i++){
                // number of ode system objects must match the number of nodes, i.e the individual odes may be multi-dimensional
                odeSystem.push_back(new AbstractChemicalOdeForCoupledPdeSystem(chemicalReactionSystem));
            }

            // pde system
            PdeSchnackenbergCoupledPdeOdeSystem<elementDim, spaceDim, probDim> pde(odeSystem[0],1e-4, 1e-2 );
            //EulerIvpOdeSolver euler_solver; 
            boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);

            //double t_end = 100;
            std::cout<<"Solver"<<std::endl;
            // solver
            LinearParabolicPdeSystemWithCoupledOdeSystemSolver<elementDim,spaceDim,probDim> solver(p_mesh, &pde, &bcc,odeSystem,p_solver);

            // solver properties
        
            solver.SetTimes(0, 30);
            solver.SetTimeStep(1e-2);
            solver.SetSamplingTimeStep(1e-2);
            solver.SetOutputDirectory("TestReactionSystemWithoutCells_ver2");
            solver.SetInitialCondition(initial_condition);
            // solve
            std::cout<<"Solve"<<std::endl;
            solver.SolveAndWriteResultsToFile();
            std::cout<<"Clean"<<std::endl;

        }
    }

    void TestFisher()
    {
        if(control.Fisher)
        {
            // system properties
            const unsigned probDim =1;
            const unsigned spaceDim=2;
            const unsigned elementDim=2;
            double MeshStepSize = 1.0;
            std::vector<unsigned> MeshDimensions = {100,100,0};
            std::vector<double> initValuesHigh = {0.5};
            std::vector<double> initValuesLow = {0.0};
            std::vector<double> bcValues = {0.0};
            std::vector<double> areBCsNeumann = {true};
            std::vector<double> diffusionRates = {1.0};
            std::vector<double> growthRates = {1.0};
            std::vector<double> carryingCapacities = {1.0};


            // mesh
            TetrahedralMesh<elementDim,spaceDim>* p_mesh = new TetrahedralMesh<elementDim,spaceDim>();

            switch (spaceDim)
            {
                case 1:
                    p_mesh->ConstructRegularSlabMesh(MeshStepSize, MeshDimensions[0]);
                    break;
                case 2:
                    p_mesh->ConstructRegularSlabMesh(MeshStepSize, MeshDimensions[0], MeshDimensions[1]);
                    break;
                case 3:
                    p_mesh->ConstructRegularSlabMesh(MeshStepSize, MeshDimensions[0], MeshDimensions[1], MeshDimensions[2]);
                    break;
                default:
                    NEVER_REACHED;
            }
            

            // Process Boundary Conditions
            std::cout<<"Process Boundary Conditions"<<std::endl;
            BoundaryConditionsContainer<elementDim,spaceDim,probDim> bcc;
            std::vector<bool> areNeumannBoundaryConditions(probDim, true);
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

            // initial conditions
            std::vector<double> init_conds(probDim*p_mesh->GetNumNodes(),0.0);
            unsigned columnNum = 0;
            unsigned rowNum = 0;
            for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
            {   // set as being a random perturbation about the boundary values
                std::cout<<"node number: "<<i<<std::endl;
                if(spaceDim==2)
                {
                    columnNum = 0;
                    rowNum = 0;
                    //std::cout<<i-rowNum*(MeshDimensions[0]+1)<<std::endl;
                    while(i >= rowNum*(MeshDimensions[0]+1))
                    {
                        rowNum = rowNum + 1;
                    
                    }
                    
                    columnNum = i - (rowNum-1)*(MeshDimensions[0]+1);
                    std::cout<<"Column number: "<<columnNum<<" rowNum: "<<rowNum<<" mesh dim: "<<MeshDimensions[0]<<std::endl;
                    if(columnNum<10)
                    {
                        for(unsigned pdeDim=0; pdeDim<probDim; pdeDim++)
                        {   // serialised for nodes
                            init_conds[probDim*i + pdeDim] = fabs(initValuesHigh[pdeDim]);// + RandomNumberGenerator::Instance()->ranf());
                        }
                    }
                    else
                    {
                        for(unsigned pdeDim=0; pdeDim<probDim; pdeDim++)
                        {   // serialised for nodes
                            init_conds[probDim*i + pdeDim] = fabs(initValuesLow[pdeDim]);// + RandomNumberGenerator::Instance()->ranf());
                        }
                    }
                    

                }
                else
                {
                    for(unsigned pdeDim=0; pdeDim<probDim; pdeDim++)
                    {   // serialised for nodes
                        init_conds[probDim*i + pdeDim] = fabs(initValuesLow[pdeDim] + RandomNumberGenerator::Instance()->ranf());
                    }
                }

            }
            // PETSc Vec
            Vec initial_condition = PetscTools::CreateVec(init_conds);

            std::vector<AbstractInhomogenousOdeSystemForCoupledPdeSystem*> odeSystem;
            std::vector<boost::shared_ptr<AbstractIvpOdeSolver> > solverSystem;
            for (unsigned i=0; i<p_mesh->GetNumNodes(); i++){
                // number of ode system objects must match the number of nodes, i.e the individual odes may be multi-dimensional
            //    odeSystem.push_back(new AbstractInhomogenousChemicalOdeSystemForCoupledPdeSystem(chemicalReactionSystem));
            //    boost::shared_ptr<AbstractIvpOdeSolver> p_solver(new EulerIvpOdeSolver);
            //    solverSystem.push_back(p_solver);
            }
    
            // pde system
            //InhomogenousParabolicPdeForCoupledOdeSystemTemplated<elementDim, spaceDim, probDim> pde(chemical_structure->rGetPtrChemicalDomain());

            InhomogenousFisherPde<elementDim, spaceDim, probDim> pde(diffusionRates,growthRates,carryingCapacities);

            std::cout<<"Solver"<<std::endl;
            // solver
            InhomogenousCoupledPdeOdeSolverTemplated<elementDim,spaceDim,probDim> solver(p_mesh, &pde, &bcc);

            // solver properties
        
            solver.SetTimes(0, 60);
            solver.SetTimeStep(1);
            solver.SetSamplingTimeStep(1);
            solver.SetOutputDirectory("TestFisherPaper");
            solver.SetInitialCondition(initial_condition);
            // solve
            std::cout<<"Solve"<<std::endl;
            solver.SolveAndWriteResultsToFile();
            std::cout<<"Clean"<<std::endl;
        }
    }

    void TestFisherDiffusiveInhibition()
    {
        if(control.FisherDiffusiveInhibition)
        {
            // system properties
            const unsigned probDim =1;
            const unsigned spaceDim=2;
            const unsigned elementDim=2;
            double MeshStepSize = 1.0;
            std::vector<unsigned> MeshDimensions = {100,100,0};
            std::vector<double> initValuesHigh = {0.5};
            std::vector<double> initValuesLow = {0.0};
            std::vector<double> bcValues = {0.0};
            std::vector<double> areBCsNeumann = {true};
            std::vector<double> diffusionRates = {1.0};
            std::vector<double> growthRates = {1.0};
            std::vector<double> carryingCapacities = {1.0};


            // mesh
            TetrahedralMesh<elementDim,spaceDim>* p_mesh = new TetrahedralMesh<elementDim,spaceDim>();

            switch (spaceDim)
            {
                case 1:
                    p_mesh->ConstructRegularSlabMesh(MeshStepSize, MeshDimensions[0]);
                    break;
                case 2:
                    p_mesh->ConstructRegularSlabMesh(MeshStepSize, MeshDimensions[0], MeshDimensions[1]);
                    break;
                case 3:
                    p_mesh->ConstructRegularSlabMesh(MeshStepSize, MeshDimensions[0], MeshDimensions[1], MeshDimensions[2]);
                    break;
                default:
                    NEVER_REACHED;
            }
            

            // Process Boundary Conditions
            std::cout<<"Process Boundary Conditions"<<std::endl;
            BoundaryConditionsContainer<elementDim,spaceDim,probDim> bcc;
            std::vector<bool> areNeumannBoundaryConditions(probDim, true);
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

            // initial conditions
            std::vector<double> init_conds(probDim*p_mesh->GetNumNodes(),0.0);
            unsigned columnNum = 0;
            unsigned rowNum = 0;
            for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
            {   // set as being a random perturbation about the boundary values
                std::cout<<"node number: "<<i<<std::endl;
                if(spaceDim==2)
                {
                    columnNum = 0;
                    rowNum = 0;
                    //std::cout<<i-rowNum*(MeshDimensions[0]+1)<<std::endl;
                    while(i >= rowNum*(MeshDimensions[0]+1))
                    {
                        rowNum = rowNum + 1;
                    
                    }
                    
                    columnNum = i - (rowNum-1)*(MeshDimensions[0]+1);
                    std::cout<<"Column number: "<<columnNum<<" rowNum: "<<rowNum<<" mesh dim: "<<MeshDimensions[0]<<std::endl;
                    if(columnNum<10)
                    {
                        for(unsigned pdeDim=0; pdeDim<probDim; pdeDim++)
                        {   // serialised for nodes
                            init_conds[probDim*i + pdeDim] = fabs(initValuesHigh[pdeDim]);// + RandomNumberGenerator::Instance()->ranf());
                        }
                    }
                    else
                    {
                        for(unsigned pdeDim=0; pdeDim<probDim; pdeDim++)
                        {   // serialised for nodes
                            init_conds[probDim*i + pdeDim] = fabs(initValuesLow[pdeDim]);// + RandomNumberGenerator::Instance()->ranf());
                        }
                    }
                    

                }
                else
                {
                    for(unsigned pdeDim=0; pdeDim<probDim; pdeDim++)
                    {   // serialised for nodes
                        init_conds[probDim*i + pdeDim] = fabs(initValuesLow[pdeDim] + RandomNumberGenerator::Instance()->ranf());
                    }
                }

            }
            // PETSc Vec
            Vec initial_condition = PetscTools::CreateVec(init_conds);

            std::vector<AbstractInhomogenousOdeSystemForCoupledPdeSystem*> odeSystem;
            std::vector<boost::shared_ptr<AbstractIvpOdeSolver> > solverSystem;

            // pde system
            InhomogenousFisherDiffusiveInhibitionPde<elementDim, spaceDim, probDim> pde(diffusionRates,growthRates,carryingCapacities);

            std::cout<<"Solver"<<std::endl;
            // solver
            InhomogenousCoupledPdeOdeSolverTemplated<elementDim,spaceDim,probDim> solver(p_mesh, &pde, &bcc);

            // solver properties
        
            solver.SetTimes(0, 60);
            solver.SetTimeStep(1);
            solver.SetSamplingTimeStep(1);
            solver.SetOutputDirectory("TestFisherPaperDiffusionInhibition_2_close_regions");
            solver.SetInitialCondition(initial_condition);
            // solve
            std::cout<<"Solve"<<std::endl;
            solver.SolveAndWriteResultsToFile();
            std::cout<<"Clean"<<std::endl;
        }
    }

    void TestReactionSystemWithoutCellsInhomogenousSolver()
    {

        if(control.ReactionSystemWithoutCellsUsingInhomogenousSolver)
        {
            std::cout<<"Mass action kinetics of Schnackenberg using the ChemicalDomain Class"<<std::endl;
            std::cout<<"-----------------------------"<<std::endl;
            
            ChemicalStructuresForTests* chemical_structure = new ChemicalStructuresForTests();
            chemical_structure->SetUpPdeChemicalReactionSystem();
            chemical_structure->SetUpChemicalDomainField();
            AbstractReactionSystem* chemicalReactionSystem = chemical_structure -> rGetPtrPdeChemicalReactionSystem();

            // form ode system
            std::cout<<"SchnackenbergCoupledPdeOdeSystem  -As AbstractChemicalODeSystem"<<std::endl;
            // system properties
            const unsigned probDim =2;
            const unsigned spaceDim=2;
            const unsigned elementDim=2;
            double MeshStepSize = 1.0;
            std::vector<unsigned> MeshDimensions = {100,100,0};
            std::vector<double> initValues = {2.0, 0.75};
            std::vector<double> bcValues = {0.0, 0.0};
            std::vector<double> diffusionRates = {0.05, 0.05};

            // mesh
            TetrahedralMesh<elementDim,spaceDim>* p_mesh = new TetrahedralMesh<elementDim,spaceDim>();

            switch (spaceDim)
            {
                case 1:
                    p_mesh->ConstructRegularSlabMesh(MeshStepSize, MeshDimensions[0]);
                    break;
                case 2:
                    p_mesh->ConstructRegularSlabMesh(MeshStepSize, MeshDimensions[0], MeshDimensions[1]);
                    break;
                case 3:
                    p_mesh->ConstructRegularSlabMesh(MeshStepSize, MeshDimensions[0], MeshDimensions[1], MeshDimensions[2]);
                    break;
                default:
                    NEVER_REACHED;
            }

            // Process Boundary Conditions
            std::cout<<"Process Boundary Conditions"<<std::endl;
            BoundaryConditionsContainer<elementDim,spaceDim,probDim> bcc;
            std::vector<bool> areNeumannBoundaryConditions(probDim, true);
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
            Vec initial_condition = PetscTools::CreateVec(init_conds);
            std::cout<<"Here"<<std::endl;
            std::vector<AbstractInhomogenousOdeSystemForCoupledPdeSystem*> odeSystem;
            std::vector<boost::shared_ptr<AbstractIvpOdeSolver> > solverSystem;
            for (unsigned i=0; i<p_mesh->GetNumNodes(); i++){
                // number of ode system objects must match the number of nodes, i.e the individual odes may be multi-dimensional
                odeSystem.push_back(new AbstractInhomogenousChemicalOdeSystemForCoupledPdeSystem(chemicalReactionSystem));
                boost::shared_ptr<AbstractIvpOdeSolver> p_solver(new EulerIvpOdeSolver);
                solverSystem.push_back(p_solver);//std::dynamic_pointer_cast<AbstractIvpOdeSolver>(p_solver));
            }

            // pde system
            InhomogenousParabolicPdeForCoupledOdeSystemTemplated<elementDim, spaceDim, probDim> pde(chemical_structure->rGetPtrChemicalDomain());


            std::cout<<"Solver"<<std::endl;
            // solver
            InhomogenousCoupledPdeOdeSolverTemplated<elementDim,spaceDim,probDim> solver(p_mesh, &pde, &bcc,odeSystem,solverSystem);

            // solver properties
        
            solver.SetTimes(0, 20);
            solver.SetTimeStep(1e-2);
            solver.SetSamplingTimeStep(1e-2);
            solver.SetOutputDirectory("TestReactionSystemWithoutCellsUsingInhomogenousSolver_single_larger_only_reaction_within");
            solver.SetInitialCondition(initial_condition);
            // solve
            std::cout<<"Solve"<<std::endl;
            solver.SolveAndWriteResultsToFile();
            std::cout<<"Clean"<<std::endl;

        }

    }

    void TestReactionSystemInhibitedDiffusionWithoutCellsUsingInhomogenousSolver()
    {

        if(control.ReactionSystemInhibitedDiffusionWithoutCellsUsingInhomogenousSolver)
        {
            std::cout<<"Mass action kinetics of Schnackenberg using the ChemicalDomain Class"<<std::endl;
            std::cout<<"-----------------------------"<<std::endl;
            
            ChemicalStructuresForTests* chemical_structure = new ChemicalStructuresForTests();
            chemical_structure->SetUpPdeChemicalReactionSystem();
            chemical_structure->SetUpChemicalDomainField();
            AbstractReactionSystem* chemicalReactionSystem = chemical_structure -> rGetPtrPdeChemicalReactionSystem();

            // form ode system
            std::cout<<"SchnackenbergCoupledPdeOdeSystem  -As AbstractChemicalODeSystem"<<std::endl;
            // system properties
            const unsigned probDim =2;
            const unsigned spaceDim=2;
            const unsigned elementDim=2;
            double MeshStepSize = 1.0;
            std::vector<unsigned> MeshDimensions = {100,100,0};
            std::vector<double> initValues = {2.0, 0.75};
            std::vector<double> bcValues = {0.0, 0.0};
            std::vector<double> diffusionRates = {0.05, 0.05};

            // mesh
            TetrahedralMesh<elementDim,spaceDim>* p_mesh = new TetrahedralMesh<elementDim,spaceDim>();

            switch (spaceDim)
            {
                case 1:
                    p_mesh->ConstructRegularSlabMesh(MeshStepSize, MeshDimensions[0]);
                    break;
                case 2:
                    p_mesh->ConstructRegularSlabMesh(MeshStepSize, MeshDimensions[0], MeshDimensions[1]);
                    break;
                case 3:
                    p_mesh->ConstructRegularSlabMesh(MeshStepSize, MeshDimensions[0], MeshDimensions[1], MeshDimensions[2]);
                    break;
                default:
                    NEVER_REACHED;
            }
            // Process Boundary Conditions
            std::cout<<"Process Boundary Conditions"<<std::endl;
            BoundaryConditionsContainer<elementDim,spaceDim,probDim> bcc;
            std::vector<bool> areNeumannBoundaryConditions(probDim, true);
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
            Vec initial_condition = PetscTools::CreateVec(init_conds);

            std::vector<AbstractInhomogenousOdeSystemForCoupledPdeSystem*> odeSystem;
            std::vector<boost::shared_ptr<AbstractIvpOdeSolver> > solverSystem;
            for (unsigned i=0; i<p_mesh->GetNumNodes(); i++){
                // number of ode system objects must match the number of nodes, i.e the individual odes may be multi-dimensional
                odeSystem.push_back(new AbstractInhomogenousChemicalOdeSystemForCoupledPdeSystem(chemicalReactionSystem));
                boost::shared_ptr<AbstractIvpOdeSolver> p_solver(new EulerIvpOdeSolver);
                solverSystem.push_back(p_solver);//std::dynamic_pointer_cast<AbstractIvpOdeSolver>(p_solver));
            }

            // pde system
            InhomogenousParabolicPdeForCoupledOdeSystemInhibitedDiffusionTemplated<elementDim, spaceDim, probDim> pde(chemical_structure->rGetPtrChemicalDomain());


            std::cout<<"Solver"<<std::endl;
            // solver
            InhomogenousCoupledPdeOdeSolverTemplated<elementDim,spaceDim,probDim> solver(p_mesh, &pde, &bcc,odeSystem,solverSystem);

            // solver properties
        
            solver.SetTimes(0, 20);
            solver.SetTimeStep(1e-2);
            solver.SetSamplingTimeStep(1e-2);
            solver.SetOutputDirectory("TestReactionSystemInhibitedDiffusionWithoutCellsUsingInhomogenousSolver_single_larger_only_reaction_within");
            solver.SetInitialCondition(initial_condition);
            // solve
            std::cout<<"Solve"<<std::endl;
            solver.SolveAndWriteResultsToFile();
            std::cout<<"Clean"<<std::endl;

        }

    }


    void TestReactionSystemWithCells()
    {
    
        if(control.ReactionSystemWithCells)
        {
            std::cout<<"=============================================="<<std::endl;
            std::cout<<"TestReactionSystemWithCells()"<<std::endl;
            std::cout<<"=============================================="<<std::endl;

            // make mesh of cell with associated mesh based cell population
            HoneycombMeshGenerator generator(2, 2);    // Parameters are: cells across, cells up
            MutableMesh<2,2>* p_mesh = generator.GetMesh();

            std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

            std::vector<CellPtr> cells;
            // assume cell at each node in cell layer mesh
            for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
            {
                // provide each cell with a transport cell property and membrane property, cell cycle, wild type states
                ChemicalStructuresForTests* chemical_structure = new ChemicalStructuresForTests();
                
                if(i==0)
                {
                    CellPtr p_cell = chemical_structure->SetUpChemicalCellObjectA();
                    cells.push_back(p_cell);
                }
                //else if(i==8)
                //{
                //    CellPtr p_cell = chemical_structure->SetUpChemicalCellObjectB();
                //    cells.push_back(p_cell);
                //}
                else
                {
                    //CellPtr p_cell = chemical_structure->SetUpNullCellObject();  
                    CellPtr p_cell = chemical_structure->SetUpChemicalCellObjectA();

                    cells.push_back(p_cell);
                }
                
            }
        
            MeshBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices);
            
            

            // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
            ChastePoint<2> lower(-4.0, -4.0);
            ChastePoint<2> upper(7.0, 7.0);
            MAKE_PTR_ARGS(ChasteCuboid<2>, p_cuboid, (lower, upper));

            // set up chemical domain field
            ChemicalStructuresForTests* chemical_structure = new ChemicalStructuresForTests();

            chemical_structure -> SetUpChemicalDomainFieldForCellCoupling();
            
            ChemicalDomainFieldForCellCoupling<2,2,2>*& p_Pde_field = chemical_structure -> rGetPtrChemicalDomainFieldForCellCoupling();

            std::cout<<"initial conditions"<<std::endl;
            std::vector<double> init_conditions = p_Pde_field -> GetInitialNodeConditions();
            for(unsigned i=0; i<init_conditions.size(); i++)
            {
                std::cout<<init_conditions[i]<<std::endl;
            }

            
            boost::shared_ptr<ParabolicBoxDomainPdeSystemModifier<2,2,2>> p_pde_modifier(new ParabolicBoxDomainPdeSystemModifier<2,2,2>(p_Pde_field, p_cuboid));
            
            boost::shared_ptr<ChemicalTrackingModifier<2,2>> p_chemical_tracking_modifier(new ChemicalTrackingModifier<2,2>()); //= chemical_structure -> rGetPtrChemicalTrackingModifier();
            
            // writers
            cell_population.SetWriteVtkAsPoints(false);
            cell_population.AddPopulationWriter<VoronoiDataWriter>();

            cell_population.AddCellWriter<CellIdWriter>();
            cell_population.AddCellWriter<CellAgesWriter>();
            cell_population.AddCellWriter<CellLocationIndexWriter>();

            std::vector<std::string> chemicalCellASpeciesNames = chemical_structure->GetChemicalCellASpeciesNames();

            for(unsigned i=0; i<chemicalCellASpeciesNames.size(); i++)
            {
                boost::shared_ptr<CellDataItemWriter<2,2>> dataWriter(new CellDataItemWriter<2,2>(chemicalCellASpeciesNames[i]));
                cell_population.AddCellWriter(dataWriter);
            }
                
        
            OffLatticeSimulation<2> simulator(cell_population);

            simulator.AddSimulationModifier(p_pde_modifier);

            simulator.AddSimulationModifier(p_chemical_tracking_modifier);
        
            simulator.SetOutputDirectory("TestReactionSystemWithCells_OscillatingCase_meshWriters");
            simulator.SetEndTime(4.0);

            MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
            p_linear_force->SetCutOffLength(1.5);
            simulator.AddForce(p_linear_force);

            std::cout<<"=============================================="<<std::endl;
            std::cout<<"OffLatticeSimulation -> AbstractCellBasedSumulation :: Solve()"<<std::endl;
            std::cout<<"=============================================="<<std::endl;
            simulator.Solve();
        }

    }

    void TestReactionSystemWithNodeBasedCells()
    {
        if(control.ReactionSystemWithNodeBasedCells)
        {
            std::cout<<"=============================================="<<std::endl;
            std::cout<<"TestReactionSystemWithCellsNodes()"<<std::endl;
            std::cout<<"=============================================="<<std::endl;

            // make mesh of cell with associated mesh based cell population
            HoneycombMeshGenerator generator(2, 2);    // Parameters are: cells across, cells up
            MutableMesh<2,2>* p_mesh = generator.GetMesh();

            NodesOnlyMesh<2> mesh;
            mesh.ConstructNodesWithoutMesh(*p_mesh, 1.5);


            std::vector<CellPtr> cells;
            // assume cell at each node in cell layer mesh
            for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
            {
                // provide each cell with a transport cell property and membrane property, cell cycle, wild type states
                ChemicalStructuresForTests* chemical_structure = new ChemicalStructuresForTests();
                
                if(i==0)
                {
                    CellPtr p_cell = chemical_structure->SetUpChemicalCellObjectA();
                    cells.push_back(p_cell);
                }
                //else if(i==8)
                //{
                //    CellPtr p_cell = chemical_structure->SetUpChemicalCellObjectB();
                //    cells.push_back(p_cell);
                //}
                else
                {
                    //CellPtr p_cell = chemical_structure->SetUpNullCellObject();  
                    CellPtr p_cell = chemical_structure->SetUpChemicalCellObjectA();

                    cells.push_back(p_cell);
                }
                
            }

            NodeBasedCellPopulation<2> cell_population(mesh, cells);
        
            // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
            ChastePoint<2> lower(-4.0, -4.0);
            ChastePoint<2> upper(7.0, 7.0);
            MAKE_PTR_ARGS(ChasteCuboid<2>, p_cuboid, (lower, upper));

            // set up chemical domain field
            ChemicalStructuresForTests* chemical_structure = new ChemicalStructuresForTests();

            chemical_structure -> SetUpChemicalDomainFieldForCellCoupling();
            
            ChemicalDomainFieldForCellCoupling<2,2,2>*& p_Pde_field = chemical_structure -> rGetPtrChemicalDomainFieldForCellCoupling();

            std::cout<<"initial conditions"<<std::endl;
            std::vector<double> init_conditions = p_Pde_field -> GetInitialNodeConditions();
            for(unsigned i=0; i<init_conditions.size(); i++)
            {
                std::cout<<init_conditions[i]<<std::endl;
            }

            
            boost::shared_ptr<ParabolicBoxDomainPdeSystemModifier<2,2,2>> p_pde_modifier(new ParabolicBoxDomainPdeSystemModifier<2,2,2>(p_Pde_field, p_cuboid));
            
            boost::shared_ptr<ChemicalTrackingModifier<2,2>> p_chemical_tracking_modifier(new ChemicalTrackingModifier<2,2>()); //= chemical_structure -> rGetPtrChemicalTrackingModifier();
            
        
        
            OffLatticeSimulation<2> simulator(cell_population);

            simulator.AddSimulationModifier(p_pde_modifier);

            simulator.AddSimulationModifier(p_chemical_tracking_modifier);
        
            simulator.SetOutputDirectory("TestReactionSystemWithCells_OscillatingCase_meshWriters_nodeBased");
            simulator.SetEndTime(4.0);

            MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
            p_linear_force->SetCutOffLength(1.5);
            simulator.AddForce(p_linear_force);

            std::cout<<"=============================================="<<std::endl;
            std::cout<<"OffLatticeSimulation -> AbstractCellBasedSumulation :: Solve()"<<std::endl;
            std::cout<<"=============================================="<<std::endl;
            simulator.Solve();

        }

    }


    void TestCellLayerRead()
    {
        std::string dataFileRoot = "/home/chaste/projects/ChemChaste/DataInput/Data/MulticellCase/DomainField/";
        std::string cellFileRoot = "/home/chaste/projects/ChemChaste/DataInput/Data/MulticellCase/Cell/";
        std::string cellLabelFilename = "CellLayerTopology.csv";
        std::string cellKeyFilename = "CellLayerKey.csv";
        std::string domainFilename = "Domain.csv";
        std::string domainKeyFilename = "DomainKey.csv";
        std::string odeLabelFilename = "NodeSelector.csv";
        std::string odeKeyFilename = "OdeReactionFileKey.csv";
        std::string diffusionFilename = "DiffusionDatabaseFile.csv";
        std::string initialConditionsFilename = "InitialConditionFile.csv";
        std::string boundaryConditionsFilename = "BoundaryConditionFile.csv";
        
        // System properties
        const unsigned probDim =4; // need to set manually to the number of diffusive variables for the pde solver to solve
        const unsigned spaceDim=2;
        const unsigned elementDim=2;

        //TetrahedralMesh<elementDim,spaceDim>* p_field = new TetrahedralMesh<elementDim,spaceDim>();
        // generate domain
        // run the domain field set up and parse files
        ChemicalDomainFieldForCellCoupling<elementDim,spaceDim,probDim>* p_Pde_field = new ChemicalDomainFieldForCellCoupling<elementDim,spaceDim,probDim>(dataFileRoot,cellFileRoot+cellLabelFilename,cellFileRoot+cellKeyFilename,dataFileRoot+domainFilename, dataFileRoot+domainKeyFilename, dataFileRoot+odeLabelFilename, dataFileRoot+odeKeyFilename, dataFileRoot+diffusionFilename, dataFileRoot+initialConditionsFilename, dataFileRoot+boundaryConditionsFilename);

        //TetrahedralMesh<elementDim,spaceDim>* p_cell_mesh = p_Pde_field->rGetCellMesh();
    std::cout<<"here"<<std::endl;
        TetrahedralMesh<elementDim,spaceDim>* p_cell_mesh = p_Pde_field->rGetCellMesh();
        //std::vector<unsigned> location_indices = generator.GetCellLocationIndices();
  std::cout<<"here"<<std::endl;

        NodesOnlyMesh<spaceDim> mesh;
        mesh.ConstructNodesWithoutMesh(*p_cell_mesh, 1.5);

        std::cout<<"here"<<std::endl;
        std::vector<CellPtr> cells;
        // assume cell at each node in cell layer mesh
        std::cout<<"here"<<std::endl;
        std::string cell_label;
        std::string cell_key;
        std::string given_cell_root;
        for (unsigned i=0; i<p_cell_mesh->GetNumNodes(); i++)
        {
            cell_label = p_Pde_field->GetCellLabelByIndex(i);
            cell_key = p_Pde_field->ReturnCellKeyFromCellLabel(cell_label);
            std::cout<<"Cell i: "<<i<<" label: "<<cell_label<<std::endl;
            std::cout<<"Cell i: "<<i<<" label: "<<cell_key<<std::endl;
            given_cell_root = cellFileRoot+cell_key+"/";

            ChemicalCellFromFile* p_cell_reader = new ChemicalCellFromFile(
                                given_cell_root+"SpeciesThreshold.csv", 
                                given_cell_root+"Srn.txt",
                                given_cell_root+"InitialCellConcentrations.csv",
                                given_cell_root+"TransportReactions.txt",
                                given_cell_root+"MembraneReactions.txt"
                                );

            cells.push_back(p_cell_reader -> GetCellPtr());
        }   


        NodeBasedCellPopulation<2> cell_population(mesh, cells);
        
        // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
        ChastePoint<2> lower(-4.0, -4.0);
        ChastePoint<2> upper(7.0, 7.0);
        MAKE_PTR_ARGS(ChasteCuboid<2>, p_cuboid, (lower, upper));

        std::cout<<"initial conditions"<<std::endl;
        std::vector<double> init_conditions = p_Pde_field -> GetInitialNodeConditions();
        for(unsigned i=0; i<init_conditions.size(); i++)
        {
            std::cout<<init_conditions[i]<<std::endl;
        }

        
        boost::shared_ptr<ParabolicBoxDomainPdeSystemModifier<elementDim,spaceDim,probDim>> p_pde_modifier(new ParabolicBoxDomainPdeSystemModifier<elementDim,spaceDim,probDim>(p_Pde_field, p_cuboid));
        
        boost::shared_ptr<ChemicalTrackingModifier<2,2>> p_chemical_tracking_modifier(new ChemicalTrackingModifier<2,2>()); //= chemical_structure -> rGetPtrChemicalTrackingModifier();
        
        OffLatticeSimulation<2> simulator(cell_population);

        simulator.AddSimulationModifier(p_pde_modifier);

        simulator.AddSimulationModifier(p_chemical_tracking_modifier);
    
        simulator.SetOutputDirectory("TestNodeBasedCellsFromFile");
        simulator.SetEndTime(10.0);

        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        p_linear_force->SetCutOffLength(1.5);
        simulator.AddForce(p_linear_force);

        std::cout<<"=============================================="<<std::endl;
        std::cout<<"OffLatticeSimulation -> AbstractCellBasedSumulation :: Solve()"<<std::endl;
        std::cout<<"=============================================="<<std::endl;
        simulator.Solve();


    }
};

#endif

