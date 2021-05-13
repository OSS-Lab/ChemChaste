#ifndef TESTCHASTEPDELEAKAGE_HPP_
#define TESTCHASTEPDELEAKAGE_HPP_

#include <cxxtest/TestSuite.h>
#include "UblasIncludes.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "RandomNumberGenerator.hpp"

#include "TetrahedralMesh.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "LinearParabolicPdeSystemWithCoupledOdeSystemSolver.hpp"
#include "SimpleLinearParabolicSolver.hpp"

#include "LinearParabolicHeatEquationPde.hpp"
#include "SimpleHeatEquation.hpp"


struct ParameterStruct {
    // standard values
    double MeshStepSize = 1.0;
    std::vector<unsigned> MeshDimensions = {100,10};
    std::vector<double> initValuesHigh = {1.0};
    std::vector<double> initValuesLow = {0.0};
    std::vector<double> bcValuesDir = {0.0};
    std::vector<double> bcValuesNeu = {0.0};
    std::vector<double> diffusionRates = {100.0};

    // solver properties
    double startTime= 0.0;
    double endTime = 10.0;
    double timestep = 0.01;
    double samplingTimestep = 0.1;
    std::string outputDirName = "ChasteLeakage/"; //+ specific test solver

} params;


class TestChastePdeLeakage : public AbstractCellBasedTestSuite
{
public:


    void TestSimpleHeatDiffusionWithoutSource()
    {

        // mesh
        TetrahedralMesh<2,2>* p_mesh = new TetrahedralMesh<2,2>();

        p_mesh->ConstructRegularSlabMesh(params.MeshStepSize, params.MeshDimensions[0], params.MeshDimensions[1]);

        
        // Process Boundary Conditions
        BoundaryConditionsContainer<2,2,1> bcc;
    
        std::vector<ConstBoundaryCondition<2>*> vectorConstBCs;
        
        for(unsigned pdeDim=0; pdeDim<1; pdeDim++){
            vectorConstBCs.push_back(new ConstBoundaryCondition<2>(params.bcValuesNeu[pdeDim]));
        }
        
        for(unsigned pdeDim=0; pdeDim<1; pdeDim++)
        {
            for (TetrahedralMesh<2,2>::BoundaryElementIterator boundary_iter = p_mesh->GetBoundaryElementIteratorBegin();
            boundary_iter != p_mesh->GetBoundaryElementIteratorEnd();
            boundary_iter++)
            {
                bcc.AddNeumannBoundaryCondition(*boundary_iter, vectorConstBCs[pdeDim], pdeDim);
            }
        }


        // initial conditions
        std::vector<double> init_conds(1*p_mesh->GetNumNodes(),0.0);
        unsigned columnNum = 0;
        unsigned rowNum = 0;
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {   // set as being a random perturbation about the boundary values
            
            columnNum = 0;
            rowNum = 0;

            while(i >= rowNum*(params.MeshDimensions[0]+1))
            {
                rowNum = rowNum + 1;
            
            }
            
            columnNum = i - (rowNum-1)*(params.MeshDimensions[0]+1);
            if((columnNum==3 || columnNum==4 || columnNum==5 || columnNum==6) && (rowNum ==3 || rowNum ==4 || rowNum ==5 || rowNum ==6 ))
            {
                for(unsigned pdeDim=0; pdeDim<1; pdeDim++)
                {   // serialised for nodes
                    init_conds[1*i + pdeDim] = fabs(params.initValuesHigh[pdeDim]);// + RandomNumberGenerator::Instance()->ranf());
                }
            }
            else
            {
                for(unsigned pdeDim=0; pdeDim<1; pdeDim++)
                {   // serialised for nodes
                    init_conds[1*i + pdeDim] = fabs(params.initValuesLow[pdeDim]);// + RandomNumberGenerator::Instance()->ranf());
                }
            }

        }
        // PETSc Vec
        Vec initial_condition = PetscTools::CreateVec(init_conds);

        SimpleHeatEquation<2> pde;
        // solver
        SimpleLinearParabolicSolver<2,2> solver(p_mesh, &pde, &bcc);


        solver.SetTimes(params.startTime, params.endTime);

        solver.SetTimeStep(params.timestep);
        solver.SetInitialCondition(initial_condition);
   
        solver.SetOutputDirectoryAndPrefix(params.outputDirName+"ChasteLeakageSimpleLinearParabolicSolver","results");
        
        solver.SetOutputToVtk(true);
        solver.SetPrintingTimestepMultiple(10);
        // solve
        Vec solution = solver.Solve();
        ReplicatableVector solution_repl(solution);

    }

    void TestHeatDiffusionWithoutSource()
    {
std::cout<<"Here0"<<std::endl;
        // mesh
        TetrahedralMesh<2,2>* p_mesh = new TetrahedralMesh<2,2>();

        p_mesh->ConstructRegularSlabMesh(params.MeshStepSize, params.MeshDimensions[0], params.MeshDimensions[1]);
    
        
        // Process Boundary Conditions
        BoundaryConditionsContainer<2,2,1> bcc;
    
        std::vector<ConstBoundaryCondition<2>*> vectorConstBCs;
        
        vectorConstBCs.push_back(new ConstBoundaryCondition<2>(params.bcValuesDir[0]));
        vectorConstBCs.push_back(new ConstBoundaryCondition<2>(params.bcValuesNeu[0]));
        
     std::cout<<"Here1"<<std::endl;   

        for (TetrahedralMesh<2,2>::BoundaryNodeIterator node_iter = p_mesh->GetBoundaryNodeIteratorBegin();
        node_iter != p_mesh->GetBoundaryNodeIteratorEnd();
        ++node_iter)
        {

            bcc.AddDirichletBoundaryCondition(*node_iter, vectorConstBCs[0]);
        }

        for (TetrahedralMesh<2,2>::BoundaryElementIterator boundary_iter = p_mesh->GetBoundaryElementIteratorBegin();
        boundary_iter != p_mesh->GetBoundaryElementIteratorEnd();
        boundary_iter++)
        {
            bcc.AddNeumannBoundaryCondition(*boundary_iter, vectorConstBCs[1]);
        }
    

std::cout<<"Here2"<<std::endl;
        // initial conditions
        std::vector<double> init_conds(1*p_mesh->GetNumNodes(),0.0);
        unsigned columnNum = 0;
        unsigned rowNum = 0;
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {   // set as being a random perturbation about the boundary values
            
            columnNum = 0;
            rowNum = 0;

            while(i >= rowNum*(params.MeshDimensions[0]+1))
            {
                rowNum = rowNum + 1;
            
            }
            
            columnNum = i - (rowNum-1)*(params.MeshDimensions[0]+1);
            if((columnNum==3 || columnNum==4 || columnNum==5 || columnNum==6) && (rowNum ==3 || rowNum ==4 || rowNum ==5 || rowNum ==6 ))
            {
                for(unsigned pdeDim=0; pdeDim<1; pdeDim++)
                {   // serialised for nodes
                    init_conds[1*i + pdeDim] = fabs(params.initValuesHigh[pdeDim]);// + RandomNumberGenerator::Instance()->ranf());
                }
            }
            else
            {
                for(unsigned pdeDim=0; pdeDim<1; pdeDim++)
                {   // serialised for nodes
                    init_conds[1*i + pdeDim] = fabs(params.initValuesLow[pdeDim]);// + RandomNumberGenerator::Instance()->ranf());
                }
            }

        }
        // PETSc Vec
        Vec initial_condition = PetscTools::CreateVec(init_conds);
std::cout<<"Here3"<<std::endl;

        LinearParabolicHeatEquationPde<2, 2, 1> pde(params.diffusionRates);
        std::cout<<"Here4"<<std::endl;
        // solver
        LinearParabolicPdeSystemWithCoupledOdeSystemSolver<2,2,1> solver(p_mesh, &pde, &bcc);
std::cout<<"Here5"<<std::endl;

        solver.SetTimes(params.startTime, params.endTime);

        solver.SetTimeStep(params.timestep);
        if(params.timestep>params.samplingTimestep)
        {
            solver.SetSamplingTimeStep(params.timestep);   
        }
        else
        {
            solver.SetSamplingTimeStep(params.samplingTimestep);
        }
   std::cout<<"Here6"<<std::endl;
        solver.SetOutputDirectory(params.outputDirName+"LinearParabolicPdeSystemWithCoupledOdeSystemSolver/");
        solver.SetInitialCondition(initial_condition);
        // solve
        solver.SolveAndWriteResultsToFile();
        

    }



};

#endif