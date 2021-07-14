#ifndef TESTTOTALMASS_HPP_
#define TESTTOTALMASS_HPP_


#include "ChemChasteFeAssemblerCommon.hpp"
#include "ChemChasteVolumeAssembler.hpp"
#include "ChemChasteSurfaceAssembler.hpp"

#include <cxxtest/TestSuite.h>
#include "UblasIncludes.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "RandomNumberGenerator.hpp"


#include "TetrahedralMesh.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "InhomogenousCoupledPdeOdeSolver_templated.hpp"
#include "LinearParabolicPdeSystemWithCoupledOdeSystemSolver.hpp"
#include "SimpleLinearParabolicSolver.hpp"

#include "InhomogenousHeatEquationPde.hpp"
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
    std::vector<double> diffusionRates10 = {10.0};
    std::vector<double> diffusionRates100 = {1.0};

    // solver properties
    double startTime= 0.0;
    double endTime = 1.0;
    double timestep = 0.1;
    double samplingTimestep = 0.1;
    std::string outputDirName = "TotalMass"; //+ specific test solver

} params;

/*
const c_vector<double, SPACE_DIM> node_loc = rSurfaceElement.GetNode(i)->rGetLocation();
            x.rGetLocation() += phi(i)*node_loc;

            // Allow the concrete version of the assembler to interpolate any desired quantities
            this->IncrementInterpolatedQuantities(phi(i), rSurfaceElement.GetNode(i));

*/


class TestTotalMass : public AbstractCellBasedTestSuite
{
public:

    void TestTotalMassHeatDiffusionWithoutSource()
    {

        ofstream myfile("/home/chaste/testoutput/totalMassExample.txt");
        if (myfile.is_open())
        {
            myfile << "Writing this to a file.\n";
            myfile.close();
        }else{
            std::cout<<"not open"<<std::endl;
        }

      

        // mesh
        TetrahedralMesh<2,2>* p_mesh = new TetrahedralMesh<2,2>();

        p_mesh->ConstructRegularSlabMesh(params.MeshStepSize, params.MeshDimensions[0], params.MeshDimensions[1]);
    
        
        // Process Boundary Conditions
        BoundaryConditionsContainer<2,2,1> bcc;
    
        std::vector<ConstBoundaryCondition<2>*> vectorConstBCs;
        
        vectorConstBCs.push_back(new ConstBoundaryCondition<2>(params.bcValuesDir[0]));
        vectorConstBCs.push_back(new ConstBoundaryCondition<2>(params.bcValuesNeu[0]));
        

        for (TetrahedralMesh<2,2>::BoundaryElementIterator boundary_iter = p_mesh->GetBoundaryElementIteratorBegin();
        boundary_iter != p_mesh->GetBoundaryElementIteratorEnd();
        boundary_iter++)
        {
            bcc.AddNeumannBoundaryCondition(*boundary_iter, vectorConstBCs[1]);
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
            if((columnNum==3 || columnNum==4 || columnNum==5 || columnNum==6))
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


        InhomogenousHeatEquationPde<2, 2, 1> pde(params.diffusionRates);

        // solver
        InhomogenousCoupledPdeOdeSolverTemplated<2,2,1> solver(p_mesh, &pde, &bcc);

        solver.SetTimes(params.startTime, params.endTime);

        solver.SetTimeStep(params.timestep/10);
        if(params.timestep>params.samplingTimestep)
        {
            solver.SetSamplingTimeStep(params.timestep);   
        }
        else
        {
            solver.SetSamplingTimeStep(params.samplingTimestep);
        }

        solver.SetOutputDirectory(params.outputDirName+"LinearParabolicPdeSystemWithCoupledOdeSystemSolverScaleBy1dt001");
        solver.SetInitialCondition(initial_condition);
        // solve
        solver.SolveAndWriteResultsToFile();
        

    }


    void TestTotalMassHeatDiffusionWithoutSourceScaled1000dt0001()
    {


        // mesh
        TetrahedralMesh<2,2>* p_mesh = new TetrahedralMesh<2,2>();

        p_mesh->ConstructRegularSlabMesh(params.MeshStepSize, params.MeshDimensions[0], params.MeshDimensions[1]);
        p_mesh->Scale(0.001,0.001);
        
        // Process Boundary Conditions
        BoundaryConditionsContainer<2,2,1> bcc;
    
        std::vector<ConstBoundaryCondition<2>*> vectorConstBCs;
        
        vectorConstBCs.push_back(new ConstBoundaryCondition<2>(params.bcValuesDir[0]));
        vectorConstBCs.push_back(new ConstBoundaryCondition<2>(params.bcValuesNeu[0]));
        
     std::cout<<"Here1"<<std::endl;   
/*
        for (TetrahedralMesh<2,2>::BoundaryNodeIterator node_iter = p_mesh->GetBoundaryNodeIteratorBegin();
        node_iter != p_mesh->GetBoundaryNodeIteratorEnd();
        ++node_iter)
        {

            bcc.AddDirichletBoundaryCondition(*node_iter, vectorConstBCs[0]);
        }
*/
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
            if((columnNum==3 || columnNum==4 || columnNum==5 || columnNum==6))
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


        InhomogenousHeatEquationPde<2, 2, 1> pde(params.diffusionRates10);
        std::cout<<"Here4"<<std::endl;
        // solver
        InhomogenousCoupledPdeOdeSolverTemplated<2,2,1> solver(p_mesh, &pde, &bcc);
std::cout<<"Here5"<<std::endl;

        solver.SetTimes(params.startTime, params.endTime);

        solver.SetTimeStep(params.timestep/100);
        if(params.timestep>params.samplingTimestep)
        {
            solver.SetSamplingTimeStep(params.timestep);   
        }
        else
        {
            solver.SetSamplingTimeStep(params.samplingTimestep);
        }
   std::cout<<"Here6"<<std::endl;
        solver.SetOutputDirectory(params.outputDirName+"LinearParabolicPdeSystemWithCoupledOdeSystemSolverScaledBy1dt0001");
        solver.SetInitialCondition(initial_condition);
        // solve
        solver.SolveAndWriteResultsToFile();
        

    }

    void TestTotalMassHeatDiffusionWithoutSourceScaled1000dt00001()
    {


        // mesh
        TetrahedralMesh<2,2>* p_mesh = new TetrahedralMesh<2,2>();

        p_mesh->ConstructRegularSlabMesh(params.MeshStepSize, params.MeshDimensions[0], params.MeshDimensions[1]);
        p_mesh->Scale(0.001,0.001);
        
        // Process Boundary Conditions
        BoundaryConditionsContainer<2,2,1> bcc;
    
        std::vector<ConstBoundaryCondition<2>*> vectorConstBCs;
        
        vectorConstBCs.push_back(new ConstBoundaryCondition<2>(params.bcValuesDir[0]));
        vectorConstBCs.push_back(new ConstBoundaryCondition<2>(params.bcValuesNeu[0]));
        
     std::cout<<"Here1"<<std::endl;   
/*
        for (TetrahedralMesh<2,2>::BoundaryNodeIterator node_iter = p_mesh->GetBoundaryNodeIteratorBegin();
        node_iter != p_mesh->GetBoundaryNodeIteratorEnd();
        ++node_iter)
        {

            bcc.AddDirichletBoundaryCondition(*node_iter, vectorConstBCs[0]);
        }
*/
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
            if((columnNum==3 || columnNum==4 || columnNum==5 || columnNum==6))
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


        InhomogenousHeatEquationPde<2, 2, 1> pde(params.diffusionRates10);
        std::cout<<"Here4"<<std::endl;
        // solver
        InhomogenousCoupledPdeOdeSolverTemplated<2,2,1> solver(p_mesh, &pde, &bcc);
std::cout<<"Here5"<<std::endl;

        solver.SetTimes(params.startTime, params.endTime);

        solver.SetTimeStep(params.timestep/10);
        if(params.timestep>params.samplingTimestep)
        {
            solver.SetSamplingTimeStep(params.timestep);   
        }
        else
        {
            solver.SetSamplingTimeStep(params.samplingTimestep);
        }
   std::cout<<"Here6"<<std::endl;
        solver.SetOutputDirectory(params.outputDirName+"LinearParabolicPdeSystemWithCoupledOdeSystemSolverScaledBy10dt00001");
        solver.SetInitialCondition(initial_condition);
        // solve
        solver.SolveAndWriteResultsToFile();
        

    }

    void TestTotalMassHeatDiffusionWithoutSourceScaled10dt0000001()
    {

      /*

        // mesh
        TetrahedralMesh<2,2>* p_mesh = new TetrahedralMesh<2,2>();

        p_mesh->ConstructRegularSlabMesh(params.MeshStepSize, params.MeshDimensions[0], params.MeshDimensions[1]);
        p_mesh->Scale(0.1,0.1);
        
        // Process Boundary Conditions
        BoundaryConditionsContainer<2,2,1> bcc;
    
        std::vector<ConstBoundaryCondition<2>*> vectorConstBCs;
        
        vectorConstBCs.push_back(new ConstBoundaryCondition<2>(params.bcValuesDir[0]));
        vectorConstBCs.push_back(new ConstBoundaryCondition<2>(params.bcValuesNeu[0]));
        


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
            if((columnNum==3 || columnNum==4 || columnNum==5 || columnNum==6))
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


        InhomogenousHeatEquationPde<2, 2, 1> pde(params.diffusionRates10);
        std::cout<<"Here4"<<std::endl;
        // solver
        InhomogenousCoupledPdeOdeSolverTemplated<2,2,1> solver(p_mesh, &pde, &bcc);
std::cout<<"Here5"<<std::endl;

        solver.SetTimes(params.startTime, params.endTime);

        solver.SetTimeStep(params.timestep/10000);
        if(params.timestep>params.samplingTimestep)
        {
            solver.SetSamplingTimeStep(params.timestep);   
        }
        else
        {
            solver.SetSamplingTimeStep(params.samplingTimestep);
        }
   std::cout<<"Here6"<<std::endl;
        solver.SetOutputDirectory(params.outputDirName+"LinearParabolicPdeSystemWithCoupledOdeSystemSolverScaledBy10dt0000001");
        solver.SetInitialCondition(initial_condition);
        // solve
        solver.SolveAndWriteResultsToFile();
        
*/
    }



    void TestTotalMassHeatDiffusionWithoutSourceScaled1000dt001()
    {
/*

        // mesh
        TetrahedralMesh<2,2>* p_mesh = new TetrahedralMesh<2,2>();

        p_mesh->ConstructRegularSlabMesh(params.MeshStepSize, params.MeshDimensions[0], params.MeshDimensions[1]);
        p_mesh->Scale(0.001,0.001);
        
        // Process Boundary Conditions
        BoundaryConditionsContainer<2,2,1> bcc;
    
        std::vector<ConstBoundaryCondition<2>*> vectorConstBCs;
        
        vectorConstBCs.push_back(new ConstBoundaryCondition<2>(params.bcValuesDir[0]));
        vectorConstBCs.push_back(new ConstBoundaryCondition<2>(params.bcValuesNeu[0]));

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
            if((columnNum==3 || columnNum==4 || columnNum==5 || columnNum==6))
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


        InhomogenousHeatEquationPde<2, 2, 1> pde(params.diffusionRates10);
        std::cout<<"Here4"<<std::endl;
        // solver
        InhomogenousCoupledPdeOdeSolverTemplated<2,2,1> solver(p_mesh, &pde, &bcc);
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
        solver.SetOutputDirectory(params.outputDirName+"LinearParabolicPdeSystemWithCoupledOdeSystemSolverScaledBy1000dt001");
        solver.SetInitialCondition(initial_condition);
        // solve
        solver.SolveAndWriteResultsToFile();
        */

    }

    void TestTotalMassHeatDiffusionWithoutSourceScaled1000dt0001()
    {


        // mesh
        TetrahedralMesh<2,2>* p_mesh = new TetrahedralMesh<2,2>();

        p_mesh->ConstructRegularSlabMesh(params.MeshStepSize, params.MeshDimensions[0], params.MeshDimensions[1]);
        p_mesh->Scale(0.001,0.001);
        
        // Process Boundary Conditions
        BoundaryConditionsContainer<2,2,1> bcc;
    
        std::vector<ConstBoundaryCondition<2>*> vectorConstBCs;
        
        vectorConstBCs.push_back(new ConstBoundaryCondition<2>(params.bcValuesDir[0]));
        vectorConstBCs.push_back(new ConstBoundaryCondition<2>(params.bcValuesNeu[0]));
        
     std::cout<<"Here1"<<std::endl;   
/*
        for (TetrahedralMesh<2,2>::BoundaryNodeIterator node_iter = p_mesh->GetBoundaryNodeIteratorBegin();
        node_iter != p_mesh->GetBoundaryNodeIteratorEnd();
        ++node_iter)
        {

            bcc.AddDirichletBoundaryCondition(*node_iter, vectorConstBCs[0]);
        }
*/
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
            if((columnNum==3 || columnNum==4 || columnNum==5 || columnNum==6))
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


        InhomogenousHeatEquationPde<2, 2, 1> pde(params.diffusionRates10);
        std::cout<<"Here4"<<std::endl;
        // solver
        InhomogenousCoupledPdeOdeSolverTemplated<2,2,1> solver(p_mesh, &pde, &bcc);
std::cout<<"Here5"<<std::endl;

        solver.SetTimes(params.startTime, params.endTime);

        solver.SetTimeStep(params.timestep/10);
        if(params.timestep>params.samplingTimestep)
        {
            solver.SetSamplingTimeStep(params.timestep);   
        }
        else
        {
            solver.SetSamplingTimeStep(params.samplingTimestep);
        }
   std::cout<<"Here6"<<std::endl;
        solver.SetOutputDirectory(params.outputDirName+"LinearParabolicPdeSystemWithCoupledOdeSystemSolverScaledBy1000dt0001");
        solver.SetInitialCondition(initial_condition);
        // solve
        solver.SolveAndWriteResultsToFile();
        

    }

    void TestTotalMassHeatDiffusionWithoutSourceScaled1000dt00001()
    {

      

        // mesh
        TetrahedralMesh<2,2>* p_mesh = new TetrahedralMesh<2,2>();

        p_mesh->ConstructRegularSlabMesh(params.MeshStepSize, params.MeshDimensions[0], params.MeshDimensions[1]);
        p_mesh->Scale(0.001,0.001);
        
        // Process Boundary Conditions
        BoundaryConditionsContainer<2,2,1> bcc;
    
        std::vector<ConstBoundaryCondition<2>*> vectorConstBCs;
        
        vectorConstBCs.push_back(new ConstBoundaryCondition<2>(params.bcValuesDir[0]));
        vectorConstBCs.push_back(new ConstBoundaryCondition<2>(params.bcValuesNeu[0]));
        
     std::cout<<"Here1"<<std::endl;   
/*
        for (TetrahedralMesh<2,2>::BoundaryNodeIterator node_iter = p_mesh->GetBoundaryNodeIteratorBegin();
        node_iter != p_mesh->GetBoundaryNodeIteratorEnd();
        ++node_iter)
        {

            bcc.AddDirichletBoundaryCondition(*node_iter, vectorConstBCs[0]);
        }
*/
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
            if((columnNum==3 || columnNum==4 || columnNum==5 || columnNum==6))
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


        InhomogenousHeatEquationPde<2, 2, 1> pde(params.diffusionRates10);
        std::cout<<"Here4"<<std::endl;
        // solver
        InhomogenousCoupledPdeOdeSolverTemplated<2,2,1> solver(p_mesh, &pde, &bcc);
std::cout<<"Here5"<<std::endl;

        solver.SetTimes(params.startTime, params.endTime);

        solver.SetTimeStep(params.timestep/1000);
        if(params.timestep>params.samplingTimestep)
        {
            solver.SetSamplingTimeStep(params.timestep);   
        }
        else
        {
            solver.SetSamplingTimeStep(params.samplingTimestep);
        }
   std::cout<<"Here6"<<std::endl;
        solver.SetOutputDirectory(params.outputDirName+"LinearParabolicPdeSystemWithCoupledOdeSystemSolverScaledBy1000dt00001");
        solver.SetInitialCondition(initial_condition);
        // solve
        solver.SolveAndWriteResultsToFile();
        

    }







};

#endif