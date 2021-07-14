#ifndef TESTSPATIALPATTERNING_HPP_
#define TESTSPATIALPATTERNING_HPP_

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

#include "SchnackenbergCoupledPdeSystem.hpp"
#include "InhomogenousSchnackenbergNonDimPde.hpp"



class TestSpatialPatterning : public AbstractCellBasedTestSuite
{
public:


    void TestSchnackenbergSystemOnSquareMeshChemChaste() throw (Exception)
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4096_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        double original_width_x = mesh.GetWidth(0);
        double original_width_y = mesh.GetWidth(1);
        double L = 100;
        mesh.Scale(L/original_width_x, L/original_width_y);

        SchnackenbergCoupledPdeSystem<2> pde(1, 40, 0.1, 1, 0.9, 1);

        BoundaryConditionsContainer<2,2,2> bcc;
        ConstBoundaryCondition<2>* p_zero_bc = new ConstBoundaryCondition<2>(0.0);
        for (TetrahedralMesh<2,2>::BoundaryElementIterator elem_iter = mesh.GetBoundaryElementIteratorBegin();
            elem_iter != mesh.GetBoundaryElementIteratorEnd();
            elem_iter++)
        {
            bcc.AddNeumannBoundaryCondition(*elem_iter, p_zero_bc, 0);
            bcc.AddNeumannBoundaryCondition(*elem_iter, p_zero_bc, 1);
        }

        LinearParabolicPdeSystemWithCoupledOdeSystemSolver<2,2,2> solver(&mesh, &pde, &bcc);

        double t_end = 100;
        solver.SetTimes(0, t_end);
        solver.SetTimeStep(1e-1);
        solver.SetSamplingTimeStep(1);
        solver.SetOutputDirectory("TestSchnackenbergSystemOnSquareMeshChemChaste");

        std::vector<double> init_conds(2*mesh.GetNumNodes());
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            init_conds[2*i] = fabs(1.0 + RandomNumberGenerator::Instance()->ranf());
            init_conds[2*i + 1] = fabs(0.9 + RandomNumberGenerator::Instance()->ranf());
        }
        Vec initial_condition = PetscTools::CreateVec(init_conds);
        solver.SetInitialCondition(initial_condition);

        solver.SolveAndWriteResultsToFile();
        PetscTools::Destroy(initial_condition);
    }


    void TestSchnackenbergSystemOnSquareMeshChemChasteInhomogenous() throw (Exception)
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4096_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        double original_width_x = mesh.GetWidth(0);
        double original_width_y = mesh.GetWidth(1);
        double L = 100;
        mesh.Scale(L/original_width_x, L/original_width_y);

        std::vector<double> diffusionRates = {1, 40};
        InhomogenousSchnackenbergNonDimPde<2,2,2> pde(diffusionRates, 0.1, 1.0, 0.9, 1.0);

        BoundaryConditionsContainer<2,2,2> bcc;
        ConstBoundaryCondition<2>* p_zero_bc = new ConstBoundaryCondition<2>(0.0);
        for (TetrahedralMesh<2,2>::BoundaryElementIterator elem_iter = mesh.GetBoundaryElementIteratorBegin();
            elem_iter != mesh.GetBoundaryElementIteratorEnd();
            elem_iter++)
        {
            bcc.AddNeumannBoundaryCondition(*elem_iter, p_zero_bc, 0);
            bcc.AddNeumannBoundaryCondition(*elem_iter, p_zero_bc, 1);
        }

        InhomogenousCoupledPdeOdeSolverTemplated<2,2,2> solver(&mesh, &pde, &bcc);



        double t_end = 100;
        solver.SetTimes(0, t_end);
        solver.SetTimeStep(1e-1);
        solver.SetSamplingTimeStep(1);
        solver.SetOutputDirectory("TestSchnackenbergSystemOnSquareMeshChemChasteInhomogenous");

        std::vector<double> init_conds(2*mesh.GetNumNodes());
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            init_conds[2*i] = fabs(1.0 + RandomNumberGenerator::Instance()->ranf());
            init_conds[2*i + 1] = fabs(0.9 + RandomNumberGenerator::Instance()->ranf());
        }
        Vec initial_condition = PetscTools::CreateVec(init_conds);
        solver.SetInitialCondition(initial_condition);

        solver.SolveAndWriteResultsToFile();
        PetscTools::Destroy(initial_condition);
    }

};

#endif