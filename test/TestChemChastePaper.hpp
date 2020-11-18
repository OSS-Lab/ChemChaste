#ifndef TESTCHEMCHASTEPAPER_HPP_
#define TESTCHEMCHASTEPAPER_HPP_

// chaste includes
#include "BoundaryConditionsContainer_extended.hpp"
#include "ChasteHeaders.hpp"

// general includes
#include "GeneralHeaders.hpp"

//ChemChaste includes
#include "ChemChasteHeaders.hpp"

#include "ChemicalCellFromFile.hpp"

struct ControlStruct {
    bool VariableDiffusionConstantCase = false;
    bool VariableDiffusionStepChangeCase = false;
    bool ExtracellularReaction = false;
    bool MolenaarSystem = true;
} control;

struct VariableDiffusionConstantCaseStruct {
    
    // Variables for the user modify

    // simulation
    std::string outputFilename = "TestChemChastePaperHetergenousCellTest2";
    double simultationEndTime = 10.0;


    // bulk
    std::string bulkDataFileRoot = "/home/chaste/projects/ChemChaste/src/Data/VariableDiffusionConstantCase/DomainField/";
    std::string domainFilename = "Domain.csv";
    std::string domainKeyFilename = "DomainKey.csv";
    std::string odeLabelFilename = "NodeSelector.csv";
    std::string odeKeyFilename = "OdeReactionFileKey.csv";
    std::string diffusionFilename = "DiffusionDatabaseFile.csv";
    std::string initialConditionsFilename = "InitialConditionFile.csv";
    std::string boundaryConditionsFilename = "BoundaryConditionFile.csv";

    // cell
    std::string cellDataFileRoot = "/home/chaste/projects/ChemChaste/src/Data/VariableDiffusionConstantCase/Cell/";
    std::string initialConditionFilename = "InitialCellConcentrations.csv";
    std::string membranePropertyFilename = "MembraneReactions.txt";
    std::string transportPropertyFilename = "TransportReactions.txt";
    std::string srnFilename = "Srn.txt";
    std::string cellCycleFilename = "SpeciesThreshold.csv";

    // mesh
    unsigned numberOfCellsAcross = 1;
    unsigned numberOfCellsUpwards = 1;
    bool nodeBasedSimulation = true; // default to mesh based simulation, controls the interactions between the cells
    double nodesOnlyCutOff = 1.5; // for node based simulation case
    double originOfCellMeshRelativePdeFEMesh = -4.0;
    double linearForceCutOffLength = 1.5;

    


}VDCC;

struct VariableDiffusionStepChangeCaseStruct {
        
    // Variables for the user modify

    // simulation
    std::string outputFilename = "TestChemChastePaperVariableDiffusionStepChangeCase";
    double simultationEndTime = 100.0;


    // bulk
    std::string bulkDataFileRoot = "/home/chaste/projects/ChemChaste/src/Data/VariableDiffusionStepChangeCase/DomainField/";
    std::string domainFilename = "Domain.csv";
    std::string domainKeyFilename = "DomainKey.csv";
    std::string odeLabelFilename = "NodeSelector.csv";
    std::string odeKeyFilename = "OdeReactionFileKey.csv";
    std::string diffusionFilename = "DiffusionDatabaseFile.csv";
    std::string initialConditionsFilename = "InitialConditionFile.csv";
    std::string boundaryConditionsFilename = "BoundaryConditionFile.csv";

    // cell
    std::string cellDataFileRoot = "/home/chaste/projects/ChemChaste/src/Data/VariableDiffusionStepChangeCase/Cell/";
    std::string initialConditionFilename = "InitialCellConcentrations.csv";
    std::string membranePropertyFilename = "MembraneReactions.txt";
    std::string transportPropertyFilename = "TransportReactions.txt";
    std::string srnFilename = "Srn.txt";
    std::string cellCycleFilename = "SpeciesThreshold.csv";

    // mesh
    unsigned numberOfCellsAcross = 2;
    unsigned numberOfCellsUpwards = 2;
    bool nodeBasedSimulation = true; // default to mesh based simulation, controls the interactions between the cells
    double nodesOnlyCutOff = 1.5; // for node based simulation case
    double originOfCellMeshRelativePdeFEMesh = -4.0;
    double linearForceCutOffLength = 1.5;

    // System properties
    const unsigned probDim =2; // need to set manually to the number of diffusive variables for the pde solver to solve
    const unsigned spaceDim=2;
    const unsigned elementDim=2;

}VDSC;

struct ExtracellularReactionCaseStruct {
        
    // Variables for the user modify

    // simulation
    std::string outputFilename = "TestChemChastePaperExtracellularReactionCaseTest";
    double simultationEndTime = 10.0;


    // bulk
    std::string bulkDataFileRoot = "/home/chaste/projects/ChemChaste/src/Data/ExtracellularReaction/DomainField/";
    std::string domainFilename = "Domain.csv";
    std::string domainKeyFilename = "DomainKey.csv";
    std::string odeLabelFilename = "NodeSelector.csv";
    std::string odeKeyFilename = "OdeReactionFileKey.csv";
    std::string diffusionFilename = "DiffusionDatabaseFile.csv";
    std::string initialConditionsFilename = "InitialConditionFile.csv";
    std::string boundaryConditionsFilename = "BoundaryConditionFile.csv";

    // cell
    std::string cellDataFileRoot = "/home/chaste/projects/ChemChaste/src/Data/ExtracellularReaction/Cell/";
    std::string initialConditionFilename = "InitialCellConcentrations.csv";
    std::string membranePropertyFilename = "";
    std::string transportPropertyFilename = "TransportReactions.txt";
    std::string srnFilename = "Srn.txt";
    std::string cellCycleFilename = "SpeciesThreshold.csv";

    // mesh
    unsigned numberOfCellsAcross = 2;
    unsigned numberOfCellsUpwards = 2;
    bool nodeBasedSimulation = true; // default to mesh based simulation, controls the interactions between the cells
    double nodesOnlyCutOff = 1.5; // for node based simulation case
    double originOfCellMeshRelativePdeFEMesh = -4.0;
    double linearForceCutOffLength = 1.5;

    // System properties
    const unsigned probDim =3; // need to set manually to the number of diffusive variables for the pde solver to solve
    const unsigned spaceDim=2;
    const unsigned elementDim=2;

}ER;

struct MolenaarSystemCaseStruct {
        
    // Variables for the user modify

    // simulation
    std::string outputFilename = "TestChemChastePaperMolenaarSystem";
    double simultationEndTime = 10.0;


    // bulk
    std::string bulkDataFileRoot = "/home/chaste/projects/ChemChaste/src/Data/MolenaarPaper/DomainField/";
    std::string domainFilename = "Domain.csv";
    std::string domainKeyFilename = "DomainKey.csv";
    std::string odeLabelFilename = "NodeSelector.csv";
    std::string odeKeyFilename = "OdeReactionFileKey.csv";
    std::string diffusionFilename = "DiffusionDatabaseFile.csv";
    std::string initialConditionsFilename = "InitialConditionFile.csv";
    std::string boundaryConditionsFilename = "BoundaryConditionFile.csv";

    // cell
    std::string cellDataFileRoot = "/home/chaste/projects/ChemChaste/src/Data/MolenaarPaper/Cell/";
    std::string initialConditionFilename = "InitialCellConcentrations.csv";
    std::string membranePropertyFilename = "MembraneReactions.txt";
    std::string transportPropertyFilename = "TransportReactions.txt";
    std::string srnFilename = "Srn.txt";
    std::string cellCycleFilename = "SpeciesThreshold.csv";

    // mesh
    unsigned numberOfCellsAcross = 1;
    unsigned numberOfCellsUpwards = 1;
    bool nodeBasedSimulation = true; // default to mesh based simulation, controls the interactions between the cells
    double nodesOnlyCutOff = 1.5; // for node based simulation case
    double originOfCellMeshRelativePdeFEMesh = -4.0;
    double linearForceCutOffLength = 1.5;

    // System properties
    const unsigned probDim =1; // need to set manually to the number of diffusive variables for the pde solver to solve
    const unsigned spaceDim=2;
    const unsigned elementDim=2;

}ML;



class TestChemChastePaper : public AbstractCellBasedTestSuite
{
public:

    void TestVariableDiffusionConstantCase()
    {
        
        if(control.VariableDiffusionConstantCase)
        {

            // System properties
            const unsigned probDim =2; // need to set manually to the number of diffusive variables for the pde solver to solve
            const unsigned spaceDim=2;
            const unsigned elementDim=2;

            // bulk
            std::string bulkDataFileRoot = VDCC.bulkDataFileRoot;
            std::string cellLabelFilename = "";
            std::string cellKeyFilename = "";
            std::string domainFilename = VDCC.domainFilename;
            std::string domainKeyFilename = VDCC.domainKeyFilename;
            std::string odeLabelFilename = VDCC.odeLabelFilename;
            std::string odeKeyFilename = VDCC.odeKeyFilename;
            std::string diffusionFilename = VDCC.diffusionFilename;
            std::string initialConditionsFilename = VDCC.initialConditionsFilename;
            std::string boundaryConditionsFilename = VDCC.boundaryConditionsFilename;

            // cell
            std::string cellDataFileRoot = VDCC.cellDataFileRoot;
            std::string initialConditionFilename = VDCC.initialConditionFilename;
            std::string membranePropertyFilename = VDCC.membranePropertyFilename;
            std::string transportPropertyFilename = VDCC.transportPropertyFilename;
            std::string srnFilename = VDCC.srnFilename;
            std::string cellCycleFilename = VDCC.cellCycleFilename;

        
            // make mesh of cell with associated mesh based cell population
            HoneycombMeshGenerator generator(VDCC.numberOfCellsAcross,VDCC.numberOfCellsUpwards);    // Parameters are: cells across, cells up
            MutableMesh<elementDim,spaceDim>* p_mesh = generator.GetMesh();


            std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

            std::vector<CellPtr> cells;

            std::vector<std::string> chemicalCellSpeciesNames;

            // assume cell at each node in cell layer mesh
            bool IsFirstCell = true;
            for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
            {
                // provide each cell with a transport cell property and membrane property, cell cycle, wild type states

                ChemicalCellFromFile* p_cell_reader = new ChemicalCellFromFile(
                                    VDCC.cellDataFileRoot+VDCC.cellCycleFilename, 
                                    VDCC.cellDataFileRoot+VDCC.srnFilename,
                                    VDCC.cellDataFileRoot+VDCC.initialConditionFilename,
                                    VDCC.cellDataFileRoot+VDCC.transportPropertyFilename,
                                    VDCC.cellDataFileRoot+VDCC.membranePropertyFilename);
                
                if(IsFirstCell)
                {
                    chemicalCellSpeciesNames =  p_cell_reader->GetFullChemicalNamesVector();
                    IsFirstCell = false;
                }
                cells.push_back(p_cell_reader -> GetCellPtr());
            }
    

            // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
            ChastePoint<spaceDim> lower(VDCC.originOfCellMeshRelativePdeFEMesh,VDCC.originOfCellMeshRelativePdeFEMesh);
            ChastePoint<spaceDim> upper(10.0, 10.0);
            MAKE_PTR_ARGS(ChasteCuboid<spaceDim>, p_cuboid, (lower, upper));

            // set up chemical domain field
            ChemicalDomainFieldForCellCoupling<elementDim,spaceDim,probDim>* p_Pde_field = new ChemicalDomainFieldForCellCoupling<elementDim,spaceDim,probDim>
                                    (   VDCC.bulkDataFileRoot,
                                        VDCC.bulkDataFileRoot+cellLabelFilename,
                                        VDCC.bulkDataFileRoot+cellKeyFilename,
                                        VDCC.bulkDataFileRoot+VDCC.domainFilename, 
                                        VDCC.bulkDataFileRoot+VDCC.domainKeyFilename, 
                                        VDCC.bulkDataFileRoot+VDCC.odeLabelFilename, 
                                        VDCC.bulkDataFileRoot+VDCC.odeKeyFilename, 
                                        VDCC.bulkDataFileRoot+VDCC.diffusionFilename, 
                                        VDCC.bulkDataFileRoot+VDCC.initialConditionsFilename, 
                                        VDCC.bulkDataFileRoot+VDCC.boundaryConditionsFilename);

            
            boost::shared_ptr<ParabolicBoxDomainPdeSystemModifier<elementDim,spaceDim,probDim>> p_pde_modifier(new ParabolicBoxDomainPdeSystemModifier<elementDim,spaceDim,probDim>(p_Pde_field, p_cuboid));
            
            boost::shared_ptr<ChemicalTrackingModifier<elementDim,spaceDim>> p_chemical_tracking_modifier(new ChemicalTrackingModifier<elementDim,spaceDim>());
            
            if(VDCC.nodeBasedSimulation)
            {
                NodesOnlyMesh<spaceDim> mesh;
           
                mesh.ConstructNodesWithoutMesh(*p_mesh, VDCC.nodesOnlyCutOff);
            
                NodeBasedCellPopulation<spaceDim> cell_population(mesh, cells);   

                // writers

                cell_population.AddCellWriter<CellIdWriter>();
                cell_population.AddCellWriter<CellAgesWriter>();
                cell_population.AddCellWriter<CellLocationIndexWriter>();

                

                for(unsigned i=0; i<chemicalCellSpeciesNames.size(); i++)
                {
                    boost::shared_ptr<CellDataItemWriter<elementDim,spaceDim>> dataWriter(new CellDataItemWriter<elementDim,spaceDim>(chemicalCellSpeciesNames[i]));
                    cell_population.AddCellWriter(dataWriter);
                }
                    
            
                OffLatticeSimulation<spaceDim> simulator(cell_population);

                simulator.AddSimulationModifier(p_pde_modifier);

                simulator.AddSimulationModifier(p_chemical_tracking_modifier);
            
                simulator.SetOutputDirectory(VDCC.outputFilename);
                simulator.SetEndTime(VDCC.simultationEndTime);

                MAKE_PTR(GeneralisedLinearSpringForce<spaceDim>, p_linear_force);
                p_linear_force->SetCutOffLength(VDCC.linearForceCutOffLength);
                simulator.AddForce(p_linear_force);

                std::cout<<"=============================================="<<std::endl;
                std::cout<<"OffLatticeSimulation -> AbstractCellBasedSimulation :: Solve()"<<std::endl;
                std::cout<<"=============================================="<<std::endl;
                simulator.Solve();

            }
            else
            {
                
                MeshBasedCellPopulation<spaceDim> cell_population(*p_mesh, cells, location_indices);
                // writers

                cell_population.AddCellWriter<CellIdWriter>();
                cell_population.AddCellWriter<CellAgesWriter>();
                cell_population.AddCellWriter<CellLocationIndexWriter>();

                

                for(unsigned i=0; i<chemicalCellSpeciesNames.size(); i++)
                {
                    boost::shared_ptr<CellDataItemWriter<elementDim,spaceDim>> dataWriter(new CellDataItemWriter<elementDim,spaceDim>(chemicalCellSpeciesNames[i]));
                    cell_population.AddCellWriter(dataWriter);
                }
                    
            
                OffLatticeSimulation<spaceDim> simulator(cell_population);

                cell_population.SetWriteVtkAsPoints(false);
                cell_population.AddPopulationWriter<VoronoiDataWriter>();

                simulator.AddSimulationModifier(p_pde_modifier);

                simulator.AddSimulationModifier(p_chemical_tracking_modifier);
            
                simulator.SetOutputDirectory(VDCC.outputFilename);
                simulator.SetEndTime(VDCC.simultationEndTime);

                MAKE_PTR(GeneralisedLinearSpringForce<spaceDim>, p_linear_force);
                p_linear_force->SetCutOffLength(VDCC.linearForceCutOffLength);
                simulator.AddForce(p_linear_force);

                std::cout<<"=============================================="<<std::endl;
                std::cout<<"OffLatticeSimulation -> AbstractCellBasedSimulation :: Solve()"<<std::endl;
                std::cout<<"=============================================="<<std::endl;
                simulator.Solve();
            }

        }
    
    }

    void TestVariableDiffusionStepChangeCase()
    {
        if(control.VariableDiffusionStepChangeCase)
        {
            // System properties
            const unsigned probDim =2; // need to set manually to the number of diffusive variables for the pde solver to solve
            const unsigned spaceDim=2;
            const unsigned elementDim=2;

            // bulk
            std::string bulkDataFileRoot = VDSC.bulkDataFileRoot;
            std::string cellLabelFilename = "";
            std::string cellKeyFilename = "";
            std::string domainFilename = VDSC.domainFilename;
            std::string domainKeyFilename = VDSC.domainKeyFilename;
            std::string odeLabelFilename = VDSC.odeLabelFilename;
            std::string odeKeyFilename = VDSC.odeKeyFilename;
            std::string diffusionFilename = VDSC.diffusionFilename;
            std::string initialConditionsFilename = VDSC.initialConditionsFilename;
            std::string boundaryConditionsFilename = VDSC.boundaryConditionsFilename;

            // cell
            std::string cellDataFileRoot = VDSC.cellDataFileRoot;
            std::string initialConditionFilename = VDSC.initialConditionFilename;
            std::string membranePropertyFilename = VDSC.membranePropertyFilename;
            std::string transportPropertyFilename = VDSC.transportPropertyFilename;
            std::string srnFilename = VDSC.srnFilename;
            std::string cellCycleFilename = VDSC.cellCycleFilename;

        
            // make mesh of cell with associated mesh based cell population
            HoneycombMeshGenerator generator(VDSC.numberOfCellsAcross,VDSC.numberOfCellsUpwards);    // Parameters are: cells across, cells up
            MutableMesh<elementDim,spaceDim>* p_mesh = generator.GetMesh();


            std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

            std::vector<CellPtr> cells;

            std::vector<std::string> chemicalCellSpeciesNames;

            // assume cell at each node in cell layer mesh
            bool IsFirstCell = true;
            for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
            {
                // provide each cell with a transport cell property and membrane property, cell cycle, wild type states

                ChemicalCellFromFile* p_cell_reader = new ChemicalCellFromFile(
                                    VDSC.cellDataFileRoot+VDSC.cellCycleFilename, 
                                    VDSC.cellDataFileRoot+VDSC.srnFilename,
                                    VDSC.cellDataFileRoot+VDSC.initialConditionFilename,
                                    VDSC.cellDataFileRoot+VDSC.transportPropertyFilename,
                                    VDSC.cellDataFileRoot+VDSC.membranePropertyFilename);
                
                if(IsFirstCell)
                {
                    chemicalCellSpeciesNames =  p_cell_reader->GetFullChemicalNamesVector();
                    IsFirstCell = false;
                }
                cells.push_back(p_cell_reader -> GetCellPtr());
            }
    

            // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
            ChastePoint<spaceDim> lower(VDSC.originOfCellMeshRelativePdeFEMesh,VDSC.originOfCellMeshRelativePdeFEMesh);
            ChastePoint<spaceDim> upper(10.0, 10.0);
            MAKE_PTR_ARGS(ChasteCuboid<spaceDim>, p_cuboid, (lower, upper));

            // set up chemical domain field
            ChemicalDomainFieldForCellCoupling<elementDim,spaceDim,probDim>* p_Pde_field = new ChemicalDomainFieldForCellCoupling<elementDim,spaceDim,probDim>
                                    (   VDSC.bulkDataFileRoot,
                                        VDSC.bulkDataFileRoot+cellLabelFilename,
                                        VDSC.bulkDataFileRoot+cellKeyFilename,
                                        VDSC.bulkDataFileRoot+VDSC.domainFilename, 
                                        VDSC.bulkDataFileRoot+VDSC.domainKeyFilename, 
                                        VDSC.bulkDataFileRoot+VDSC.odeLabelFilename, 
                                        VDSC.bulkDataFileRoot+VDSC.odeKeyFilename, 
                                        VDSC.bulkDataFileRoot+VDSC.diffusionFilename, 
                                        VDSC.bulkDataFileRoot+VDSC.initialConditionsFilename, 
                                        VDSC.bulkDataFileRoot+VDSC.boundaryConditionsFilename);

            
            boost::shared_ptr<ParabolicBoxDomainPdeSystemModifier<elementDim,spaceDim,probDim>> p_pde_modifier(new ParabolicBoxDomainPdeSystemModifier<elementDim,spaceDim,probDim>(p_Pde_field, p_cuboid));
            
            boost::shared_ptr<ChemicalTrackingModifier<elementDim,spaceDim>> p_chemical_tracking_modifier(new ChemicalTrackingModifier<elementDim,spaceDim>());
            
            if(VDSC.nodeBasedSimulation)
            {
                NodesOnlyMesh<spaceDim> mesh;
           
                mesh.ConstructNodesWithoutMesh(*p_mesh, VDSC.nodesOnlyCutOff);
            
                NodeBasedCellPopulation<spaceDim> cell_population(mesh, cells);   

                // writers

                cell_population.AddCellWriter<CellIdWriter>();
                cell_population.AddCellWriter<CellAgesWriter>();
                cell_population.AddCellWriter<CellLocationIndexWriter>();

                

                for(unsigned i=0; i<chemicalCellSpeciesNames.size(); i++)
                {
                    boost::shared_ptr<CellDataItemWriter<elementDim,spaceDim>> dataWriter(new CellDataItemWriter<elementDim,spaceDim>(chemicalCellSpeciesNames[i]));
                    cell_population.AddCellWriter(dataWriter);
                }
                    
            
                OffLatticeSimulation<spaceDim> simulator(cell_population);

                simulator.AddSimulationModifier(p_pde_modifier);

                simulator.AddSimulationModifier(p_chemical_tracking_modifier);
            
                simulator.SetOutputDirectory(VDSC.outputFilename);
                simulator.SetEndTime(VDSC.simultationEndTime);

                MAKE_PTR(GeneralisedLinearSpringForce<spaceDim>, p_linear_force);
                p_linear_force->SetCutOffLength(VDSC.linearForceCutOffLength);
                simulator.AddForce(p_linear_force);

                std::cout<<"=============================================="<<std::endl;
                std::cout<<"OffLatticeSimulation -> AbstractCellBasedSimulation :: Solve()"<<std::endl;
                std::cout<<"=============================================="<<std::endl;
                simulator.Solve();

            }
            else
            {
                
                MeshBasedCellPopulation<spaceDim> cell_population(*p_mesh, cells, location_indices);
                // writers

                cell_population.AddCellWriter<CellIdWriter>();
                cell_population.AddCellWriter<CellAgesWriter>();
                cell_population.AddCellWriter<CellLocationIndexWriter>();

                

                for(unsigned i=0; i<chemicalCellSpeciesNames.size(); i++)
                {
                    boost::shared_ptr<CellDataItemWriter<elementDim,spaceDim>> dataWriter(new CellDataItemWriter<elementDim,spaceDim>(chemicalCellSpeciesNames[i]));
                    cell_population.AddCellWriter(dataWriter);
                }
                    
            
                OffLatticeSimulation<spaceDim> simulator(cell_population);

                cell_population.SetWriteVtkAsPoints(false);
                cell_population.AddPopulationWriter<VoronoiDataWriter>();

                simulator.AddSimulationModifier(p_pde_modifier);

                simulator.AddSimulationModifier(p_chemical_tracking_modifier);
            
                simulator.SetOutputDirectory(VDSC.outputFilename);
                simulator.SetEndTime(VDSC.simultationEndTime);

                MAKE_PTR(GeneralisedLinearSpringForce<spaceDim>, p_linear_force);
                p_linear_force->SetCutOffLength(VDSC.linearForceCutOffLength);
                simulator.AddForce(p_linear_force);

                std::cout<<"=============================================="<<std::endl;
                std::cout<<"OffLatticeSimulation -> AbstractCellBasedSimulation :: Solve()"<<std::endl;
                std::cout<<"=============================================="<<std::endl;
                simulator.Solve();
            }

        }
    


    }

    void TestExtracellularReaction()
    {
        if(control.ExtracellularReaction)
        {
            // System properties
            const unsigned probDim =4; // need to set manually to the number of diffusive variables for the pde solver to solve
            const unsigned spaceDim=2;
            const unsigned elementDim=2;

            // bulk
            std::string bulkDataFileRoot = ER.bulkDataFileRoot;
            std::string cellLabelFilename = "";
            std::string cellKeyFilename = "";
            std::string domainFilename = ER.domainFilename;
            std::string domainKeyFilename = ER.domainKeyFilename;
            std::string odeLabelFilename = ER.odeLabelFilename;
            std::string odeKeyFilename = ER.odeKeyFilename;
            std::string diffusionFilename = ER.diffusionFilename;
            std::string initialConditionsFilename = ER.initialConditionsFilename;
            std::string boundaryConditionsFilename = ER.boundaryConditionsFilename;

            // cell
            std::string cellDataFileRoot = ER.cellDataFileRoot;
            std::string initialConditionFilename = ER.initialConditionFilename;
            std::string membranePropertyFilename = ER.membranePropertyFilename;
            std::string transportPropertyFilename = ER.transportPropertyFilename;
            std::string srnFilename = ER.srnFilename;
            std::string cellCycleFilename = ER.cellCycleFilename;

        
            // make mesh of cell with associated mesh based cell population
            HoneycombMeshGenerator generator(ER.numberOfCellsAcross,ER.numberOfCellsUpwards);    // Parameters are: cells across, cells up
            MutableMesh<elementDim,spaceDim>* p_mesh = generator.GetMesh();


            std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

            std::vector<CellPtr> cells;

            std::vector<std::string> chemicalCellSpeciesNames;

            // assume cell at each node in cell layer mesh
            bool IsFirstCell = true;
            for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
            {
                // provide each cell with a transport cell property and membrane property, cell cycle, wild type states

                ChemicalCellFromFile* p_cell_reader = new ChemicalCellFromFile(
                                    ER.cellDataFileRoot+ER.cellCycleFilename, 
                                    ER.cellDataFileRoot+ER.srnFilename,
                                    ER.cellDataFileRoot+ER.initialConditionFilename,
                                    ER.cellDataFileRoot+ER.transportPropertyFilename,
                                    ER.cellDataFileRoot+ER.membranePropertyFilename);
                
                if(IsFirstCell)
                {
                    chemicalCellSpeciesNames =  p_cell_reader->GetFullChemicalNamesVector();
                    IsFirstCell = false;
                }
                cells.push_back(p_cell_reader -> GetCellPtr());
            }
    

            // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
            ChastePoint<spaceDim> lower(ER.originOfCellMeshRelativePdeFEMesh,ER.originOfCellMeshRelativePdeFEMesh);
            ChastePoint<spaceDim> upper(10.0, 10.0);
            MAKE_PTR_ARGS(ChasteCuboid<spaceDim>, p_cuboid, (lower, upper));

            // set up chemical domain field
            ChemicalDomainFieldForCellCoupling<elementDim,spaceDim,probDim>* p_Pde_field = new ChemicalDomainFieldForCellCoupling<elementDim,spaceDim,probDim>
                                    (   ER.bulkDataFileRoot,
                                        ER.bulkDataFileRoot+cellLabelFilename,
                                        ER.bulkDataFileRoot+cellKeyFilename,
                                        ER.bulkDataFileRoot+ER.domainFilename, 
                                        ER.bulkDataFileRoot+ER.domainKeyFilename, 
                                        ER.bulkDataFileRoot+ER.odeLabelFilename, 
                                        ER.bulkDataFileRoot+ER.odeKeyFilename, 
                                        ER.bulkDataFileRoot+ER.diffusionFilename, 
                                        ER.bulkDataFileRoot+ER.initialConditionsFilename, 
                                        ER.bulkDataFileRoot+ER.boundaryConditionsFilename);

            
            boost::shared_ptr<ParabolicBoxDomainPdeSystemModifier<elementDim,spaceDim,probDim>> p_pde_modifier(new ParabolicBoxDomainPdeSystemModifier<elementDim,spaceDim,probDim>(p_Pde_field, p_cuboid));
            
            boost::shared_ptr<ChemicalTrackingModifier<elementDim,spaceDim>> p_chemical_tracking_modifier(new ChemicalTrackingModifier<elementDim,spaceDim>());
            
            if(ER.nodeBasedSimulation)
            {
                NodesOnlyMesh<spaceDim> mesh;
           
                mesh.ConstructNodesWithoutMesh(*p_mesh, ER.nodesOnlyCutOff);
            
                NodeBasedCellPopulation<spaceDim> cell_population(mesh, cells);   

                // writers

                cell_population.AddCellWriter<CellIdWriter>();
                cell_population.AddCellWriter<CellAgesWriter>();
                cell_population.AddCellWriter<CellLocationIndexWriter>();

                

                for(unsigned i=0; i<chemicalCellSpeciesNames.size(); i++)
                {
                    boost::shared_ptr<CellDataItemWriter<elementDim,spaceDim>> dataWriter(new CellDataItemWriter<elementDim,spaceDim>(chemicalCellSpeciesNames[i]));
                    cell_population.AddCellWriter(dataWriter);
                }
                    
            
                OffLatticeSimulation<spaceDim> simulator(cell_population);

                simulator.AddSimulationModifier(p_pde_modifier);

                simulator.AddSimulationModifier(p_chemical_tracking_modifier);
            
                simulator.SetOutputDirectory(ER.outputFilename);
                simulator.SetEndTime(ER.simultationEndTime);

                MAKE_PTR(GeneralisedLinearSpringForce<spaceDim>, p_linear_force);
                p_linear_force->SetCutOffLength(ER.linearForceCutOffLength);
                simulator.AddForce(p_linear_force);

                std::cout<<"=============================================="<<std::endl;
                std::cout<<"OffLatticeSimulation -> AbstractCellBasedSimulation :: Solve()"<<std::endl;
                std::cout<<"=============================================="<<std::endl;
                simulator.Solve();

            }
            else
            {
                
                MeshBasedCellPopulation<spaceDim> cell_population(*p_mesh, cells, location_indices);
                // writers

                cell_population.SetWriteVtkAsPoints(false);
            cell_population.AddPopulationWriter<VoronoiDataWriter>();

                cell_population.AddCellWriter<CellIdWriter>();
                cell_population.AddCellWriter<CellAgesWriter>();
                cell_population.AddCellWriter<CellLocationIndexWriter>();

                

                for(unsigned i=0; i<chemicalCellSpeciesNames.size(); i++)
                {
                    boost::shared_ptr<CellDataItemWriter<elementDim,spaceDim>> dataWriter(new CellDataItemWriter<elementDim,spaceDim>(chemicalCellSpeciesNames[i]));
                    cell_population.AddCellWriter(dataWriter);
                }
                    
            
                OffLatticeSimulation<spaceDim> simulator(cell_population);

                simulator.AddSimulationModifier(p_pde_modifier);

                simulator.AddSimulationModifier(p_chemical_tracking_modifier);
            
                simulator.SetOutputDirectory(ER.outputFilename);
                simulator.SetEndTime(ER.simultationEndTime);

                MAKE_PTR(GeneralisedLinearSpringForce<spaceDim>, p_linear_force);
                p_linear_force->SetCutOffLength(ER.linearForceCutOffLength);
                simulator.AddForce(p_linear_force);

                std::cout<<"=============================================="<<std::endl;
                std::cout<<"OffLatticeSimulation -> AbstractCellBasedSimulation :: Solve()"<<std::endl;
                std::cout<<"=============================================="<<std::endl;
                simulator.Solve();
            }

        }


    }

    void TestMolenaarSystem()
    {
        if(control.MolenaarSystem)
        {
            // System properties
            const unsigned probDim =1; // need to set manually to the number of diffusive variables for the pde solver to solve
            const unsigned spaceDim=2;
            const unsigned elementDim=2;

            // bulk
            std::string bulkDataFileRoot = ML.bulkDataFileRoot;
            std::string cellLabelFilename = "";
            std::string cellKeyFilename = "";
            std::string domainFilename = ML.domainFilename;
            std::string domainKeyFilename = ML.domainKeyFilename;
            std::string odeLabelFilename = ML.odeLabelFilename;
            std::string odeKeyFilename = ML.odeKeyFilename;
            std::string diffusionFilename = ML.diffusionFilename;
            std::string initialConditionsFilename = ML.initialConditionsFilename;
            std::string boundaryConditionsFilename = ML.boundaryConditionsFilename;

            // cell
            std::string cellDataFileRoot = ML.cellDataFileRoot;
            std::string initialConditionFilename = ML.initialConditionFilename;
            std::string membranePropertyFilename = ML.membranePropertyFilename;
            std::string transportPropertyFilename = ML.transportPropertyFilename;
            std::string srnFilename = ML.srnFilename;
            std::string cellCycleFilename = ML.cellCycleFilename;

        
            // make mesh of cell with associated mesh based cell population
            HoneycombMeshGenerator generator(ML.numberOfCellsAcross,ML.numberOfCellsUpwards);    // Parameters are: cells across, cells up
            MutableMesh<elementDim,spaceDim>* p_mesh = generator.GetMesh();


            std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

            std::vector<CellPtr> cells;

            std::vector<std::string> chemicalCellSpeciesNames;

            // assume cell at each node in cell layer mesh
            bool IsFirstCell = true;
            for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
            {
                // provide each cell with a transport cell property and membrane property, cell cycle, wild type states

                ChemicalCellFromFile* p_cell_reader = new ChemicalCellFromFile(
                                    ML.cellDataFileRoot+ML.cellCycleFilename, 
                                    ML.cellDataFileRoot+ML.srnFilename,
                                    ML.cellDataFileRoot+ML.initialConditionFilename,
                                    ML.cellDataFileRoot+ML.transportPropertyFilename,
                                    ML.cellDataFileRoot+ML.membranePropertyFilename);
                
                if(IsFirstCell)
                {
                    chemicalCellSpeciesNames =  p_cell_reader->GetFullChemicalNamesVector();
                    IsFirstCell = false;
                }
                cells.push_back(p_cell_reader -> GetCellPtr());
            }
    

            // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
            ChastePoint<spaceDim> lower(ML.originOfCellMeshRelativePdeFEMesh,ML.originOfCellMeshRelativePdeFEMesh);
            ChastePoint<spaceDim> upper(10.0, 10.0);
            MAKE_PTR_ARGS(ChasteCuboid<spaceDim>, p_cuboid, (lower, upper));

            // set up chemical domain field
            ChemicalDomainFieldForCellCoupling<elementDim,spaceDim,probDim>* p_Pde_field = new ChemicalDomainFieldForCellCoupling<elementDim,spaceDim,probDim>
                                    (   ML.bulkDataFileRoot,
                                        ML.bulkDataFileRoot+cellLabelFilename,
                                        ML.bulkDataFileRoot+cellKeyFilename,
                                        ML.bulkDataFileRoot+ML.domainFilename, 
                                        ML.bulkDataFileRoot+ML.domainKeyFilename, 
                                        ML.bulkDataFileRoot+ML.odeLabelFilename, 
                                        ML.bulkDataFileRoot+ML.odeKeyFilename, 
                                        ML.bulkDataFileRoot+ML.diffusionFilename, 
                                        ML.bulkDataFileRoot+ML.initialConditionsFilename, 
                                        ML.bulkDataFileRoot+ML.boundaryConditionsFilename);

            
            boost::shared_ptr<ParabolicBoxDomainPdeSystemModifier<elementDim,spaceDim,probDim>> p_pde_modifier(new ParabolicBoxDomainPdeSystemModifier<elementDim,spaceDim,probDim>(p_Pde_field, p_cuboid));
            
            boost::shared_ptr<ChemicalTrackingModifier<elementDim,spaceDim>> p_chemical_tracking_modifier(new ChemicalTrackingModifier<elementDim,spaceDim>());
            
            if(ML.nodeBasedSimulation)
            {
                NodesOnlyMesh<spaceDim> mesh;
           
                mesh.ConstructNodesWithoutMesh(*p_mesh, ML.nodesOnlyCutOff);
            
                NodeBasedCellPopulation<spaceDim> cell_population(mesh, cells);   

                // writers

                cell_population.AddCellWriter<CellIdWriter>();
                cell_population.AddCellWriter<CellAgesWriter>();
                cell_population.AddCellWriter<CellLocationIndexWriter>();

                

                for(unsigned i=0; i<chemicalCellSpeciesNames.size(); i++)
                {
                    boost::shared_ptr<CellDataItemWriter<elementDim,spaceDim>> dataWriter(new CellDataItemWriter<elementDim,spaceDim>(chemicalCellSpeciesNames[i]));
                    cell_population.AddCellWriter(dataWriter);
                }
                    
            
                OffLatticeSimulation<spaceDim> simulator(cell_population);

                simulator.AddSimulationModifier(p_pde_modifier);

                simulator.AddSimulationModifier(p_chemical_tracking_modifier);
            
                simulator.SetOutputDirectory(ML.outputFilename);
                simulator.SetEndTime(ML.simultationEndTime);

                MAKE_PTR(GeneralisedLinearSpringForce<spaceDim>, p_linear_force);
                p_linear_force->SetCutOffLength(ML.linearForceCutOffLength);
                simulator.AddForce(p_linear_force);

                std::cout<<"=============================================="<<std::endl;
                std::cout<<"OffLatticeSimulation -> AbstractCellBasedSimulation :: Solve()"<<std::endl;
                std::cout<<"=============================================="<<std::endl;
                simulator.Solve();

            }
            else
            {
                
                MeshBasedCellPopulation<spaceDim> cell_population(*p_mesh, cells, location_indices);
                // writers
                cell_population.SetWriteVtkAsPoints(false);
                cell_population.AddPopulationWriter<VoronoiDataWriter>();

                cell_population.AddCellWriter<CellIdWriter>();
                cell_population.AddCellWriter<CellAgesWriter>();
                cell_population.AddCellWriter<CellLocationIndexWriter>();

                

                for(unsigned i=0; i<chemicalCellSpeciesNames.size(); i++)
                {
                    boost::shared_ptr<CellDataItemWriter<elementDim,spaceDim>> dataWriter(new CellDataItemWriter<elementDim,spaceDim>(chemicalCellSpeciesNames[i]));
                    cell_population.AddCellWriter(dataWriter);
                }
                    
            
                OffLatticeSimulation<spaceDim> simulator(cell_population);

                simulator.AddSimulationModifier(p_pde_modifier);

                simulator.AddSimulationModifier(p_chemical_tracking_modifier);
            
                simulator.SetOutputDirectory(ML.outputFilename);
                simulator.SetEndTime(ML.simultationEndTime);

                MAKE_PTR(GeneralisedLinearSpringForce<spaceDim>, p_linear_force);
                p_linear_force->SetCutOffLength(ML.linearForceCutOffLength);
                simulator.AddForce(p_linear_force);

                std::cout<<"=============================================="<<std::endl;
                std::cout<<"OffLatticeSimulation -> AbstractCellBasedSimulation :: Solve()"<<std::endl;
                std::cout<<"=============================================="<<std::endl;
                simulator.Solve();
            }

        }


    }

};

#endif