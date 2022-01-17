//#include "Cell_virtual.hpp"
//#include "AbstractCellBasedTestSuite.hpp"
#include "BoundaryConditionsContainer_extended.hpp"
#include "ChemChasteFeAssemblerCommon.hpp"
#include "ChemChasteVolumeAssembler.hpp"
#include "ChemChasteSurfaceAssembler.hpp"


#include <cxxtest/TestSuite.h>

//#include <cxxtest/GlobalFixture.h>
#include "CheckpointArchiveTypes.hpp"
#include "SmartPointers.hpp"

#include <string>
#include <fstream>
#include <iostream>
#include <chrono>



//ChemChaste includes
#include "ChemChasteExecutableHeaders.hpp"

#include "SchnackenbergCoupledPdeSystem.hpp"

#include <boost/lexical_cast.hpp>

#include "ExecutableSupport.hpp"

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>


void SetupSingletons();
void DestroySingletons();
void SetupAndRunSimulation(unsigned simulation_id,boost::program_options::variables_map& variables_map);
void WriteSimulationToChemChasteLog(std::string,unsigned);

int main(int argc, char *argv[])
{
    // This sets up PETSc and prints out copyright information, etc.
    ExecutableSupport::StartupWithoutShowingCopyright(&argc, &argv);

    

    int exit_code = ExecutableSupport::EXIT_OK;

    // You should put all the main code within a try-catch, to ensure that
    // you clean up PETSc before quitting.
    try
    {
        // begin time recording
        auto start = std::chrono::high_resolution_clock::now();

        // Declare a group of options that will be 
        // allowed only on command line
        boost::program_options::options_description general("General options");
        general.add_options()
        ("help", "produce help message")
        ("config", boost::program_options::value<std::string>()->default_value("/home/chaste/projects/ChemChaste/DataInput/ChemChasteConfig.txt"), "Config file")
        ("ID", boost::program_options::value<unsigned>()->default_value(0),"ID of the simulation (for output)")   
        ("simulation_type", boost::program_options::value<std::string>()->default_value("coupled_cell"),"Define the type of simualtion to run: domain_only, coupled_cell")
        ("simulation_timestep", boost::program_options::value<double>()->default_value(1/120),"Timestep length of the simulation (default: 1/120)")
        ("sampling_timestep", boost::program_options::value<double>()->default_value(1/120),"Timestep length of theresult sampling in the simulation (default: 1/120)")
          
        ;
        
        // Declare a group of options that will be 
        // allowed both on command line and in
        // config file
        boost::program_options::options_description config("Configuration");
        config.add_options()
        ("output_filename,F", boost::program_options::value<std::string>()->default_value("ChemChasteExecutable"),"Base filename to output the executable results (default: ChemChasteExecutable)")
        ("simulation_end_time", boost::program_options::value<double>()->default_value(10),"Temporal length of the simulation (default: 10)")
        ("domain_file_root",boost::program_options::value<std::string>()->default_value(""),"Absolute root of the domain folder (default: "")")
        ("domain_file",boost::program_options::value<std::string>()->default_value(""),"CSV file containing the domain topology (default: "")")
        ("domain_key_file",boost::program_options::value<std::string>()->default_value(""),"Key file relating labels in domain_file topology to subdomains (default: "")")
        ("domain_width",boost::program_options::value<double>()->default_value(0),"Width of the domain, to be used in domain scaling (default: 0)")
        ("domain_height",boost::program_options::value<double>()->default_value(0),"Height of the domain, to be used in domain scaling (default: 0)")
        ("scale_domain_width",boost::program_options::value<double>()->default_value(1),"Scale factor for the Width of the domain (default: 1)")
        ("scale_domain_height",boost::program_options::value<double>()->default_value(1),"Scale factor for the height of the domain (default: 1)")
        ("FE_mesh_step_size",boost::program_options::value<double>()->default_value(1),"Mesh step size for the FE (PDE) mesh in a coupled cell simulation (default: 1)")
        ("mesh_location",boost::program_options::value<std::string>()->default_value(""),"Where to overlap the FE (PDE) mesh and cell mesh in a coupled cell simualtion (default: "")")
        ("ode_file",boost::program_options::value<std::string>()->default_value(""),"CSV file containing the ode topology (default: "")")
        ("ode_key_file",boost::program_options::value<std::string>()->default_value(""),"Key file relating labels in ode_file topology to ode domain systems (default: "")")
        ("diffusion_database",boost::program_options::value<std::string>()->default_value(""),"CSV containing the standard diffusion rates of the species in each of the subdomains (default: "")")
        ("initial_conditions",boost::program_options::value<std::string>()->default_value(""),"CSV file containing the initial conditions of each species in each of the subdomains specified in domain_file (default: "")")
        ("boundary_conditions",boost::program_options::value<std::string>()->default_value(""),"CSV file containing the outer domain boundary conditions for each species (default: "")")
        ("cell_file_root",boost::program_options::value<std::string>()->default_value(""),"Absolute root of the cell folder (default: "")")
        ("cell_file",boost::program_options::value<std::string>()->default_value(""),"CSV file containing the cell layer topology (default: "")")
        ("cell_key_file",boost::program_options::value<std::string>()->default_value(""),"Key file relating labels cell_file topology to cell populations (default: "")")
        ("cell_initial_conditions",boost::program_options::value<std::string>()->default_value(""),"CSV file containing the initial cell internal concentrations (default: "")")
        ("membrane_property",boost::program_options::value<std::string>()->default_value(""),"Membrane file, reactions bound at the membrane (default: "")")
        ("cell_configuration",boost::program_options::value<std::string>()->default_value(""),"Cell configuration file (default: "")")
        ("transport_property",boost::program_options::value<std::string>()->default_value(""),"Transport file, reactions across the cell membrane (default: "")")
        ("environment_property",boost::program_options::value<std::string>()->default_value(""),"Environemnt file, preferred environment of the cell (default: "")")
        ("state_switching_property",boost::program_options::value<std::string>()->default_value(""),"State switching file (default: "")")
        ("srn_file",boost::program_options::value<std::string>()->default_value(""),"Subcellular reaction network file (default: "")")
        ("cell_cycle_file",boost::program_options::value<std::string>()->default_value(""),"Species threshold (default: "")")
        ("cell_division_rules",boost::program_options::value<std::string>()->default_value(""),"Cell division rules (default: "")")
        ("number_cells_across",boost::program_options::value<unsigned>()->default_value(1),"Width of honeycomb cell mesh (default: 1)")
        ("number_cells_high",boost::program_options::value<unsigned>()->default_value(1),"Height of honeycomb cell mesh (default: 1)")
        ("number_of_reaction_pdes",boost::program_options::value<unsigned>()->default_value(1),"Dimension of the pde system to solve, the number of unique chemical species diffusing in the bulk (default: 1)")
        ("spatial_dimensions",boost::program_options::value<unsigned>()->default_value(1),"The spatial dimension of the domain (default: 1)")
        ("FE_element_dimension",boost::program_options::value<unsigned>()->default_value(1),"The dimension of the finite Element mesh to be created and for whom the pdes are solved over (default: 1)")
        ("node_cutoff_length",boost::program_options::value<double>()->default_value(1.5),"Connectivity length of nodes only mesh (default: 1.5)")
        ("cell_mesh_origin",boost::program_options::value<double>()->default_value(-4.0),"Mesh overlapping point (default: -4.0)")
        ("linear_force_cutoff",boost::program_options::value<double>()->default_value(1.5),"Cut off length of linear spring force between cells (default: 1.5)")
        ("parameter_filename",boost::program_options::value<std::string>()->default_value(""),"Absolute root of the parameter file (default: "")")
        ;

        boost::program_options::options_description cmdline_options;
        cmdline_options.add(general);

        boost::program_options::options_description config_options;
        config_options.add(config).add(general);


       
        // define parse command line into variables_map
        boost::program_options::variables_map variables_map;
        boost::program_options::store(parse_command_line(argc, argv, cmdline_options), variables_map);
        
        // parse config file

        if (variables_map.count("config"))
        {
            std::cout<<"parse config file: "<<variables_map["config"].as<std::string>()<<" for run: "<<variables_map["ID"].as<unsigned>()<<std::endl;
            
            std::ifstream Config_File(variables_map["config"].as<std::string>());
            if(Config_File)
            {
                boost::program_options::store(parse_config_file(Config_File, config_options), variables_map);
            }
            else
            {
                std::cout<<"Error: config file not found: "<<variables_map["config"].as<std::string>()<<std::endl;
            }
        }
        boost::program_options::notify(variables_map);    

        // print help message if wanted
        if (variables_map.count("help"))
        {
            std::cout << setprecision(3) << cmdline_options << "\n";
            return 1;
        }

        SetupSingletons();

        SetupAndRunSimulation(variables_map["ID"].as<unsigned>(),variables_map);
        DestroySingletons();

        auto finish = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = finish - start;


        WriteSimulationToChemChasteLog(variables_map["config"].as<std::string>(),variables_map["ID"].as<unsigned>());


        std::cout<<"Simulation end: "<<variables_map["config"].as<std::string>()<<" for run: "<<variables_map["ID"].as<unsigned>()<<std::endl;
        std::cout << "Elapsed time: " << elapsed.count() << " s\n";
    }
    catch (const Exception& e)
    {
        ExecutableSupport::PrintError(e.GetMessage());
        exit_code = ExecutableSupport::EXIT_ERROR;
    }

}
//(a subfolder of tmp/[USERNAME]/testoutput)
void SetupSingletons()
{
    // Set up what the test suite would do
    SimulationTime::Instance()->SetStartTime(0.0);

    // we want every realisation of the simulation to be different!
    RandomNumberGenerator::Instance()->Reseed(time(NULL));
    CellPropertyRegistry::Instance()->Clear();
    CellId::ResetMaxCellId();
}

void DestroySingletons()
{
    /*
     * Finally, we must tidy up by destroying the {{{WntConcentration}}}
     * singleton object. This avoids memory leaks occurring.
     */

    // this is from the tearDown method of the test suite
    SimulationTime::Destroy();
    RandomNumberGenerator::Destroy();
}

void SetupAndRunSimulation(unsigned simulation_id, boost::program_options::variables_map& variables_map)
{
    // We would like the simulation id to be the random seed of the simulation
    // in order to ensure that each simulation is different.
    RandomNumberGenerator::Instance()->Reseed(simulation_id);
    
    // System properties
    const unsigned probDim =2; // need to set manually to the number of diffusive variables for the pde solver to solve
    const unsigned spaceDim=2;
    const unsigned elementDim=2;


    if(variables_map["simulation_type"].as<std::string>()=="coupled_cell")
    {

        if(variables_map["cell_file"].as<std::string>()!="")
        {
            // cell population topology is defined so read cell layer from files
   
            // run the domain field set up and parse files
            ChemicalDomainFieldForCellCoupling<elementDim,spaceDim,probDim>* p_Pde_field = new ChemicalDomainFieldForCellCoupling<elementDim,spaceDim,probDim>
                                    (   variables_map["domain_file_root"].as<std::string>(),
                                        variables_map["cell_file_root"].as<std::string>()+variables_map["cell_file"].as<std::string>(),
                                        variables_map["cell_file_root"].as<std::string>()+variables_map["cell_key_file"].as<std::string>(),
                                        variables_map["domain_file_root"].as<std::string>()+variables_map["domain_file"].as<std::string>(), 
                                        variables_map["domain_file_root"].as<std::string>()+variables_map["domain_key_file"].as<std::string>(), 
                                        variables_map["domain_file_root"].as<std::string>()+variables_map["ode_file"].as<std::string>(), 
                                        variables_map["domain_file_root"].as<std::string>()+variables_map["ode_key_file"].as<std::string>(), 
                                        variables_map["domain_file_root"].as<std::string>()+variables_map["diffusion_database"].as<std::string>(), 
                                        variables_map["domain_file_root"].as<std::string>()+variables_map["initial_conditions"].as<std::string>(), 
                                        variables_map["domain_file_root"].as<std::string>()+variables_map["boundary_conditions"].as<std::string>());
            
            TetrahedralMesh<elementDim,spaceDim>* p_cell_mesh = p_Pde_field->rGetCellMesh();

            NodesOnlyMesh<spaceDim> mesh;
            mesh.ConstructNodesWithoutMesh(*p_cell_mesh, variables_map["node_cutoff_length"].as<double>());

            std::vector<CellPtr> cells;
            // assume cell at each node in cell layer mesh
            std::string cell_label;
            unsigned numericalCellID;
            std::string cell_key;
            std::string given_cell_root;
            for (unsigned i=0; i<p_cell_mesh->GetNumNodes(); i++)
            {
                cell_label = p_Pde_field->GetCellLabelByIndex(i);
                cell_key = p_Pde_field->ReturnCellKeyFromCellLabel(cell_label);
                numericalCellID = p_Pde_field->ReturnUnsignedIDFromCellKeyString(cell_key);
                given_cell_root = variables_map["cell_file_root"].as<std::string>()+cell_key+"/";

                ChemicalCellFromFile* p_cell_reader = new ChemicalCellFromFile(
                                    given_cell_root+"SpeciesThreshold.csv", 
                                    given_cell_root+"Srn.txt",
                                    given_cell_root+"InitialCellConcentrations.csv",
                                    given_cell_root+"TransportReactions.txt",
                                    given_cell_root+"MembraneReactions.txt",
                                    numericalCellID,
                                    true                                                                    
                                    );
                p_cell_reader->SetCellTypeName(cell_key);
                cells.push_back(p_cell_reader -> GetCellPtr());
            }   


            NodeBasedCellPopulation<spaceDim> cell_population(mesh, cells);
            
            // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
            ChastePoint<spaceDim> lower(variables_map["cell_mesh_origin"].as<double>(),variables_map["cell_mesh_origin"].as<double>());
            ChastePoint<spaceDim> upper(10.0, 10.0);
            MAKE_PTR_ARGS(ChasteCuboid<spaceDim>, p_cuboid, (lower, upper));

            double feMeshStepSize = variables_map["FE_mesh_step_size"].as<double>();
            
            bool isCenterMesh = false;
            if(variables_map["mesh_location"].as<std::string>()=="Center")
            {
                isCenterMesh = true;
            }
            
            boost::shared_ptr<ParabolicBoxDomainPdeSystemModifier<elementDim,spaceDim,probDim>> p_pde_modifier(new ParabolicBoxDomainPdeSystemModifier<elementDim,spaceDim,probDim>(p_Pde_field, p_cuboid, feMeshStepSize, isCenterMesh));
            
            boost::shared_ptr<ChemicalTrackingModifier<spaceDim,elementDim>> p_chemical_tracking_modifier(new ChemicalTrackingModifier<spaceDim,elementDim>()); //= chemical_structure -> rGetPtrChemicalTrackingModifier();
            
            cell_population.AddCellWriter<CellIdWriter>();
            cell_population.AddCellWriter<CellAgesWriter>();
            cell_population.AddCellWriter<CellLocationIndexWriter>();
            cell_population.AddCellWriter<CellAnalyticsWriter>();
            //cell_population.AddCellWriter<CellStateWriter>();



            OffLatticeSimulation<spaceDim> simulator(cell_population);

            simulator.AddSimulationModifier(p_pde_modifier);

            simulator.AddSimulationModifier(p_chemical_tracking_modifier);
        
            simulator.SetOutputDirectory(variables_map["output_filename"].as<std::string>()+"/"+boost::lexical_cast<std::string>(simulation_id));
            simulator.SetEndTime(variables_map["simulation_end_time"].as<double>());
            simulator.SetDt(variables_map["simulation_timestep"].as<double>());

            MAKE_PTR(GeneralisedLinearSpringForce<spaceDim>, p_linear_force);
            p_linear_force->SetCutOffLength(variables_map["linear_force_cutoff"].as<double>());
            simulator.AddForce(p_linear_force);

            simulator.Solve();

        }
        else
        {
            // set up chemical domain field
            ChemicalDomainFieldForCellCoupling<elementDim,spaceDim,probDim>* p_Pde_field = new ChemicalDomainFieldForCellCoupling<elementDim,spaceDim,probDim>
                                    (   variables_map["domain_file_root"].as<std::string>(),
                                        variables_map["domain_file_root"].as<std::string>()+"",
                                        variables_map["domain_file_root"].as<std::string>()+"",
                                        variables_map["domain_file_root"].as<std::string>()+variables_map["domain_file"].as<std::string>(), 
                                        variables_map["domain_file_root"].as<std::string>()+variables_map["domain_key_file"].as<std::string>(), 
                                        variables_map["domain_file_root"].as<std::string>()+variables_map["ode_file"].as<std::string>(), 
                                        variables_map["domain_file_root"].as<std::string>()+variables_map["ode_key_file"].as<std::string>(), 
                                        variables_map["domain_file_root"].as<std::string>()+variables_map["diffusion_database"].as<std::string>(), 
                                        variables_map["domain_file_root"].as<std::string>()+variables_map["initial_conditions"].as<std::string>(), 
                                        variables_map["domain_file_root"].as<std::string>()+variables_map["boundary_conditions"].as<std::string>());
            


            // make mesh of cell with associated mesh based cell population
            HoneycombMeshGenerator generator(variables_map["number_cells_across"].as<unsigned>(),variables_map["number_cells_high"].as<unsigned>());    // Parameters are: cells across, cells up
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
                                    variables_map["cell_file_root"].as<std::string>()+variables_map["cell_cycle_file"].as<std::string>(), 
                                    variables_map["cell_file_root"].as<std::string>()+variables_map["srn_file"].as<std::string>(),
                                    variables_map["cell_file_root"].as<std::string>()+variables_map["cell_initial_conditions"].as<std::string>(),
                                    variables_map["cell_file_root"].as<std::string>()+variables_map["transport_property"].as<std::string>(),
                                    variables_map["cell_file_root"].as<std::string>()+variables_map["membrane_property"].as<std::string>()
                                    );

                if(IsFirstCell)
                {
                    chemicalCellSpeciesNames =  p_cell_reader->GetFullChemicalNamesVector();
                    IsFirstCell = false;
                }
                cells.push_back(p_cell_reader -> GetCellPtr());
            }

            // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
            ChastePoint<spaceDim> lower(variables_map["cell_mesh_origin"].as<double>(),variables_map["cell_mesh_origin"].as<double>());
            ChastePoint<spaceDim> upper(10.0, 10.0);
            MAKE_PTR_ARGS(ChasteCuboid<spaceDim>, p_cuboid, (lower, upper));

            double feMeshStepSize = variables_map["FE_mesh_step_size"].as<double>();
            
            bool isCenterMesh = false;
            if(variables_map["mesh_location"].as<std::string>()=="Center")
            {
                isCenterMesh = true;
            }
            
            boost::shared_ptr<ParabolicBoxDomainPdeSystemModifier<elementDim,spaceDim,probDim>> p_pde_modifier(new ParabolicBoxDomainPdeSystemModifier<elementDim,spaceDim,probDim>(p_Pde_field, p_cuboid, feMeshStepSize, isCenterMesh));
            
            boost::shared_ptr<ChemicalTrackingModifier<elementDim,spaceDim>> p_chemical_tracking_modifier(new ChemicalTrackingModifier<elementDim,spaceDim>());
            

            NodesOnlyMesh<spaceDim> mesh;

            mesh.ConstructNodesWithoutMesh(*p_mesh, variables_map["node_cutoff_length"].as<double>());

            NodeBasedCellPopulation<spaceDim> cell_population(mesh, cells);   

            // writers

            cell_population.AddCellWriter<CellIdWriter>();
            cell_population.AddCellWriter<CellAgesWriter>();
            cell_population.AddCellWriter<CellLocationIndexWriter>();

            OffLatticeSimulation<spaceDim> simulator(cell_population);

            simulator.AddSimulationModifier(p_pde_modifier);

            simulator.AddSimulationModifier(p_chemical_tracking_modifier);
            simulator.SetOutputDirectory(variables_map["output_filename"].as<std::string>()+"/"+boost::lexical_cast<std::string>(simulation_id));
            simulator.SetEndTime(variables_map["simulation_end_time"].as<double>());
            simulator.SetDt(variables_map["simulation_timestep"].as<double>());

            MAKE_PTR(GeneralisedLinearSpringForce<spaceDim>, p_linear_force);
            p_linear_force->SetCutOffLength(variables_map["linear_force_cutoff"].as<double>());
            simulator.AddForce(p_linear_force);

    
            simulator.Solve();
        }

    }
    else if(variables_map["simulation_type"].as<std::string>()=="complex_cell")
    {

        if(variables_map["cell_file"].as<std::string>()!="")
        {
            // cell population topology is defined so read cell layer from files
            std::cout<<"---------------------------here"<<std::endl;
            // run the domain field set up and parse files
            ChemicalDomainFieldForCellCoupling<elementDim,spaceDim,probDim>* p_Pde_field = new ChemicalDomainFieldForCellCoupling<elementDim,spaceDim,probDim>
                                    (   variables_map["domain_file_root"].as<std::string>(),
                                        variables_map["cell_file_root"].as<std::string>()+variables_map["cell_file"].as<std::string>(),
                                        variables_map["cell_file_root"].as<std::string>()+variables_map["cell_key_file"].as<std::string>(),
                                        variables_map["domain_file_root"].as<std::string>()+variables_map["domain_file"].as<std::string>(), 
                                        variables_map["domain_file_root"].as<std::string>()+variables_map["domain_key_file"].as<std::string>(), 
                                        variables_map["domain_file_root"].as<std::string>()+variables_map["ode_file"].as<std::string>(), 
                                        variables_map["domain_file_root"].as<std::string>()+variables_map["ode_key_file"].as<std::string>(), 
                                        variables_map["domain_file_root"].as<std::string>()+variables_map["diffusion_database"].as<std::string>(), 
                                        variables_map["domain_file_root"].as<std::string>()+variables_map["initial_conditions"].as<std::string>(), 
                                        variables_map["domain_file_root"].as<std::string>()+variables_map["boundary_conditions"].as<std::string>());
            
            TetrahedralMesh<elementDim,spaceDim>* p_cell_mesh = p_Pde_field->rGetCellMesh();

            NodesOnlyMesh<spaceDim> mesh;
            mesh.ConstructNodesWithoutMesh(*p_cell_mesh, variables_map["node_cutoff_length"].as<double>());


            std::vector<CellPtr> cells;
            // assume cell at each node in cell layer mesh
            std::string cell_label;
            unsigned numericalCellID;
            std::string cell_key;
            std::string given_cell_root;
            for (unsigned i=0; i<p_cell_mesh->GetNumNodes(); i++)
            {
                cell_label = p_Pde_field->GetCellLabelByIndex(i);
                cell_key = p_Pde_field->ReturnCellKeyFromCellLabel(cell_label);
                numericalCellID = p_Pde_field->ReturnUnsignedIDFromCellKeyString(cell_key);
                given_cell_root = variables_map["cell_file_root"].as<std::string>()+cell_key+"/";

                ComplexCellFromFile* p_cell_reader = new ComplexCellFromFile(
                                    given_cell_root,
                                    "SpeciesThreshold.csv", 
                                    "SpeciesDivisionRules.csv", 
                                    "Srn.txt",
                                    "InitialCellConcentrations.csv",
                                    "TransportReactions.txt",
                                    "MembraneReactions.txt",
                                    "",
                                    numericalCellID,
                                    true                                                                    
                                    );
                p_cell_reader->SetCellTypeName(cell_key);
                p_cell_reader->SetUp();

                                    
                std::cout<<"Division file root: "<<given_cell_root+"SpeciesDivisionRules.csv"<<std::endl;
                std::cout<<"test for cell properties"<<std::endl;
                std::cout<<"has membrane: "<<p_cell_reader -> GetCellPtr() -> rGetCellPropertyCollection().HasProperty<MembraneCellProperty>()<<std::endl;
std::cout<<"has transport: "<<p_cell_reader -> GetCellPtr() -> rGetCellPropertyCollection().HasProperty<TransportCellProperty>()<<std::endl;
         
                boost::static_pointer_cast<CellAnalyticsProperty>(p_cell_reader -> GetCellPtr() -> rGetCellPropertyCollection().GetPropertiesType<CellAnalyticsProperty>().GetProperty()) -> SetPopulationCellTypeNames(p_Pde_field -> GetCellKeyVector());

                cells.push_back(p_cell_reader -> GetCellPtr());
            }   


            NodeBasedCellPopulation<spaceDim> cell_population(mesh, cells);
            
            // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
            ChastePoint<spaceDim> lower(variables_map["cell_mesh_origin"].as<double>(),variables_map["cell_mesh_origin"].as<double>());
            ChastePoint<spaceDim> upper(10.0, 10.0);
            MAKE_PTR_ARGS(ChasteCuboid<spaceDim>, p_cuboid, (lower, upper));


            
            double feMeshStepSize = variables_map["FE_mesh_step_size"].as<double>();
            
            bool isCenterMesh = false;
            //if(variables_map["mesh_location"].as<std::string>()=="Center")
            //{
                isCenterMesh = true;
            //}
            
            boost::shared_ptr<ParabolicBoxDomainPdeSystemModifier<elementDim,spaceDim,probDim>> p_pde_modifier(new ParabolicBoxDomainPdeSystemModifier<elementDim,spaceDim,probDim>(p_Pde_field, p_cuboid, feMeshStepSize, isCenterMesh));
            
            boost::shared_ptr<ChemicalTrackingModifier<spaceDim,elementDim>> p_chemical_tracking_modifier(new ChemicalTrackingModifier<spaceDim,elementDim>()); //= chemical_structure -> rGetPtrChemicalTrackingModifier();
            
            boost::shared_ptr<CellNeighbourTrackingModifier<spaceDim,elementDim>> p_neighbour_tracking_modifier(new CellNeighbourTrackingModifier<spaceDim,elementDim>(p_Pde_field -> GetCellKeyVector()));

            mesh.Translate(p_pde_modifier ->GetMeshOffset());

            for(unsigned i=0; i<p_Pde_field -> GetCellKeyVector().size(); i++)
            {
                std::cout<<"Key: "<<p_Pde_field -> GetCellKeyVector()[i]<<std::endl;
            }

            cell_population.AddCellWriter<CellIdWriter>();
            cell_population.AddCellWriter<CellAgesWriter>();
            cell_population.AddCellWriter<CellLocationIndexWriter>();
            cell_population.AddCellWriter<CellAnalyticsWriter>();

            cell_population.AddCellWriter<CellNeighbourDiversityWriter>();
            cell_population.AddCellWriter<CellNeighbourLocalMoranWriter>();
            cell_population.AddCellWriter<CellNeighbourGetisOrdWriter>();
            cell_population.AddCellPopulationCountWriter<PopulationDiversityWriter>();

            OffLatticeSimulation<spaceDim> simulator(cell_population);

            simulator.AddSimulationModifier(p_pde_modifier);

            simulator.AddSimulationModifier(p_chemical_tracking_modifier);

            simulator.AddSimulationModifier(p_neighbour_tracking_modifier);
        
            simulator.SetOutputDirectory(variables_map["output_filename"].as<std::string>()+"/"+boost::lexical_cast<std::string>(simulation_id));
            simulator.SetEndTime(variables_map["simulation_end_time"].as<double>());
            simulator.SetDt(variables_map["simulation_timestep"].as<double>());

            MAKE_PTR(GeneralisedLinearSpringForce<spaceDim>, p_linear_force);
            p_linear_force->SetCutOffLength(variables_map["linear_force_cutoff"].as<double>());
            simulator.AddForce(p_linear_force);

            std::vector<unsigned>  domainDimensions = p_Pde_field->GetMeshDimensions();
            std::vector<double> domainMax;
            std::vector<double> domainMin;

            for(unsigned dim=0; dim<domainDimensions.size(); dim++)
            {
                domainMax.push_back(domainDimensions[dim]);
                domainMin.push_back(0.0);
            }

            boost::shared_ptr<FixedCellBoundaryCondition<elementDim,spaceDim>> p_cell_BC(new FixedCellBoundaryCondition<elementDim,spaceDim>(&cell_population, domainMax, domainMin));

            simulator.AddCellPopulationBoundaryCondition(p_cell_BC);


            simulator.Solve();

        }
        else
        {
            
            // set up chemical domain field
            ChemicalDomainFieldForCellCoupling<elementDim,spaceDim,probDim>* p_Pde_field = new ChemicalDomainFieldForCellCoupling<elementDim,spaceDim,probDim>
                                    (   variables_map["domain_file_root"].as<std::string>(),
                                        variables_map["domain_file_root"].as<std::string>()+"",
                                        variables_map["domain_file_root"].as<std::string>()+"",
                                        variables_map["domain_file_root"].as<std::string>()+variables_map["domain_file"].as<std::string>(), 
                                        variables_map["domain_file_root"].as<std::string>()+variables_map["domain_key_file"].as<std::string>(), 
                                        variables_map["domain_file_root"].as<std::string>()+variables_map["ode_file"].as<std::string>(), 
                                        variables_map["domain_file_root"].as<std::string>()+variables_map["ode_key_file"].as<std::string>(), 
                                        variables_map["domain_file_root"].as<std::string>()+variables_map["diffusion_database"].as<std::string>(), 
                                        variables_map["domain_file_root"].as<std::string>()+variables_map["initial_conditions"].as<std::string>(), 
                                        variables_map["domain_file_root"].as<std::string>()+variables_map["boundary_conditions"].as<std::string>());
            


            // make mesh of cell with associated mesh based cell population
            HoneycombMeshGenerator generator(variables_map["number_cells_across"].as<unsigned>(),variables_map["number_cells_high"].as<unsigned>());    // Parameters are: cells across, cells up
            MutableMesh<elementDim,spaceDim>* p_mesh = generator.GetMesh();


            std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

            std::vector<CellPtr> cells;

            std::vector<std::string> chemicalCellSpeciesNames;

            // assume cell at each node in cell layer mesh
            bool IsFirstCell = true;
            for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
            {
                // provide each cell with a transport cell property and membrane property, cell cycle, wild type states

                ComplexCellFromFile* p_cell_reader = new ComplexCellFromFile(
                                    variables_map["cell_file_root"].as<std::string>(),
                                    variables_map["cell_cycle_file"].as<std::string>(), 
                                    variables_map["cell_division_rules"].as<std::string>(),
                                    variables_map["srn_file"].as<std::string>(),
                                    variables_map["cell_initial_conditions"].as<std::string>(),
                                    variables_map["transport_property"].as<std::string>(),
                                    variables_map["membrane_property"].as<std::string>(),
                                    variables_map["environment_property"].as<std::string>()
                                    );
                p_cell_reader->SetUp();

                if(IsFirstCell)
                {
                    chemicalCellSpeciesNames =  p_cell_reader->GetFullChemicalNamesVector();
                    IsFirstCell = false;
                }
                cells.push_back(p_cell_reader -> GetCellPtr());
            }

            // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
            ChastePoint<spaceDim> lower(variables_map["cell_mesh_origin"].as<double>(),variables_map["cell_mesh_origin"].as<double>());
            ChastePoint<spaceDim> upper(10.0, 10.0);
            MAKE_PTR_ARGS(ChasteCuboid<spaceDim>, p_cuboid, (lower, upper));

            
            
            double feMeshStepSize = variables_map["FE_mesh_step_size"].as<double>();
            
            bool isCenterMesh = false;
            if(variables_map["mesh_location"].as<std::string>()=="Center")
            {
                isCenterMesh = true;
            }
            
            boost::shared_ptr<ParabolicBoxDomainPdeSystemModifier<elementDim,spaceDim,probDim>> p_pde_modifier(new ParabolicBoxDomainPdeSystemModifier<elementDim,spaceDim,probDim>(p_Pde_field, p_cuboid, feMeshStepSize, isCenterMesh));
            
         
            boost::shared_ptr<ChemicalTrackingModifier<elementDim,spaceDim>> p_chemical_tracking_modifier(new ChemicalTrackingModifier<elementDim,spaceDim>());
            

            NodesOnlyMesh<spaceDim> mesh;

            mesh.ConstructNodesWithoutMesh(*p_mesh, variables_map["node_cutoff_length"].as<double>());

            NodeBasedCellPopulation<spaceDim> cell_population(mesh, cells);   

            // writers

            cell_population.AddCellWriter<CellIdWriter>();
            cell_population.AddCellWriter<CellAgesWriter>();
            cell_population.AddCellWriter<CellLocationIndexWriter>();

            OffLatticeSimulation<spaceDim> simulator(cell_population);

            simulator.AddSimulationModifier(p_pde_modifier);

            simulator.AddSimulationModifier(p_chemical_tracking_modifier);
            simulator.SetOutputDirectory(variables_map["output_filename"].as<std::string>()+"/"+boost::lexical_cast<std::string>(simulation_id));
            simulator.SetEndTime(variables_map["simulation_end_time"].as<double>());
            simulator.SetDt(variables_map["simulation_timestep"].as<double>());

            MAKE_PTR(GeneralisedLinearSpringForce<spaceDim>, p_linear_force);
            p_linear_force->SetCutOffLength(variables_map["linear_force_cutoff"].as<double>());
            simulator.AddForce(p_linear_force);

            std::vector<unsigned>  domainDimensions = p_Pde_field->GetMeshDimensions();
            std::vector<double> domainMax;
            std::vector<double> domainMin;

            for(unsigned dim=0; dim<domainDimensions.size(); dim++)
            {
                domainMax.push_back(domainDimensions[dim]);
                domainMin.push_back(0.0);
            }

            boost::shared_ptr<FixedCellBoundaryCondition<elementDim,spaceDim>> p_cell_BC(new FixedCellBoundaryCondition<elementDim,spaceDim>(&cell_population, domainMax, domainMin));

            simulator.AddCellPopulationBoundaryCondition(p_cell_BC);
    
            simulator.Solve();
        }

    }
    else if(variables_map["simulation_type"].as<std::string>()=="environment_cell")
    {

        if(variables_map["cell_file"].as<std::string>()!="")
        {
            // cell population topology is defined so read cell layer from files
            std::cout<<"---------------------------here"<<std::endl;
            // run the domain field set up and parse files
            ChemicalDomainFieldForCellCoupling<elementDim,spaceDim,probDim>* p_Pde_field = new ChemicalDomainFieldForCellCoupling<elementDim,spaceDim,probDim>
                                    (   variables_map["domain_file_root"].as<std::string>(),
                                        variables_map["cell_file_root"].as<std::string>()+variables_map["cell_file"].as<std::string>(),
                                        variables_map["cell_file_root"].as<std::string>()+variables_map["cell_key_file"].as<std::string>(),
                                        variables_map["domain_file_root"].as<std::string>()+variables_map["domain_file"].as<std::string>(), 
                                        variables_map["domain_file_root"].as<std::string>()+variables_map["domain_key_file"].as<std::string>(), 
                                        variables_map["domain_file_root"].as<std::string>()+variables_map["ode_file"].as<std::string>(), 
                                        variables_map["domain_file_root"].as<std::string>()+variables_map["ode_key_file"].as<std::string>(), 
                                        variables_map["domain_file_root"].as<std::string>()+variables_map["diffusion_database"].as<std::string>(), 
                                        variables_map["domain_file_root"].as<std::string>()+variables_map["initial_conditions"].as<std::string>(), 
                                        variables_map["domain_file_root"].as<std::string>()+variables_map["boundary_conditions"].as<std::string>());
            
            TetrahedralMesh<elementDim,spaceDim>* p_cell_mesh = p_Pde_field->rGetCellMesh();

            NodesOnlyMesh<spaceDim> mesh;
            mesh.ConstructNodesWithoutMesh(*p_cell_mesh, variables_map["node_cutoff_length"].as<double>());

            std::vector<CellPtr> cells;
            // assume cell at each node in cell layer mesh
            std::string cell_label;
            unsigned numericalCellID;
            std::string cell_key;
            std::string given_cell_root;
            for (unsigned i=0; i<p_cell_mesh->GetNumNodes(); i++)
            {
                cell_label = p_Pde_field->GetCellLabelByIndex(i);
                cell_key = p_Pde_field->ReturnCellKeyFromCellLabel(cell_label);
                numericalCellID = p_Pde_field->ReturnUnsignedIDFromCellKeyString(cell_key);
                given_cell_root = variables_map["cell_file_root"].as<std::string>()+cell_key+"/";

                ComplexCellFromFile* p_cell_reader = new ComplexCellFromFile(
                                    given_cell_root,
                                    "SpeciesThreshold.csv", 
                                    "SpeciesDivisionRules.csv", 
                                    "Srn.txt",
                                    "InitialCellConcentrations.csv",
                                    "TransportReactions.txt",
                                    "MembraneReactions.txt",
                                    "Environment.csv",
                                    numericalCellID,
                                    true                                                                    
                                    );
                p_cell_reader->SetCellTypeName(cell_key);
                p_cell_reader->SetUp();

                                    
                std::cout<<"Division file root: "<<given_cell_root+"SpeciesDivisionRules.csv"<<std::endl;
                std::cout<<"test for cell properties"<<std::endl;
                std::cout<<"has membrane: "<<p_cell_reader -> GetCellPtr() -> rGetCellPropertyCollection().HasProperty<MembraneCellProperty>()<<std::endl;
std::cout<<"has transport: "<<p_cell_reader -> GetCellPtr() -> rGetCellPropertyCollection().HasProperty<TransportCellProperty>()<<std::endl;

                boost::static_pointer_cast<CellAnalyticsProperty>(p_cell_reader -> GetCellPtr() -> rGetCellPropertyCollection().GetPropertiesType<CellAnalyticsProperty>().GetProperty()) -> SetPopulationCellTypeNames(p_Pde_field -> GetCellKeyVector());

                cells.push_back(p_cell_reader -> GetCellPtr());
            }   


            NodeBasedCellPopulation<spaceDim> cell_population(mesh, cells);
            
            // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
            ChastePoint<spaceDim> lower(variables_map["cell_mesh_origin"].as<double>(),variables_map["cell_mesh_origin"].as<double>());
            ChastePoint<spaceDim> upper(10.0, 10.0);
            MAKE_PTR_ARGS(ChasteCuboid<spaceDim>, p_cuboid, (lower, upper));


            
            double feMeshStepSize = variables_map["FE_mesh_step_size"].as<double>();
            
            bool isCenterMesh = false;
            //if(variables_map["mesh_location"].as<std::string>()=="Center")
            //{
                isCenterMesh = true;
            //}
            
            boost::shared_ptr<ParabolicBoxDomainPdeSystemModifier<elementDim,spaceDim,probDim>> p_pde_modifier(new ParabolicBoxDomainPdeSystemModifier<elementDim,spaceDim,probDim>(p_Pde_field, p_cuboid, feMeshStepSize, isCenterMesh));
            
            boost::shared_ptr<ChemicalTrackingModifier<spaceDim,elementDim>> p_chemical_tracking_modifier(new ChemicalTrackingModifier<spaceDim,elementDim>()); //= chemical_structure -> rGetPtrChemicalTrackingModifier();
            
            boost::shared_ptr<CellNeighbourTrackingModifier<spaceDim,elementDim>> p_neighbour_tracking_modifier(new CellNeighbourTrackingModifier<spaceDim,elementDim>(p_Pde_field -> GetCellKeyVector()));

            mesh.Translate(p_pde_modifier ->GetMeshOffset());

            for(unsigned i=0; i<p_Pde_field -> GetCellKeyVector().size(); i++)
            {
                std::cout<<"Key: "<<p_Pde_field -> GetCellKeyVector()[i]<<std::endl;
            }

            cell_population.AddCellWriter<CellIdWriter>();
            cell_population.AddCellWriter<CellAgesWriter>();
            cell_population.AddCellWriter<CellLocationIndexWriter>();
            cell_population.AddCellWriter<CellAnalyticsWriter>();

            cell_population.AddCellWriter<CellNeighbourDiversityWriter>();
            cell_population.AddCellWriter<CellNeighbourLocalMoranWriter>();
            cell_population.AddCellWriter<CellNeighbourGetisOrdWriter>();
            cell_population.AddCellPopulationCountWriter<PopulationDiversityWriter>();

            OffLatticeSimulation<spaceDim> simulator(cell_population);

            simulator.AddSimulationModifier(p_pde_modifier);

            simulator.AddSimulationModifier(p_chemical_tracking_modifier);

            simulator.AddSimulationModifier(p_neighbour_tracking_modifier);
        
            simulator.SetOutputDirectory(variables_map["output_filename"].as<std::string>()+"/"+boost::lexical_cast<std::string>(simulation_id));
            simulator.SetEndTime(variables_map["simulation_end_time"].as<double>());
            simulator.SetDt(variables_map["simulation_timestep"].as<double>());

            MAKE_PTR(GeneralisedLinearSpringForce<spaceDim>, p_linear_force);
            p_linear_force->SetCutOffLength(variables_map["linear_force_cutoff"].as<double>());
            simulator.AddForce(p_linear_force);

            std::vector<unsigned>  domainDimensions = p_Pde_field->GetMeshDimensions();
            std::vector<double> domainMax;
            std::vector<double> domainMin;

            for(unsigned dim=0; dim<domainDimensions.size(); dim++)
            {
                domainMax.push_back(domainDimensions[dim]);
                domainMin.push_back(0.0);
            }

            boost::shared_ptr<FixedCellBoundaryCondition<elementDim,spaceDim>> p_cell_BC(new FixedCellBoundaryCondition<elementDim,spaceDim>(&cell_population, domainMax, domainMin));

            simulator.AddCellPopulationBoundaryCondition(p_cell_BC);

            simulator.Solve();

        }
        else
        {
            
            // set up chemical domain field
            ChemicalDomainFieldForCellCoupling<elementDim,spaceDim,probDim>* p_Pde_field = new ChemicalDomainFieldForCellCoupling<elementDim,spaceDim,probDim>
                                    (   variables_map["domain_file_root"].as<std::string>(),
                                        variables_map["domain_file_root"].as<std::string>()+"",
                                        variables_map["domain_file_root"].as<std::string>()+"",
                                        variables_map["domain_file_root"].as<std::string>()+variables_map["domain_file"].as<std::string>(), 
                                        variables_map["domain_file_root"].as<std::string>()+variables_map["domain_key_file"].as<std::string>(), 
                                        variables_map["domain_file_root"].as<std::string>()+variables_map["ode_file"].as<std::string>(), 
                                        variables_map["domain_file_root"].as<std::string>()+variables_map["ode_key_file"].as<std::string>(), 
                                        variables_map["domain_file_root"].as<std::string>()+variables_map["diffusion_database"].as<std::string>(), 
                                        variables_map["domain_file_root"].as<std::string>()+variables_map["initial_conditions"].as<std::string>(), 
                                        variables_map["domain_file_root"].as<std::string>()+variables_map["boundary_conditions"].as<std::string>());
            


            // make mesh of cell with associated mesh based cell population
            HoneycombMeshGenerator generator(variables_map["number_cells_across"].as<unsigned>(),variables_map["number_cells_high"].as<unsigned>());    // Parameters are: cells across, cells up
            MutableMesh<elementDim,spaceDim>* p_mesh = generator.GetMesh();


            std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

            std::vector<CellPtr> cells;

            std::vector<std::string> chemicalCellSpeciesNames;

            // assume cell at each node in cell layer mesh
            bool IsFirstCell = true;
            for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
            {
                // provide each cell with a transport cell property and membrane property, cell cycle, wild type states

                ComplexCellFromFile* p_cell_reader = new ComplexCellFromFile(
                                    variables_map["cell_file_root"].as<std::string>(),
                                    variables_map["cell_cycle_file"].as<std::string>(), 
                                    variables_map["cell_division_rules"].as<std::string>(),
                                    variables_map["srn_file"].as<std::string>(),
                                    variables_map["cell_initial_conditions"].as<std::string>(),
                                    variables_map["transport_property"].as<std::string>(),
                                    variables_map["membrane_property"].as<std::string>(),
                                    variables_map["environment_property"].as<std::string>()
                                    );
                p_cell_reader->SetUp();

                if(IsFirstCell)
                {
                    chemicalCellSpeciesNames =  p_cell_reader->GetFullChemicalNamesVector();
                    IsFirstCell = false;
                }
                cells.push_back(p_cell_reader -> GetCellPtr());
            }

            NodesOnlyMesh<spaceDim> mesh;

            mesh.ConstructNodesWithoutMesh(*p_mesh, variables_map["node_cutoff_length"].as<double>());

            NodeBasedCellPopulation<spaceDim> cell_population(mesh, cells);

            // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
            ChastePoint<spaceDim> lower(variables_map["cell_mesh_origin"].as<double>(),variables_map["cell_mesh_origin"].as<double>());
            ChastePoint<spaceDim> upper(10.0, 10.0);
            MAKE_PTR_ARGS(ChasteCuboid<spaceDim>, p_cuboid, (lower, upper));
            
            double feMeshStepSize = variables_map["FE_mesh_step_size"].as<double>();
            
            bool isCenterMesh = false;
            //if(variables_map["mesh_location"].as<std::string>()=="Center")
            //{
                isCenterMesh = true;
            //}
            
            boost::shared_ptr<ParabolicBoxDomainPdeSystemModifier<elementDim,spaceDim,probDim>> p_pde_modifier(new ParabolicBoxDomainPdeSystemModifier<elementDim,spaceDim,probDim>(p_Pde_field, p_cuboid, feMeshStepSize, isCenterMesh));
            
            boost::shared_ptr<ChemicalTrackingModifier<spaceDim,elementDim>> p_chemical_tracking_modifier(new ChemicalTrackingModifier<spaceDim,elementDim>()); //= chemical_structure -> rGetPtrChemicalTrackingModifier();
            
            boost::shared_ptr<CellNeighbourTrackingModifier<spaceDim,elementDim>> p_neighbour_tracking_modifier(new CellNeighbourTrackingModifier<spaceDim,elementDim>(p_Pde_field -> GetCellKeyVector()));

            mesh.Translate(p_pde_modifier ->GetMeshOffset());

            for(unsigned i=0; i<p_Pde_field -> GetCellKeyVector().size(); i++)
            {
                std::cout<<"Key: "<<p_Pde_field -> GetCellKeyVector()[i]<<std::endl;
            }

            cell_population.AddCellWriter<CellIdWriter>();
            cell_population.AddCellWriter<CellAgesWriter>();
            cell_population.AddCellWriter<CellLocationIndexWriter>();
            cell_population.AddCellWriter<CellAnalyticsWriter>();

            cell_population.AddCellWriter<CellNeighbourDiversityWriter>();
            cell_population.AddCellWriter<CellNeighbourLocalMoranWriter>();
            cell_population.AddCellWriter<CellNeighbourGetisOrdWriter>();
            cell_population.AddCellPopulationCountWriter<PopulationDiversityWriter>();

            OffLatticeSimulation<spaceDim> simulator(cell_population);

            simulator.AddSimulationModifier(p_pde_modifier);

            simulator.AddSimulationModifier(p_chemical_tracking_modifier);

            simulator.AddSimulationModifier(p_neighbour_tracking_modifier);


            simulator.AddSimulationModifier(p_chemical_tracking_modifier);
            simulator.SetOutputDirectory(variables_map["output_filename"].as<std::string>()+"/"+boost::lexical_cast<std::string>(simulation_id));
            simulator.SetEndTime(variables_map["simulation_end_time"].as<double>());
            simulator.SetDt(variables_map["simulation_timestep"].as<double>());

            MAKE_PTR(GeneralisedLinearSpringForce<spaceDim>, p_linear_force);
            p_linear_force->SetCutOffLength(variables_map["linear_force_cutoff"].as<double>());
            simulator.AddForce(p_linear_force);

            std::vector<unsigned>  domainDimensions = p_Pde_field->GetMeshDimensions();
            std::vector<double> domainMax;
            std::vector<double> domainMin;

            for(unsigned dim=0; dim<domainDimensions.size(); dim++)
            {
                domainMax.push_back(domainDimensions[dim]);
                domainMin.push_back(0.0);
            }

            boost::shared_ptr<FixedCellBoundaryCondition<elementDim,spaceDim>> p_cell_BC(new FixedCellBoundaryCondition<elementDim,spaceDim>(&cell_population, domainMax, domainMin));

            simulator.AddCellPopulationBoundaryCondition(p_cell_BC);
    
            simulator.Solve();
        }

    }
    /*
    else if(variables_map["simulation_type"].as<std::string>()=="custom_cell")
    {

        if(variables_map["cell_file"].as<std::string>()!="")
        {
            // cell population topology is defined so read cell layer from files
   
            // run the domain field set up and parse files
            ChemicalDomainFieldForCellCoupling<elementDim,spaceDim,probDim>* p_Pde_field = new ChemicalDomainFieldForCellCoupling<elementDim,spaceDim,probDim>
                                    (   variables_map["domain_file_root"].as<std::string>(),
                                        variables_map["cell_file_root"].as<std::string>()+variables_map["cell_file"].as<std::string>(),
                                        variables_map["cell_file_root"].as<std::string>()+variables_map["cell_key_file"].as<std::string>(),
                                        variables_map["domain_file_root"].as<std::string>()+variables_map["domain_file"].as<std::string>(), 
                                        variables_map["domain_file_root"].as<std::string>()+variables_map["domain_key_file"].as<std::string>(), 
                                        variables_map["domain_file_root"].as<std::string>()+variables_map["ode_file"].as<std::string>(), 
                                        variables_map["domain_file_root"].as<std::string>()+variables_map["ode_key_file"].as<std::string>(), 
                                        variables_map["domain_file_root"].as<std::string>()+variables_map["diffusion_database"].as<std::string>(), 
                                        variables_map["domain_file_root"].as<std::string>()+variables_map["initial_conditions"].as<std::string>(), 
                                        variables_map["domain_file_root"].as<std::string>()+variables_map["boundary_conditions"].as<std::string>());
            
            TetrahedralMesh<elementDim,spaceDim>* p_cell_mesh = p_Pde_field->rGetCellMesh();

            NodesOnlyMesh<spaceDim> mesh;
            mesh.ConstructNodesWithoutMesh(*p_cell_mesh, variables_map["node_cutoff_length"].as<double>());

            std::vector<CellPtr> cells;
            // assume cell at each node in cell layer mesh
            std::string cell_label;
            unsigned numericalCellID;
            std::string cell_key;
            std::string given_cell_root;
            for (unsigned i=0; i<p_cell_mesh->GetNumNodes(); i++)
            {
                cell_label = p_Pde_field->GetCellLabelByIndex(i);
                cell_key = p_Pde_field->ReturnCellKeyFromCellLabel(cell_label);
                numericalCellID = p_Pde_field->ReturnUnsignedIDFromCellKeyString(cell_key);
                given_cell_root = variables_map["cell_file_root"].as<std::string>()+cell_key+"/";

                CustomCellFromFile* p_cell_reader = new CustomCellFromFile(
                                    given_cell_root,
                                    "SpeciesThreshold.csv", 
                                    "SpeciesDivisionRules.csv", 
                                    "Srn.txt",
                                    "InitialCellConcentrations.csv",
                                    "TransportReactions.txt",
                                    "MembraneReactions.txt",
                                    "Environment.csv",
                                    "StateSwitches.csv",
                                    "cell_configuration.txt",
                                    numericalCellID,
                                    true                                                                    
                                    );
                p_cell_reader->SetCellTypeName(cell_key);
                p_cell_reader->SetUp();

                cells.push_back(p_cell_reader -> GetCellPtr());
            }   


            NodeBasedCellPopulation<spaceDim> cell_population(mesh, cells);
            
            // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
            ChastePoint<spaceDim> lower(variables_map["cell_mesh_origin"].as<double>(),variables_map["cell_mesh_origin"].as<double>());
            ChastePoint<spaceDim> upper(10.0, 10.0);
            MAKE_PTR_ARGS(ChasteCuboid<spaceDim>, p_cuboid, (lower, upper));

            p_Pde_field ->FeMeshScaling(    variables_map["domain_width"].as<double>(),
                                            variables_map["domain_height"].as<double>(),
                                            variables_map["scale_domain_width"].as<double>(),
                                            variables_map["scale_domain_height"].as<double>()
                                            );
            
            double feMeshStepSize = variables_map["FE_mesh_step_size"].as<double>();
            
            bool isCenterMesh = false;
            if(variables_map["mesh_location"].as<std::string>()=="Center")
            {
                isCenterMesh = true;
            }
            
            boost::shared_ptr<ParabolicBoxDomainPdeSystemModifier<elementDim,spaceDim,probDim>> p_pde_modifier(new ParabolicBoxDomainPdeSystemModifier<elementDim,spaceDim,probDim>(p_Pde_field, p_cuboid, feMeshStepSize, isCenterMesh));
            
         
            boost::shared_ptr<ChemicalTrackingPreSwitchModifier<elementDim,spaceDim>> p_chemical_tracking_pre_switch_modifier(new ChemicalTrackingPreSwitchModifier<elementDim,spaceDim>());

            boost::shared_ptr<CellStateTrackingModifier<elementDim,spaceDim>> p_cell_tracking_modifier(new CellStateTrackingModifier<elementDim,spaceDim>());
            
            boost::shared_ptr<ChemicalTrackingPostSwitchModifier<elementDim,spaceDim>> p_chemical_tracking_post_switch_modifier(new ChemicalTrackingPostSwitchModifier<elementDim,spaceDim>());


            NodesOnlyMesh<spaceDim> mesh;

            mesh.ConstructNodesWithoutMesh(*p_mesh, variables_map["node_cutoff_length"].as<double>());

            NodeBasedCellPopulation<spaceDim> cell_population(mesh, cells);   

            // writers

            cell_population.AddCellWriter<CellIdWriter>();
            cell_population.AddCellWriter<CellAgesWriter>();
            cell_population.AddCellWriter<CellLocationIndexWriter>();
            cell_population.AddCellWriter<CellAnalyticsWriter>();
            cell_population.AddCellWriter<CellStateWriter>();

            OffLatticeSimulation<spaceDim> simulator(cell_population);

            simulator.AddSimulationModifier(p_pde_modifier);

            simulator.AddSimulationModifier(p_chemical_tracking_pre_switch_modifier);
            simulator.AddSimulationModifier(p_cell_tracking_modifier);
            simulator.AddSimulationModifier(p_chemical_tracking_post_switch_modifier);

            simulator.SetOutputDirectory(variables_map["output_filename"].as<std::string>()+"/"+boost::lexical_cast<std::string>(simulation_id));
            simulator.SetEndTime(variables_map["simulation_end_time"].as<double>());
            simulator.SetDt(variables_map["simulation_timestep"].as<double>());


            MAKE_PTR(GeneralisedLinearSpringForce<spaceDim>, p_linear_force);
            p_linear_force->SetCutOffLength(variables_map["linear_force_cutoff"].as<double>());
            simulator.AddForce(p_linear_force);

            std::vector<unsigned>  domainDimensions = p_Pde_field->GetMeshDimensions();
            std::vector<double> domainMax;
            std::vector<double> domainMin;

            for(unsigned dim=0; dim<domainDimensions.size(); dim++)
            {
                domainMax.push_back(domainDimensions[dim]);
                domainMin.push_back(0.0);
            }

            boost::shared_ptr<FixedCellBoundaryCondition<elementDim,spaceDim>> p_cell_BC(new FixedCellBoundaryCondition<elementDim,spaceDim>(&cell_population, domainMax, domainMin));

            simulator.AddCellPopulationBoundaryCondition(p_cell_BC);

            simulator.Solve();

        }
        else
        {
            // set up chemical domain field
            ChemicalDomainFieldForCellCoupling<elementDim,spaceDim,probDim>* p_Pde_field = new ChemicalDomainFieldForCellCoupling<elementDim,spaceDim,probDim>
                                    (   variables_map["domain_file_root"].as<std::string>(),
                                        variables_map["domain_file_root"].as<std::string>()+"",
                                        variables_map["domain_file_root"].as<std::string>()+"",
                                        variables_map["domain_file_root"].as<std::string>()+variables_map["domain_file"].as<std::string>(), 
                                        variables_map["domain_file_root"].as<std::string>()+variables_map["domain_key_file"].as<std::string>(), 
                                        variables_map["domain_file_root"].as<std::string>()+variables_map["ode_file"].as<std::string>(), 
                                        variables_map["domain_file_root"].as<std::string>()+variables_map["ode_key_file"].as<std::string>(), 
                                        variables_map["domain_file_root"].as<std::string>()+variables_map["diffusion_database"].as<std::string>(), 
                                        variables_map["domain_file_root"].as<std::string>()+variables_map["initial_conditions"].as<std::string>(), 
                                        variables_map["domain_file_root"].as<std::string>()+variables_map["boundary_conditions"].as<std::string>());
            


            // make mesh of cell with associated mesh based cell population
            HoneycombMeshGenerator generator(variables_map["number_cells_across"].as<unsigned>(),variables_map["number_cells_high"].as<unsigned>());    // Parameters are: cells across, cells up
            MutableMesh<elementDim,spaceDim>* p_mesh = generator.GetMesh();


            std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

            std::vector<CellPtr> cells;

            std::vector<std::string> chemicalCellSpeciesNames;

            // assume cell at each node in cell layer mesh
            bool IsFirstCell = true;
            for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
            {
                // provide each cell with a transport cell property and membrane property, cell cycle, wild type states

                CustomCellFromFile* p_cell_reader = new CustomCellFromFile(
                                    variables_map["cell_file_root"].as<std::string>(),
                                    variables_map["cell_cycle_file"].as<std::string>(), 
                                    variables_map["cell_division_rules"].as<std::string>(),
                                    variables_map["srn_file"].as<std::string>(),
                                    variables_map["cell_initial_conditions"].as<std::string>(),
                                    variables_map["transport_property"].as<std::string>(),
                                    variables_map["membrane_property"].as<std::string>(),
                                    variables_map["environment_property"].as<std::string>(),
                                    variables_map["state_switching_property"].as<std::string>()
                                    );

                p_cell_reader->SetUp();

                if(IsFirstCell)
                {
                    chemicalCellSpeciesNames =  p_cell_reader->GetFullChemicalNamesVector();
                    IsFirstCell = false;
                }
                cells.push_back(p_cell_reader -> GetCellPtr());
            }

            // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
            ChastePoint<spaceDim> lower(variables_map["cell_mesh_origin"].as<double>(),variables_map["cell_mesh_origin"].as<double>());
            ChastePoint<spaceDim> upper(10.0, 10.0);
            MAKE_PTR_ARGS(ChasteCuboid<spaceDim>, p_cuboid, (lower, upper));

            
            
            double feMeshStepSize = variables_map["FE_mesh_step_size"].as<double>();
            
            bool isCenterMesh = false;
            if(variables_map["mesh_location"].as<std::string>()=="Center")
            {
                isCenterMesh = true;
            }
            
            boost::shared_ptr<ParabolicBoxDomainPdeSystemModifier<elementDim,spaceDim,probDim>> p_pde_modifier(new ParabolicBoxDomainPdeSystemModifier<elementDim,spaceDim,probDim>(p_Pde_field, p_cuboid, feMeshStepSize, isCenterMesh));
            
         
            boost::shared_ptr<ChemicalTrackingPreSwitchModifier<elementDim,spaceDim>> p_chemical_tracking_pre_switch_modifier(new ChemicalTrackingPreSwitchModifier<elementDim,spaceDim>());

            boost::shared_ptr<CellStateTrackingModifier<elementDim,spaceDim>> p_cell_tracking_modifier(new CellStateTrackingModifier<elementDim,spaceDim>());
            
            boost::shared_ptr<ChemicalTrackingPostSwitchModifier<elementDim,spaceDim>> p_chemical_tracking_post_switch_modifier(new ChemicalTrackingPostSwitchModifier<elementDim,spaceDim>());


            NodesOnlyMesh<spaceDim> mesh;

            mesh.ConstructNodesWithoutMesh(*p_mesh, variables_map["node_cutoff_length"].as<double>());

            NodeBasedCellPopulation<spaceDim> cell_population(mesh, cells);   

            // writers

            cell_population.AddCellWriter<CellIdWriter>();
            cell_population.AddCellWriter<CellAgesWriter>();
            cell_population.AddCellWriter<CellLocationIndexWriter>();
            cell_population.AddCellWriter<CellAnalyticsWriter>();
            cell_population.AddCellWriter<CellStateWriter>();

            OffLatticeSimulation<spaceDim> simulator(cell_population);

            simulator.AddSimulationModifier(p_pde_modifier);

            simulator.AddSimulationModifier(p_chemical_tracking_pre_switch_modifier);
            simulator.AddSimulationModifier(p_cell_tracking_modifier);
            simulator.AddSimulationModifier(p_chemical_tracking_post_switch_modifier);

            simulator.SetOutputDirectory(variables_map["output_filename"].as<std::string>()+"/"+boost::lexical_cast<std::string>(simulation_id));
            simulator.SetEndTime(variables_map["simulation_end_time"].as<double>());
            simulator.SetDt(variables_map["simulation_timestep"].as<double>());

            MAKE_PTR(GeneralisedLinearSpringForce<spaceDim>, p_linear_force);
            p_linear_force->SetCutOffLength(variables_map["linear_force_cutoff"].as<double>());
            simulator.AddForce(p_linear_force);

            std::vector<unsigned>  domainDimensions = p_Pde_field->GetMeshDimensions();
            std::vector<double> domainMax;
            std::vector<double> domainMin;

            for(unsigned dim=0; dim<domainDimensions.size(); dim++)
            {
                domainMax.push_back(domainDimensions[dim]);
                domainMin.push_back(0.0);
            }

            boost::shared_ptr<FixedCellBoundaryCondition<elementDim,spaceDim>> p_cell_BC(new FixedCellBoundaryCondition<elementDim,spaceDim>(&cell_population, domainMax, domainMin));

            simulator.AddCellPopulationBoundaryCondition(p_cell_BC);

    
            simulator.Solve();
        }

    }

    */
    else if(variables_map["simulation_type"].as<std::string>()=="domain_only")
    {

        // run the domain field set up and parse files
        ChemicalDomainFieldTemplated<elementDim,spaceDim,probDim>* p_field = new ChemicalDomainFieldTemplated<elementDim,spaceDim,probDim>(
                    variables_map["domain_file_root"].as<std::string>(),
                    variables_map["domain_file_root"].as<std::string>()+variables_map["domain_file"].as<std::string>(), 
                    variables_map["domain_file_root"].as<std::string>()+variables_map["domain_key_file"].as<std::string>(), 
                    variables_map["domain_file_root"].as<std::string>()+variables_map["ode_file"].as<std::string>(), 
                    variables_map["domain_file_root"].as<std::string>()+variables_map["ode_key_file"].as<std::string>(), 
                    variables_map["domain_file_root"].as<std::string>()+variables_map["diffusion_database"].as<std::string>() 
                    );

        TetrahedralMesh<elementDim,spaceDim>* p_mesh = p_field ->rGetDomainFeMesh();

        std::cout<<"Domain width: "<<p_mesh->GetWidth(0)<<" Domain height: "<<p_mesh->GetWidth(1)<<std::endl;
        
        if(std::abs(variables_map["domain_width"].as<double>()-0)>1e-10)
        {
            p_mesh -> Scale(variables_map["domain_width"].as<double>()/p_mesh->GetWidth(0), variables_map["domain_height"].as<double>()/p_mesh ->GetWidth(1));

        }
        
        p_mesh -> Scale(variables_map["scale_domain_width"].as<double>(), variables_map["scale_domain_height"].as<double>());

        std::cout<<"Scaled width: "<<p_mesh->GetWidth(0)<<" Scaled height: "<<p_mesh->GetWidth(1)<<std::endl;

        p_field -> ParseInitialConditionsFromFile(variables_map["domain_file_root"].as<std::string>()+variables_map["initial_conditions"].as<std::string>());

        Vec initial_condition = PetscTools::CreateVec(p_field -> GetInitialNodeConditions());

        // process boundary conditions
        p_field -> ParseBoundaryConditionsFromFile(variables_map["domain_file_root"].as<std::string>()+variables_map["boundary_conditions"].as<std::string>());

        std::vector<std::string> boundaryConditionTypes = p_field -> GetBoundaryConditionTypes();

        std::vector<double> boundaryConditionValues = p_field -> GetBoundaryConditionValues();

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

        std::vector<AbstractInhomogenousOdeSystemForCoupledPdeSystem*> odeSystem;
        std::vector<boost::shared_ptr<AbstractIvpOdeSolver> > solverSystem;
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++){
            // number of ode system objects must match the number of nodes, i.e the individual odes may be multi-dimensional
            odeSystem.push_back(p_field->GetOdeSystem()[i]);
            boost::shared_ptr<AbstractIvpOdeSolver> p_solver(new EulerIvpOdeSolver);
            solverSystem.push_back(p_solver);
        }


        // pde system
        InhomogenousParabolicPdeForCoupledOdeSystemTemplated<elementDim, spaceDim, probDim> pde(p_field);

        // solver
        InhomogenousCoupledPdeOdeSolverTemplated<elementDim,spaceDim,probDim> solver(p_mesh, &pde, &bcc,odeSystem,solverSystem);

        // solver properties
    
        solver.SetTimes(0.0, variables_map["simulation_end_time"].as<double>());
        solver.SetTimeStep(variables_map["simulation_timestep"].as<double>());
        if(variables_map["simulation_timestep"].as<double>()>variables_map["sampling_timestep"].as<double>())
        {
            solver.SetSamplingTimeStep(variables_map["simulation_timestep"].as<double>());
        }
        else
        {
            solver.SetSamplingTimeStep(variables_map["sampling_timestep"].as<double>());
        }
        solver.SetOutputDirectory(variables_map["output_filename"].as<std::string>()+"_dt_"+boost::lexical_cast<std::string>(variables_map["simulation_timestep"].as<double>()));
        solver.SetInitialCondition(initial_condition);
 
        solver.SolveAndWriteResultsToFile();

    }
    else if(variables_map["simulation_type"].as<std::string>()=="fisher_equation")
    {

        double MeshStepSize = 1.0;
        std::vector<unsigned> MeshDimensions = {100,10,0};
        std::vector<double> initValuesHigh = {1.0};
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
        

        std::cout<<"Domain width: "<<p_mesh->GetWidth(0)<<" Domain height: "<<p_mesh->GetWidth(1)<<std::endl;
        
        if(std::abs(variables_map["domain_width"].as<double>()-0)>1e-10)
        {
            p_mesh -> Scale(variables_map["domain_width"].as<double>()/p_mesh->GetWidth(0), variables_map["domain_height"].as<double>()/p_mesh ->GetWidth(1));

        }
        
        p_mesh -> Scale(variables_map["scale_domain_width"].as<double>(), variables_map["scale_domain_height"].as<double>());

        std::cout<<"Scaled width: "<<p_mesh->GetWidth(0)<<" Scaled height: "<<p_mesh->GetWidth(1)<<std::endl;

        // Process Boundary Conditions
        BoundaryConditionsContainer<elementDim,spaceDim,1> bcc;
        std::vector<bool> areNeumannBoundaryConditions(1, true);
        std::vector<ConstBoundaryCondition<spaceDim>*> vectorConstBCs;
        
        for(unsigned pdeDim=0; pdeDim<1; pdeDim++){
            vectorConstBCs.push_back(new ConstBoundaryCondition<spaceDim>(bcValues[pdeDim]));
        }
        
        for(unsigned pdeDim=0; pdeDim<1; pdeDim++)
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
        std::vector<double> init_conds(1*p_mesh->GetNumNodes(),0.0);
        unsigned columnNum = 0;
        unsigned rowNum = 0;
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {   // set as being a random perturbation about the boundary values
            if(spaceDim==2)
            {
                columnNum = 0;
                rowNum = 0;

                while(i >= rowNum*(MeshDimensions[0]+1))
                {
                    rowNum = rowNum + 1;
                
                }
                
                columnNum = i - (rowNum-1)*(MeshDimensions[0]+1);
                if(columnNum<1)
                {
                    for(unsigned pdeDim=0; pdeDim<1; pdeDim++)
                    {   // serialised for nodes
                        init_conds[1*i + pdeDim] = fabs(initValuesHigh[pdeDim]);// + RandomNumberGenerator::Instance()->ranf());
                    }
                }
                else
                {
                    for(unsigned pdeDim=0; pdeDim<1; pdeDim++)
                    {   // serialised for nodes
                        init_conds[1*i + pdeDim] = fabs(initValuesLow[pdeDim]);// + RandomNumberGenerator::Instance()->ranf());
                    }
                }
            }
            else
            {
                for(unsigned pdeDim=0; pdeDim<1; pdeDim++)
                {   // serialised for nodes
                    init_conds[1*i + pdeDim] = fabs(initValuesLow[pdeDim] + RandomNumberGenerator::Instance()->ranf());
                }
            }

        }
        // PETSc Vec
        Vec initial_condition = PetscTools::CreateVec(init_conds);

        std::vector<AbstractInhomogenousOdeSystemForCoupledPdeSystem*> odeSystem;
        std::vector<boost::shared_ptr<AbstractIvpOdeSolver> > solverSystem;

        InhomogenousFisherPde<elementDim, spaceDim, 1> pde(diffusionRates,growthRates,carryingCapacities);

        // solver
        InhomogenousCoupledPdeOdeSolverTemplated<elementDim,spaceDim,1> solver(p_mesh, &pde, &bcc);

        // solver properties
        double startTime= 0.0;
        double endTime = variables_map["simulation_end_time"].as<double>();

        solver.SetTimes(startTime, endTime);

        solver.SetTimeStep(variables_map["simulation_timestep"].as<double>());
        if(variables_map["simulation_timestep"].as<double>()>variables_map["sampling_timestep"].as<double>())
        {
            solver.SetSamplingTimeStep(variables_map["simulation_timestep"].as<double>());
        }
        else
        {
            solver.SetSamplingTimeStep(variables_map["sampling_timestep"].as<double>());
        }
        
        solver.SetOutputDirectory(variables_map["output_filename"].as<std::string>()+"_dt_"+boost::lexical_cast<std::string>(variables_map["simulation_timestep"].as<double>()));
        solver.SetInitialCondition(initial_condition);
        // solve
        solver.SolveAndWriteResultsToFile();
    }
    else if(variables_map["simulation_type"].as<std::string>()=="fisher_equation_long")
    {

        double MeshStepSize = variables_map["scale_domain_width"].as<double>();
        std::vector<unsigned> MeshDimensions = {100,10,0};
        std::vector<double> initValuesHigh = {1.0};
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
        
       // std::cout<<"Domain width: "<<p_mesh->GetWidth(0)<<" Domain height: "<<p_mesh->GetWidth(1)<<std::endl;
        
    //    if(std::abs(variables_map["domain_width"].as<double>()-0)>1e-10)
    //    {
    //        p_mesh -> Scale(variables_map["domain_width"].as<double>()/p_mesh->GetWidth(0), variables_map["domain_height"].as<double>()/p_mesh ->GetWidth(1));

    //    }
        
    //    p_mesh -> Scale(variables_map["scale_domain_width"].as<double>(), variables_map["scale_domain_height"].as<double>());

    //    std::cout<<"Scaled width: "<<p_mesh->GetWidth(0)<<" Scaled height: "<<p_mesh->GetWidth(1)<<std::endl;

        // Process Boundary Conditions
        BoundaryConditionsContainer<elementDim,spaceDim,1> bcc;
        std::vector<bool> areNeumannBoundaryConditions(1, true);
        std::vector<ConstBoundaryCondition<spaceDim>*> vectorConstBCs;
        
        for(unsigned pdeDim=0; pdeDim<1; pdeDim++){
            vectorConstBCs.push_back(new ConstBoundaryCondition<spaceDim>(bcValues[pdeDim]));
        }
        
        for(unsigned pdeDim=0; pdeDim<1; pdeDim++)
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
        std::vector<double> init_conds(1*p_mesh->GetNumNodes(),0.0);
        unsigned columnNum = 0;
        unsigned rowNum = 0;
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {   // set as being a random perturbation about the boundary values
            
            double x = p_mesh->GetNode(i)->rGetLocation()[0];
            
            if(x<=1.0)
            {
                for(unsigned pdeDim=0; pdeDim<1; pdeDim++)
                {   // serialised for nodes
                    init_conds[1*i + pdeDim] = fabs(initValuesHigh[pdeDim]);// + RandomNumberGenerator::Instance()->ranf());
                }

            }
            else
            {
                for(unsigned pdeDim=0; pdeDim<1; pdeDim++)
                {   // serialised for nodes
                    init_conds[1*i + pdeDim] = fabs(initValuesLow[pdeDim]);// + RandomNumberGenerator::Instance()->ranf());
                }
            }
        }
        // PETSc Vec
        Vec initial_condition = PetscTools::CreateVec(init_conds);

        std::vector<AbstractInhomogenousOdeSystemForCoupledPdeSystem*> odeSystem;
        std::vector<boost::shared_ptr<AbstractIvpOdeSolver> > solverSystem;

        InhomogenousFisherPde<elementDim, spaceDim, 1> pde(diffusionRates,growthRates,carryingCapacities);

        // solver
        InhomogenousCoupledPdeOdeSolverTemplated<elementDim,spaceDim,1> solver(p_mesh, &pde, &bcc);

        // solver properties
        double startTime= 0.0;
        double endTime = variables_map["simulation_end_time"].as<double>();

        solver.SetTimes(startTime, endTime);

        solver.SetTimeStep(variables_map["simulation_timestep"].as<double>());
        if(variables_map["simulation_timestep"].as<double>()>variables_map["sampling_timestep"].as<double>())
        {
            solver.SetSamplingTimeStep(variables_map["simulation_timestep"].as<double>());
        }
        else
        {
            solver.SetSamplingTimeStep(variables_map["sampling_timestep"].as<double>());
        }
        
        solver.SetOutputDirectory(variables_map["output_filename"].as<std::string>()+"_dt_"+boost::lexical_cast<std::string>(variables_map["simulation_timestep"].as<double>()));
        solver.SetInitialCondition(initial_condition);
        // solve
        solver.SolveAndWriteResultsToFile();
    }
    else if(variables_map["simulation_type"].as<std::string>()=="schnackenberg_equation_D")
    {


        //double a =2;
        //double b = 3;
        //double gamma = 1000;

        double a =0.2;
        double b = 1.3;
        double gamma = 920;

        std::vector<double> diffusionRates = {40,1};

        //std::vector<double> diffusionRates = {1.0,1000};
        std::vector<double> initValues = {a+b,b/((a+b)*(a+b))};
        std::vector<double> bcValues = initValues;
        //std::vector<double> diffusionRates = {1e-4,1e-2};
        //double k1 = 0.1;
        //double k_1 = 0.2;
        //double k2 = 0.3;
        //double k3 = 0.1;

        double MeshStepSize = 1.0;
        std::vector<unsigned> MeshDimensions = {100,100,0};
        //std::vector<double> initValues = {2.0,0.75};
        //std::vector<double> bcValues = {2.0,0.75};
        std::vector<bool> areBCsNeumann = {false,false};

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
        BoundaryConditionsContainer<elementDim,spaceDim,probDim> bcc;
        //std::vector<bool> areNeumannBoundaryConditions(probDim, true);
        std::vector<ConstBoundaryCondition<spaceDim>*> vectorConstBCs;
        
        for(unsigned pdeDim=0; pdeDim<probDim; pdeDim++){
            vectorConstBCs.push_back(new ConstBoundaryCondition<spaceDim>(bcValues[pdeDim]));
        }
        
        for(unsigned pdeDim=0; pdeDim<probDim; pdeDim++)
        {
            if(areBCsNeumann[pdeDim]==false)
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
            
            for(unsigned pdeDim=0; pdeDim<probDim; pdeDim++)
            {   // serialised for nodes
                init_conds[probDim*i + pdeDim] = fabs(initValues[pdeDim] + RandomNumberGenerator::Instance()->ranf());
            }
        
        }
        // PETSc Vec
        Vec initial_condition = PetscTools::CreateVec(init_conds);

        std::vector<AbstractInhomogenousOdeSystemForCoupledPdeSystem*> odeSystem;
        std::vector<boost::shared_ptr<AbstractIvpOdeSolver> > solverSystem;

        InhomogenousSchnackenbergPde<elementDim, spaceDim, probDim> pde(diffusionRates,a,b,gamma);

        // solver
        InhomogenousCoupledPdeOdeSolverTemplated<elementDim,spaceDim,probDim> solver(p_mesh, &pde, &bcc);

        // solver properties
        double startTime= 0.0;
        double endTime = variables_map["simulation_end_time"].as<double>();

        solver.SetTimes(startTime, endTime);

        solver.SetTimeStep(variables_map["simulation_timestep"].as<double>());
        if(variables_map["simulation_timestep"].as<double>()>variables_map["sampling_timestep"].as<double>())
        {
            solver.SetSamplingTimeStep(variables_map["simulation_timestep"].as<double>());
        }
        else
        {
            solver.SetSamplingTimeStep(variables_map["sampling_timestep"].as<double>());
        }
        
        solver.SetOutputDirectory(variables_map["output_filename"].as<std::string>()+"_dt_"+boost::lexical_cast<std::string>(variables_map["simulation_timestep"].as<double>()));
        solver.SetInitialCondition(initial_condition);
        // solve
        solver.SolveAndWriteResultsToFile();
    }
    else if(variables_map["simulation_type"].as<std::string>()=="schnackenberg_equation_E")
    {

        double a =0.2;
        double b = 1.3;
        double gamma = 920;

        std::vector<double> diffusionRates = {1e-4,8.4e-7};

        

        double MeshStepSize = 1.0;
        std::vector<unsigned> MeshDimensions = {100,100,0};
        std::vector<double> initValues = {a+b,b/((a+b)*(a+b))};
        std::vector<double> bcValues = {0.0};
        std::vector<double> areBCsNeumann = {true,true};

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
        BoundaryConditionsContainer<elementDim,spaceDim,probDim> bcc;
        
        std::vector<ConstBoundaryCondition<spaceDim>*> vectorConstBCs;
        
        for(unsigned pdeDim=0; pdeDim<probDim; pdeDim++){
            vectorConstBCs.push_back(new ConstBoundaryCondition<spaceDim>(bcValues[pdeDim]));
        }
        
        for(unsigned pdeDim=0; pdeDim<probDim; pdeDim++)
        {
            if(areBCsNeumann[pdeDim]==false)
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
            
            for(unsigned pdeDim=0; pdeDim<probDim; pdeDim++)
            {   // serialised for nodes
                init_conds[probDim*i + pdeDim] = fabs(initValues[pdeDim] + RandomNumberGenerator::Instance()->ranf());
            }
        
        }
        // PETSc Vec
        Vec initial_condition = PetscTools::CreateVec(init_conds);

        std::vector<AbstractInhomogenousOdeSystemForCoupledPdeSystem*> odeSystem;
        std::vector<boost::shared_ptr<AbstractIvpOdeSolver> > solverSystem;

        InhomogenousSchnackenbergPde<elementDim, spaceDim, probDim> pde(diffusionRates,a,b,gamma);

        // solver
        InhomogenousCoupledPdeOdeSolverTemplated<elementDim,spaceDim,probDim> solver(p_mesh, &pde, &bcc);

        // solver properties
        double startTime= 0.0;
        double endTime = variables_map["simulation_end_time"].as<double>();

        solver.SetTimes(startTime, endTime);

        solver.SetTimeStep(variables_map["simulation_timestep"].as<double>());
        if(variables_map["simulation_timestep"].as<double>()>variables_map["sampling_timestep"].as<double>())
        {
            solver.SetSamplingTimeStep(variables_map["simulation_timestep"].as<double>());
        }
        else
        {
            solver.SetSamplingTimeStep(variables_map["sampling_timestep"].as<double>());
        }
        
        solver.SetOutputDirectory(variables_map["output_filename"].as<std::string>()+"_dt_"+boost::lexical_cast<std::string>(variables_map["simulation_timestep"].as<double>()));
        solver.SetInitialCondition(initial_condition);
        // solve
        solver.SolveAndWriteResultsToFile();
    }
    else if(variables_map["simulation_type"].as<std::string>()=="schnackenberg_equation_F")
    {


        double a =0.2;
        double b = 1.3;
        double gamma = 220;

        std::vector<double> diffusionRates = {1e-4,119e-4};

        double MeshStepSize = 1.0;
        std::vector<unsigned> MeshDimensions = {100,100,0};
        std::vector<double> initValues = {a+b,b/((a+b)*(a+b))};
        std::vector<double> bcValues = {0.0};
        std::vector<double> areBCsNeumann = {true,true};

     

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
        BoundaryConditionsContainer<elementDim,spaceDim,probDim> bcc;
  
        std::vector<ConstBoundaryCondition<spaceDim>*> vectorConstBCs;
        
        for(unsigned pdeDim=0; pdeDim<probDim; pdeDim++){
            vectorConstBCs.push_back(new ConstBoundaryCondition<spaceDim>(bcValues[pdeDim]));
        }
        
        for(unsigned pdeDim=0; pdeDim<probDim; pdeDim++)
        {
            if(areBCsNeumann[pdeDim]==false)
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
            
            for(unsigned pdeDim=0; pdeDim<probDim; pdeDim++)
            {   // serialised for nodes
                init_conds[probDim*i + pdeDim] = fabs(initValues[pdeDim] + RandomNumberGenerator::Instance()->ranf());
            }
        
        }
        // PETSc Vec
        Vec initial_condition = PetscTools::CreateVec(init_conds);

        std::vector<AbstractInhomogenousOdeSystemForCoupledPdeSystem*> odeSystem;
        std::vector<boost::shared_ptr<AbstractIvpOdeSolver> > solverSystem;

        InhomogenousSchnackenbergPde<elementDim, spaceDim, probDim> pde(diffusionRates,a,b,gamma);

        // solver
        InhomogenousCoupledPdeOdeSolverTemplated<elementDim,spaceDim,probDim> solver(p_mesh, &pde, &bcc);

        // solver properties
        double startTime= 0.0;
        double endTime = variables_map["simulation_end_time"].as<double>();

        solver.SetTimes(startTime, endTime);

        solver.SetTimeStep(variables_map["simulation_timestep"].as<double>());
        if(variables_map["simulation_timestep"].as<double>()>variables_map["sampling_timestep"].as<double>())
        {
            solver.SetSamplingTimeStep(variables_map["simulation_timestep"].as<double>());
        }
        else
        {
            solver.SetSamplingTimeStep(variables_map["sampling_timestep"].as<double>());
        }
        
        solver.SetOutputDirectory(variables_map["output_filename"].as<std::string>()+"_dt_"+boost::lexical_cast<std::string>(variables_map["simulation_timestep"].as<double>()));
        solver.SetInitialCondition(initial_condition);
        // solve
        solver.SolveAndWriteResultsToFile();
    }
    else if(variables_map["simulation_type"].as<std::string>()=="heat_equation")
    {

        double MeshStepSize = 1.0;
        std::vector<unsigned> MeshDimensions = {100,10,0};
        std::vector<double> initValuesHigh = {1.0};
        std::vector<double> initValuesLow = {0.0};
        std::vector<double> bcValues = {0.0};
        std::vector<double> areBCsNeumann = {true};
        std::vector<double> diffusionRates = {100.0};


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
        BoundaryConditionsContainer<elementDim,spaceDim,1> bcc;
        std::vector<bool> areNeumannBoundaryConditions(1, true);
        std::vector<ConstBoundaryCondition<spaceDim>*> vectorConstBCs;
        
        for(unsigned pdeDim=0; pdeDim<1; pdeDim++){
            vectorConstBCs.push_back(new ConstBoundaryCondition<spaceDim>(bcValues[pdeDim]));
        }
        
        for(unsigned pdeDim=0; pdeDim<1; pdeDim++)
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
        std::vector<double> init_conds(1*p_mesh->GetNumNodes(),0.0);
        unsigned columnNum = 0;
        unsigned rowNum = 0;
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {   // set as being a random perturbation about the boundary values
            if(spaceDim==2)
            {
                columnNum = 0;
                rowNum = 0;

                while(i >= rowNum*(MeshDimensions[0]+1))
                {
                    rowNum = rowNum + 1;
                
                }
                
                columnNum = i - (rowNum-1)*(MeshDimensions[0]+1);
                if(columnNum<1)
                {
                    for(unsigned pdeDim=0; pdeDim<1; pdeDim++)
                    {   // serialised for nodes
                        init_conds[1*i + pdeDim] = fabs(initValuesHigh[pdeDim]);// + RandomNumberGenerator::Instance()->ranf());
                    }
                }
                else
                {
                    for(unsigned pdeDim=0; pdeDim<1; pdeDim++)
                    {   // serialised for nodes
                        init_conds[1*i + pdeDim] = fabs(initValuesLow[pdeDim]);// + RandomNumberGenerator::Instance()->ranf());
                    }
                }
            }
            else
            {
                for(unsigned pdeDim=0; pdeDim<1; pdeDim++)
                {   // serialised for nodes
                    init_conds[1*i + pdeDim] = fabs(initValuesLow[pdeDim] + RandomNumberGenerator::Instance()->ranf());
                }
            }

        }
        // PETSc Vec
        Vec initial_condition = PetscTools::CreateVec(init_conds);

        std::vector<AbstractInhomogenousOdeSystemForCoupledPdeSystem*> odeSystem;
        std::vector<boost::shared_ptr<AbstractIvpOdeSolver> > solverSystem;

        InhomogenousHeatEquationPde<elementDim, spaceDim, 1> pde(diffusionRates);
        // solver
        InhomogenousCoupledPdeOdeSolverTemplated<elementDim,spaceDim,1> solver(p_mesh, &pde, &bcc);

        // solver properties
        double startTime= 0.0;
        double endTime = variables_map["simulation_end_time"].as<double>();

        solver.SetTimes(startTime, endTime);

        solver.SetTimeStep(variables_map["simulation_timestep"].as<double>());
        if(variables_map["simulation_timestep"].as<double>()>variables_map["sampling_timestep"].as<double>())
        {
            solver.SetSamplingTimeStep(variables_map["simulation_timestep"].as<double>());
        }
        else
        {
            solver.SetSamplingTimeStep(variables_map["sampling_timestep"].as<double>());
        }
        
        solver.SetOutputDirectory(variables_map["output_filename"].as<std::string>()+"_dt_"+boost::lexical_cast<std::string>(variables_map["simulation_timestep"].as<double>()));
        solver.SetInitialCondition(initial_condition);
        // solve
        solver.SolveAndWriteResultsToFile();
    }
    else if(variables_map["simulation_type"].as<std::string>()=="heat_equation_long")
    {

        double MeshStepSize = 1.0;
        std::vector<unsigned> MeshDimensions = {100,10,0};
        std::vector<double> initValuesHigh = {1.0};
        std::vector<double> initValuesLow = {0.0};
        std::vector<double> bcValues = {0.0};
        std::vector<double> areBCsNeumann = {true};
        std::vector<double> diffusionRates = {100.0};



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
        BoundaryConditionsContainer<elementDim,spaceDim,1> bcc;
        std::vector<bool> areNeumannBoundaryConditions(1, true);
        std::vector<ConstBoundaryCondition<spaceDim>*> vectorConstBCs;
        
        for(unsigned pdeDim=0; pdeDim<1; pdeDim++){
            vectorConstBCs.push_back(new ConstBoundaryCondition<spaceDim>(bcValues[pdeDim]));
        }
        
        for(unsigned pdeDim=0; pdeDim<1; pdeDim++)
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
        std::vector<double> init_conds(1*p_mesh->GetNumNodes(),0.0);
        unsigned columnNum = 0;
        unsigned rowNum = 0;
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {   // set as being a random perturbation about the boundary values
            if(spaceDim==2)
            {
                columnNum = 0;
                rowNum = 0;

                while(i >= rowNum*(MeshDimensions[0]+1))
                {
                    rowNum = rowNum + 1;
                
                }
                
                columnNum = i - (rowNum-1)*(MeshDimensions[0]+1);
                if(columnNum<1)
                {
                    for(unsigned pdeDim=0; pdeDim<1; pdeDim++)
                    {   // serialised for nodes
                        init_conds[1*i + pdeDim] = fabs(initValuesHigh[pdeDim]);// + RandomNumberGenerator::Instance()->ranf());
                    }
                }
                else
                {
                    for(unsigned pdeDim=0; pdeDim<1; pdeDim++)
                    {   // serialised for nodes
                        init_conds[1*i + pdeDim] = fabs(initValuesLow[pdeDim]);// + RandomNumberGenerator::Instance()->ranf());
                    }
                }
            }
            else
            {
                for(unsigned pdeDim=0; pdeDim<1; pdeDim++)
                {   // serialised for nodes
                    init_conds[1*i + pdeDim] = fabs(initValuesLow[pdeDim] + RandomNumberGenerator::Instance()->ranf());
                }
            }

        }
        // PETSc Vec
        Vec initial_condition = PetscTools::CreateVec(init_conds);

        std::vector<AbstractInhomogenousOdeSystemForCoupledPdeSystem*> odeSystem;
        std::vector<boost::shared_ptr<AbstractIvpOdeSolver> > solverSystem;

        InhomogenousHeatEquationPde<elementDim, spaceDim, 1> pde(diffusionRates);

        // solver
        InhomogenousCoupledPdeOdeSolverTemplated<elementDim,spaceDim,1> solver(p_mesh, &pde, &bcc);

        // solver properties
        double startTime= 0.0;
        double endTime = variables_map["simulation_end_time"].as<double>();

        solver.SetTimes(startTime, endTime);

        solver.SetTimeStep(variables_map["simulation_timestep"].as<double>());
        if(variables_map["simulation_timestep"].as<double>()>variables_map["sampling_timestep"].as<double>())
        {
            solver.SetSamplingTimeStep(variables_map["simulation_timestep"].as<double>());
        }
        else
        {
            solver.SetSamplingTimeStep(variables_map["sampling_timestep"].as<double>());
        }
        
        solver.SetOutputDirectory(variables_map["output_filename"].as<std::string>()+"_dt_"+boost::lexical_cast<std::string>(variables_map["simulation_timestep"].as<double>()));
        solver.SetInitialCondition(initial_condition);
        // solve
        solver.SolveAndWriteResultsToFile();
    }
    else if(variables_map["simulation_type"].as<std::string>()=="heat_equation_source")
    {

        double MeshStepSize = 1.0;
        std::vector<unsigned> MeshDimensions = {100,10,0};
        std::vector<double> initValuesHigh = {1.0};
        std::vector<double> initValuesLow = {0.0};
        std::vector<double> bcValues = {0.0};
        std::vector<double> areBCsNeumann = {true};
        std::vector<double> diffusionRates = {100.0};
        double source_term = 0.1;


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
        BoundaryConditionsContainer<elementDim,spaceDim,1> bcc;
        std::vector<bool> areNeumannBoundaryConditions(1, true);
        std::vector<ConstBoundaryCondition<spaceDim>*> vectorConstBCs;
        
        for(unsigned pdeDim=0; pdeDim<1; pdeDim++){
            vectorConstBCs.push_back(new ConstBoundaryCondition<spaceDim>(bcValues[pdeDim]));
        }
        
        for(unsigned pdeDim=0; pdeDim<1; pdeDim++)
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
        std::vector<double> init_conds(1*p_mesh->GetNumNodes(),0.0);
        unsigned columnNum = 0;
        unsigned rowNum = 0;
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {   // set as being a random perturbation about the boundary values
            if(spaceDim==2)
            {
                columnNum = 0;
                rowNum = 0;

                while(i >= rowNum*(MeshDimensions[0]+1))
                {
                    rowNum = rowNum + 1;
                
                }
                
                columnNum = i - (rowNum-1)*(MeshDimensions[0]+1);
                if(columnNum<1)
                {
                    for(unsigned pdeDim=0; pdeDim<1; pdeDim++)
                    {   // serialised for nodes
                        init_conds[1*i + pdeDim] = fabs(initValuesHigh[pdeDim]);// + RandomNumberGenerator::Instance()->ranf());
                    }
                }
                else
                {
                    for(unsigned pdeDim=0; pdeDim<1; pdeDim++)
                    {   // serialised for nodes
                        init_conds[1*i + pdeDim] = fabs(initValuesLow[pdeDim]);// + RandomNumberGenerator::Instance()->ranf());
                    }
                }
            }
            else
            {
                for(unsigned pdeDim=0; pdeDim<1; pdeDim++)
                {   // serialised for nodes
                    init_conds[1*i + pdeDim] = fabs(initValuesLow[pdeDim] + RandomNumberGenerator::Instance()->ranf());
                }
            }

        }
        // PETSc Vec
        Vec initial_condition = PetscTools::CreateVec(init_conds);

        std::vector<AbstractInhomogenousOdeSystemForCoupledPdeSystem*> odeSystem;
        std::vector<boost::shared_ptr<AbstractIvpOdeSolver> > solverSystem;

        InhomogenousHeatEquationWithSourcePde<elementDim, spaceDim, 1> pde(diffusionRates,source_term);

        // solver
        InhomogenousCoupledPdeOdeSolverTemplated<elementDim,spaceDim,1> solver(p_mesh, &pde, &bcc);

        // solver properties
        double startTime= 0.0;
        double endTime = variables_map["simulation_end_time"].as<double>();

        solver.SetTimes(startTime, endTime);

        solver.SetTimeStep(variables_map["simulation_timestep"].as<double>());
        if(variables_map["simulation_timestep"].as<double>()>variables_map["sampling_timestep"].as<double>())
        {
            solver.SetSamplingTimeStep(variables_map["simulation_timestep"].as<double>());
        }
        else
        {
            solver.SetSamplingTimeStep(variables_map["sampling_timestep"].as<double>());
        }
        
        solver.SetOutputDirectory(variables_map["output_filename"].as<std::string>()+"_dt_"+boost::lexical_cast<std::string>(variables_map["simulation_timestep"].as<double>()));
        solver.SetInitialCondition(initial_condition);
        // solve
        solver.SolveAndWriteResultsToFile();
    }
    else if(variables_map["simulation_type"].as<std::string>()=="heat_equation_source_long")
    {

        double MeshStepSize = 1.0;
        std::vector<unsigned> MeshDimensions = {100,10,0};
        std::vector<double> initValuesHigh = {1.0};
        std::vector<double> initValuesLow = {0.0};
        std::vector<double> bcValues = {0.0};
        std::vector<double> areBCsNeumann = {true};
        std::vector<double> diffusionRates = {100.0};
        double source_term = 0.1;


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
        BoundaryConditionsContainer<elementDim,spaceDim,1> bcc;
        std::vector<bool> areNeumannBoundaryConditions(1, true);
        std::vector<ConstBoundaryCondition<spaceDim>*> vectorConstBCs;
        
        for(unsigned pdeDim=0; pdeDim<1; pdeDim++){
            vectorConstBCs.push_back(new ConstBoundaryCondition<spaceDim>(bcValues[pdeDim]));
        }
        
        for(unsigned pdeDim=0; pdeDim<1; pdeDim++)
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
        std::vector<double> init_conds(1*p_mesh->GetNumNodes(),0.0);
        unsigned columnNum = 0;
        unsigned rowNum = 0;
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {   // set as being a random perturbation about the boundary values
            if(spaceDim==2)
            {
                columnNum = 0;
                rowNum = 0;

                while(i >= rowNum*(MeshDimensions[0]+1))
                {
                    rowNum = rowNum + 1;
                
                }
                
                columnNum = i - (rowNum-1)*(MeshDimensions[0]+1);
                if(columnNum<1)
                {
                    for(unsigned pdeDim=0; pdeDim<1; pdeDim++)
                    {   // serialised for nodes
                        init_conds[1*i + pdeDim] = fabs(initValuesHigh[pdeDim]);// + RandomNumberGenerator::Instance()->ranf());
                    }
                }
                else
                {
                    for(unsigned pdeDim=0; pdeDim<1; pdeDim++)
                    {   // serialised for nodes
                        init_conds[1*i + pdeDim] = fabs(initValuesLow[pdeDim]);// + RandomNumberGenerator::Instance()->ranf());
                    }
                }
            }
            else
            {
                for(unsigned pdeDim=0; pdeDim<1; pdeDim++)
                {   // serialised for nodes
                    init_conds[1*i + pdeDim] = fabs(initValuesLow[pdeDim] + RandomNumberGenerator::Instance()->ranf());
                }
            }

        }
        // PETSc Vec
        Vec initial_condition = PetscTools::CreateVec(init_conds);

        std::vector<AbstractInhomogenousOdeSystemForCoupledPdeSystem*> odeSystem;
        std::vector<boost::shared_ptr<AbstractIvpOdeSolver> > solverSystem;

        InhomogenousHeatEquationWithSourcePde<elementDim, spaceDim, 1> pde(diffusionRates,source_term);

        // solver
        InhomogenousCoupledPdeOdeSolverTemplated<elementDim,spaceDim,1> solver(p_mesh, &pde, &bcc);

        // solver properties
        double startTime= 0.0;
        double endTime = variables_map["simulation_end_time"].as<double>();

        solver.SetTimes(startTime, endTime);

        solver.SetTimeStep(variables_map["simulation_timestep"].as<double>());
        if(variables_map["simulation_timestep"].as<double>()>variables_map["sampling_timestep"].as<double>())
        {
            solver.SetSamplingTimeStep(variables_map["simulation_timestep"].as<double>());
        }
        else
        {
            solver.SetSamplingTimeStep(variables_map["sampling_timestep"].as<double>());
        }
        
        solver.SetOutputDirectory(variables_map["output_filename"].as<std::string>()+"_dt_"+boost::lexical_cast<std::string>(variables_map["simulation_timestep"].as<double>()));
        solver.SetInitialCondition(initial_condition);
        // solve
        solver.SolveAndWriteResultsToFile();
    }
    else if(variables_map["simulation_type"].as<std::string>()=="chaste_fisher_equation")
    {

        double MeshStepSize = 1.0;
        std::vector<unsigned> MeshDimensions = {100,10,0};
        std::vector<double> initValuesHigh = {1.0};
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
        BoundaryConditionsContainer<elementDim,spaceDim,1> bcc;
        std::vector<bool> areNeumannBoundaryConditions(1, true);
        std::vector<ConstBoundaryCondition<spaceDim>*> vectorConstBCs;
        
        for(unsigned pdeDim=0; pdeDim<1; pdeDim++){
            vectorConstBCs.push_back(new ConstBoundaryCondition<spaceDim>(bcValues[pdeDim]));
        }
        
        for(unsigned pdeDim=0; pdeDim<1; pdeDim++)
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
        std::vector<double> init_conds(1*p_mesh->GetNumNodes(),0.0);
        unsigned columnNum = 0;
        unsigned rowNum = 0;
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {   // set as being a random perturbation about the boundary values
            if(spaceDim==2)
            {
                columnNum = 0;
                rowNum = 0;

                while(i >= rowNum*(MeshDimensions[0]+1))
                {
                    rowNum = rowNum + 1;
                
                }
                
                columnNum = i - (rowNum-1)*(MeshDimensions[0]+1);
                if(columnNum<1)
                {
                    for(unsigned pdeDim=0; pdeDim<1; pdeDim++)
                    {   // serialised for nodes
                        init_conds[1*i + pdeDim] = fabs(initValuesHigh[pdeDim]);// + RandomNumberGenerator::Instance()->ranf());
                    }
                }
                else
                {
                    for(unsigned pdeDim=0; pdeDim<1; pdeDim++)
                    {   // serialised for nodes
                        init_conds[1*i + pdeDim] = fabs(initValuesLow[pdeDim]);// + RandomNumberGenerator::Instance()->ranf());
                    }
                }
            }
            else
            {
                for(unsigned pdeDim=0; pdeDim<1; pdeDim++)
                {   // serialised for nodes
                    init_conds[1*i + pdeDim] = fabs(initValuesLow[pdeDim] + RandomNumberGenerator::Instance()->ranf());
                }
            }

        }
        // PETSc Vec
        Vec initial_condition = PetscTools::CreateVec(init_conds);

        std::vector<AbstractInhomogenousOdeSystemForCoupledPdeSystem*> odeSystem;
        std::vector<boost::shared_ptr<AbstractIvpOdeSolver> > solverSystem;

        LinearParabolicFisherPde<elementDim, spaceDim, 1> pde(diffusionRates,growthRates,carryingCapacities);

        // solver
        LinearParabolicPdeSystemWithCoupledOdeSystemSolver<elementDim,spaceDim,1> solver(p_mesh, &pde, &bcc);

        // solver properties
        double startTime= 0.0;
        double endTime = variables_map["simulation_end_time"].as<double>();

        solver.SetTimes(startTime, endTime);

        solver.SetTimeStep(variables_map["simulation_timestep"].as<double>());
        if(variables_map["simulation_timestep"].as<double>()>variables_map["sampling_timestep"].as<double>())
        {
            solver.SetSamplingTimeStep(variables_map["simulation_timestep"].as<double>());
        }
        else
        {
            solver.SetSamplingTimeStep(variables_map["sampling_timestep"].as<double>());
        }
        
        solver.SetOutputDirectory(variables_map["output_filename"].as<std::string>()+"_chaste"+"_dt_"+boost::lexical_cast<std::string>(variables_map["simulation_timestep"].as<double>()));
        solver.SetInitialCondition(initial_condition);
        // solve
        solver.SolveAndWriteResultsToFile();
    }
    else if(variables_map["simulation_type"].as<std::string>()=="chaste_fisher_equation_long")
    {

        double MeshStepSize = 1.0;
        std::vector<unsigned> MeshDimensions = {1000,10,0};
        std::vector<double> initValuesHigh = {1.0};
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
        BoundaryConditionsContainer<elementDim,spaceDim,1> bcc;
        std::vector<bool> areNeumannBoundaryConditions(1, true);
        std::vector<ConstBoundaryCondition<spaceDim>*> vectorConstBCs;
        
        for(unsigned pdeDim=0; pdeDim<1; pdeDim++){
            vectorConstBCs.push_back(new ConstBoundaryCondition<spaceDim>(bcValues[pdeDim]));
        }
        
        for(unsigned pdeDim=0; pdeDim<1; pdeDim++)
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
        std::vector<double> init_conds(1*p_mesh->GetNumNodes(),0.0);
        unsigned columnNum = 0;
        unsigned rowNum = 0;
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {   // set as being a random perturbation about the boundary values
            if(spaceDim==2)
            {
                columnNum = 0;
                rowNum = 0;

                while(i >= rowNum*(MeshDimensions[0]+1))
                {
                    rowNum = rowNum + 1;
                
                }
                
                columnNum = i - (rowNum-1)*(MeshDimensions[0]+1);
                if(columnNum<1)
                {
                    for(unsigned pdeDim=0; pdeDim<1; pdeDim++)
                    {   // serialised for nodes
                        init_conds[1*i + pdeDim] = fabs(initValuesHigh[pdeDim]);// + RandomNumberGenerator::Instance()->ranf());
                    }
                }
                else
                {
                    for(unsigned pdeDim=0; pdeDim<1; pdeDim++)
                    {   // serialised for nodes
                        init_conds[1*i + pdeDim] = fabs(initValuesLow[pdeDim]);// + RandomNumberGenerator::Instance()->ranf());
                    }
                }
            }
            else
            {
                for(unsigned pdeDim=0; pdeDim<1; pdeDim++)
                {   // serialised for nodes
                    init_conds[1*i + pdeDim] = fabs(initValuesLow[pdeDim] + RandomNumberGenerator::Instance()->ranf());
                }
            }

        }
        // PETSc Vec
        Vec initial_condition = PetscTools::CreateVec(init_conds);

        std::vector<AbstractInhomogenousOdeSystemForCoupledPdeSystem*> odeSystem;
        std::vector<boost::shared_ptr<AbstractIvpOdeSolver> > solverSystem;

        LinearParabolicFisherPde<elementDim, spaceDim, 1> pde(diffusionRates,growthRates,carryingCapacities);

        // solver
        LinearParabolicPdeSystemWithCoupledOdeSystemSolver<elementDim,spaceDim,1> solver(p_mesh, &pde, &bcc);

        // solver properties
        double startTime= 0.0;
        double endTime = variables_map["simulation_end_time"].as<double>();

        solver.SetTimes(startTime, endTime);

        solver.SetTimeStep(variables_map["simulation_timestep"].as<double>());
        if(variables_map["simulation_timestep"].as<double>()>variables_map["sampling_timestep"].as<double>())
        {
            solver.SetSamplingTimeStep(variables_map["simulation_timestep"].as<double>());
        }
        else
        {
            solver.SetSamplingTimeStep(variables_map["sampling_timestep"].as<double>());
        }
        
        solver.SetOutputDirectory(variables_map["output_filename"].as<std::string>()+"_chaste"+"_dt_"+boost::lexical_cast<std::string>(variables_map["simulation_timestep"].as<double>()));
        solver.SetInitialCondition(initial_condition);
        // solve
        solver.SolveAndWriteResultsToFile();
    }
    else if(variables_map["simulation_type"].as<std::string>()=="chaste_schnackenberg_equation_D")
    {


        //double a =2;
        //double b = 3;
        //double gamma = 1000;

        double a =0.2;
        double b = 1.3;
        double gamma = 920;

        std::vector<double> diffusionRates = {40,1};

        //std::vector<double> diffusionRates = {1.0,1000};
        std::vector<double> initValues = {a+b,b/((a+b)*(a+b))};
        std::vector<double> bcValues = initValues;
        //std::vector<double> diffusionRates = {1e-4,1e-2};
        //double k1 = 0.1;
        //double k_1 = 0.2;
        //double k2 = 0.3;
        //double k3 = 0.1;

        double MeshStepSize = 1.0;
        std::vector<unsigned> MeshDimensions = {100,100,0};
        //std::vector<double> initValues = {2.0,0.75};
        //std::vector<double> bcValues = {2.0,0.75};
        std::vector<bool> areBCsNeumann = {false,false};

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
        BoundaryConditionsContainer<elementDim,spaceDim,probDim> bcc;
        //std::vector<bool> areNeumannBoundaryConditions(probDim, true);
        std::vector<ConstBoundaryCondition<spaceDim>*> vectorConstBCs;
        
        for(unsigned pdeDim=0; pdeDim<probDim; pdeDim++){
            vectorConstBCs.push_back(new ConstBoundaryCondition<spaceDim>(bcValues[pdeDim]));
        }
        
        for(unsigned pdeDim=0; pdeDim<probDim; pdeDim++)
        {
            if(areBCsNeumann[pdeDim]==false)
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
            
            for(unsigned pdeDim=0; pdeDim<probDim; pdeDim++)
            {   // serialised for nodes
                init_conds[probDim*i + pdeDim] = fabs(initValues[pdeDim] + RandomNumberGenerator::Instance()->ranf());
            }
        
        }
        // PETSc Vec
        Vec initial_condition = PetscTools::CreateVec(init_conds);

        std::vector<AbstractInhomogenousOdeSystemForCoupledPdeSystem*> odeSystem;
        std::vector<boost::shared_ptr<AbstractIvpOdeSolver> > solverSystem;

        LinearParabolicSchnackenbergPde<elementDim, spaceDim, probDim> pde(diffusionRates,a,b,gamma);

        // solver
        LinearParabolicPdeSystemWithCoupledOdeSystemSolver<elementDim,spaceDim,probDim> solver(p_mesh, &pde, &bcc);

        // solver properties
        double startTime= 0.0;
        double endTime = variables_map["simulation_end_time"].as<double>();

        solver.SetTimes(startTime, endTime);

        solver.SetTimeStep(variables_map["simulation_timestep"].as<double>());
        if(variables_map["simulation_timestep"].as<double>()>variables_map["sampling_timestep"].as<double>())
        {
            solver.SetSamplingTimeStep(variables_map["simulation_timestep"].as<double>());
        }
        else
        {
            solver.SetSamplingTimeStep(variables_map["sampling_timestep"].as<double>());
        }
        
        solver.SetOutputDirectory(variables_map["output_filename"].as<std::string>()+"_chaste"+"_dt_"+boost::lexical_cast<std::string>(variables_map["simulation_timestep"].as<double>()));
        solver.SetInitialCondition(initial_condition);
        // solve
        solver.SolveAndWriteResultsToFile();
    }
    else if(variables_map["simulation_type"].as<std::string>()=="chaste_schnackenberg_equation_E")
    {

        double a =0.2;
        double b = 1.3;
        double gamma = 920;

        std::vector<double> diffusionRates = {1e-4,8.4e-7};

        

        double MeshStepSize = 1.0;
        std::vector<unsigned> MeshDimensions = {100,100,0};
        std::vector<double> initValues = {a+b,b/((a+b)*(a+b))};
        std::vector<double> bcValues = {0.0};
        std::vector<double> areBCsNeumann = {true,true};

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
        BoundaryConditionsContainer<elementDim,spaceDim,probDim> bcc;
        
        std::vector<ConstBoundaryCondition<spaceDim>*> vectorConstBCs;
        
        for(unsigned pdeDim=0; pdeDim<probDim; pdeDim++){
            vectorConstBCs.push_back(new ConstBoundaryCondition<spaceDim>(bcValues[pdeDim]));
        }
        
        for(unsigned pdeDim=0; pdeDim<probDim; pdeDim++)
        {
            if(areBCsNeumann[pdeDim]==false)
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
            
            for(unsigned pdeDim=0; pdeDim<probDim; pdeDim++)
            {   // serialised for nodes
                init_conds[probDim*i + pdeDim] = fabs(initValues[pdeDim] + RandomNumberGenerator::Instance()->ranf());
            }
        
        }
        // PETSc Vec
        Vec initial_condition = PetscTools::CreateVec(init_conds);

        std::vector<AbstractInhomogenousOdeSystemForCoupledPdeSystem*> odeSystem;
        std::vector<boost::shared_ptr<AbstractIvpOdeSolver> > solverSystem;

        LinearParabolicSchnackenbergPde<elementDim, spaceDim, probDim> pde(diffusionRates,a,b,gamma);

        // solver
        LinearParabolicPdeSystemWithCoupledOdeSystemSolver<elementDim,spaceDim,probDim> solver(p_mesh, &pde, &bcc);

        // solver properties
        double startTime= 0.0;
        double endTime = variables_map["simulation_end_time"].as<double>();

        solver.SetTimes(startTime, endTime);

        solver.SetTimeStep(variables_map["simulation_timestep"].as<double>());
        if(variables_map["simulation_timestep"].as<double>()>variables_map["sampling_timestep"].as<double>())
        {
            solver.SetSamplingTimeStep(variables_map["simulation_timestep"].as<double>());
        }
        else
        {
            solver.SetSamplingTimeStep(variables_map["sampling_timestep"].as<double>());
        }
        
        solver.SetOutputDirectory(variables_map["output_filename"].as<std::string>()+"_chaste"+"_dt_"+boost::lexical_cast<std::string>(variables_map["simulation_timestep"].as<double>()));
        solver.SetInitialCondition(initial_condition);
        // solve
        solver.SolveAndWriteResultsToFile();
    }
    else if(variables_map["simulation_type"].as<std::string>()=="chaste_schnackenberg_equation_F")
    {


        double a =0.2;
        double b = 1.3;
        double gamma = 220;

        std::vector<double> diffusionRates = {1e-4,119e-4};

        double MeshStepSize = 1.0;
        std::vector<unsigned> MeshDimensions = {100,100,0};
        std::vector<double> initValues = {a+b,b/((a+b)*(a+b))};
        std::vector<double> bcValues = {0.0};
        std::vector<double> areBCsNeumann = {true,true};

     

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
        BoundaryConditionsContainer<elementDim,spaceDim,probDim> bcc;
  
        std::vector<ConstBoundaryCondition<spaceDim>*> vectorConstBCs;
        
        for(unsigned pdeDim=0; pdeDim<probDim; pdeDim++){
            vectorConstBCs.push_back(new ConstBoundaryCondition<spaceDim>(bcValues[pdeDim]));
        }
        
        for(unsigned pdeDim=0; pdeDim<probDim; pdeDim++)
        {
            if(areBCsNeumann[pdeDim]==false)
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
            
            for(unsigned pdeDim=0; pdeDim<probDim; pdeDim++)
            {   // serialised for nodes
                init_conds[probDim*i + pdeDim] = fabs(initValues[pdeDim] + RandomNumberGenerator::Instance()->ranf());
            }
        
        }
        // PETSc Vec
        Vec initial_condition = PetscTools::CreateVec(init_conds);

        std::vector<AbstractInhomogenousOdeSystemForCoupledPdeSystem*> odeSystem;
        std::vector<boost::shared_ptr<AbstractIvpOdeSolver> > solverSystem;

        LinearParabolicSchnackenbergPde<elementDim, spaceDim, probDim> pde(diffusionRates,a,b,gamma);

        // solver
        LinearParabolicPdeSystemWithCoupledOdeSystemSolver<elementDim,spaceDim,probDim> solver(p_mesh, &pde, &bcc);

        // solver properties
        double startTime= 0.0;
        double endTime = variables_map["simulation_end_time"].as<double>();

        solver.SetTimes(startTime, endTime);

        solver.SetTimeStep(variables_map["simulation_timestep"].as<double>());
        if(variables_map["simulation_timestep"].as<double>()>variables_map["sampling_timestep"].as<double>())
        {
            solver.SetSamplingTimeStep(variables_map["simulation_timestep"].as<double>());
        }
        else
        {
            solver.SetSamplingTimeStep(variables_map["sampling_timestep"].as<double>());
        }
        
        solver.SetOutputDirectory(variables_map["output_filename"].as<std::string>()+"_chaste"+"_dt_"+boost::lexical_cast<std::string>(variables_map["simulation_timestep"].as<double>()));
        solver.SetInitialCondition(initial_condition);
        // solve
        solver.SolveAndWriteResultsToFile();
    }
    else if(variables_map["simulation_type"].as<std::string>()=="chaste_kinetic_schnakenberg")
    {


        double D1 = 1e-4;
        double D2 = 1e-2;
        double k1 = 0.1;
        double k_1= 0.2;
        double k2 = 0.3;
        double k3 = 0.1;

        double MeshStepSize = 1.0;
        std::vector<unsigned> MeshDimensions = {100,100,0};
        std::vector<double> initValues = {(k1+k2)/k_1,k2*k_1*k_1/k3*(k1+k2)*(k1+k2)};
        std::vector<double> bcValues = {0.0};
        std::vector<double> areBCsNeumann = {true,true};

     

        // mesh
        TetrahedralMesh<2,2>* p_mesh = new TetrahedralMesh<2,2>();

        switch (2)
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
        BoundaryConditionsContainer<2,2,2> bcc;
  
        std::vector<ConstBoundaryCondition<2>*> vectorConstBCs;
        
        for(unsigned pdeDim=0; pdeDim<2; pdeDim++){
            vectorConstBCs.push_back(new ConstBoundaryCondition<2>(bcValues[pdeDim]));
        }
        
        for(unsigned pdeDim=0; pdeDim<2; pdeDim++)
        {
            if(areBCsNeumann[pdeDim]==false)
            {
                for (TetrahedralMesh<2,2>::BoundaryNodeIterator node_iter = p_mesh->GetBoundaryNodeIteratorBegin();
                node_iter != p_mesh->GetBoundaryNodeIteratorEnd();
                ++node_iter)
                {

                    bcc.AddDirichletBoundaryCondition(*node_iter, vectorConstBCs[pdeDim], pdeDim);
                }
            }else{
                for (TetrahedralMesh<2,2>::BoundaryElementIterator boundary_iter = p_mesh->GetBoundaryElementIteratorBegin();
                boundary_iter != p_mesh->GetBoundaryElementIteratorEnd();
                boundary_iter++)
                {
                    bcc.AddNeumannBoundaryCondition(*boundary_iter, vectorConstBCs[pdeDim], pdeDim);
                }
            }
        }

        // initial conditions
        std::vector<double> init_conds(2*p_mesh->GetNumNodes(),0.0);
        unsigned columnNum = 0;
        unsigned rowNum = 0;
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {   // set as being a random perturbation about the boundary values
            
            for(unsigned pdeDim=0; pdeDim<2; pdeDim++)
            {   // serialised for nodes
                init_conds[2*i + pdeDim] = fabs(initValues[pdeDim] + RandomNumberGenerator::Instance()->ranf());
            }
        
        }
        // PETSc Vec
        Vec initial_condition = PetscTools::CreateVec(init_conds);

        std::vector<AbstractInhomogenousOdeSystemForCoupledPdeSystem*> odeSystem;
        std::vector<boost::shared_ptr<AbstractIvpOdeSolver> > solverSystem;

        //LinearParabolicSchnackenbergPde<elementDim, spaceDim, probDim> pde(diffusionRates,a,b,gamma);
        SchnackenbergCoupledPdeSystem<2> pde(D1, D2, k1,k_1,k2,k3);
        // solver
        LinearParabolicPdeSystemWithCoupledOdeSystemSolver<2,2,2> solver(p_mesh, &pde, &bcc);

        // solver properties
        double startTime= 0.0;
        double endTime = variables_map["simulation_end_time"].as<double>();

        solver.SetTimes(startTime, endTime);

        solver.SetTimeStep(variables_map["simulation_timestep"].as<double>());
        if(variables_map["simulation_timestep"].as<double>()>variables_map["sampling_timestep"].as<double>())
        {
            solver.SetSamplingTimeStep(variables_map["simulation_timestep"].as<double>());
        }
        else
        {
            solver.SetSamplingTimeStep(variables_map["sampling_timestep"].as<double>());
        }
        
        solver.SetOutputDirectory(variables_map["output_filename"].as<std::string>()+"_chaste"+"_dt_"+boost::lexical_cast<std::string>(variables_map["simulation_timestep"].as<double>()));
        solver.SetInitialCondition(initial_condition);
        // solve
        solver.SolveAndWriteResultsToFile();
    }
    else if(variables_map["simulation_type"].as<std::string>()=="a")
    {
 
        double MeshStepSize = 1.0;
        std::vector<unsigned> MeshDimensions = {100,10,0};
        std::vector<double> initValuesHigh = {1.0};
        std::vector<double> initValuesLow = {0.0};
        std::vector<double> bcValues = {0.0};
        std::vector<double> areBCsNeumann = {true};
        std::vector<double> diffusionRates = {100.0};

        bool bothBCs = false;
        bool scaleByDiffusion=false;
        
        unsigned test_id = variables_map["ID"].as<unsigned>();

        switch(test_id) {
        case 0:
            // don't change anything
            break;
        case 1:
            diffusionRates = {0.0};
            break;
        case 2:
            diffusionRates = {1.0};
            break;
        case 3:
            diffusionRates = {10.0};
            break;
        case 4:
            diffusionRates = {100.0};
            break;
        case 5:
            diffusionRates = {1000.0};
            break;
        case 6:
            initValuesHigh = {0.0};
            break;
        case 7:
            initValuesHigh = {1.0};
            break;
        case 8:
            initValuesHigh = {10.0};
            break;
        case 9:
            initValuesHigh = {100.0};
            break;
        case 10:
            bcValues = {0.0};
            areBCsNeumann = {true};
            break;
        case 11:
            bcValues = {1.0};
            areBCsNeumann = {true};
            break;
        case 12:
            bcValues = {0.0};
            areBCsNeumann = {false};
            break;
        case 13:
            bcValues = {0.0};
            areBCsNeumann = {true};
            bothBCs = true;
            break;
        case 14:
            MeshDimensions = {10,10,0};
            break;
        case 15:
            MeshDimensions = {100,10,0};
            break;
        case 16:
            MeshDimensions = {1000,10,0};
            break;

        case 17:
            diffusionRates = {1.0};
            scaleByDiffusion=true;
            break;
        case 18:
            diffusionRates = {10.0};
            scaleByDiffusion=true;
            break;
        case 19:
            diffusionRates = {100.0};
            scaleByDiffusion=true;
            break;
        case 20:
            diffusionRates = {1000.0};
            scaleByDiffusion=true;
            break;

        case 21:
            diffusionRates = {2.0};
            break;
        case 22:
            diffusionRates = {5.0};
            break;
        case 23:
            diffusionRates = {50.0};
            break;

        default:
            // don't change anything
            bothBCs = false;
        }






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
        BoundaryConditionsContainer<elementDim,spaceDim,1> bcc;
        //std::vector<bool> areNeumannBoundaryConditions(1, true);
        std::vector<ConstBoundaryCondition<spaceDim>*> vectorConstBCs;
        
        for(unsigned pdeDim=0; pdeDim<1; pdeDim++){
            vectorConstBCs.push_back(new ConstBoundaryCondition<spaceDim>(bcValues[pdeDim]));
        }
        
        for(unsigned pdeDim=0; pdeDim<1; pdeDim++)
        {
            if(bothBCs==false)
            {
                if(areBCsNeumann[pdeDim]==false)
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
            else
            {
                for (TetrahedralMesh<elementDim,spaceDim>::BoundaryNodeIterator node_iter = p_mesh->GetBoundaryNodeIteratorBegin();
                node_iter != p_mesh->GetBoundaryNodeIteratorEnd();
                ++node_iter)
                {

                    bcc.AddDirichletBoundaryCondition(*node_iter, vectorConstBCs[pdeDim], pdeDim);
                }

                for (TetrahedralMesh<elementDim,spaceDim>::BoundaryElementIterator boundary_iter = p_mesh->GetBoundaryElementIteratorBegin();
                boundary_iter != p_mesh->GetBoundaryElementIteratorEnd();
                boundary_iter++)
                {
                    bcc.AddNeumannBoundaryCondition(*boundary_iter, vectorConstBCs[pdeDim], pdeDim);
                }
            }

        }


        // initial conditions
        std::vector<double> init_conds(1*p_mesh->GetNumNodes(),0.0);
        unsigned columnNum = 0;
        unsigned rowNum = 0;
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {   // set as being a random perturbation about the boundary values
            if(spaceDim==2)
            {
                columnNum = 0;
                rowNum = 0;

                while(i >= rowNum*(MeshDimensions[0]+1))
                {
                    rowNum = rowNum + 1;
                
                }
                
                columnNum = i - (rowNum-1)*(MeshDimensions[0]+1);
                if(columnNum==2 && (rowNum ==3 || rowNum ==4 || rowNum ==5 || rowNum ==6 ))
                {
                    for(unsigned pdeDim=0; pdeDim<1; pdeDim++)
                    {   // serialised for nodes
                        init_conds[1*i + pdeDim] = fabs(initValuesHigh[pdeDim]);// + RandomNumberGenerator::Instance()->ranf());
                    }
                }
                else
                {
                    for(unsigned pdeDim=0; pdeDim<1; pdeDim++)
                    {   // serialised for nodes
                        init_conds[1*i + pdeDim] = fabs(initValuesLow[pdeDim]);// + RandomNumberGenerator::Instance()->ranf());
                    }
                }
            }
            else
            {
                for(unsigned pdeDim=0; pdeDim<1; pdeDim++)
                {   // serialised for nodes
                    init_conds[1*i + pdeDim] = fabs(initValuesLow[pdeDim] + RandomNumberGenerator::Instance()->ranf());
                }
            }

        }
        // PETSc Vec
        Vec initial_condition = PetscTools::CreateVec(init_conds);

        std::vector<AbstractInhomogenousOdeSystemForCoupledPdeSystem*> odeSystem;
        std::vector<boost::shared_ptr<AbstractIvpOdeSolver> > solverSystem;

        LinearParabolicHeatEquationPde<elementDim, spaceDim, 1> pde(diffusionRates);
        // solver
        LinearParabolicPdeSystemWithCoupledOdeSystemSolver<elementDim,spaceDim,1> solver(p_mesh, &pde, &bcc);

        // solver properties
        double startTime= 0.0;
        double endTime = variables_map["simulation_end_time"].as<double>();

        if(scaleByDiffusion==true)
        {
            endTime = endTime/diffusionRates[0];
        }


        solver.SetTimes(startTime, endTime);

        solver.SetTimeStep(variables_map["simulation_timestep"].as<double>());
        if(variables_map["simulation_timestep"].as<double>()>variables_map["sampling_timestep"].as<double>())
        {
            if(scaleByDiffusion==true)
            {
                solver.SetSamplingTimeStep(variables_map["simulation_timestep"].as<double>()/diffusionRates[0]);
            }
            else
            {
                solver.SetSamplingTimeStep(variables_map["simulation_timestep"].as<double>());
            }
       
        }
        else
        {
          
            if(scaleByDiffusion==true)
            {
                solver.SetSamplingTimeStep(variables_map["sampling_timestep"].as<double>()/diffusionRates[0]);
            }
            else
            {
                solver.SetSamplingTimeStep(variables_map["sampling_timestep"].as<double>());
            }
        }
        
        solver.SetOutputDirectory(variables_map["output_filename"].as<std::string>()+"/chaste"+"_test_"+std::to_string(test_id)+"/dt_"+boost::lexical_cast<std::string>(variables_map["simulation_timestep"].as<double>()));
        
        //solver.SetOutputDirectory("ChemChaste/chaste_test/Heat/5"+"/");
        solver.SetInitialCondition(initial_condition);
        // solve
        solver.SolveAndWriteResultsToFile();
        
    }
    else if(variables_map["simulation_type"].as<std::string>()=="chaste_heat_equation_long")
    {

        double MeshStepSize = 1.0;
        std::vector<unsigned> MeshDimensions = {100,10,0};
        std::vector<double> initValuesHigh = {1.0};
        std::vector<double> initValuesLow = {0.0};
        std::vector<double> bcValues = {0.0};
        std::vector<double> areBCsNeumann = {true};
        std::vector<double> diffusionRates = {100.0};



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
        BoundaryConditionsContainer<elementDim,spaceDim,1> bcc;
        std::vector<bool> areNeumannBoundaryConditions(1, true);
        std::vector<ConstBoundaryCondition<spaceDim>*> vectorConstBCs;
        
        for(unsigned pdeDim=0; pdeDim<1; pdeDim++){
            vectorConstBCs.push_back(new ConstBoundaryCondition<spaceDim>(bcValues[pdeDim]));
        }
        
        for(unsigned pdeDim=0; pdeDim<1; pdeDim++)
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
        std::vector<double> init_conds(1*p_mesh->GetNumNodes(),0.0);
        unsigned columnNum = 0;
        unsigned rowNum = 0;
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {   // set as being a random perturbation about the boundary values
            if(spaceDim==2)
            {
                columnNum = 0;
                rowNum = 0;

                while(i >= rowNum*(MeshDimensions[0]+1))
                {
                    rowNum = rowNum + 1;
                
                }
                
                columnNum = i - (rowNum-1)*(MeshDimensions[0]+1);
                if(columnNum<1)
                {
                    for(unsigned pdeDim=0; pdeDim<1; pdeDim++)
                    {   // serialised for nodes
                        init_conds[1*i + pdeDim] = fabs(initValuesHigh[pdeDim]);// + RandomNumberGenerator::Instance()->ranf());
                    }
                }
                else
                {
                    for(unsigned pdeDim=0; pdeDim<1; pdeDim++)
                    {   // serialised for nodes
                        init_conds[1*i + pdeDim] = fabs(initValuesLow[pdeDim]);// + RandomNumberGenerator::Instance()->ranf());
                    }
                }
            }
            else
            {
                for(unsigned pdeDim=0; pdeDim<1; pdeDim++)
                {   // serialised for nodes
                    init_conds[1*i + pdeDim] = fabs(initValuesLow[pdeDim] + RandomNumberGenerator::Instance()->ranf());
                }
            }

        }
        // PETSc Vec
        Vec initial_condition = PetscTools::CreateVec(init_conds);

        std::vector<AbstractInhomogenousOdeSystemForCoupledPdeSystem*> odeSystem;
        std::vector<boost::shared_ptr<AbstractIvpOdeSolver> > solverSystem;

        LinearParabolicHeatEquationPde<elementDim, spaceDim, 1> pde(diffusionRates);

        // solver
        LinearParabolicPdeSystemWithCoupledOdeSystemSolver<elementDim,spaceDim,1> solver(p_mesh, &pde, &bcc);

        // solver properties
        double startTime= 0.0;
        double endTime = variables_map["simulation_end_time"].as<double>();

        solver.SetTimes(startTime, endTime);

        solver.SetTimeStep(variables_map["simulation_timestep"].as<double>());
        if(variables_map["simulation_timestep"].as<double>()>variables_map["sampling_timestep"].as<double>())
        {
            solver.SetSamplingTimeStep(variables_map["simulation_timestep"].as<double>());
        }
        else
        {
            solver.SetSamplingTimeStep(variables_map["sampling_timestep"].as<double>());
        }
        
        solver.SetOutputDirectory(variables_map["output_filename"].as<std::string>()+"_chaste"+"_dt_"+boost::lexical_cast<std::string>(variables_map["simulation_timestep"].as<double>()));
        solver.SetInitialCondition(initial_condition);
        // solve
        solver.SolveAndWriteResultsToFile();
    }
    else if(variables_map["simulation_type"].as<std::string>()=="chaste_heat_equation_source")
    {

        double MeshStepSize = 1.0;
        std::vector<unsigned> MeshDimensions = {100,10,0};
        std::vector<double> initValuesHigh = {1.0};
        std::vector<double> initValuesLow = {0.0};
        std::vector<double> bcValues = {0.0};
        std::vector<double> areBCsNeumann = {true};
        std::vector<double> diffusionRates = {100.0};
        double source_term = 0.1;



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
        BoundaryConditionsContainer<elementDim,spaceDim,1> bcc;
        std::vector<bool> areNeumannBoundaryConditions(1, true);
        std::vector<ConstBoundaryCondition<spaceDim>*> vectorConstBCs;
        
        for(unsigned pdeDim=0; pdeDim<1; pdeDim++){
            vectorConstBCs.push_back(new ConstBoundaryCondition<spaceDim>(bcValues[pdeDim]));
        }
        
        for(unsigned pdeDim=0; pdeDim<1; pdeDim++)
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
        std::vector<double> init_conds(1*p_mesh->GetNumNodes(),0.0);
        unsigned columnNum = 0;
        unsigned rowNum = 0;
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {   // set as being a random perturbation about the boundary values
            if(spaceDim==2)
            {
                columnNum = 0;
                rowNum = 0;

                while(i >= rowNum*(MeshDimensions[0]+1))
                {
                    rowNum = rowNum + 1;
                
                }
                
                columnNum = i - (rowNum-1)*(MeshDimensions[0]+1);
                if(columnNum<1)
                {
                    for(unsigned pdeDim=0; pdeDim<1; pdeDim++)
                    {   // serialised for nodes
                        init_conds[1*i + pdeDim] = fabs(initValuesHigh[pdeDim]);// + RandomNumberGenerator::Instance()->ranf());
                    }
                }
                else
                {
                    for(unsigned pdeDim=0; pdeDim<1; pdeDim++)
                    {   // serialised for nodes
                        init_conds[1*i + pdeDim] = fabs(initValuesLow[pdeDim]);// + RandomNumberGenerator::Instance()->ranf());
                    }
                }
            }
            else
            {
                for(unsigned pdeDim=0; pdeDim<1; pdeDim++)
                {   // serialised for nodes
                    init_conds[1*i + pdeDim] = fabs(initValuesLow[pdeDim] + RandomNumberGenerator::Instance()->ranf());
                }
            }

        }
        // PETSc Vec
        Vec initial_condition = PetscTools::CreateVec(init_conds);

        std::vector<AbstractInhomogenousOdeSystemForCoupledPdeSystem*> odeSystem;
        std::vector<boost::shared_ptr<AbstractIvpOdeSolver> > solverSystem;

        LinearParabolicHeatEquationWithSourcePde<elementDim, spaceDim, 1> pde(diffusionRates,source_term);

        // solver
        LinearParabolicPdeSystemWithCoupledOdeSystemSolver<elementDim,spaceDim,1> solver(p_mesh, &pde, &bcc);

        // solver properties
        double startTime= 0.0;
        double endTime = variables_map["simulation_end_time"].as<double>();

        solver.SetTimes(startTime, endTime);

        solver.SetTimeStep(variables_map["simulation_timestep"].as<double>());
        if(variables_map["simulation_timestep"].as<double>()>variables_map["sampling_timestep"].as<double>())
        {
            solver.SetSamplingTimeStep(variables_map["simulation_timestep"].as<double>());
        }
        else
        {
            solver.SetSamplingTimeStep(variables_map["sampling_timestep"].as<double>());
        }
        
        solver.SetOutputDirectory(variables_map["output_filename"].as<std::string>()+"_chaste"+"_dt_"+boost::lexical_cast<std::string>(variables_map["simulation_timestep"].as<double>()));
        solver.SetInitialCondition(initial_condition);
        // solve
        solver.SolveAndWriteResultsToFile();
    }
    else if(variables_map["simulation_type"].as<std::string>()=="chaste_heat_equation_source_long")
    {

        double MeshStepSize = 1.0;
        std::vector<unsigned> MeshDimensions = {100,10,0};
        std::vector<double> initValuesHigh = {1.0};
        std::vector<double> initValuesLow = {0.0};
        std::vector<double> bcValues = {0.0};
        std::vector<double> areBCsNeumann = {true};
        std::vector<double> diffusionRates = {100.0};
        double source_term = 0.1;


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
        BoundaryConditionsContainer<elementDim,spaceDim,1> bcc;
        std::vector<bool> areNeumannBoundaryConditions(1, true);
        std::vector<ConstBoundaryCondition<spaceDim>*> vectorConstBCs;
        
        for(unsigned pdeDim=0; pdeDim<1; pdeDim++){
            vectorConstBCs.push_back(new ConstBoundaryCondition<spaceDim>(bcValues[pdeDim]));
        }
        
        for(unsigned pdeDim=0; pdeDim<1; pdeDim++)
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
        std::vector<double> init_conds(1*p_mesh->GetNumNodes(),0.0);
        unsigned columnNum = 0;
        unsigned rowNum = 0;
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {   // set as being a random perturbation about the boundary values
            if(spaceDim==2)
            {
                columnNum = 0;
                rowNum = 0;

                while(i >= rowNum*(MeshDimensions[0]+1))
                {
                    rowNum = rowNum + 1;
                
                }
                
                columnNum = i - (rowNum-1)*(MeshDimensions[0]+1);
                if(columnNum<1)
                {
                    for(unsigned pdeDim=0; pdeDim<1; pdeDim++)
                    {   // serialised for nodes
                        init_conds[1*i + pdeDim] = fabs(initValuesHigh[pdeDim]);// + RandomNumberGenerator::Instance()->ranf());
                    }
                }
                else
                {
                    for(unsigned pdeDim=0; pdeDim<1; pdeDim++)
                    {   // serialised for nodes
                        init_conds[1*i + pdeDim] = fabs(initValuesLow[pdeDim]);// + RandomNumberGenerator::Instance()->ranf());
                    }
                }
            }
            else
            {
                for(unsigned pdeDim=0; pdeDim<1; pdeDim++)
                {   // serialised for nodes
                    init_conds[1*i + pdeDim] = fabs(initValuesLow[pdeDim] + RandomNumberGenerator::Instance()->ranf());
                }
            }

        }
        // PETSc Vec
        Vec initial_condition = PetscTools::CreateVec(init_conds);

        std::vector<AbstractInhomogenousOdeSystemForCoupledPdeSystem*> odeSystem;
        std::vector<boost::shared_ptr<AbstractIvpOdeSolver> > solverSystem;

        LinearParabolicHeatEquationWithSourcePde<elementDim, spaceDim, 1> pde(diffusionRates,source_term);

        // solver
        LinearParabolicPdeSystemWithCoupledOdeSystemSolver<elementDim,spaceDim,1> solver(p_mesh, &pde, &bcc);

        // solver properties
        double startTime= 0.0;
        double endTime = variables_map["simulation_end_time"].as<double>();

        solver.SetTimes(startTime, endTime);

        solver.SetTimeStep(variables_map["simulation_timestep"].as<double>());
        if(variables_map["simulation_timestep"].as<double>()>variables_map["sampling_timestep"].as<double>())
        {
            solver.SetSamplingTimeStep(variables_map["simulation_timestep"].as<double>());
        }
        else
        {
            solver.SetSamplingTimeStep(variables_map["sampling_timestep"].as<double>());
        }
        
        solver.SetOutputDirectory(variables_map["output_filename"].as<std::string>()+"_chaste"+"_dt_"+boost::lexical_cast<std::string>(variables_map["simulation_timestep"].as<double>()));
        solver.SetInitialCondition(initial_condition);
        // solve
        solver.SolveAndWriteResultsToFile();
    }
    
    else
    {
        std::cout<<"Simulation type ("<<variables_map["simulation_type"].as<std::string>()<<") not recognised"<<std::endl;
    }
    
}

void WriteSimulationToChemChasteLog(std::string config_filename,unsigned run_id)
{

    std::ofstream outfile;

    outfile.open("ChemChaste_log.txt", std::ios_base::app); // append instead of overwrite
    outfile << config_filename << ", " << run_id <<"\n";  
    return; 
}
