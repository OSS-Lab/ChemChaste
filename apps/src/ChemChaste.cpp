
//#include "AbstractCellBasedTestSuite.hpp"

#include <cxxtest/TestSuite.h>

//#include <cxxtest/GlobalFixture.h>
#include "CheckpointArchiveTypes.hpp"
#include "SmartPointers.hpp"

#include <string>
#include <fstream>
#include <iostream>

//ChemChaste includes
#include "ChemChasteExecutableHeaders.hpp"

#include <boost/lexical_cast.hpp>

#include "ExecutableSupport.hpp"

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>


void SetupSingletons();
void DestroySingletons();
void SetupAndRunSimulation(unsigned simulation_id,boost::program_options::variables_map& variables_map);

int main(int argc, char *argv[])
{
    std::cout<<"ChemChaste Executable"<<std::endl;

    //std::string relative_config_path = "../../DataInput/";
    //std::string relative_config_path = "/projects/DataInput/";
    // This sets up PETSc and prints out copyright information, etc.
    ExecutableSupport::StartupWithoutShowingCopyright(&argc, &argv);

    int exit_code = ExecutableSupport::EXIT_OK;

    // You should put all the main code within a try-catch, to ensure that
    // you clean up PETSc before quitting.
    try
    {

         // Declare a group of options that will be 
        // allowed only on command line
        boost::program_options::options_description general("General options");
        general.add_options()
        ("help", "produce help message")
        ("config", boost::program_options::value<std::string>()->default_value("/home/chaste/projects/ChemChaste/DataInput/ChemChasteConfig.txt"), "Config file")
        ("ID", boost::program_options::value<unsigned>()->default_value(0),"ID of the simulation (for output)")   
        ("simulation_type,S", boost::program_options::value<std::string>()->default_value("domain_only"),"Define the type of simualtion to run: domain_only, coupled_cell")
        ;
        
        // Declare a group of options that will be 
        // allowed both on command line and in
        // config file
        boost::program_options::options_description config("Configuration");
        config.add_options()
        ("output_filename,F", boost::program_options::value<std::string>()->default_value("ChemChasteExecutable"),"Base filename to output the executable results")
        ("simulation_end_time,T", boost::program_options::value<double>()->default_value(10),"Temporal length of the simualtion")
        ("domain_file_root",boost::program_options::value<std::string>()->default_value(""),"Absolute root of the domain folder")
        ("domain_file",boost::program_options::value<std::string>()->default_value(""),"CSV file containing the domain topology")
        ("domain_key_file",boost::program_options::value<std::string>()->default_value(""),"Key file relating labels in domain_file topology to subdomains")
        ("ode_file",boost::program_options::value<std::string>()->default_value(""),"CSV file containing the ode topology")
        ("ode_key_file",boost::program_options::value<std::string>()->default_value(""),"Key file relating labels in ode_file topology to ode domain systems")
        ("diffusion_database",boost::program_options::value<std::string>()->default_value(""),"CSV containing the standard diffusion rates of the species in each of the subdomains")
        ("initial_conditions",boost::program_options::value<std::string>()->default_value(""),"CSV file containing the initial conditions of each species in each of the subdomains specified in domain_file")
        ("boundary_conditions",boost::program_options::value<std::string>()->default_value(""),"CSV file containing the outer domain boundary conditions for each species")
        ("cell_file_root",boost::program_options::value<std::string>()->default_value(""),"Absolute root of the cell folder")
        ("cell_initial_conditions",boost::program_options::value<std::string>()->default_value(""),"CSV file containing the initial cell internal concentrations")
        ("membrane_property",boost::program_options::value<std::string>()->default_value(""),"Membrane file, reactions bound at the membrane")
        ("transport_property",boost::program_options::value<std::string>()->default_value(""),"Transport file, reactions across the cell membrane")
        ("srn_file",boost::program_options::value<std::string>()->default_value(""),"Subcellular reaction network file")
        ("cell_cycle_file",boost::program_options::value<std::string>()->default_value(""),"Species threshold")
        ("number_cells_across",boost::program_options::value<unsigned>()->default_value(1),"Width of honeycomb cell mesh")
        ("number_cells_high",boost::program_options::value<unsigned>()->default_value(1),"Height of honeycomb cell mesh")
        ("node_cutoff_length",boost::program_options::value<double>()->default_value(1.5),"Connectivity length of nodes only mesh")
        ("cell_mesh_origin",boost::program_options::value<double>()->default_value(-4.0),"Mesh overlapping point")
        ("linear_force_cutoff",boost::program_options::value<double>()->default_value(1.5),"Cut off length of linear spring force between cells")
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

        std::cout<<variables_map["cell_cycle_file"].as<std::string>()<<std::endl;

        SetupSingletons();
        SetupAndRunSimulation(variables_map["ID"].as<unsigned>(),variables_map);
        DestroySingletons();

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
                            variables_map["cell_file_root"].as<std::string>()+variables_map["membrane_property"].as<std::string>());

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
    
    
    boost::shared_ptr<ParabolicBoxDomainPdeSystemModifier<elementDim,spaceDim,probDim>> p_pde_modifier(new ParabolicBoxDomainPdeSystemModifier<elementDim,spaceDim,probDim>(p_Pde_field, p_cuboid));
    
    boost::shared_ptr<ChemicalTrackingModifier<elementDim,spaceDim>> p_chemical_tracking_modifier(new ChemicalTrackingModifier<elementDim,spaceDim>());
    

    NodesOnlyMesh<spaceDim> mesh;

    mesh.ConstructNodesWithoutMesh(*p_mesh, variables_map["node_cutoff_length"].as<double>());

    NodeBasedCellPopulation<spaceDim> cell_population(mesh, cells);   

    // writers

    cell_population.AddCellWriter<CellIdWriter>();
    cell_population.AddCellWriter<CellAgesWriter>();
    cell_population.AddCellWriter<CellLocationIndexWriter>();

    
    /*
    for(unsigned i=0; i<chemicalCellSpeciesNames.size(); i++)
    {
        boost::shared_ptr<CellDataItemWriter<elementDim,spaceDim>> dataWriter(new CellDataItemWriter<elementDim,spaceDim>(chemicalCellSpeciesNames[i]));
        cell_population.AddCellWriter(dataWriter);
    }
    */

    OffLatticeSimulation<spaceDim> simulator(cell_population);

    simulator.AddSimulationModifier(p_pde_modifier);

    simulator.AddSimulationModifier(p_chemical_tracking_modifier);
    //std::cout<<"Ouput directory: "<<variables_map["output_filename"].as<std::string>()+boost::lexical_cast<std::string>(simulation_id)<<std::endl;
    simulator.SetOutputDirectory("TestChemChastePaperExecutable/"+boost::lexical_cast<std::string>(simulation_id));
    simulator.SetEndTime(variables_map["simulation_end_time"].as<double>());

    MAKE_PTR(GeneralisedLinearSpringForce<spaceDim>, p_linear_force);
    p_linear_force->SetCutOffLength(variables_map["linear_force_cutoff"].as<double>());
    simulator.AddForce(p_linear_force);

    std::cout<<"Running simulation: "<<simulation_id<<std::endl;
    simulator.Solve();
}


