#ifndef TESTCHEMICALCELLOBJECT_HPP_
#define TESTCHEMICALCELLOBJECT_HPP_

// chaste includes
#include "ChasteHeaders.hpp"

// general includes
#include "GeneralHeaders.hpp"

//ChemChaste includes
#include "ChemChasteHeaders.hpp"

// test specific include
#include "ChemicalStructuresForTests.hpp"

#include "ChemicalCellFromFile.hpp"

struct ControlStruct {
    bool ChemicalCellObjectReadFromFile = false;
    bool ChemicalCellObjectInPde = true;
} control;

class TestChemicalCellObject : public AbstractCellBasedTestSuite
{
public:

    void TestChemicalCellObjectReadFromFile()
    {
        if(control.ChemicalCellObjectReadFromFile)
        {
            // Variables for the user modify
            std::string dataFileRoot = "/home/chaste/projects/ChemChaste/src/Data/CoupledCellPdeTest/Cell/";
            std::string initialConditionFilename = "InitialCellConcentrations.csv";
            std::string membranePropertyFilename = "MembraneReactions.txt";
            std::string transportPropertyFilename = "TransportReactions.txt";
            std::string srnFilename = "Srn.txt";
            std::string cellCycleFilename = "SpeciesThreshold.csv";

            ChemicalCellFromFile* p_cell_reader = new ChemicalCellFromFile(
                                    dataFileRoot+cellCycleFilename, 
                                    dataFileRoot+srnFilename,
                                    dataFileRoot+initialConditionFilename,
                                    dataFileRoot+transportPropertyFilename,
                                    dataFileRoot+membranePropertyFilename);

            std::cout<<"===================================================="<<std::endl;

            CellPtr p_cell = p_cell_reader -> GetCellPtr();

            std::cout<<"Chemical cell property"<<std::endl;
        
            boost::shared_ptr<ChemicalCellProperty> p_cell_chemical = p_cell_reader -> GetChemicalCellProperty();

            std::cout<<"Transport cell property"<<std::endl;

            boost::shared_ptr<TransportCellProperty> p_cell_transport = p_cell_reader -> GetTransportCellProperty();

            StateVariableRegister* p_transport_bulk_register = p_cell_transport -> GetBulkStateVariableRegister();
            std::cout<<"size register: "<< p_transport_bulk_register -> GetNumberOfStateVariables()<<std::endl;
            for(unsigned i=0; i<p_transport_bulk_register -> GetNumberOfStateVariables(); i++)
            {
                std::cout<<"bulk: "<< p_transport_bulk_register -> RetrieveStateVariableName(i)<<std::endl;
            }
            
            StateVariableRegister* p_transport_cell_register = p_cell_transport -> GetCellStateVariableRegister();
            std::cout<<"size register: "<< p_transport_cell_register -> GetNumberOfStateVariables()<<std::endl;
            for(unsigned i=0; i<p_transport_cell_register -> GetNumberOfStateVariables(); i++)
            {
                std::cout<<"cell: "<< p_transport_cell_register -> RetrieveStateVariableName(i)<<std::endl;
            }



            std::vector<double> change_environment_concentration_vector= {0.0,0.0,0.0};
            std::vector<double>  change_cell_concentration_vector= {0.0,0.0,0.0};

            std::vector<double> environment_concentration_at_point= {1.0,2.0,3.0};
            std::vector<double> cell_concentration= p_cell_chemical -> GetCellConcentrationVector();

            p_cell_transport -> PerformTransportSystem(environment_concentration_at_point,cell_concentration,change_environment_concentration_vector,change_cell_concentration_vector);


            std::cout<<"Membrane cell property"<<std::endl;

            boost::shared_ptr<MembraneCellProperty> p_cell_membrane = p_cell_reader -> GetMembraneCellProperty();

            StateVariableRegister* p_membrane_bulk_register = p_cell_membrane -> GetBulkStateVariableRegister();
            std::cout<<"size register: "<< p_membrane_bulk_register -> GetNumberOfStateVariables()<<std::endl;
            for(unsigned i=0; i<p_membrane_bulk_register -> GetNumberOfStateVariables(); i++)
            {
                std::cout<<"bulk: "<< p_membrane_bulk_register -> RetrieveStateVariableName(i)<<std::endl;
            }

            StateVariableRegister* p_membrane_cell_register = p_cell_membrane -> GetCellStateVariableRegister();
            std::cout<<"size register: "<< p_membrane_cell_register -> GetNumberOfStateVariables()<<std::endl;
            for(unsigned i=0; i<p_membrane_cell_register -> GetNumberOfStateVariables(); i++)
            {
                std::cout<<"cell: "<< p_membrane_cell_register -> RetrieveStateVariableName(i)<<std::endl;
            }
            


            
            std::cout<<"Membrane thickness :"<<p_cell_membrane -> GetMembraneThickness()<<std::endl;
        /*
            AbstractTransportOdeSystem transportOde(p_cell_transport -> GetTransportReactionSystem());
            EulerIvpOdeSolver euler_solver;
            std::vector<double> initial_condition = {1.0, 1.0};

            OdeSolution solutions = euler_solver.Solve(&transportOde, initial_condition, 0, 1, 0.01, 0.1);
            for (unsigned i=0; i<solutions.rGetTimes().size(); i++)
            {
                std::cout << solutions.rGetTimes()[i] << " " << solutions.rGetSolutions()[i][0] << " " << solutions.rGetSolutions()[i][1]<< "\n";
            }

*/
            ChemicalSrnModel* p_model = static_cast<ChemicalSrnModel*>(p_cell->GetSrnModel());
          
            AbstractChemistry* this_cell_chemistry = p_model->GetCellChemistry();
            unsigned numberOfChemicals = this_cell_chemistry->GetNumberChemicals();
            std::cout<<"number of chemicals: "<<numberOfChemicals<<std::endl;





        }
    }

    void TestChemicalCellObjectInPde()
    {
        if(control.ChemicalCellObjectInPde)
        {
            // Variables for the user modify
            std::string dataFileRoot = "/home/chaste/projects/ChemChaste/src/Data/CoupledCellPdeTest/Cell/";
            std::string initialConditionFilename = "InitialCellConcentrations.csv";
            std::string membranePropertyFilename = "MembraneReactions.txt";
            std::string transportPropertyFilename = "TransportReactions.txt";
            std::string srnFilename = "Srn.txt";
            std::string cellCycleFilename = "SpeciesThreshold.csv";


            // make mesh of cell with associated mesh based cell population
            HoneycombMeshGenerator generator(2, 2);    // Parameters are: cells across, cells up
            MutableMesh<2,2>* p_mesh = generator.GetMesh();

            std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

            std::vector<CellPtr> cells;
            // assume cell at each node in cell layer mesh
            for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
            {
                // provide each cell with a transport cell property and membrane property, cell cycle, wild type states

                ChemicalCellFromFile* p_cell_reader = new ChemicalCellFromFile(
                                    dataFileRoot+cellCycleFilename, 
                                    dataFileRoot+srnFilename,
                                    dataFileRoot+initialConditionFilename,
                                    dataFileRoot+transportPropertyFilename,
                                    dataFileRoot+membranePropertyFilename);


                //CellPtr p_cell = ;

                cells.push_back(p_cell_reader -> GetCellPtr());
                
                
            }
        
            MeshBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices);
            
            // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
            ChastePoint<2> lower(-4.0, -4.0);
            ChastePoint<2> upper(10.0, 10.0);
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
        
            simulator.SetOutputDirectory("TestChemicalCellObjectInPde");
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

};

#endif