#ifndef CHEMICALSTRUCTURESFORTESTS_HPP_
#define CHEMICALSTRUCTURESFORTESTS_HPP_

// chaste includes
#include "ChasteHeaders.hpp"

// general includes
#include "GeneralHeaders.hpp"

//ChemChaste includes
#include "ChemChasteHeaders.hpp"

#include "SimpleChemicalThresholdCellCycleModel.hpp"

class ChemicalStructuresForTests
{
    private:

        // r1: 2Bulk_U + Bulk_V -> 3Bulk_U              forwardRate = 0.1
        // r2: 0 <-> Bulk_U                             forwardRate = 0.1 reverseRate = 0.2
        // r3: 0 -> Bulk_V                              forwardRate = 0.3
        std::vector<double> pde_reaction_rate_vector = {1.0,0.1,0.2,0.3};

        // r1: Bulk_U <-> Cell_U                        forwardRate = 0.0 reverseRate = 0.0
        // r2: Bulk_V <-> Cell_V                        forwardRate = 0.0 reverseRate = 0.0
        std::vector<double> transport_reaction_rate_vector = {1.0,0.0,1.0,0.0};

        // r1: Bulk_V <-> 2Bulk_V | 2Cell_U + Cell_V <-> Cell_U  forwardRate = 0.0 reverseRate = 0.0
        std::vector<double> membrane_reaction_rate_vector = {0.0,0.0};

        // r1: Cell_U <-> Cell_Biomass                  forwardRate = 0.0 reverseRate = 0.0
        std::vector<double> srn_reaction_rate_vector = {0.1,0.0};

        // SimpleChemicalCellCycleModel {U,Biomass,V}
        std::vector<double> MaximumSpeciesThreshold = {2.0,2.0,0};
        std::vector<double> MinimumSpeciesThreshold = {0,0.1,0};
        std::vector<bool> MaximumThresholdCheck = {false,true,false};
        std::vector<bool> MinimumThresholdCheck = {false,true,false};

        // cells
        // cellObjectA : {U,V} initial cell/srn conditions, SchnackenbergSrnModel
        std::vector<double> cellObjectAConcentrationVector = {1.0,0.0};

        // cellObjectB : {U,V} initial cell/srn conditions, SchnackenbergSrnModel
        std::vector<double> cellObjectBConcentrationVector = {1.0,0.0};

        // cellObjectB : {U,V} initial cell, no srn model
        std::vector<double> nullCellObjectConcentrationVector = {0.0,0.0};

        // chemicalCellA, inital cell/srn conditions, custom srn model
        std::vector<std::string> chemicalCellASpeciesNames = {"U","Biomass","V"}; // needs reflecting Srn
        std::vector<double> chemicalCellAConcentrationVector = {5.0,1.9,1.0};

        // chemicalCellB, inital cell/srn conditions, custom srn model
        std::vector<std::string> chemicalCellBSpeciesNames = {"U","V"}; // needs reflecting Srn
        std::vector<double> chemicalCellBConcentrationVector = {0.0,1.0};

        // nullChemicalCell, inital conditions, no srn model
        std::vector<std::string> nullChemicalCellSpeciesNames = {"U","V"};
        std::vector<double> nullChemicalCellConcentrationVector = {0.0,0.0};

        AbstractReactionSystem* mp_pde_chemical_reaction_system;

        AbstractChemicalOdeSystem* mp_pde_chemical_ode;

        AbstractTransportReactionSystem* mp_transport_system;

        AbstractTransportReactionSystem* mp_null_transport_system;

        AbstractTransportOdeSystem* mp_transport_ode;

        AbstractTransportOdeSystem* mp_null_transport_ode;

        AbstractMembraneReactionSystem* mp_membrane_system;

        AbstractMembraneReactionSystem* mp_null_membrane_system;

        AbstractMembraneOdeSystem* mp_membrane_ode;

        AbstractMembraneOdeSystem* mp_null_membrane_ode;

        ChemicalDomainFieldTemplated<2,2,2>* mpChemicalDomain;

        ChemicalDomainFieldForCellCoupling<2,2,2>* mpChemicalDomainForCellCoupling;

        AbstractReactionSystem* mp_srn_chemical_reaction_system;

        AbstractChemicalOdeSystem* mp_srn_chemical_ode;

        ChemicalSrnModel*   mpChemicalSrnModel;

        std::vector<std::string> mThresholdNames;

        std::vector<double> mMaxThresholdVector;

        std::vector<double> mMinThresholdVector;

        SimpleChemicalThresholdCellCycleModel*  mpSimpleChemicalThresholdCellCycleModel;

        boost::shared_ptr<ChemicalTrackingModifier<2,2>> mpChemicalTrackingModifier;


    public:

        ChemicalStructuresForTests()
        {
        }

        ~ChemicalStructuresForTests()
        {
        }

        void SetUpPdeChemicalReactionSystem()
        {
            // r1: 2U + V -> 3U     forwardRate = 0.1
            // r2: 0 <-> U          forwardRate = 0.1 reverseRate = 0.2
            // r3: 0 -> V           forwardRate = 0.3

            AbstractChemistry* p_system_chemistry = new AbstractChemistry();

            std::vector<AbstractChemical*> p_substrates_1 = std::vector<AbstractChemical*>();
            std::vector<AbstractChemical*> p_substrates_2 = std::vector<AbstractChemical*>();
            std::vector<AbstractChemical*> p_substrates_3 = std::vector<AbstractChemical*>();
            std::vector<AbstractChemical*> p_products_1 = std::vector<AbstractChemical*>();
            std::vector<AbstractChemical*> p_products_2 = std::vector<AbstractChemical*>();
            std::vector<AbstractChemical*> p_products_3 = std::vector<AbstractChemical*>();

            std::vector<unsigned> stoich_substrates_1 = std::vector<unsigned>();
            std::vector<unsigned> stoich_products_1 = std::vector<unsigned>();
            std::vector<unsigned> stoich_substrates_2 = std::vector<unsigned>();
            std::vector<unsigned> stoich_products_2 = std::vector<unsigned>();
            std::vector<unsigned> stoich_substrates_3 = std::vector<unsigned>();
            std::vector<unsigned> stoich_products_3 = std::vector<unsigned>();

            AbstractChemical *p_chemical_U = new AbstractChemical("U");
            p_system_chemistry -> AddChemical(p_chemical_U);
            // add U to reactions
            p_substrates_1.push_back(p_chemical_U);
            stoich_substrates_1.push_back(2);
            p_products_1.push_back(p_chemical_U);
            stoich_products_1.push_back(3);
            p_products_2.push_back(p_chemical_U);
            stoich_products_2.push_back(1);

            AbstractChemical *p_chemical_V = new AbstractChemical("V");
            p_system_chemistry -> AddChemical(p_chemical_V);
            // add V to reactions
            p_substrates_1.push_back(p_chemical_V);
            stoich_substrates_1.push_back(1);
            p_products_3.push_back(p_chemical_V);
            stoich_products_3.push_back(1);

            double reaction_1_rate = pde_reaction_rate_vector[0];
            double reaction_2_forward_rate = pde_reaction_rate_vector[1];
            double reaction_2_reverse_rate = pde_reaction_rate_vector[2];
            double reaction_3_rate = pde_reaction_rate_vector[3];

            MassActionReaction* p_mass_action_reaction_1 = new MassActionReaction(p_substrates_1, p_products_1, stoich_substrates_1, stoich_products_1,false,false,reaction_1_rate);
            MassActionReaction* p_mass_action_reaction_2 = new MassActionReaction(p_substrates_2, p_products_2, stoich_substrates_2, stoich_products_2,false, true, reaction_2_forward_rate, reaction_2_reverse_rate);
            MassActionReaction* p_mass_action_reaction_3 = new MassActionReaction(p_substrates_3, p_products_3, stoich_substrates_3, stoich_products_3,false, false,reaction_3_rate);

            std::vector<AbstractReaction*> p_mass_action_reaction_vector;
            p_mass_action_reaction_vector.push_back(p_mass_action_reaction_1);
            p_mass_action_reaction_vector.push_back(p_mass_action_reaction_2);
            p_mass_action_reaction_vector.push_back(p_mass_action_reaction_3);

            AbstractReactionSystem* p_mass_action_reaction_system = new AbstractReactionSystem(p_system_chemistry, p_mass_action_reaction_vector);
            
            SetPtrPdeChemicalReactionSystem(p_mass_action_reaction_system);
        }

        void SetUpTransportReactionSystem()
        {
            AbstractChemistry* p_bulk_system_chemistry = new AbstractChemistry();

            AbstractChemistry* p_cell_system_chemistry = new AbstractChemistry();

            std::vector<AbstractChemical*> p_bulk_2 = std::vector<AbstractChemical*>();
            std::vector<AbstractChemical*> p_bulk_3 = std::vector<AbstractChemical*>();
            std::vector<AbstractChemical*> p_cell_2 = std::vector<AbstractChemical*>();
            std::vector<AbstractChemical*> p_cell_3 = std::vector<AbstractChemical*>();


            std::vector<unsigned> stoich_bulk_2 = std::vector<unsigned>();
            std::vector<unsigned> stoich_cell_2 = std::vector<unsigned>();
            std::vector<unsigned> stoich_bulk_3 = std::vector<unsigned>();
            std::vector<unsigned> stoich_cell_3 = std::vector<unsigned>();


            // r1: 2U + V -> 3U     forwardRate = 0.1
            // r2: U <-> U          forwardRate = 0.1 reverseRate = 0.2
            // r3: V <- V           forwardRate = 0.3

            AbstractChemical *p_chemical_U = new AbstractChemical("U");
            p_bulk_system_chemistry -> AddChemical(p_chemical_U);
            p_cell_system_chemistry -> AddChemical(p_chemical_U);
            // add U to reactions

            p_bulk_2.push_back(p_chemical_U);
            stoich_bulk_2.push_back(1);
            p_cell_2.push_back(p_chemical_U);
            stoich_cell_2.push_back(1);

            AbstractChemical *p_chemical_V = new AbstractChemical("V");
            p_bulk_system_chemistry -> AddChemical(p_chemical_V);
            p_cell_system_chemistry -> AddChemical(p_chemical_V);
            // add U to reactions

            p_bulk_3.push_back(p_chemical_V);
            stoich_bulk_3.push_back(1);
            p_cell_3.push_back(p_chemical_V);
            stoich_cell_3.push_back(1);

            // for the sake of testing multiple species in the bulk which do not take part in the transport process
            //AbstractChemical *p_chemical_C = new AbstractChemical("C");
            //p_bulk_system_chemistry -> AddChemical(p_chemical_C);

            //double reaction_2_forward_rate = 1.0;
           // double reaction_2_reverse_rate = 2.0;
            //double reaction_3_rate = 5.0;

            double reaction_2_forward_rate = transport_reaction_rate_vector[0];
            double reaction_2_reverse_rate = transport_reaction_rate_vector[1];
            double reaction_3_forward_rate = transport_reaction_rate_vector[2];
            double reaction_3_reverse_rate = transport_reaction_rate_vector[3];

            //AbstractTransportReaction* p_reaction_1 = new AbstractTransportReaction(p_bulk_1, p_cell_1, stoich_bulk_1, stoich_cell_1,reaction_1_rate);
            MassActionTransportReaction* p_reaction_2 = new MassActionTransportReaction(p_bulk_2, p_cell_2, stoich_bulk_2, stoich_cell_2,false,true,reaction_2_forward_rate, reaction_2_reverse_rate);
            MassActionTransportReaction* p_reaction_3 = new MassActionTransportReaction(p_bulk_3, p_cell_3, stoich_bulk_3, stoich_cell_3,false,true,reaction_3_forward_rate, reaction_3_reverse_rate);

            std::vector<AbstractTransportReaction*> p_transport_reaction_vector;
            //p_transport_reaction_vector.push_back(p_reaction_1);
            p_transport_reaction_vector.push_back(p_reaction_2);
            p_transport_reaction_vector.push_back(p_reaction_3);

            AbstractTransportReactionSystem* p_transport_reaction_system = new AbstractTransportReactionSystem(p_bulk_system_chemistry, p_cell_system_chemistry, p_transport_reaction_vector);

            SetPtrTransportReactionSystem(p_transport_reaction_system);
        }

        void SetUpNullTransportReactionSystem()
        {
            AbstractChemistry* p_bulk_system_chemistry = new AbstractChemistry();

            AbstractChemistry* p_cell_system_chemistry = new AbstractChemistry();

            std::vector<AbstractChemical*> p_bulk_1 = std::vector<AbstractChemical*>();
            std::vector<AbstractChemical*> p_bulk_2 = std::vector<AbstractChemical*>();
            std::vector<AbstractChemical*> p_bulk_3 = std::vector<AbstractChemical*>();
            std::vector<AbstractChemical*> p_cell_1 = std::vector<AbstractChemical*>();
            std::vector<AbstractChemical*> p_cell_2 = std::vector<AbstractChemical*>();
            std::vector<AbstractChemical*> p_cell_3 = std::vector<AbstractChemical*>();

            std::vector<unsigned> stoich_bulk_1 = std::vector<unsigned>();
            std::vector<unsigned> stoich_cell_1 = std::vector<unsigned>();
            std::vector<unsigned> stoich_bulk_2 = std::vector<unsigned>();
            std::vector<unsigned> stoich_cell_2 = std::vector<unsigned>();
            std::vector<unsigned> stoich_bulk_3 = std::vector<unsigned>();
            std::vector<unsigned> stoich_cell_3 = std::vector<unsigned>();


            // r1: 2U + V -> 3U     forwardRate = 0.1
            // r2: U <-> U          forwardRate = 0.1 reverseRate = 0.2
            // r3: V <- V           forwardRate = 0.3

            AbstractChemical *p_chemical_U = new AbstractChemical("U");
            p_bulk_system_chemistry -> AddChemical(p_chemical_U);
            p_cell_system_chemistry -> AddChemical(p_chemical_U);
            // add U to reactions
            p_bulk_1.push_back(p_chemical_U);
            stoich_bulk_1.push_back(2);
            p_cell_1.push_back(p_chemical_U);
            stoich_cell_1.push_back(3);
            p_bulk_2.push_back(p_chemical_U);
            stoich_bulk_2.push_back(1);
            p_cell_2.push_back(p_chemical_U);
            stoich_cell_2.push_back(1);

            AbstractChemical *p_chemical_V = new AbstractChemical("V");
            p_bulk_system_chemistry -> AddChemical(p_chemical_V);
            p_cell_system_chemistry -> AddChemical(p_chemical_V);
            // add U to reactions
            p_bulk_1.push_back(p_chemical_V);
            stoich_bulk_1.push_back(1);
            p_bulk_3.push_back(p_chemical_V);
            stoich_bulk_3.push_back(1);
            p_cell_3.push_back(p_chemical_V);
            stoich_cell_3.push_back(1);

            // for the sake of testing multiple species in the bulk which do not take part in the transport process
            //AbstractChemical *p_chemical_C = new AbstractChemical("C");
            //p_bulk_system_chemistry -> AddChemical(p_chemical_C);

            double reaction_1_rate = 0.0;
            //double reaction_2_forward_rate = 0.0;
            //double reaction_2_reverse_rate = 0.0;
            //double reaction_3_rate = 0.0;

            AbstractTransportReaction* p_reaction_1 = new AbstractTransportReaction(p_bulk_1, p_cell_1, stoich_bulk_1, stoich_cell_1,reaction_1_rate);
            //AbstractReversibleTransportReaction* p_reaction_2 = new AbstractReversibleTransportReaction(p_bulk_2, p_cell_2, stoich_bulk_2, stoich_cell_2,reaction_2_forward_rate, reaction_2_reverse_rate);
            //AbstractTransportOutReaction* p_reaction_3 = new AbstractTransportOutReaction(p_bulk_3, p_cell_3, stoich_bulk_3, stoich_cell_3,reaction_3_rate);

            std::vector<AbstractTransportReaction*> p_transport_reaction_vector;
            p_transport_reaction_vector.push_back(p_reaction_1);
            //p_transport_reaction_vector.push_back(p_reaction_2);
            //p_transport_reaction_vector.push_back(p_reaction_3);

            AbstractTransportReactionSystem* p_transport_reaction_system = new AbstractTransportReactionSystem(p_bulk_system_chemistry, p_cell_system_chemistry, p_transport_reaction_vector);

            SetPtrNullTransportReactionSystem(p_transport_reaction_system);
        }

        void SetUpMembraneReactionSystem()
        {
            AbstractChemistry* p_bulk_system_chemistry = new AbstractChemistry();

            AbstractChemistry* p_cell_system_chemistry = new AbstractChemistry();

            std::vector<AbstractChemical*> p_bulk_substrates_1 = std::vector<AbstractChemical*>();
            std::vector<AbstractChemical*> p_cell_substrates_1 = std::vector<AbstractChemical*>();
            std::vector<AbstractChemical*> p_bulk_products_1 = std::vector<AbstractChemical*>();
            std::vector<AbstractChemical*> p_cell_products_1 = std::vector<AbstractChemical*>();


            std::vector<unsigned> stoich_bulk_substrates_1 = std::vector<unsigned>();
            std::vector<unsigned> stoich_cell_substrates_1 = std::vector<unsigned>();
            std::vector<unsigned> stoich_bulk_products_1 = std::vector<unsigned>();
            std::vector<unsigned> stoich_cell_products_1 = std::vector<unsigned>();


            AbstractChemical *p_chemical_U = new AbstractChemical("U");
            p_cell_system_chemistry -> AddChemical(p_chemical_U);
            // add U to reactions
            p_cell_substrates_1.push_back(p_chemical_U);
            p_cell_products_1.push_back(p_chemical_U);
            stoich_cell_substrates_1.push_back(2);
            stoich_cell_products_1.push_back(1);

            AbstractChemical *p_chemical_V = new AbstractChemical("V");
            p_cell_system_chemistry -> AddChemical(p_chemical_V);

            p_cell_substrates_1.push_back(p_chemical_V);
            stoich_cell_substrates_1.push_back(1);

            //p_bulk_system_chemistry -> AddChemical(p_chemical_U);
            p_bulk_system_chemistry -> AddChemical(p_chemical_V);
            // add U to reactions
            p_bulk_substrates_1.push_back(p_chemical_V);
            p_bulk_products_1.push_back(p_chemical_V);
            stoich_bulk_substrates_1.push_back(1);
            stoich_bulk_products_1.push_back(2);

            //AbstractChemical *p_chemical_C = new AbstractChemical("C");
            //p_bulk_system_chemistry -> AddChemical(p_chemical_C);

            //double reaction_1_forward_rate = 2.0;
            //double reaction_1_reverse_rate = 1.0;

            

            double reaction_1_forward_rate = membrane_reaction_rate_vector[0];
            double reaction_1_reverse_rate = membrane_reaction_rate_vector[1];
   
            MassActionCoupledMembraneReaction* p_reaction_1 = new MassActionCoupledMembraneReaction(p_bulk_substrates_1,p_bulk_products_1,p_cell_substrates_1,p_cell_products_1,stoich_bulk_substrates_1,stoich_bulk_products_1,stoich_cell_substrates_1,stoich_cell_products_1,false,true,reaction_1_forward_rate, reaction_1_reverse_rate);

            std::vector<AbstractMembraneReaction*> p_membrane_reaction_vector;

            p_membrane_reaction_vector.push_back(p_reaction_1);

            AbstractMembraneReactionSystem* p_membrane_reaction_system = new AbstractMembraneReactionSystem(p_bulk_system_chemistry, p_cell_system_chemistry, p_membrane_reaction_vector);

            SetPtrMembraneReactionSystem(p_membrane_reaction_system);
        }

        void SetUpNullMembraneReactionSystem()
        {
            AbstractChemistry* p_bulk_system_chemistry = new AbstractChemistry();

            AbstractChemistry* p_cell_system_chemistry = new AbstractChemistry();

            std::vector<AbstractChemical*> p_bulk_substrates_1 = std::vector<AbstractChemical*>();
            std::vector<AbstractChemical*> p_cell_substrates_1 = std::vector<AbstractChemical*>();
            std::vector<AbstractChemical*> p_bulk_products_1 = std::vector<AbstractChemical*>();
            std::vector<AbstractChemical*> p_cell_products_1 = std::vector<AbstractChemical*>();


            std::vector<unsigned> stoich_bulk_substrates_1 = std::vector<unsigned>();
            std::vector<unsigned> stoich_cell_substrates_1 = std::vector<unsigned>();
            std::vector<unsigned> stoich_bulk_products_1 = std::vector<unsigned>();
            std::vector<unsigned> stoich_cell_products_1 = std::vector<unsigned>();


            AbstractChemical *p_chemical_U = new AbstractChemical("U");
            p_cell_system_chemistry -> AddChemical(p_chemical_U);
            // add U to reactions
            p_cell_substrates_1.push_back(p_chemical_U);
            p_cell_products_1.push_back(p_chemical_U);
            stoich_cell_substrates_1.push_back(1);
            stoich_cell_products_1.push_back(1);


            AbstractChemical *p_chemical_V = new AbstractChemical("V");
            p_bulk_system_chemistry -> AddChemical(p_chemical_V);
            // add U to reactions
            p_cell_substrates_1.push_back(p_chemical_V);
            p_cell_products_1.push_back(p_chemical_V);
            stoich_bulk_substrates_1.push_back(1);
            stoich_bulk_products_1.push_back(1);

            // for the sake of testing multiple species in the bulk which do not take part in the transport process
            //AbstractChemical *p_chemical_C = new AbstractChemical("C");
            //p_bulk_system_chemistry -> AddChemical(p_chemical_C);

            double reaction_1_rate = 0.0;
  

            AbstractMembraneReaction* p_reaction_1 = new AbstractMembraneReaction(p_bulk_substrates_1,p_bulk_products_1,p_cell_substrates_1,p_cell_products_1,stoich_bulk_substrates_1,stoich_bulk_products_1,stoich_cell_substrates_1,stoich_cell_products_1,reaction_1_rate);
        
            std::vector<AbstractMembraneReaction*> p_membrane_reaction_vector;
            p_membrane_reaction_vector.push_back(p_reaction_1);
   
            AbstractMembraneReactionSystem* p_membrane_reaction_system = new AbstractMembraneReactionSystem(p_bulk_system_chemistry, p_cell_system_chemistry, p_membrane_reaction_vector);

            SetPtrNullMembraneReactionSystem(p_membrane_reaction_system);
        }

        void SetUpPdeChemicalOdeSystem()
        {
            SetUpPdeChemicalReactionSystem();

            AbstractChemicalOdeSystem* p_chemical_ode = new AbstractChemicalOdeSystem(mp_pde_chemical_reaction_system);

            SetPtrPdeChemicalOdeSystem(p_chemical_ode);
        }

        void SetUpTransportOdeSystem()
        {
            SetUpTransportReactionSystem();

            AbstractTransportOdeSystem* p_transport_ode = new AbstractTransportOdeSystem(mp_transport_system);
   
            SetPtrTransportOdeSystem(p_transport_ode);
        }

        void SetUpNullTransportOdeSystem()
        {
            SetUpNullTransportReactionSystem();

            AbstractTransportOdeSystem* p_transport_ode = new AbstractTransportOdeSystem(mp_null_transport_system);
   
            SetPtrNullTransportOdeSystem(p_transport_ode);
        }

        void SetUpMembraneOdeSystem()
        {
            SetUpMembraneReactionSystem();

            AbstractMembraneOdeSystem* p_membrane_ode = new AbstractMembraneOdeSystem(mp_membrane_system);
   
            SetPtrMembraneOdeSystem(p_membrane_ode);
        }

        void SetUpNullMembraneOdeSystem()
        {
            SetUpNullMembraneReactionSystem();

            AbstractMembraneOdeSystem* p_membrane_ode = new AbstractMembraneOdeSystem(mp_null_membrane_system);
   
            SetPtrNullMembraneOdeSystem(p_membrane_ode);
        }

        void SetUpChemicalDomainField()
        {
                   
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
            const unsigned probDim =2; // need to set manually to the number of diffusive variables for the pde solver to solve
            const unsigned spaceDim=2;
            const unsigned elementDim=2;
    
            //TetrahedralMesh<elementDim,spaceDim>* p_field = new TetrahedralMesh<elementDim,spaceDim>();
            // generate domain
            // run the domain field set up and parse files
            ChemicalDomainFieldTemplated<elementDim,spaceDim,probDim>* p_Pde_field = new ChemicalDomainFieldTemplated<elementDim,spaceDim,probDim>(dataFileRoot,dataFileRoot+domainFilename, dataFileRoot+domainKeyFilename, dataFileRoot+odeLabelFilename, dataFileRoot+odeKeyFilename, dataFileRoot+diffusionFilename);
            //std::cout<<"chemical structure number of mesh nodes: "<<p_Pde_field ->rGetDomainFeMesh()->GetNumNodes()<<std::endl;
            
            SetPtrChemicalDomain(p_Pde_field);
        }

        void SetUpChemicalDomainFieldForCellCoupling()
        {
            //std::cout<<"SetUpChemicalDomainFieldForCellCoupling()"<<std::endl; 
            // Variables for the user modify
            std::string dataFileRoot = "/home/chaste/projects/ChemChaste/src/Data/TemplateChemicalSimulation/";
            std::string cellLabelFilename = "";
            std::string cellKeyFilename = "";
            std::string domainFilename = "Domain.csv";
            std::string domainKeyFilename = "DomainKey.csv";
            std::string odeLabelFilename = "NodeSelector.csv";
            std::string odeKeyFilename = "OdeReactionFileKey.csv";
            std::string diffusionFilename = "DiffusionDatabaseFile.csv";
            std::string initialConditionsFilename = "InitialConditionFile.csv";
            std::string boundaryConditionsFilename = "BoundaryConditionFile.csv";
            
            // System properties
            const unsigned probDim =2; // need to set manually to the number of diffusive variables for the pde solver to solve
            const unsigned spaceDim=2;
            const unsigned elementDim=2;
    
            //TetrahedralMesh<elementDim,spaceDim>* p_field = new TetrahedralMesh<elementDim,spaceDim>();
            // generate domain
            // run the domain field set up and parse files
            ChemicalDomainFieldForCellCoupling<elementDim,spaceDim,probDim>* p_Pde_field = new ChemicalDomainFieldForCellCoupling<elementDim,spaceDim,probDim>(dataFileRoot,dataFileRoot+cellLabelFilename,dataFileRoot+cellKeyFilename,dataFileRoot+domainFilename, dataFileRoot+domainKeyFilename, dataFileRoot+odeLabelFilename, dataFileRoot+odeKeyFilename, dataFileRoot+diffusionFilename, dataFileRoot+initialConditionsFilename, dataFileRoot+boundaryConditionsFilename);

            SetPtrChemicalDomainFieldForCellCoupling(p_Pde_field);
            //std::cout<<"SetUpChemicalDomainFieldForCellCoupling() - end"<<std::endl; 
        }

        void SetUpSrnChemicalReactionSystem()
        {

            AbstractChemistry* p_cell_system_chemistry = new AbstractChemistry();

            std::vector<AbstractChemical*> p_substrates_1 = std::vector<AbstractChemical*>();
            std::vector<AbstractChemical*> p_products_1 = std::vector<AbstractChemical*>();

            std::vector<unsigned> stoich_substrates_1 = std::vector<unsigned>();
            std::vector<unsigned> stoich_products_1 = std::vector<unsigned>();
            
            AbstractChemical *p_chemical_U = new AbstractChemical("U");
            p_cell_system_chemistry -> AddChemical(p_chemical_U);
            // add U to reactions
            p_substrates_1.push_back(p_chemical_U);
            stoich_substrates_1.push_back(1);


            AbstractChemical *p_chemical_biomass = new AbstractChemical("Biomass");
            p_cell_system_chemistry -> AddChemical(p_chemical_biomass);
            // add V to reactions
            p_products_1.push_back(p_chemical_biomass);
            stoich_products_1.push_back(1);
  
            double reaction_1_forward_rate = srn_reaction_rate_vector[0];
            double reaction_1_reverse_rate = srn_reaction_rate_vector[1];

            MassActionReaction* p_mass_action_reaction_1 = new MassActionReaction(p_substrates_1, p_products_1, stoich_substrates_1, stoich_products_1,false,false,reaction_1_forward_rate,reaction_1_reverse_rate);

            std::vector<AbstractReaction*> p_mass_action_reaction_vector;
            p_mass_action_reaction_vector.push_back(p_mass_action_reaction_1);

            AbstractReactionSystem* p_mass_action_reaction_system = new AbstractReactionSystem(p_cell_system_chemistry, p_mass_action_reaction_vector);

            SetPtrSrnChemicalReactionSystem(p_mass_action_reaction_system);
        }

        void SetUpSrnChemicalOdeSystem()
        {

            SetUpSrnChemicalReactionSystem();

            AbstractChemicalOdeSystem* p_chemical_ode = new AbstractChemicalOdeSystem(mp_srn_chemical_reaction_system);

            SetPtrSrnChemicalOdeSystem(p_chemical_ode);
        }

        void SetUpChemicalSRNModel()
        {

            SetUpSrnChemicalOdeSystem();

            boost::shared_ptr<AbstractCellCycleModelOdeSolver> pCellOdeSolver = boost::shared_ptr<AbstractCellCycleModelOdeSolver>();

            ChemicalSrnModel* pSRNModel = new ChemicalSrnModel(mp_srn_chemical_reaction_system,pCellOdeSolver);

            pSRNModel -> Initialise();

            SetPtrChemicalSrnModel(pSRNModel);
        }

        void SetUpSimpleChemicalCellCycleModel()
        {

            SetUpSrnChemicalReactionSystem();

            AbstractChemistry* pCellChemistry = new AbstractChemistry();

            pCellChemistry -> AddChemistry(mp_srn_chemical_reaction_system -> GetSystemChemistry());
            AbstractChemical* pNewChemical = new AbstractChemical("V");
            pCellChemistry -> AddChemical(pNewChemical);

            SimpleChemicalThresholdCellCycleModel* pCellCycle = new SimpleChemicalThresholdCellCycleModel();
            pCellCycle->SetUp(pCellChemistry);
            // same species order as in the cell reaction system chemistry; mp_srn_chemical_reaction_system -> GetSystemChemistry()
        
            pCellCycle -> SetMaximumSpeciesThreshold(MaximumSpeciesThreshold);

            pCellCycle -> SetMinimumSpeciesThreshold(MinimumSpeciesThreshold);

            pCellCycle -> SetNumberThresholdSpecies(MaximumSpeciesThreshold.size());

            pCellCycle -> SetMaximumThresholdCheck(MaximumThresholdCheck);

            pCellCycle -> SetMinimumThresholdCheck(MinimumThresholdCheck);

            SetPtrSimpleChemicalThresholdCellCycleModel(pCellCycle);
        }

        void SetUpChemicalTrackingModifier()
        {

            boost::shared_ptr<ChemicalTrackingModifier<2,2>> p_chemical_tracking_modifier(new ChemicalTrackingModifier<2,2>());

            SetPtrChemicalTrackingModifier(p_chemical_tracking_modifier);
        }

        void SetUpCellInitialConditions(CellPtr p_cell, std::vector<std::string> speciesNames, std::vector<double> initValue)
        {
            for(unsigned i=0; i<speciesNames.size();i++)
            {
                p_cell->GetCellData()->SetItem(speciesNames[i],initValue[i]);
            }
        }

        CellPtr SetUpCellObjectA()
        {
            boost::shared_ptr<ChemicalCellProperty> p_cell_chemical(new ChemicalCellProperty());

            StateVariableRegister* p_state_register = new StateVariableRegister({"U","V"});

            p_cell_chemical -> InitialiseCell(p_state_register,cellObjectAConcentrationVector);

            boost::shared_ptr<MembraneCellProperty> p_cell_membrane(new MembraneCellProperty());
            
            p_cell_membrane -> SetMembraneThickness(1.0);

            boost::shared_ptr<TransportCellProperty> p_cell_transport(new TransportCellProperty());
            CellPropertyCollection collection;
            MAKE_PTR(WildTypeCellMutationState, p_state);
            
            collection.AddProperty(p_cell_chemical);
            collection.AddProperty(p_cell_membrane);
            collection.AddProperty(p_cell_transport);
            
            NoCellCycleModel* p_cell_cycle = new NoCellCycleModel();
            SchnackenbergSrnModel* p_srn_model = new SchnackenbergSrnModel();
            p_srn_model->SetInitialConditions(cellObjectAConcentrationVector);

            CellPtr p_cell(new Cell(p_state, p_cell_cycle, p_srn_model, false, collection));
            MAKE_PTR(StemCellProliferativeType, p_stem_type);
            p_cell->SetCellProliferativeType(p_stem_type);
            
            boost::shared_ptr<TransportCellProperty> p_transport = boost::static_pointer_cast<TransportCellProperty>(p_cell->rGetCellPropertyCollection().GetPropertiesType<TransportCellProperty>().GetProperty());

            SetUpTransportReactionSystem();

            AbstractTransportReactionSystem* p_transport_reaction_system = rGetPtrTransportReactionSystem();

            p_transport -> SetUp(p_transport_reaction_system,p_cell);

            boost::shared_ptr<MembraneCellProperty> p_membrane = boost::static_pointer_cast<MembraneCellProperty>(p_cell->rGetCellPropertyCollection().GetPropertiesType<MembraneCellProperty>().GetProperty());

            SetUpMembraneReactionSystem();

            AbstractMembraneReactionSystem* p_membrane_reaction_system = rGetPtrMembraneReactionSystem();

            p_membrane -> SetUp(p_membrane_reaction_system,p_cell);
            
            return p_cell;
        }

        CellPtr SetUpCellObjectB()
        {
            boost::shared_ptr<ChemicalCellProperty> p_cell_chemical(new ChemicalCellProperty());

            StateVariableRegister* p_state_register = new StateVariableRegister({"U","V"});

            p_cell_chemical -> InitialiseCell(p_state_register,cellObjectBConcentrationVector);

            boost::shared_ptr<MembraneCellProperty> p_cell_membrane(new MembraneCellProperty());
            
            p_cell_membrane -> SetMembraneThickness(1.0);

            boost::shared_ptr<TransportCellProperty> p_cell_transport(new TransportCellProperty());
            CellPropertyCollection collection;
            MAKE_PTR(WildTypeCellMutationState, p_state);
            collection.AddProperty(p_cell_chemical);
            collection.AddProperty(p_cell_membrane);
            collection.AddProperty(p_cell_transport);
            NoCellCycleModel* p_cell_cycle = new NoCellCycleModel();
            SchnackenbergSrnModel* p_srn_model = new SchnackenbergSrnModel();
            p_srn_model->SetInitialConditions(cellObjectBConcentrationVector);

            CellPtr p_cell(new Cell(p_state, p_cell_cycle, p_srn_model, false, collection));
            MAKE_PTR(StemCellProliferativeType, p_stem_type);
            p_cell->SetCellProliferativeType(p_stem_type);

            boost::shared_ptr<TransportCellProperty> p_transport = boost::static_pointer_cast<TransportCellProperty>(p_cell->rGetCellPropertyCollection().GetPropertiesType<TransportCellProperty>().GetProperty());

            SetUpTransportReactionSystem();

            AbstractTransportReactionSystem* p_transport_reaction_system = rGetPtrTransportReactionSystem();

            p_transport -> SetUp(p_transport_reaction_system,p_cell);

            boost::shared_ptr<MembraneCellProperty> p_membrane = boost::static_pointer_cast<MembraneCellProperty>(p_cell->rGetCellPropertyCollection().GetPropertiesType<MembraneCellProperty>().GetProperty());

            SetUpMembraneReactionSystem();

            AbstractMembraneReactionSystem* p_membrane_reaction_system = rGetPtrMembraneReactionSystem();

            p_membrane -> SetUp(p_membrane_reaction_system,p_cell);
            
            return p_cell;
        }

        CellPtr SetUpNullCellObject()
        {
            boost::shared_ptr<ChemicalCellProperty> p_cell_chemical(new ChemicalCellProperty());

            StateVariableRegister* p_state_register = new StateVariableRegister({"U","V"});

            p_cell_chemical -> InitialiseCell(p_state_register,nullCellObjectConcentrationVector);

            boost::shared_ptr<MembraneCellProperty> p_cell_membrane(new MembraneCellProperty());
            
            p_cell_membrane -> SetMembraneThickness(1.0);

            boost::shared_ptr<TransportCellProperty> p_cell_transport(new TransportCellProperty());
            CellPropertyCollection collection;
            MAKE_PTR(WildTypeCellMutationState, p_state);
            collection.AddProperty(p_cell_chemical);
            collection.AddProperty(p_cell_membrane);
            collection.AddProperty(p_cell_transport);
            NoCellCycleModel* p_cell_cycle = new NoCellCycleModel();
            CellPtr p_cell(new Cell(p_state, p_cell_cycle, NULL, false, collection));
            
            MAKE_PTR(StemCellProliferativeType, p_stem_type);
            p_cell->SetCellProliferativeType(p_stem_type);

            boost::shared_ptr<TransportCellProperty> p_transport = boost::static_pointer_cast<TransportCellProperty>(p_cell->rGetCellPropertyCollection().GetPropertiesType<TransportCellProperty>().GetProperty());
            boost::shared_ptr<MembraneCellProperty> p_membrane = boost::static_pointer_cast<MembraneCellProperty>(p_cell->rGetCellPropertyCollection().GetPropertiesType<MembraneCellProperty>().GetProperty());

            SetUpNullTransportReactionSystem();

            SetUpNullMembraneReactionSystem();

            AbstractTransportReactionSystem* p_transport_reaction_system = rGetPtrNullTransportReactionSystem();

            p_transport -> SetUp(p_transport_reaction_system,p_cell);

            AbstractMembraneReactionSystem* p_membrane_reaction_system = rGetPtrNullMembraneReactionSystem();

            p_membrane -> SetUp(p_membrane_reaction_system,p_cell);
            
            
            return p_cell;
        }

        CellPtr SetUpChemicalCellObjectA()
        {
            boost::shared_ptr<ChemicalCellProperty> p_cell_chemical(new ChemicalCellProperty());
            boost::shared_ptr<MembraneCellProperty> p_cell_membrane(new MembraneCellProperty());
            boost::shared_ptr<TransportCellProperty> p_cell_transport(new TransportCellProperty());

            // set up cell cycle model
            SetUpSimpleChemicalCellCycleModel();
            SimpleChemicalThresholdCellCycleModel* p_cell_cycle = rGetPtrSimpleChemicalThresholdCellCycleModel();

            // set up SRN model
            SetUpChemicalSRNModel();
            ChemicalSrnModel* p_srn_model = rGetPtrChemicalSrnModel();

            // set up cell properties of the Chemical Cell Property(not needed?)
            StateVariableRegister* p_state_register = new StateVariableRegister(chemicalCellASpeciesNames);
            //p_srn_model->SetInitialConditions(initialConcentrationVector);

            p_cell_chemical -> InitialiseCell(p_state_register,chemicalCellAConcentrationVector);

            p_cell_membrane -> SetMembraneThickness(1.0);

            // form cell
            CellPropertyCollection collection;
            MAKE_PTR(WildTypeCellMutationState, p_state);
            collection.AddProperty(p_cell_chemical);
            collection.AddProperty(p_cell_membrane);
            collection.AddProperty(p_cell_transport);
            
            CellPtr p_cell(new ChemicalCell(p_state, p_cell_cycle, p_srn_model, false, collection));

            SetUpCellInitialConditions(p_cell, chemicalCellASpeciesNames, chemicalCellAConcentrationVector);

            MAKE_PTR(StemCellProliferativeType, p_stem_type);
            p_cell->SetCellProliferativeType(p_stem_type);

            // set up the transport property
            boost::shared_ptr<TransportCellProperty> p_transport = boost::static_pointer_cast<TransportCellProperty>(p_cell->rGetCellPropertyCollection().GetPropertiesType<TransportCellProperty>().GetProperty());
            SetUpTransportReactionSystem();
            AbstractTransportReactionSystem* p_transport_reaction_system = rGetPtrTransportReactionSystem();
            p_transport -> SetUp(p_transport_reaction_system,p_cell);

            // set up the membrane property
            boost::shared_ptr<MembraneCellProperty> p_membrane = boost::static_pointer_cast<MembraneCellProperty>(p_cell->rGetCellPropertyCollection().GetPropertiesType<MembraneCellProperty>().GetProperty());
            SetUpMembraneReactionSystem();
            AbstractMembraneReactionSystem* p_membrane_reaction_system = rGetPtrMembraneReactionSystem();
            p_membrane -> SetUp(p_membrane_reaction_system,p_cell);

            
            return p_cell;
        }

        CellPtr SetUpChemicalCellObjectB()
        {
            boost::shared_ptr<ChemicalCellProperty> p_cell_chemical(new ChemicalCellProperty());
            boost::shared_ptr<MembraneCellProperty> p_cell_membrane(new MembraneCellProperty());
            boost::shared_ptr<TransportCellProperty> p_cell_transport(new TransportCellProperty());

            // set up cell cycle model
            SetUpSimpleChemicalCellCycleModel();
            SimpleChemicalThresholdCellCycleModel* p_cell_cycle = rGetPtrSimpleChemicalThresholdCellCycleModel();

            // set up SRN model
            NullSrnModel* p_srn_model = new NullSrnModel();

            StateVariableRegister* p_state_register = new StateVariableRegister(chemicalCellBSpeciesNames);

            p_cell_chemical -> InitialiseCell(p_state_register,chemicalCellBConcentrationVector);
            
            p_cell_membrane -> SetMembraneThickness(1.0);

        
            CellPropertyCollection collection;
            MAKE_PTR(WildTypeCellMutationState, p_state);
            collection.AddProperty(p_cell_chemical);
            collection.AddProperty(p_cell_membrane);
            collection.AddProperty(p_cell_transport);

            CellPtr p_cell(new ChemicalCell(p_state, p_cell_cycle, p_srn_model, false, collection));
            
            MAKE_PTR(StemCellProliferativeType, p_stem_type);
            p_cell->SetCellProliferativeType(p_stem_type);

            // set up the transport property
            boost::shared_ptr<TransportCellProperty> p_transport = boost::static_pointer_cast<TransportCellProperty>(p_cell->rGetCellPropertyCollection().GetPropertiesType<TransportCellProperty>().GetProperty());
            SetUpTransportReactionSystem();
            AbstractTransportReactionSystem* p_transport_reaction_system = rGetPtrTransportReactionSystem();
            p_transport -> SetUp(p_transport_reaction_system,p_cell);

            // set up the membrane property
            boost::shared_ptr<MembraneCellProperty> p_membrane = boost::static_pointer_cast<MembraneCellProperty>(p_cell->rGetCellPropertyCollection().GetPropertiesType<MembraneCellProperty>().GetProperty());
            SetUpMembraneReactionSystem();
            AbstractMembraneReactionSystem* p_membrane_reaction_system = rGetPtrMembraneReactionSystem();
            p_membrane -> SetUp(p_membrane_reaction_system,p_cell);
            
            return p_cell;
        }

        CellPtr SetUpNullChemicalCellObject()
        {
            boost::shared_ptr<ChemicalCellProperty> p_cell_chemical(new ChemicalCellProperty());

            StateVariableRegister* p_state_register = new StateVariableRegister(nullChemicalCellSpeciesNames);
        
            p_cell_chemical -> InitialiseCell(p_state_register,nullChemicalCellConcentrationVector);

            boost::shared_ptr<MembraneCellProperty> p_cell_membrane(new MembraneCellProperty());
            
            p_cell_membrane -> SetMembraneThickness(1.0);

            boost::shared_ptr<TransportCellProperty> p_cell_transport(new TransportCellProperty());
            CellPropertyCollection collection;
            MAKE_PTR(WildTypeCellMutationState, p_state);
            collection.AddProperty(p_cell_chemical);
            collection.AddProperty(p_cell_membrane);
            collection.AddProperty(p_cell_transport);
            NoCellCycleModel* p_cell_cycle = new NoCellCycleModel();
            CellPtr p_cell(new ChemicalCell(p_state, p_cell_cycle, NULL, false, collection));
            
            MAKE_PTR(StemCellProliferativeType, p_stem_type);
            p_cell->SetCellProliferativeType(p_stem_type);

            boost::shared_ptr<TransportCellProperty> p_transport = boost::static_pointer_cast<TransportCellProperty>(p_cell->rGetCellPropertyCollection().GetPropertiesType<TransportCellProperty>().GetProperty());
            boost::shared_ptr<MembraneCellProperty> p_membrane = boost::static_pointer_cast<MembraneCellProperty>(p_cell->rGetCellPropertyCollection().GetPropertiesType<MembraneCellProperty>().GetProperty());

            SetUpNullTransportReactionSystem();

            SetUpNullMembraneReactionSystem();

            AbstractTransportReactionSystem* p_transport_reaction_system = rGetPtrNullTransportReactionSystem();

            p_transport -> SetUp(p_transport_reaction_system,p_cell);

            AbstractMembraneReactionSystem* p_membrane_reaction_system = rGetPtrNullMembraneReactionSystem();

            p_membrane -> SetUp(p_membrane_reaction_system,p_cell);
            
            
            return p_cell;
        }


        void SetPtrPdeChemicalReactionSystem(AbstractReactionSystem*  p_reaction_system)
        {
            mp_pde_chemical_reaction_system = p_reaction_system;
        }

        void SetPtrTransportReactionSystem(AbstractTransportReactionSystem*  p_reaction_system)
        {
            mp_transport_system = p_reaction_system;
        }

        void SetPtrNullTransportReactionSystem(AbstractTransportReactionSystem*  p_reaction_system)
        {
            mp_null_transport_system = p_reaction_system;
        }

        void SetPtrMembraneReactionSystem(AbstractMembraneReactionSystem*  p_reaction_system)
        {
            mp_membrane_system = p_reaction_system;
        }

        void SetPtrNullMembraneReactionSystem(AbstractMembraneReactionSystem*  p_reaction_system)
        {
            mp_null_membrane_system = p_reaction_system;
        }


        void SetPtrPdeChemicalOdeSystem(AbstractChemicalOdeSystem*  p_chemical_ode)
        {
            mp_pde_chemical_ode = p_chemical_ode;
        }

        void SetPtrTransportOdeSystem(AbstractTransportOdeSystem*  p_transport_ode)
        {
            mp_transport_ode = p_transport_ode;
        }

        void SetPtrNullTransportOdeSystem(AbstractTransportOdeSystem*  p_transport_ode)
        {
            mp_null_transport_ode = p_transport_ode;
        }

        void SetPtrMembraneOdeSystem(AbstractMembraneOdeSystem*  p_membrane_ode)
        {
            mp_membrane_ode = p_membrane_ode;
        }

        void SetPtrNullMembraneOdeSystem(AbstractMembraneOdeSystem*  p_membrane_ode)
        {
            mp_null_membrane_ode = p_membrane_ode;
        }

        void SetPtrChemicalDomain(ChemicalDomainFieldTemplated<2,2,2>*  p_Pde_field)
        {
            mpChemicalDomain = p_Pde_field;
        }

        void SetPtrChemicalDomainFieldForCellCoupling(ChemicalDomainFieldForCellCoupling<2,2,2>*  p_field)
        {
            mpChemicalDomainForCellCoupling = p_field;
        }

        void SetPtrSrnChemicalReactionSystem(AbstractReactionSystem*  p_srn_chemical_reaction_system)
        {
            mp_srn_chemical_reaction_system = p_srn_chemical_reaction_system;
        }
            
        void SetPtrSrnChemicalOdeSystem(AbstractChemicalOdeSystem*  p_srn_chemical_ode)
        {
            mp_srn_chemical_ode = p_srn_chemical_ode;
        }
        
        void SetPtrChemicalSrnModel(ChemicalSrnModel*  p_srnModel)
        {
            mpChemicalSrnModel = p_srnModel;
        }

        void SetThresholdNames(std::vector<std::string>  thresholdNames)
        {   
            mThresholdNames = thresholdNames;
        }

        void SetMaxThresholdVector(std::vector<double>  maxVector)
        {
            mMaxThresholdVector = maxVector;
        }

        void SetMinThresholdVector(std::vector<double>  minVector)
        {
            mMinThresholdVector = minVector;
        }

        void SetPtrSimpleChemicalThresholdCellCycleModel(SimpleChemicalThresholdCellCycleModel*  p_cell_cycle)
        {
            mpSimpleChemicalThresholdCellCycleModel = p_cell_cycle;
        }

        void SetPtrChemicalTrackingModifier(boost::shared_ptr<ChemicalTrackingModifier<2,2>>  p_modifier)
        {
            mpChemicalTrackingModifier = p_modifier;
        }


        AbstractReactionSystem*& rGetPtrPdeChemicalReactionSystem()
        {
            return mp_pde_chemical_reaction_system;
        }

        AbstractTransportReactionSystem*& rGetPtrTransportReactionSystem()
        {
            return mp_transport_system;
        }

        AbstractTransportReactionSystem*& rGetPtrNullTransportReactionSystem()
        {
            return mp_null_transport_system;
        }

        AbstractMembraneReactionSystem*& rGetPtrMembraneReactionSystem()
        {
            return mp_membrane_system;
        }

        AbstractMembraneReactionSystem*& rGetPtrNullMembraneReactionSystem()
        {
            return mp_null_membrane_system;
        }


        AbstractChemicalOdeSystem*& rGetPtrPdeChemicalOdeSystem()
        {
            return mp_pde_chemical_ode;
        }

        AbstractTransportOdeSystem*& rGetPtrTransportOdeSystem()
        {
            return mp_transport_ode;
        }

        AbstractTransportOdeSystem*& rGetPtrNullTransportOdeSystem()
        {
            return mp_null_transport_ode;
        }

        AbstractMembraneOdeSystem*& rGetPtrMembraneOdeSystem()
        {
            return mp_membrane_ode;
        }

        AbstractMembraneOdeSystem*& rGetPtrNullMembraneOdeSystem()
        {
            return mp_null_membrane_ode;
        }

        ChemicalDomainFieldTemplated<2,2,2>*& rGetPtrChemicalDomain()
        {
            return mpChemicalDomain;
        }

        ChemicalDomainFieldForCellCoupling<2,2,2>*& rGetPtrChemicalDomainFieldForCellCoupling()
        {
            return mpChemicalDomainForCellCoupling;
        }
        
        AbstractReactionSystem*& rGetPtrSrnChemicalReactionSystem()
        {
            return mp_srn_chemical_reaction_system;
        }
            
        AbstractChemicalOdeSystem*& rGetPtrSrnChemicalOdeSystem()
        {
            return mp_srn_chemical_ode;
        }
        
        ChemicalSrnModel*& rGetPtrChemicalSrnModel()
        {
            return mpChemicalSrnModel;
        }

        std::vector<std::string> GetThresholdNames()
        {   
            return mThresholdNames;
        }

        std::vector<double> GetMaxThresholdVector()
        {
            return mMaxThresholdVector;
        }

        std::vector<double> GetMinThresholdVector()
        {
            return mMinThresholdVector;
        }

        std::vector<std::string> GetChemicalCellASpeciesNames()
        {
            return chemicalCellASpeciesNames;
        }


        SimpleChemicalThresholdCellCycleModel*& rGetPtrSimpleChemicalThresholdCellCycleModel()
        {
            return mpSimpleChemicalThresholdCellCycleModel;
        }

        boost::shared_ptr<ChemicalTrackingModifier<2,2>>& rGetPtrChemicalTrackingModifier()
        {
            return mpChemicalTrackingModifier;
        }
};

#endif