#ifndef TESTCHEMICALCELLSYSTEM_HPP_
#define TESTCHEMICALCELLSYSTEM_HPP_

// chaste includes
#include <cxxtest/TestSuite.h>
#include "UblasIncludes.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "RandomNumberGenerator.hpp"
#include "SmartPointers.hpp"

// general includes
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <tuple>
#include <cmath>



// modified tumour spheroid simulation includes
#include "HoneycombMeshGenerator.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "SimpleOxygenBasedCellCycleModel.hpp"
#include "WildTypeCellMutationState.hpp"
#include "StemCellProliferativeType.hpp"
#include "CellwiseSourceEllipticPde.hpp"
#include "ConstBoundaryCondition.hpp"
#include "EllipticGrowingDomainPdeModifier.hpp"
#include "EllipticBoxDomainPdeModifier.hpp" 
#include "OffLatticeSimulation.hpp"

#include "ApoptoticCellProperty.hpp"
#include "ArchiveOpener.hpp"
//#include "AveragedSourceParabolicPde.hpp"
#include "AveragedSourceParabolicPde_test.hpp"
#include "CaBasedCellPopulation.hpp"
#include "CellsGenerator.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "ParabolicBoxDomainPdeModifier.hpp"
#include "ReplicatableVector.hpp"
#include "UniformCellCycleModel.hpp"
#include "UniformSourceParabolicPde.hpp"

#include "NoCellCycleModel.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"

// added classes
//#include "AbstractDomainField_templated.hpp"
//#include "ParabolicBoxDomainPdeSystemModifier.hpp"
//#include "ChemicalDomainField_templated.hpp"
#include "ParabolicBoxDomainPdeSystemModifier.hpp"
#include "ChemicalDomainFieldForCellCoupling.hpp"

#include "ChemicalDomainField_templated.hpp"
//#include "ChemicalDomainFieldForCellCoupling.hpp"

//#include "InhomogenousCoupledPdeOdeSolver_templated.hpp"
#include "InhomogenousParabolicPdeForCoupledOdeSystem_templated.hpp"


// cell property includes
#include "TransportCellProperty.hpp"
#include "ChemicalCellProperty.hpp"
#include "ExtendedCellProperty.hpp"
#include "MembraneCellProperty.hpp"

#include "StemCellProliferativeType.hpp"

// chemical ode includes
#include "AbstractReactionSystem.hpp"
#include "MassActionReaction.hpp"
#include "AbstractChemicalOdeSystem.hpp"
#include "AbstractChemicalOdeForCoupledPdeSystem.hpp"

// transport includes
#include "AbstractTransportReactionSystem.hpp"
#include "AbstractTransportReaction.hpp"
#include "AbstractTransportOutReaction.hpp"
#include "MassActionTransportReaction.hpp"
//#include "AbstractTransportReactionSystemFromFile.hpp"
#include "AbstractReversibleTransportReaction.hpp"
#include "AbstractTransportOdeSystem.hpp"

// membrane includes
#include "AbstractMembraneReactionSystem.hpp"
#include "AbstractMembraneReaction.hpp"
#include "MassActionCoupledMembraneReaction.hpp"
//#include "AbstractMembraneReactionSystemFromFile.hpp"
#include "AbstractReversibleMembraneReaction.hpp"
#include "AbstractMembraneOdeSystem.hpp"


#include "EulerIvpOdeSolver.hpp"

// SRN includes
#include "SchnackenbergSrnModel.hpp"
#include "ChemicalSrnModel.hpp"
#include "NullSrnModel.hpp"

// cell cycle includes
#include "SimpleChemicalThresholdCellCycleModel.hpp"

// tracking includes
#include "ChemicalTrackingModifier.hpp"

#include "ChemicalCell.hpp"

#include "AbstractBoxDomainPdeSystemModifier.hpp"
#include "AbstractPdeSystemModifier.hpp"

//#include "AbstractPdeSystemModifier.hpp"

class ChemicalStructuresForTests
{
    private:
        AbstractReactionSystem* mp_chemical_reaction_system;

        AbstractChemicalOdeSystem* mp_chemical_ode;

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



        AbstractReactionSystem* mp_cell_chemical_reaction_system;

        AbstractChemicalOdeSystem* mp_cell_chemical_ode;

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

        void SetUpChemicalReactionSystem()
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

            //double reaction_1_rate = 0.1;
            //double reaction_2_forward_rate = 0.1;
            //double reaction_2_reverse_rate = 0.2;
            //double reaction_3_rate = 0.3;

            double reaction_1_rate = 0.0;
            double reaction_2_forward_rate = 0.0;
            double reaction_2_reverse_rate = 0.0;
            double reaction_3_rate = 0.0;


            MassActionReaction* p_mass_action_reaction_1 = new MassActionReaction(p_substrates_1, p_products_1, stoich_substrates_1, stoich_products_1,false,false,reaction_1_rate);
            MassActionReaction* p_mass_action_reaction_2 = new MassActionReaction(p_substrates_2, p_products_2, stoich_substrates_2, stoich_products_2,false, true, reaction_2_forward_rate, reaction_2_reverse_rate);
            MassActionReaction* p_mass_action_reaction_3 = new MassActionReaction(p_substrates_3, p_products_3, stoich_substrates_3, stoich_products_3,false, false,reaction_3_rate);

            std::vector<AbstractReaction*> p_mass_action_reaction_vector;
            p_mass_action_reaction_vector.push_back(p_mass_action_reaction_1);
            p_mass_action_reaction_vector.push_back(p_mass_action_reaction_2);
            p_mass_action_reaction_vector.push_back(p_mass_action_reaction_3);

            AbstractReactionSystem* p_mass_action_reaction_system = new AbstractReactionSystem(p_system_chemistry, p_mass_action_reaction_vector);
            std::cout<<"NumberOfReactions = "<< p_mass_action_reaction_system -> GetNumberOfReactions()<<std::endl;
            SetPtrChemicalReactionSystem(p_mass_action_reaction_system);
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

            double reaction_2_forward_rate = 0.3;
            double reaction_2_reverse_rate = 0.1;
            double reaction_3_forward_rate = 0.1;
            double reaction_3_reverse_rate = 0.2;


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


            p_bulk_system_chemistry -> AddChemical(p_chemical_U);
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

            double reaction_1_forward_rate = 0.0;
            double reaction_1_reverse_rate = 0.0;
   
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

        void SetUpChemicalOdeSystem()
        {
            SetUpChemicalReactionSystem();

            AbstractChemicalOdeSystem* p_chemical_ode = new AbstractChemicalOdeSystem(mp_chemical_reaction_system);

            SetPtrChemicalOdeSystem(p_chemical_ode);
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
            std::cout<<"SetUpChemicalDomainFieldForCellCoupling()"<<std::endl; 
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
            std::cout<<"chemical structure number of mesh nodes: "<<p_Pde_field ->rGetDomainFeMesh()->GetNumNodes()<<std::endl;
            
            SetPtrChemicalDomainFieldForCellCoupling(p_Pde_field);
            std::cout<<"SetUpChemicalDomainFieldForCellCoupling() - end"<<std::endl; 
        }

        void SetUpCellChemicalReactionSystem()
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
  
            double reaction_1_rate = 1.0;

            MassActionReaction* p_mass_action_reaction_1 = new MassActionReaction(p_substrates_1, p_products_1, stoich_substrates_1, stoich_products_1,false,false,reaction_1_rate);

            std::vector<AbstractReaction*> p_mass_action_reaction_vector;
            p_mass_action_reaction_vector.push_back(p_mass_action_reaction_1);

            AbstractReactionSystem* p_mass_action_reaction_system = new AbstractReactionSystem(p_cell_system_chemistry, p_mass_action_reaction_vector);

            SetPtrCellChemicalReactionSystem(p_mass_action_reaction_system);
        }

        void SetUpCellChemicalOdeSystem()
        {

            SetUpCellChemicalReactionSystem();

            AbstractChemicalOdeSystem* p_chemical_ode = new AbstractChemicalOdeSystem(mp_cell_chemical_reaction_system);

            SetPtrCellChemicalOdeSystem(p_chemical_ode);
        }

        void SetUpChemicalSRNModel()
        {

            SetUpCellChemicalOdeSystem();

            boost::shared_ptr<AbstractCellCycleModelOdeSolver> pCellOdeSolver = boost::shared_ptr<AbstractCellCycleModelOdeSolver>();

            ChemicalSrnModel* pSRNModel = new ChemicalSrnModel(mp_cell_chemical_reaction_system,pCellOdeSolver);

            pSRNModel -> Initialise();

            SetPtrChemicalSrnModel(pSRNModel);
        }

        void SetUpSimpleChemicalCellCycleModel()
        {

            SetUpCellChemicalReactionSystem();

            AbstractChemistry* pCellChemistry = new AbstractChemistry();

            pCellChemistry -> AddChemistry(mp_cell_chemical_reaction_system -> GetSystemChemistry());
            AbstractChemical* pNewChemical = new AbstractChemical("V");
            pCellChemistry -> AddChemical(pNewChemical);

            SimpleChemicalThresholdCellCycleModel* pCellCycle = new SimpleChemicalThresholdCellCycleModel();
            pCellCycle->SetUp(pCellChemistry);
            // same species order as in the cell reaction system chemistry; mp_cell_chemical_reaction_system -> GetSystemChemistry()
            // {U,Biomass}
            std::vector<double> MaximumSpeciesThreshold = {0,2.0,0};

            std::vector<double> MinimumSpeciesThreshold = {0,0.1,0};

            std::vector<bool> MaximumThresholdCheck = {false,true,false};

            std::vector<bool> MinimumThresholdCheck = {false,true,false};

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
                std::cout<<speciesNames[i]<<std::endl;
                p_cell->GetCellData()->SetItem(speciesNames[i],initValue[i]);
            }
        }

        CellPtr SetUpCellObjectA()
        {
            boost::shared_ptr<ChemicalCellProperty> p_cell_chemical(new ChemicalCellProperty());

            StateVariableRegister* p_state_register = new StateVariableRegister({"U","V"});
            std::vector<double> concentrationVector = {1.0,0.0};

            p_cell_chemical -> InitialiseCell(p_state_register,concentrationVector);

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
            p_srn_model->SetInitialConditions(concentrationVector);

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
            std::vector<double> concentrationVector = {0.0,1.0};

            p_cell_chemical -> InitialiseCell(p_state_register,concentrationVector);

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
            p_srn_model->SetInitialConditions(concentrationVector);

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
            std::vector<double> concentrationVector = {0.0,0.0};

            p_cell_chemical -> InitialiseCell(p_state_register,concentrationVector);

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
            
            std::vector<std::string> speciesNames = {"U","Biomass","V"};
            StateVariableRegister* p_state_register = new StateVariableRegister(speciesNames);
            std::vector<double> initialConcentrationVector = {5.0,1.9,1.0}; // initial conditions
            
            //p_srn_model->SetInitialConditions(initialConcentrationVector);

            p_cell_chemical -> InitialiseCell(p_state_register,initialConcentrationVector);

            p_cell_membrane -> SetMembraneThickness(1.0);

            // form cell
            CellPropertyCollection collection;
            MAKE_PTR(WildTypeCellMutationState, p_state);
            collection.AddProperty(p_cell_chemical);
            collection.AddProperty(p_cell_membrane);
            collection.AddProperty(p_cell_transport);
            
            CellPtr p_cell(new ChemicalCell(p_state, p_cell_cycle, p_srn_model, false, collection));

            SetUpCellInitialConditions(p_cell, speciesNames, initialConcentrationVector);

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

            StateVariableRegister* p_state_register = new StateVariableRegister({"U","V"});
            std::vector<double> initialConcentrationVector = {0.0,1.0};

            p_cell_chemical -> InitialiseCell(p_state_register,initialConcentrationVector);
            
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

            StateVariableRegister* p_state_register = new StateVariableRegister({"U","V"});
            std::vector<double> concentrationVector = {0.0,0.0};

            p_cell_chemical -> InitialiseCell(p_state_register,concentrationVector);

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


        void SetPtrChemicalReactionSystem(AbstractReactionSystem*  p_reaction_system)
        {
            mp_chemical_reaction_system = p_reaction_system;
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


        void SetPtrChemicalOdeSystem(AbstractChemicalOdeSystem*  p_chemical_ode)
        {
            mp_chemical_ode = p_chemical_ode;
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

        void SetPtrCellChemicalReactionSystem(AbstractReactionSystem*  p_cell_chemical_reaction_system)
        {
            mp_cell_chemical_reaction_system = p_cell_chemical_reaction_system;
        }
            
        void SetPtrCellChemicalOdeSystem(AbstractChemicalOdeSystem*  p_cell_chemical_ode)
        {
            mp_cell_chemical_ode = p_cell_chemical_ode;
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


        AbstractReactionSystem*& rGetPtrChemicalReactionSystem()
        {
            return mp_chemical_reaction_system;
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


        AbstractChemicalOdeSystem*& rGetPtrChemicalOdeSystem()
        {
            return mp_chemical_ode;
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
        
        AbstractReactionSystem*& rGetPtrCellChemicalReactionSystem()
        {
            return mp_cell_chemical_reaction_system;
        }
            
        AbstractChemicalOdeSystem*& rGetPtrCellChemicalOdeSystem()
        {
            return mp_cell_chemical_ode;
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

        SimpleChemicalThresholdCellCycleModel*& rGetPtrSimpleChemicalThresholdCellCycleModel()
        {
            return mpSimpleChemicalThresholdCellCycleModel;
        }

        boost::shared_ptr<ChemicalTrackingModifier<2,2>>& rGetPtrChemicalTrackingModifier()
        {
            return mpChemicalTrackingModifier;
        }

};



class TestCellwisePdeSystem : public AbstractCellBasedTestSuite
{
public:

    void TestChemicalStructures_chemical_ODE()
    {
        
        ChemicalStructuresForTests* chemical_structure = new ChemicalStructuresForTests();
        chemical_structure->SetUpChemicalOdeSystem();
        AbstractChemicalOdeSystem* chemicalOde = chemical_structure -> rGetPtrChemicalOdeSystem();
        
        std::cout<<"Number of species: "<< chemicalOde -> GetNumberOfSpecies()<<std::endl;

        std::cout<<"Number of reactions: "<< chemicalOde -> GetNumberOfReactions()<<std::endl;
        EulerIvpOdeSolver euler_solver;
        std::vector<double> initial_condition(chemicalOde -> GetNumberOfSpecies(), 1.0);

        OdeSolution solutions = euler_solver.Solve(chemicalOde, initial_condition, 0, 1, 0.01, 0.1);

        AbstractChemistry* systemChemistry = chemicalOde->GetReactionSystem()->GetSystemChemistry();

        std::cout<<"Time: ";
        for(unsigned i=0; i<systemChemistry->GetNumberChemicals(); i++)
        {
            std::cout<<systemChemistry->GetChemicalNamesByIndex(i)<<" ";
        }
        std::cout<<std::endl;

        for (unsigned i=0; i<solutions.rGetTimes().size(); i++)
        {
            std::cout << solutions.rGetTimes()[i] << " " << solutions.rGetSolutions()[i][0] << " " << solutions.rGetSolutions()[i][1]<< "\n";
        }
        
    }

    void TestChemicalStructures_transport_reaction()
    {
        
        ChemicalStructuresForTests* chemical_structure = new ChemicalStructuresForTests();

        chemical_structure -> SetUpTransportReactionSystem();

        AbstractTransportReactionSystem* p_transport_reaction_system = chemical_structure -> rGetPtrTransportReactionSystem();

        
        // test transport reaction

        std::vector<std::string> environment_states = {"U", "V","C"};
        std::vector<double> environment_concentration_at_point = {1.0,1.0,1.0};
        std::vector<double> change_environment_concentration_vector= {0.0,0.0,0.0};

        // states bound within the cell
        std::vector<std::string> cell_states = {"U", "V"};
        std::vector<double> cell_concentration = {1.0,1.0};
        std::vector<double> change_cell_concentration_vector= {0.0,0.0};

        p_transport_reaction_system -> ReactSystem(environment_concentration_at_point,cell_concentration,change_environment_concentration_vector,change_cell_concentration_vector);
        
        for(unsigned i=0; i<environment_concentration_at_point.size();i++)
        {
            environment_concentration_at_point[i] = environment_concentration_at_point[i] + change_environment_concentration_vector[i];
        }
        for(unsigned i=0; i<cell_concentration.size();i++)
        {
            cell_concentration[i] = cell_concentration[i] + change_cell_concentration_vector[i];
        }
        for(unsigned i=0; i<p_transport_reaction_system -> GetNumberOfReactions(); i++)
        {
            std::cout<<"Reaction type: "<<i<<" "<< p_transport_reaction_system -> GetReactionByIndex(i) -> GetReactionType()<<std::endl;
        }

        for(unsigned i=0; i<environment_concentration_at_point.size();i++)
        {
            std::cout<<"Bulk State "<<environment_states[i]<<": initial: "<<environment_concentration_at_point[i]<<" change: "<<change_environment_concentration_vector[i]<<std::endl;
        }

        for(unsigned i=0; i<cell_concentration.size();i++)
        {
            std::cout<<"Cell State "<<cell_states[i]<<": initial: "<<cell_concentration[i]<<" change: "<<change_cell_concentration_vector[i]<<std::endl;
        }
        
    }

    void TestChemicalStructures_transport_ODE()
    {
        
        ChemicalStructuresForTests* chemical_structure = new ChemicalStructuresForTests();
        chemical_structure->SetUpTransportOdeSystem();

        AbstractTransportOdeSystem* transportOde = chemical_structure -> rGetPtrTransportOdeSystem();
        
        std::cout<<"Number of species: "<< transportOde -> GetNumberOfSpecies()<<std::endl;

        std::cout<<"Number of reactions: "<< transportOde -> GetNumberOfReactions()<<std::endl;

        std::cout<<"Number of equations: "<<transportOde->GetNumberOfStateVariables()<<std::endl;
        EulerIvpOdeSolver euler_solver;
        std::vector<double> initial_condition(transportOde -> GetNumberOfSpecies(), 1.0);

        OdeSolution solutions = euler_solver.Solve(transportOde, initial_condition, 0, 0.2, 0.1, 0.1);


        AbstractChemistry* bulkChemistry = transportOde->GetReactionSystem()->GetBulkChemistry();
        AbstractChemistry* cellChemistry = transportOde->GetReactionSystem()->GetCellChemistry();
        std::cout<<"Time: ";
        for(unsigned i=0; i<bulkChemistry->GetNumberChemicals(); i++)
        {
            std::cout<<"Bulk: "<< bulkChemistry->GetChemicalNamesByIndex(i)<<" ";
        }
        for(unsigned i=0; i<cellChemistry->GetNumberChemicals(); i++)
        {
            std::cout<<"Cell: "<< cellChemistry->GetChemicalNamesByIndex(i)<<" ";
        }
        std::cout<<std::endl;
        
        for (unsigned i=0; i<solutions.rGetTimes().size(); i++)
        {
            std::cout << solutions.rGetTimes()[i] << " , " << solutions.rGetSolutions()[i][0] << " , " << solutions.rGetSolutions()[i][1]<< " , " << solutions.rGetSolutions()[i][2] << " , " << solutions.rGetSolutions()[i][3]<<" , " << solutions.rGetSolutions()[i][4]<<"\n";
        }
        
    }

    void TestChemicalStructures_membrane_reaction()
    {
        /*
        ChemicalStructuresForTests* chemical_structure = new ChemicalStructuresForTests();

        chemical_structure -> SetUpMembraneReactionSystem();

        AbstractMembraneReactionSystem* p_membrane_reaction_system = chemical_structure -> rGetPtrMembraneReactionSystem();

        
        // test membrane reaction

        std::vector<std::string> environment_states = {"U", "V","C"};
        std::vector<double> environment_concentration_at_point = {1.0,1.0,1.0};
        std::vector<double> change_environment_concentration_vector= {0.0,0.0,0.0};

        // states bound within the cell
        std::vector<std::string> cell_states = {"U", "V"};
        std::vector<double> cell_concentration = {1.0,1.0};
        std::vector<double> change_cell_concentration_vector= {0.0,0.0};

        p_membrane_reaction_system -> ReactSystem(environment_concentration_at_point,cell_concentration,change_environment_concentration_vector,change_cell_concentration_vector);
        
        for(unsigned i=0; i<environment_concentration_at_point.size();i++)
        {
            //environment_concentration_at_point[i] = environment_concentration_at_point[i] + change_environment_concentration_vector[i];
        }
        for(unsigned i=0; i<cell_concentration.size();i++)
        {
          //  cell_concentration[i] = cell_concentration[i] + change_cell_concentration_vector[i];
        }
        for(unsigned i=0; i<p_membrane_reaction_system -> GetNumberOfReactions(); i++)
        {
            std::cout<<"Reaction type: "<<i<<" "<< p_membrane_reaction_system -> GetReactionByIndex(i) -> GetReactionType()<<std::endl;
        }

        for(unsigned i=0; i<environment_concentration_at_point.size();i++)
        {
            std::cout<<"Bulk State "<<environment_states[i]<<": initial: "<<environment_concentration_at_point[i]<<" change: "<<change_environment_concentration_vector[i]<<std::endl;
        }

        for(unsigned i=0; i<cell_concentration.size();i++)
        {
            std::cout<<"Cell State "<<cell_states[i]<<": initial: "<<cell_concentration[i]<<" change: "<<change_cell_concentration_vector[i]<<std::endl;
        }
        */
    }

    void TestChemicalStructures_membrane_ODE()
    {
        /*
        ChemicalStructuresForTests* chemical_structure = new ChemicalStructuresForTests();
        chemical_structure->SetUpMembraneOdeSystem();

        AbstractMembraneOdeSystem* membraneOde = chemical_structure -> rGetPtrMembraneOdeSystem();
        
        std::cout<<"Number of species: "<< membraneOde -> GetNumberOfSpecies()<<std::endl;

        std::cout<<"Number of reactions: "<< membraneOde -> GetNumberOfReactions()<<std::endl;

        std::cout<<"Number of equations: "<<membraneOde->GetNumberOfStateVariables()<<std::endl;
        EulerIvpOdeSolver euler_solver;
        std::vector<double> initial_condition(membraneOde -> GetNumberOfSpecies(), 1.0);

        OdeSolution solutions = euler_solver.Solve(membraneOde, initial_condition, 0, 1, 0.01, 0.1);

        AbstractChemistry* bulkChemistry = membraneOde->GetReactionSystem()->GetBulkChemistry();
        AbstractChemistry* cellChemistry = membraneOde->GetReactionSystem()->GetCellChemistry();
        std::cout<<"Time: ";
        for(unsigned i=0; i<bulkChemistry->GetNumberChemicals(); i++)
        {
            std::cout<<"Bulk: "<< bulkChemistry->GetChemicalNamesByIndex(i)<<" ";
        }
        for(unsigned i=0; i<cellChemistry->GetNumberChemicals(); i++)
        {
            std::cout<<"Cell: "<< cellChemistry->GetChemicalNamesByIndex(i)<<" ";
        }
        std::cout<<std::endl;

        for (unsigned i=0; i<solutions.rGetTimes().size(); i++)
        {
            std::cout << solutions.rGetTimes()[i] << " , " << solutions.rGetSolutions()[i][0] << " , " << solutions.rGetSolutions()[i][1]<< " , " << solutions.rGetSolutions()[i][2] << " , " << solutions.rGetSolutions()[i][3]<<" , " << solutions.rGetSolutions()[i][4]<<"\n";
        }
    */
    }

    void TestChemicalSRNModel()
    {
    /*   
        ChemicalStructuresForTests* chemical_structure = new ChemicalStructuresForTests();
        chemical_structure->SetUpChemicalSRNModel();
        
        ChemicalSrnModel* p_srn_model = chemical_structure->rGetPtrChemicalSrnModel();

        // Create a vector of initial conditions
        std::vector<double> starter_conditions = {0.5,0.5};
        std::vector<std::string> species_names = {"U","Biomass"};
        p_srn_model->SetInitialConditions(starter_conditions);

        UniformG1GenerationalCellCycleModel* p_cc_model = new UniformG1GenerationalCellCycleModel();

        MAKE_PTR(WildTypeCellMutationState, p_healthy_state);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);

        CellPtr p_cell(new Cell(p_healthy_state, p_cc_model, p_srn_model, false, CellPropertyCollection()));
        p_cell->SetCellProliferativeType(p_diff_type);
        p_cell->InitialiseCellCycleModel();
        //p_cell->InitialiseSrnModel();

        chemical_structure->SetUpCellInitialConditions(p_cell, species_names, starter_conditions);

        // Now update the SRN
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        unsigned num_steps = 10;
        double end_time = 10.0;
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(end_time, num_steps);

        while (p_simulation_time->GetTime() < end_time)
        {
            std::cout<<"Time: "<<p_simulation_time->GetTime()<<std::endl;
            p_simulation_time->IncrementTimeOneStep();
            p_srn_model->SimulateToCurrentTime();
            std::cout<<"U: "<<p_srn_model->GetStateValueByName("U")<<" Biomass: "<<p_srn_model->GetStateValueByName("Biomass")<<std::endl;
        }
    */
    }

    void TestSimpleChemicalCellCycleModel()
    {
    /*
        ChemicalStructuresForTests* chemical_structure = new ChemicalStructuresForTests();
        chemical_structure->SetUpSimpleChemicalCellCycleModel();

        SimpleChemicalThresholdCellCycleModel* p_cell_cycle = chemical_structure->rGetPtrSimpleChemicalThresholdCellCycleModel();

        MAKE_PTR(WildTypeCellMutationState, p_healthy_state);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        CellPtr p_cell(new Cell(p_healthy_state, p_cell_cycle));
        p_cell->SetCellProliferativeType(p_stem_type);
        p_cell->InitialiseCellCycleModel();
        std::vector<std::string> speciesNames = {"U","Biomass","V"};
        std::vector<double> initValue = {1.0,1.0,1.0};
        chemical_structure->SetUpCellInitialConditions(p_cell, speciesNames, initValue);
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(10.0, 10);
        for (unsigned i = 0; i < 10; i++)
        {
            std::cout<<"Timestep: "<<i<<std::endl;
            p_simulation_time->IncrementTimeOneStep();

            std::cout<<"Species 0: "<<p_cell_cycle->GetSpeciesConcentrationsByIndex(0)<<"  Species 1: "<<p_cell_cycle->GetSpeciesConcentrationsByIndex(1)<<std::endl;
            std::cout<<"ReadyToDivide? "<<p_cell_cycle->ReadyToDivide()<<std::endl;
            std::cout<<"here"<<std::endl;
        }
    */ 
    }

    void TestSimpleChemicalCellCycleModelToDivide()
    {
        /*
        ChemicalStructuresForTests* chemical_structure = new ChemicalStructuresForTests();
        chemical_structure->SetUpSimpleChemicalCellCycleModel();

        SimpleChemicalThresholdCellCycleModel* p_cell_cycle = chemical_structure->rGetPtrSimpleChemicalThresholdCellCycleModel();

        MAKE_PTR(WildTypeCellMutationState, p_healthy_state);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        CellPtr p_cell(new Cell(p_healthy_state, p_cell_cycle));
        p_cell->SetCellProliferativeType(p_stem_type);
        p_cell->InitialiseCellCycleModel();
        std::vector<std::string> speciesNames = {"U","Biomass","V"};
        std::vector<double> initValue = {1.0,5.0,1.0};
        chemical_structure->SetUpCellInitialConditions(p_cell, speciesNames, initValue);
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(10.0, 10);
        for (unsigned i = 0; i < 10; i++)
        {
            std::cout<<"Timestep: "<<i<<std::endl;
            p_simulation_time->IncrementTimeOneStep();

            std::cout<<"Species 0: "<<p_cell_cycle->GetSpeciesConcentrationsByIndex(0)<<"  Species 1: "<<p_cell_cycle->GetSpeciesConcentrationsByIndex(1)<<std::endl;
            std::cout<<"ReadyToDivide? "<<p_cell_cycle->ReadyToDivide()<<std::endl;
        }
        */
    }

    void TestTemplatedChemicalDomainField()
    {
        /*
        ChemicalStructuresForTests* chemical_structure = new ChemicalStructuresForTests();

        chemical_structure -> SetUpChemicalDomainField();
        

        ChemicalDomainFieldTemplated<2,2,2>*& p_Pde_field = chemical_structure -> rGetPtrChemicalDomain();
        // check that the file input problem dimension is the same as the user defined problem dimension
        

        // print out the diffusion domain
        p_Pde_field -> PrintDiffusionDomain();
        // print out the diffusion domain mapped to the chaste cartesian grid
        p_Pde_field -> PrintMappedDiffusionDomain();
        // print out the ODE labels and keys
        p_Pde_field -> PrintDomainLabelKeys();
        // print out the ode domain
        p_Pde_field -> PrintODEDomain();
        // print out the ode domain mapped to the chaste cartesian grid
        p_Pde_field -> PrintMappedODEDomain();
        // print out the ODE labels and keys
        p_Pde_field -> PrintODELabelKeys();
        // print out the diffusion database
        p_Pde_field -> PrintDiffusionDatabase();

        TetrahedralMesh<2,2>* p_mesh = p_Pde_field -> rGetDomainFeMesh();

        std::cout<<"number of mesh nodes: "<<p_mesh ->GetNumNodes()<<std::endl;
        */
    }

    void TestChemicalDomainFieldForCellCoupling()
    {
        /*
        ChemicalStructuresForTests* chemical_structure = new ChemicalStructuresForTests();

        chemical_structure -> SetUpChemicalDomainFieldForCellCoupling();
        

        ChemicalDomainFieldForCellCoupling<2,2,2>*& p_Pde_field = chemical_structure -> rGetPtrChemicalDomainFieldForCellCoupling();
        // check that the file input problem dimension is the same as the user defined problem dimension
        
        

        // print out the diffusion domain
        p_Pde_field -> PrintDiffusionDomain();
        // print out the diffusion domain mapped to the chaste cartesian grid
        p_Pde_field -> PrintMappedDiffusionDomain();
        // print out the ODE labels and keys
        p_Pde_field -> PrintDomainLabelKeys();
        // print out the ode domain
        p_Pde_field -> PrintODEDomain();
        // print out the ode domain mapped to the chaste cartesian grid
        p_Pde_field -> PrintMappedODEDomain();
        // print out the ODE labels and keys
        p_Pde_field -> PrintODELabelKeys();
        // print out the diffusion database
        p_Pde_field -> PrintDiffusionDatabase();

        std::cout<<"initial conditions"<<std::endl;
        std::vector<double> init_conditions = p_Pde_field -> GetInitialNodeConditions();
        for(unsigned i=0; i<init_conditions.size(); i++)
        {
            std::cout<<init_conditions[i]<<std::endl;
        }


        TetrahedralMesh<2,2>* p_mesh = p_Pde_field -> rGetDomainFeMesh();

        std::cout<<"number of mesh nodes: "<<p_mesh ->GetNumNodes()<<std::endl;
        
        StateVariableRegister* testRegister = p_Pde_field -> GetStateVariableVector();

        for(unsigned index=0; index<testRegister ->GetNumberOfStateVariables(); index++)
        {
            std::cout<<"state name: "<<testRegister -> RetrieveStateVariableName(index)<<std::endl;
        }
        
        */
    }



    void TestCellProperties()
    {
        /*
        //======================================================================//
        //                      make a single cell object                       //
        //======================================================================//
        
        ChemicalStructuresForTests* chemical_structure = new ChemicalStructuresForTests();
   
        CellPtr p_cell = chemical_structure->SetUpCellObjectA();

        boost::shared_ptr<TransportCellProperty> p_transport_test = boost::static_pointer_cast<TransportCellProperty>(p_cell->rGetCellPropertyCollection().GetPropertiesType<TransportCellProperty>().GetProperty());
        boost::shared_ptr<MembraneCellProperty> p_cell_membrane_test = boost::static_pointer_cast<MembraneCellProperty>(p_cell->rGetCellPropertyCollection().GetPropertiesType<MembraneCellProperty>().GetProperty());

        //boost::shared_ptr<ChemicalCellProperty> p_cell_chemical_test = boost::static_pointer_cast<ChemicalCellProperty>(p_cell->rGetCellPropertyCollection().GetPropertiesType<ChemicalCellProperty>().GetProperty());
        //std::cout<<std::endl;
        //std::cout<<"Chemical cell property"<<std::endl;
        //std::cout<<std::endl;
        //std::cout<<"Call U from cell property: "<<p_cell_chemical -> GetCellConcentrationByName("U")<<std::endl;
        //std::cout<<"Call V from cell property: "<<p_cell_chemical -> GetCellConcentrationByName("V")<<std::endl;

        std::cout<<std::endl;
        std::cout<<"Membrane cell property"<<std::endl;
        std::cout<<std::endl;

        std::cout<<"Call thickness from membrane property: "<<p_cell_membrane_test -> GetMembraneThickness()<<std::endl;
        std::cout<<"Call isDouble from membrane property: "<<p_cell_membrane_test -> GetDoubleMembraneBool()<<std::endl;

        AbstractMembraneOdeSystem* p_membrane_ode_property = p_cell_membrane_test -> GetMembraneOdeSystem();
        boost::shared_ptr<AbstractIvpOdeSolver> p_membrane_ode_solver_property = p_cell_membrane_test -> GetMembraneOdeSolver();

        std::cout<<"Number of species: "<< p_membrane_ode_property -> GetNumberOfSpecies()<<std::endl;

        std::cout<<"Number of reactions: "<< p_membrane_ode_property -> GetNumberOfReactions()<<std::endl;

        std::cout<<"Number of equations: "<<p_membrane_ode_property->GetNumberOfStateVariables()<<std::endl;
 
        std::vector<double> membrane_initial_condition(p_membrane_ode_property -> GetNumberOfSpecies(), 1.0);

        OdeSolution membrane_solutions = p_membrane_ode_solver_property -> Solve(p_membrane_ode_property, membrane_initial_condition, 0, 1, 0.01, 0.1);


        AbstractChemistry* bulkTransportChemistry = p_membrane_ode_property->GetReactionSystem()->GetBulkChemistry();
        AbstractChemistry* cellTransportChemistry = p_membrane_ode_property->GetReactionSystem()->GetCellChemistry();
        std::cout<<"Time: ";
        for(unsigned i=0; i<bulkTransportChemistry->GetNumberChemicals(); i++)
        {
            std::cout<<"Bulk: "<< bulkTransportChemistry->GetChemicalNamesByIndex(i)<<" ";
        }
        for(unsigned i=0; i<cellTransportChemistry->GetNumberChemicals(); i++)
        {
            std::cout<<"Cell: "<< cellTransportChemistry->GetChemicalNamesByIndex(i)<<" ";
        }
        std::cout<<std::endl;

        for (unsigned i=0; i<membrane_solutions.rGetTimes().size(); i++)
        {
            std::cout << membrane_solutions.rGetTimes()[i] << " , " << membrane_solutions.rGetSolutions()[i][0] << " , " << membrane_solutions.rGetSolutions()[i][1]<< " , " << membrane_solutions.rGetSolutions()[i][2] << " , " << membrane_solutions.rGetSolutions()[i][3]<<" , " << membrane_solutions.rGetSolutions()[i][4]<<"\n";
        }


        std::cout<<std::endl;
        std::cout<<"Transport cell property"<<std::endl;
        std::cout<<std::endl;

        AbstractTransportOdeSystem* p_ode_property = p_transport_test -> GetTransportOdeSystem();
        boost::shared_ptr<AbstractIvpOdeSolver> p_ode_solver_property = p_transport_test -> GetTransportOdeSolver();

        std::cout<<"Number of species: "<< p_ode_property -> GetNumberOfSpecies()<<std::endl;

        std::cout<<"Number of reactions: "<< p_ode_property -> GetNumberOfReactions()<<std::endl;

        std::cout<<"Number of equations: "<<p_ode_property->GetNumberOfStateVariables()<<std::endl;
 
        std::vector<double> initial_condition(p_ode_property -> GetNumberOfSpecies(), 1.0);

        OdeSolution solutions = p_ode_solver_property -> Solve(p_ode_property, initial_condition, 0, 1, 0.01, 0.1);

        AbstractChemistry* bulkMembraneChemistry = p_ode_property->GetReactionSystem()->GetBulkChemistry();
        AbstractChemistry* cellMembraneChemistry = p_ode_property->GetReactionSystem()->GetCellChemistry();
        std::cout<<"Time: ";
        for(unsigned i=0; i<bulkMembraneChemistry->GetNumberChemicals(); i++)
        {
            std::cout<<"Bulk: "<< bulkMembraneChemistry->GetChemicalNamesByIndex(i)<<" ";
        }
        for(unsigned i=0; i<cellMembraneChemistry->GetNumberChemicals(); i++)
        {
            std::cout<<"Cell: "<< cellMembraneChemistry->GetChemicalNamesByIndex(i)<<" ";
        }
        std::cout<<std::endl;
        
        for (unsigned i=0; i<solutions.rGetTimes().size(); i++)
        {
            std::cout << solutions.rGetTimes()[i] << " , " << solutions.rGetSolutions()[i][0] << " , " << solutions.rGetSolutions()[i][1]<< " , " << solutions.rGetSolutions()[i][2] << " , " << solutions.rGetSolutions()[i][3]<<" , " << solutions.rGetSolutions()[i][4]<<"\n";
        }
    */
    }


    void TestChemicalCell()
    {
    /*
        //======================================================================//
        //                      make a single cell object                       //
        //======================================================================//
        
        ChemicalStructuresForTests* chemical_structure = new ChemicalStructuresForTests();
   
        CellPtr p_cell = chemical_structure->SetUpChemicalCellObjectA();

        boost::shared_ptr<TransportCellProperty> p_transport_test = boost::static_pointer_cast<TransportCellProperty>(p_cell->rGetCellPropertyCollection().GetPropertiesType<TransportCellProperty>().GetProperty());
        boost::shared_ptr<MembraneCellProperty> p_cell_membrane_test = boost::static_pointer_cast<MembraneCellProperty>(p_cell->rGetCellPropertyCollection().GetPropertiesType<MembraneCellProperty>().GetProperty());

        //boost::shared_ptr<ChemicalCellProperty> p_cell_chemical_test = boost::static_pointer_cast<ChemicalCellProperty>(p_cell->rGetCellPropertyCollection().GetPropertiesType<ChemicalCellProperty>().GetProperty());
        //std::cout<<std::endl;
        //std::cout<<"Chemical cell property"<<std::endl;
        //std::cout<<std::endl;
        //std::cout<<"Call U from cell property: "<<p_cell_chemical -> GetCellConcentrationByName("U")<<std::endl;
        //std::cout<<"Call V from cell property: "<<p_cell_chemical -> GetCellConcentrationByName("V")<<std::endl;

        std::cout<<std::endl;
        std::cout<<"Membrane cell property"<<std::endl;
        std::cout<<std::endl;

        std::cout<<"Call thickness from membrane property: "<<p_cell_membrane_test -> GetMembraneThickness()<<std::endl;
        std::cout<<"Call isDouble from membrane property: "<<p_cell_membrane_test -> GetDoubleMembraneBool()<<std::endl;

        AbstractMembraneOdeSystem* p_membrane_ode_property = p_cell_membrane_test -> GetMembraneOdeSystem();
        boost::shared_ptr<AbstractIvpOdeSolver> p_membrane_ode_solver_property = p_cell_membrane_test -> GetMembraneOdeSolver();

        std::cout<<"Number of species: "<< p_membrane_ode_property -> GetNumberOfSpecies()<<std::endl;

        std::cout<<"Number of reactions: "<< p_membrane_ode_property -> GetNumberOfReactions()<<std::endl;

        std::cout<<"Number of equations: "<<p_membrane_ode_property->GetNumberOfStateVariables()<<std::endl;
 
        std::vector<double> membrane_initial_condition(p_membrane_ode_property -> GetNumberOfSpecies(), 1.0);

        OdeSolution membrane_solutions = p_membrane_ode_solver_property -> Solve(p_membrane_ode_property, membrane_initial_condition, 0, 1, 0.01, 0.1);

        AbstractChemistry* bulkMembraneChemistry = p_membrane_ode_property->GetReactionSystem()->GetBulkChemistry();
        AbstractChemistry* cellMembraneChemistry = p_membrane_ode_property->GetReactionSystem()->GetCellChemistry();
        std::cout<<"Time: ";
        for(unsigned i=0; i<bulkMembraneChemistry->GetNumberChemicals(); i++)
        {
            std::cout<<"Bulk: "<< bulkMembraneChemistry->GetChemicalNamesByIndex(i)<<" ";
        }
        for(unsigned i=0; i<cellMembraneChemistry->GetNumberChemicals(); i++)
        {
            std::cout<<"Cell: "<< cellMembraneChemistry->GetChemicalNamesByIndex(i)<<" ";
        }
        std::cout<<std::endl;
        
        for (unsigned i=0; i<membrane_solutions.rGetTimes().size(); i++)
        {
            std::cout << membrane_solutions.rGetTimes()[i] << " , " << membrane_solutions.rGetSolutions()[i][0] << " , " << membrane_solutions.rGetSolutions()[i][1]<< " , " << membrane_solutions.rGetSolutions()[i][2] << " , " << membrane_solutions.rGetSolutions()[i][3]<<" , " << membrane_solutions.rGetSolutions()[i][4]<<"\n";
        }


        std::cout<<std::endl;
        std::cout<<"Transport cell property"<<std::endl;
        std::cout<<std::endl;

        AbstractTransportOdeSystem* p_ode_property = p_transport_test -> GetTransportOdeSystem();
        boost::shared_ptr<AbstractIvpOdeSolver> p_ode_solver_property = p_transport_test -> GetTransportOdeSolver();

        std::cout<<"Number of species: "<< p_ode_property -> GetNumberOfSpecies()<<std::endl;

        std::cout<<"Number of reactions: "<< p_ode_property -> GetNumberOfReactions()<<std::endl;

        std::cout<<"Number of equations: "<<p_ode_property->GetNumberOfStateVariables()<<std::endl;
 
        std::vector<double> initial_condition(p_ode_property -> GetNumberOfSpecies(), 1.0);

        OdeSolution solutions = p_ode_solver_property -> Solve(p_ode_property, initial_condition, 0, 1, 0.01, 0.1);

        AbstractChemistry* bulkTransportChemistry = p_ode_property->GetReactionSystem()->GetBulkChemistry();
        AbstractChemistry* cellTransportChemistry = p_ode_property->GetReactionSystem()->GetCellChemistry();
        std::cout<<"Time: ";
        for(unsigned i=0; i<bulkTransportChemistry->GetNumberChemicals(); i++)
        {
            std::cout<<"Bulk: "<< bulkTransportChemistry->GetChemicalNamesByIndex(i)<<" ";
        }
        for(unsigned i=0; i<cellTransportChemistry->GetNumberChemicals(); i++)
        {
            std::cout<<"Cell: "<< cellTransportChemistry->GetChemicalNamesByIndex(i)<<" ";
        }
        std::cout<<std::endl;

        for (unsigned i=0; i<solutions.rGetTimes().size(); i++)
        {
            std::cout << solutions.rGetTimes()[i] << " , " << solutions.rGetSolutions()[i][0] << " , " << solutions.rGetSolutions()[i][1]<< " , " << solutions.rGetSolutions()[i][2] << " , " << solutions.rGetSolutions()[i][3]<<" , " << solutions.rGetSolutions()[i][4]<<"\n";
        }
    */
    }



    void TestParabolicPdeSystemCoupledCellSystem()
    {
        /*
        std::cout<<"=============================================="<<std::endl;
        std::cout<<"TestParabolicPdeSystemCoupledCellSystem()"<<std::endl;
        std::cout<<"=============================================="<<std::endl;

        // make mesh of cell with associated mesh based cell population
        HoneycombMeshGenerator generator(3, 3);    // Parameters are: cells across, cells up
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

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
    
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);
        
        // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
        ChastePoint<2> lower(-5.0, -5.0);
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
       


        simulator.SetOutputDirectory("TestParabolicPdeSystemCoupledCellSystem_ver2");
        simulator.SetEndTime(30.0);

        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        p_linear_force->SetCutOffLength(3);
        simulator.AddForce(p_linear_force);

        std::cout<<"=============================================="<<std::endl;
        std::cout<<"OffLatticeSimulation -> AbstractCellBasedSumulation :: Solve()"<<std::endl;
        std::cout<<"=============================================="<<std::endl;
        simulator.Solve();

        */
    }

};

#endif