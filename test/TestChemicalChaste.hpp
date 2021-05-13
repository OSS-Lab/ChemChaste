#ifndef TESTCHEMICALCHASTE_HPP_
#define TESTCHEMICALCHASTE_HPP_

// chaste includes
#include <cxxtest/TestSuite.h>
#include "UblasIncludes.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "SmartPointers.hpp"


// general includes
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <tuple>
#include <cmath>

// chemical includes
#include "AbstractChemical.hpp"
#include "AbstractDiffusiveChemical.hpp"
#include "AbstractChemistry.hpp"
#include "AbstractDiffusiveChemistry.hpp" 

// reaction-diffusion includes
#include "AbstractReaction.hpp"
#include "AbstractReversibleReaction.hpp"
#include "AbstractReactionSystem.hpp"
#include "MassActionReaction.hpp"
#include "AbstractReactionSystemFromFile.hpp"
#include "ReactionTypeDatabase.hpp"

// chaste PdeOde includes
#include "HoneycombMeshGenerator.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "LinearParabolicPdeSystemWithCoupledOdeSystemSolver.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "ConstBoundaryCondition.hpp"
#include "OutputFileHandler.hpp"
#include "RandomNumberGenerator.hpp"
#include "TrianglesMeshReader.hpp"

// custom pdeOde includes
#include "PdeSchnackenbergCoupledPdeOdeSystem.hpp"
#include "OdeSchnackenbergCoupledPdeOdeSystem.hpp"
#include "PdeConsumerProducer.hpp"
#include "OdeConsumerProducer.hpp"

#include "AbstractChemicalOdeSystem.hpp"
#include "AbstractChemicalOdeForCoupledPdeSystem.hpp"

// CHASTE Spheroid includes
#include "SimpleOxygenBasedCellCycleModel.hpp"
#include "WildTypeCellMutationState.hpp"
#include "StemCellProliferativeType.hpp"
#include "CellwiseSourceEllipticPde.hpp"
#include "ConstBoundaryCondition.hpp"
#include "EllipticGrowingDomainPdeModifier.hpp"
#include "OffLatticeSimulation.hpp"
#include "GeneralisedLinearSpringForce.hpp"


#include "SchnackenbergCoupledPdeSystem.hpp"



class TestChemicalChaste : public AbstractCellBasedTestSuite
{
public:
    void TestChemicalClass()
    {
        std::cout<<"-----------------------------"<<std::endl;
        std::cout<<"Chemical"<<std::endl;

        AbstractChemical *p_chemical = new AbstractChemical("A");
        AbstractDiffusiveChemical *p_chemical_diffusive = new AbstractDiffusiveChemical("B");
        p_chemical_diffusive -> AddDiffusiveDomain("test",1.0);

        std::cout<<"Name: "<<p_chemical ->GetChemicalName() <<std::endl;
        std::cout<<"Type: "<<p_chemical ->GetChemicalType() <<std::endl;
        std::cout<<"Name: "<<p_chemical_diffusive ->GetChemicalName() <<std::endl;
        std::cout<<"Type: "<<p_chemical_diffusive ->GetChemicalType() <<std::endl;
        std::cout<<"Diffusivity: "<<p_chemical_diffusive ->GetDiffusiveDomainVector()[0] <<std::endl;

        // vectorise the pointers
        std::cout<<"Chemical vector"<<std::endl;
        std::vector<AbstractChemical*> p_chemicalVector;
        // implicit upcasting
        p_chemicalVector.push_back(p_chemical);
        p_chemicalVector.push_back(p_chemical_diffusive);

        std::cout<<p_chemicalVector[0] -> GetChemicalType() <<std::endl;
        std::cout<<p_chemicalVector[1] -> GetChemicalType() <<std::endl;

        // get the diffusive value for [1]
        std::cout<<"dynamic casting"<<std::endl;
        AbstractDiffusiveChemical *p_chemical_diffusive_2 = dynamic_cast<AbstractDiffusiveChemical*>(p_chemicalVector[1]);
        std::cout<<p_chemical_diffusive_2 -> GetChemicalDiffusivityVector()[0] <<std::endl;
    }

    void TestChemistryClass()
    {
        std::cout<<"-----------------------------"<<std::endl;
        std::cout<<"Abstract chemistry"<<std::endl;
        unsigned Number_of_species  =2;
        std::vector<std::string> ChemicalNames = {"U", "V"};
        std::vector<double> DiffusionRates = {1.0, 2.0};
        std::vector<bool> IsDiffusing = {true, true};
        std::vector<std::string> DiffusionDomains = {"Bulk", "Bulk"};

        AbstractChemistry* p_chemistry = new AbstractChemistry();

        for(unsigned species=0; species<Number_of_species; species++)
        {
            AbstractChemical *p_chemical = new AbstractChemical(ChemicalNames[species]);
            p_chemistry -> AddChemical(p_chemical);
            
        }
        std::vector<std::string> ChemNames = p_chemistry -> GetChemicalNames();
        
        std::cout<<"Chemical names: "<<std::endl;

        for(unsigned i=0; i<Number_of_species; i++)
        {
            std::cout<<ChemNames[i]<<std::endl;
        }
       
        std::cout<<"Abstract diffusive chemistry"<<std::endl;
        std::vector<std::string> ChemicalNamesDiffusive = {"Ud", "Vd"};
        AbstractDiffusiveChemistry* p_diffusive_chemistry = new AbstractDiffusiveChemistry();

        for(unsigned species=0; species<Number_of_species; species++)
        {
            AbstractDiffusiveChemical *p_chemical = new AbstractDiffusiveChemical(ChemicalNamesDiffusive[species]);
            p_chemical -> AddDiffusiveDomain(DiffusionDomains[species],DiffusionRates[species]);
            p_diffusive_chemistry -> AddChemical(p_chemical);
        }
        std::cout<<"Chemical diffusivities: "<<std::endl;

        for(unsigned i=0; i<Number_of_species; i++)
        {
            std::cout<<p_diffusive_chemistry -> GetChemicalNamesByIndex(i)<<std::endl;
            std::cout<<p_diffusive_chemistry -> GetDiffusivityVectorByIndex(i)[0]<<std::endl;
        }


        std::cout<<"Add chemistries"<<std::endl;

        // upcast AbstractDiffusiveChemistry to AbstractChemistry
        // maybe need to template this function? currently add to the highest class to avoid object splicing
        AbstractChemistry *p_newChemistry = dynamic_cast<AbstractChemistry*>(p_diffusive_chemistry);

        p_chemistry -> AddChemistry(p_newChemistry);
        std::vector<std::string> ChemistryNames = p_chemistry -> GetChemicalNames();
        for(unsigned i=0; i<Number_of_species; i++)
        {
            std::cout<<ChemistryNames[i]<<std::endl;
        }
       

        std::cout<<"Add diffusive chemistries"<<std::endl;

        // upcast AbstractDiffusiveChemistry to AbstractChemistry
        // maybe need to template this function? currently add to the highest class to avoid object splicing
        AbstractDiffusiveChemistry *p_newDiffusiveChemistry = new AbstractDiffusiveChemistry();
        std::vector<std::string> NewChemicalNamesDiffusive = {"Ud", "Y"};
        std::vector<double> NewDiffusionRates = {3,4};
        for(unsigned species=0; species<Number_of_species; species++)
        {
            AbstractDiffusiveChemical *p_chemical = new AbstractDiffusiveChemical(NewChemicalNamesDiffusive[species]);
            p_chemical -> AddDiffusiveDomain(DiffusionDomains[species],NewDiffusionRates[species]);
            p_newDiffusiveChemistry -> AddChemical(p_chemical); //shouldn't add the first species Ud as is a duplicate of name and domain
        }

        p_diffusive_chemistry -> AddChemistry(p_newDiffusiveChemistry);
        
        std::cout<<"Number of chemicals: "<<p_diffusive_chemistry ->GetNumberChemicals()<<std::endl;
        std::cout<<"Number of diffusive chemicals: "<<p_diffusive_chemistry ->GetNumberDiffusiveChemicals()<<std::endl;
        for(unsigned i=0; i<4; i++)
        {
            std::cout<<p_diffusive_chemistry -> GetChemicalNamesByIndex(i)<<std::endl;
            std::cout<<p_diffusive_chemistry -> GetDiffusivityValueByChemicalName(p_diffusive_chemistry -> GetChemicalNamesByIndex(i))<<std::endl;
        }

    }

    void TestReactionClass()
    {
        std::cout<<"-----------------------------"<<std::endl;

        std::cout<<"Chemical information"<<std::endl;
        unsigned Number_of_species  =2;
        std::vector<std::string> ChemicalNames = {"U", "V"};
        std::vector<double> DiffusionRates = {1.0, 2.0};
        std::vector<bool> IsDiffusing = {true, true};
        std::vector<std::string> DiffusionDomains = {"Bulk", "Bulk"};

        std::cout<<"Form chemical vector"<<std::endl;
        std::vector<AbstractChemical*> p_chemicalVector;
        

        for(unsigned species=0; species<Number_of_species; species++)
        {
            if(IsDiffusing[species])
            {
                // use the diffusive chemical root
                AbstractDiffusiveChemical *p_chemical = new AbstractDiffusiveChemical(ChemicalNames[species]);
                p_chemical -> AddDiffusiveDomain(DiffusionDomains[species],DiffusionRates[species]);
                p_chemicalVector.push_back(p_chemical);
            }else{
                // use the non-diffusive chemical root
                AbstractChemical *p_chemical = new AbstractChemical(ChemicalNames[species]);

                p_chemicalVector.push_back(p_chemical);
            }
        }

        
        
        // check iterating through the chemical vector works, need dynamic casting of the iterator
        for (std::vector<AbstractChemical*>::iterator chem_iter = p_chemicalVector.begin();
             chem_iter != p_chemicalVector.end();
             ++chem_iter)
        {
            std::string ChemicalName = "NA";
            double Diffusivity = 0.0;
            std::string Domain = "None";
            
            AbstractChemical *p_2 = dynamic_cast<AbstractChemical*>(*chem_iter);
            
            if( p_2 -> GetChemicalType() == "AbstractDiffusiveChemical")
            {
                AbstractDiffusiveChemical *p_chemical = dynamic_cast<AbstractDiffusiveChemical*>(*chem_iter);
                ChemicalName = p_chemical -> GetChemicalName();
                Diffusivity = p_chemical -> GetChemicalDiffusivityVector()[0];
                Domain = p_chemical -> GetDiffusiveDomainVector()[0];
            }
            
            std::cout<<"Name: "<<ChemicalName<<std::endl;
            std::cout<<"Diffusivity: "<<Diffusivity<<std::endl;
            std::cout<<"Domain: "<<Domain<<std::endl;
            delete p_2;
        }

        std::cout<<"Schnackenberg Chemistry"<<std::endl;

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


        // r1: 2U + V -> 3U     forwardRate = 0.1
        // r2: 0 <-> U          forwardRate = 0.1 reverseRate = 0.2
        // r3: 0 -> V           forwardRate = 0.3

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
        // add U to reactions
        p_substrates_1.push_back(p_chemical_V);
        stoich_substrates_1.push_back(1);
        p_products_3.push_back(p_chemical_V);
        stoich_products_3.push_back(1);

        double reaction_1_rate = 0.1;
        double reaction_2_forward_rate = 0.1;
        double reaction_2_reverse_rate = 0.2;
        double reaction_3_rate = 0.3;

        AbstractReaction* p_reaction_1 = new AbstractReaction(p_substrates_1, p_products_1, stoich_substrates_1, stoich_products_1,reaction_1_rate);
        AbstractReversibleReaction* p_reaction_2 = new AbstractReversibleReaction(p_substrates_2, p_products_2, stoich_substrates_2, stoich_products_2,reaction_2_forward_rate, reaction_2_reverse_rate);
        AbstractReaction* p_reaction_3 = new AbstractReaction(p_substrates_3, p_products_3, stoich_substrates_3, stoich_products_3,reaction_3_rate);

        // concentration vector
        std::vector<double> concentration_vector= {1.0,1.0};
        std::vector<double> change_concentration_vector= {0.0,0.0};

        std::cout<<"Starting conditions"<<std::endl;
        std::cout<<"Concentration (U,V): "<<concentration_vector[0]<<" "<<concentration_vector[1]<<std::endl;
        std::cout<<"Concentration change (U,V): "<<change_concentration_vector[0]<<" "<<change_concentration_vector[1]<<std::endl;
        std::cout<<"============================"<<std::endl;
        std::cout<<"r1: 2U + V -> 3U     forwardRate = "<<reaction_1_rate<<std::endl;
        change_concentration_vector= {0.0,0.0};
        p_reaction_1 -> React(p_system_chemistry,concentration_vector,change_concentration_vector);
        //concentration_vector[0] = concentration_vector[0] + change_concentration_vector[0];
        //concentration_vector[1] = concentration_vector[1] + change_concentration_vector[1];
        std::cout<<"Concentration (U,V): "<<concentration_vector[0]<<" "<<concentration_vector[1]<<std::endl;
        std::cout<<"Concentration change (U,V): "<<change_concentration_vector[0]<<" "<<change_concentration_vector[1]<<std::endl;
        std::cout<<"============================"<<std::endl;
        std::cout<<"r2: 0 <-> U          forwardRate = "<<reaction_2_forward_rate<<" reverseRate = "<<reaction_2_reverse_rate<<std::endl;
        change_concentration_vector= {0.0,0.0};
        p_reaction_2 -> React(p_system_chemistry,concentration_vector,change_concentration_vector);
        //concentration_vector[0] = concentration_vector[0] + change_concentration_vector[0];
        //concentration_vector[1] = concentration_vector[1] + change_concentration_vector[1];
        std::cout<<"Concentration (U,V): "<<concentration_vector[0]<<" "<<concentration_vector[1]<<std::endl;
        std::cout<<"Concentration change (U,V): "<<change_concentration_vector[0]<<" "<<change_concentration_vector[1]<<std::endl;
        std::cout<<"============================"<<std::endl;
        std::cout<<"r3: 0 -> V           forwardRate = "<<reaction_3_rate<<std::endl;
        change_concentration_vector= {0.0,0.0};
        p_reaction_3 -> React(p_system_chemistry,concentration_vector,change_concentration_vector);
        //concentration_vector[0] = concentration_vector[0] + change_concentration_vector[0];
        //concentration_vector[1] = concentration_vector[1] + change_concentration_vector[1];
        std::cout<<"Concentration (U,V): "<<concentration_vector[0]<<" "<<concentration_vector[1]<<std::endl;
        std::cout<<"Concentration change (U,V): "<<change_concentration_vector[0]<<" "<<change_concentration_vector[1]<<std::endl;


        // form reaction system
        std::cout<<"Form reaction system: "<<std::endl;
        concentration_vector= {1.0,1.0};
        change_concentration_vector= {0.0,0.0};

        std::vector<AbstractReaction*> p_reaction_vector_1;
        p_reaction_vector_1.push_back(p_reaction_1);

        AbstractReactionSystem* p_reaction_system_1 = new AbstractReactionSystem(p_system_chemistry, p_reaction_vector_1);
        p_reaction_system_1 -> ReactSystem(concentration_vector,change_concentration_vector);
        std::cout<<"============================"<<std::endl;
        std::cout<<"r1: 2U + V -> 3U     forwardRate = "<<reaction_1_rate<<std::endl;
        std::cout<<"Concentration (U,V): "<<concentration_vector[0]<<" "<<concentration_vector[1]<<std::endl;
        std::cout<<"Concentration change (U,V): "<<change_concentration_vector[0]<<" "<<change_concentration_vector[1]<<std::endl;


        // form mixed reaction system (AbstractReaction, AbstractReversibleReaction, AbstractReaction)
        concentration_vector= {1.0,1.0};
        change_concentration_vector= {0.0,0.0};
        std::vector<AbstractReaction*> p_reaction_vector_2;
        p_reaction_vector_2.push_back(p_reaction_1);
        p_reaction_vector_2.push_back(p_reaction_2);
        p_reaction_vector_2.push_back(p_reaction_3);

        AbstractReactionSystem* p_reaction_system_2 = new AbstractReactionSystem(p_system_chemistry, p_reaction_vector_2);
        p_reaction_system_2 -> ReactSystem(concentration_vector,change_concentration_vector);
        std::cout<<"Concentration (U,V): "<<concentration_vector[0]<<" "<<concentration_vector[1]<<std::endl;
        std::cout<<"Concentration change (U,V): "<<change_concentration_vector[0]<<" "<<change_concentration_vector[1]<<std::endl;
        for(unsigned i=0; i<p_reaction_system_2 -> GetNumberOfReactions(); i++)
        {
            std::cout<<"Reaction type: "<<i<<" "<< p_reaction_system_2 -> GetReactionByIndex(i) -> GetReactionType()<<std::endl;
        }
    }

    void TestMassActionReactionClass()
    {
        std::cout<<"-----------------------------"<<std::endl;

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


        // r1: 2U + V -> 3U     forwardRate = 0.1
        // r2: 0 <-> U          forwardRate = 0.1 reverseRate = 0.2
        // r3: 0 -> V           forwardRate = 0.3

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
        // add U to reactions
        p_substrates_1.push_back(p_chemical_V);
        stoich_substrates_1.push_back(1);
        p_products_3.push_back(p_chemical_V);
        stoich_products_3.push_back(1);

        double reaction_1_rate = 0.1;
        double reaction_2_forward_rate = 0.1;
        double reaction_2_reverse_rate = 0.2;
        double reaction_3_rate = 0.3;
    
        MassActionReaction* p_reaction_1 = new MassActionReaction(p_substrates_1, p_products_1, stoich_substrates_1, stoich_products_1,false,false,reaction_1_rate);
        MassActionReaction* p_reaction_2 = new MassActionReaction(p_substrates_2, p_products_2, stoich_substrates_2, stoich_products_2,false, true, reaction_2_forward_rate, reaction_2_reverse_rate);
        MassActionReaction* p_reaction_3 = new MassActionReaction(p_substrates_3, p_products_3, stoich_substrates_3, stoich_products_3,false, false,reaction_3_rate);

        // concentration vector
        std::vector<double> concentration_vector= {1.0,1.0};
        std::vector<double> change_concentration_vector= {0.0,0.0};

        std::cout<<"Mass action kinetics"<<std::endl;
        std::cout<<"Starting conditions"<<std::endl;
        std::cout<<"Concentration (U,V): "<<concentration_vector[0]<<" "<<concentration_vector[1]<<std::endl;
        std::cout<<"Concentration change (U,V): "<<change_concentration_vector[0]<<" "<<change_concentration_vector[1]<<std::endl;
        std::cout<<"============================"<<std::endl;
        std::cout<<"r1: 2U + V -> 3U     forwardRate = "<<reaction_1_rate<<std::endl;
        change_concentration_vector= {0.0,0.0};
        p_reaction_1 -> React(p_system_chemistry,concentration_vector,change_concentration_vector);
    
        std::cout<<"Concentration (U,V): "<<concentration_vector[0]<<" "<<concentration_vector[1]<<std::endl;
        std::cout<<"Concentration change (U,V): "<<change_concentration_vector[0]<<" "<<change_concentration_vector[1]<<std::endl;
        std::cout<<"============================"<<std::endl;
        std::cout<<"r2: 0 <-> U          forwardRate = "<<reaction_2_forward_rate<<" reverseRate = "<<reaction_2_reverse_rate<<std::endl;
        change_concentration_vector= {0.0,0.0};
        p_reaction_2 -> React(p_system_chemistry,concentration_vector,change_concentration_vector);
        std::cout<<"Concentration (U,V): "<<concentration_vector[0]<<" "<<concentration_vector[1]<<std::endl;
        std::cout<<"Concentration change (U,V): "<<change_concentration_vector[0]<<" "<<change_concentration_vector[1]<<std::endl;
        std::cout<<"============================"<<std::endl;
        std::cout<<"r3: 0 -> V           forwardRate = "<<reaction_3_rate<<std::endl;
        change_concentration_vector= {0.0,0.0};
        p_reaction_3 -> React(p_system_chemistry,concentration_vector,change_concentration_vector);
        std::cout<<"Concentration (U,V): "<<concentration_vector[0]<<" "<<concentration_vector[1]<<std::endl;
        std::cout<<"Concentration change (U,V): "<<change_concentration_vector[0]<<" "<<change_concentration_vector[1]<<std::endl;

        // form mixed mass action reaction system 
        std::cout<<"============================"<<std::endl;
        std::cout<<"Form mass action reaction sysstem"<<std::endl;
        std::cout<<"============================"<<std::endl;
        concentration_vector= {1.0,1.0};
        change_concentration_vector= {0.0,0.0};
        std::vector<AbstractReaction*> p_mass_action_reaction_vector;
        p_mass_action_reaction_vector.push_back(p_reaction_1);
        p_mass_action_reaction_vector.push_back(p_reaction_2);
        p_mass_action_reaction_vector.push_back(p_reaction_3);

        AbstractReactionSystem* p_mass_action_reaction_system = new AbstractReactionSystem(p_system_chemistry, p_mass_action_reaction_vector);
        p_mass_action_reaction_system -> ReactSystem(concentration_vector,change_concentration_vector);
        std::cout<<"Concentration (U,V): "<<concentration_vector[0]<<" "<<concentration_vector[1]<<std::endl;
        std::cout<<"Concentration change (U,V): "<<change_concentration_vector[0]<<" "<<change_concentration_vector[1]<<std::endl;
        for(unsigned i=0; i<p_mass_action_reaction_system -> GetNumberOfReactions(); i++)
        {
            std::cout<<"Reaction type: "<<i<<" "<< p_mass_action_reaction_system -> GetReactionByIndex(i) -> GetReactionType()<<std::endl;
        }
    }

    void TestSpatialPdeOdeSolver()
    {
        /*
        std::cout<<"SchnackenbergCoupledPdeOdeSystem"<<std::endl;
        // system properties
        const unsigned probDim =2;
        const unsigned spaceDim=2;
        const unsigned elementDim=2;
        std::vector<double> initValues = {2.0, 0.75};
        std::vector<double> bcValues = {2.0, 0.75};

        // mesh
        HoneycombMeshGenerator generator(10, 10, 0);
        MutableMesh<elementDim,spaceDim>* p_mesh = generator.GetMesh();

        // Process Boundary Conditions
        std::cout<<"Process Boundary Conditions"<<std::endl;
        BoundaryConditionsContainer<elementDim,spaceDim,probDim> bcc;
        std::vector<bool> areNeumannBoundaryConditions(probDim, false);
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
                init_conds[probDim*i + pdeDim] = fabs(initValues[pdeDim] + RandomNumberGenerator::Instance()->ranf());
            }
        }
        // PETSc Vec
        std::cout<<"PETSc Vec"<<std::endl;
        Vec initial_condition = PetscTools::CreateVec(init_conds);

        // pde system
        std::cout<<"Pde"<<std::endl;
        PdeSchnackenbergCoupledPdeOdeSystem<elementDim, spaceDim, probDim> pde(1e-4, 1e-2);
  
        // coupled ode system
        std::cout<<"Ode loop"<<std::endl;
        std::vector<AbstractOdeSystemForCoupledPdeSystem*> odeSystem;
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++){
            // number of ode system objects must match the number of nodes, i.e the individual odes may be multi-dimensional
            if(i==22 || i==77)
            {
                odeSystem.push_back(new OdeSchnackenbergCoupledPdeOdeSystem(0.1, 0.2, 0.3, 0.1));

            }else if(i==27 || i==72)
            {
                odeSystem.push_back(new OdeSchnackenbergCoupledPdeOdeSystem(0.1, 0.2, 0.3, 0.1));  
            }else{
                odeSystem.push_back(new OdeSchnackenbergCoupledPdeOdeSystem(0.1, 0.2, 0.3, 0.1));
            }
        }
        std::cout<<"Number odesAtNodes: "<<p_mesh->GetNumNodes()<<std::endl;

        // used the explicitly defined EulerIvpOdeSolver rather than the default defined BackwardEuler method

        //EulerIvpOdeSolver euler_solver; 
        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);

        std::cout<<"Solver"<<std::endl;
        // solver
        LinearParabolicPdeSystemWithCoupledOdeSystemSolver<elementDim,spaceDim,probDim> solver(p_mesh, &pde, &bcc, odeSystem, p_solver);

        // solver properties
        double t_end = 10;
        solver.SetTimes(0, t_end);
        solver.SetTimeStep(1e-2);
        solver.SetSamplingTimeStep(1e-2);
        solver.SetOutputDirectory("TestVectorisedSchnackenbergOutput_correct_ks");
        solver.SetInitialCondition(initial_condition);
        // solve
        std::cout<<"Solve"<<std::endl;
        //solver.SolveAndWriteResultsToFile();
        solver.SolveAndWriteResultsToFile();
        std::cout<<"Clean"<<std::endl;

        // clean
        PetscTools::Destroy(initial_condition);
        */
    }

    void TestSpatialConsumerProducerSolver()
    {
        /*
        std::cout<<"ConsumerProducer"<<std::endl;

        //  reaction system involving two species, A and B
        // 0 -> A   rateConstant = k1
        // B -> 0   rateConstant = k2
        // A <-> B  rateConstantForward = k3    rateConstantReverse = k_3
        // A diffuses at rate Da
        // B diffuses at rate Db

        // nodal ode of the form:   OdeConsumerProducer(k1, k2, k3, k4)
        // with nodes 22, 27, 72, 77 selected forming corners of a square domain offset from the boundary

        double Da = 1e-1;
        double Db = 5e-2;

        // system properties
        const unsigned probDim =2;
        const unsigned spaceDim=2;
        const unsigned elementDim=2;
        std::vector<double> initValues = {0, 0};
        std::vector<double> bcValues = {0, 0};

        // mesh
        HoneycombMeshGenerator generator(10, 10, 0);
        MutableMesh<elementDim,spaceDim>* p_mesh = generator.GetMesh();

        // Process Boundary Conditions
        std::cout<<"Process Boundary Conditions"<<std::endl;
        BoundaryConditionsContainer<elementDim,spaceDim,probDim> bcc;
        std::vector<bool> areNeumannBoundaryConditions(probDim, false);
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
                init_conds[probDim*i + pdeDim] =    initValues[pdeDim] ;//fabs(initValues[pdeDim] + RandomNumberGenerator::Instance()->ranf());
            }
        }
        // PETSc Vec
        std::cout<<"PETSc Vec"<<std::endl;
        Vec initial_condition = PetscTools::CreateVec(init_conds);

        // pde system
        std::cout<<"Pde"<<std::endl;
        PdeConsumerProducer<elementDim, spaceDim, probDim> pde(Da, Db);
  
        // coupled ode system
        std::cout<<"Ode loop"<<std::endl;
        std::vector<AbstractOdeSystemForCoupledPdeSystem*> odeSystem;
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++){
            // number of ode system objects must match the number of nodes, i.e the individual odes may be multi-dimensional
            if(i==22)
            {
                odeSystem.push_back(new OdeConsumerProducer(1, 0, 0, 0));

            }else if(i==77)
            {
                odeSystem.push_back(new OdeConsumerProducer(0, 1, 0, 0));  
            }else{
                odeSystem.push_back(new OdeConsumerProducer(0, 0, 0.1, 0));
            }
        }

        std::cout<<"Number odesAtNodes: "<<p_mesh->GetNumNodes()<<std::endl;

        // used the explicitly defined EulerIvpOdeSolver rather than the default defined BackwardEuler method

        //EulerIvpOdeSolver euler_solver; 
        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);

        std::cout<<"Solver"<<std::endl;
        // solver
        LinearParabolicPdeSystemWithCoupledOdeSystemSolver<elementDim,spaceDim,probDim> solver(p_mesh, &pde, &bcc, odeSystem, p_solver);

        // solver properties
        double t_end = 10;
        solver.SetTimes(0, t_end);
        solver.SetTimeStep(1e-2);
        solver.SetSamplingTimeStep(1e-2);
        solver.SetOutputDirectory("TestConsumerProducerOutput_1");
        solver.SetInitialCondition(initial_condition);
        // solve
        std::cout<<"Solve"<<std::endl;
        //solver.SolveAndWriteResultsToFile();
        solver.SolveAndWriteResultsToFile();
        std::cout<<"Clean"<<std::endl;

        // clean
        PetscTools::Destroy(initial_condition);
        */
    }

    void TestChemicalOde()
    {

        /*
        std::cout<<"Mass action kinetics of Schnackenberg using the AbstractChemicalOdeSystem Class"<<std::endl;
        std::cout<<"-----------------------------"<<std::endl;
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
        // add U to reactions
        p_substrates_1.push_back(p_chemical_V);
        stoich_substrates_1.push_back(1);
        p_products_3.push_back(p_chemical_V);
        stoich_products_3.push_back(1);

        double reaction_1_rate = 0.1;
        double reaction_2_forward_rate = 0.1;
        double reaction_2_reverse_rate = 2;
        double reaction_3_rate = 0.3;

        AbstractReaction* p_reaction_1 = new AbstractReaction(p_substrates_1, p_products_1, stoich_substrates_1, stoich_products_1,reaction_1_rate);
        AbstractReversibleReaction* p_reaction_2 = new AbstractReversibleReaction(p_substrates_2, p_products_2, stoich_substrates_2, stoich_products_2,reaction_2_forward_rate, reaction_2_reverse_rate);
        AbstractReaction* p_reaction_3 = new AbstractReaction(p_substrates_3, p_products_3, stoich_substrates_3, stoich_products_3,reaction_3_rate);

        std::vector<AbstractReaction*> p_reaction_vector;
        p_reaction_vector.push_back(p_reaction_1);
        p_reaction_vector.push_back(p_reaction_2);
        p_reaction_vector.push_back(p_reaction_3);
    


        MassActionReaction* p_Mass_action_reaction_1 = new MassActionReaction(p_substrates_1, p_products_1, stoich_substrates_1, stoich_products_1,false,false,reaction_1_rate);
        MassActionReaction* p_Mass_action_reaction_2 = new MassActionReaction(p_substrates_2, p_products_2, stoich_substrates_2, stoich_products_2,false, true, reaction_2_forward_rate, reaction_2_reverse_rate);
        MassActionReaction* p_Mass_action_reaction_3 = new MassActionReaction(p_substrates_3, p_products_3, stoich_substrates_3, stoich_products_3,false, false,reaction_3_rate);

        std::vector<AbstractReaction*> p_mass_action_reaction_vector;
        p_mass_action_reaction_vector.push_back(p_Mass_action_reaction_1);
        p_mass_action_reaction_vector.push_back(p_Mass_action_reaction_2);
        p_mass_action_reaction_vector.push_back(p_Mass_action_reaction_3);

        AbstractReactionSystem* p_mass_action_reaction_system = new AbstractReactionSystem(p_system_chemistry, p_mass_action_reaction_vector);



        AbstractChemicalOdeSystem chemicalOde(p_mass_action_reaction_system);
        EulerIvpOdeSolver euler_solver;
        std::vector<double> initial_condition = {1.0, 1.0};

        OdeSolution solutions = euler_solver.Solve(&chemicalOde, initial_condition, 0, 1, 0.01, 0.1);
        for (unsigned i=0; i<solutions.rGetTimes().size(); i++)
        {
            std::cout << solutions.rGetTimes()[i] << " " << solutions.rGetSolutions()[i][0] << " " << solutions.rGetSolutions()[i][1]<< "\n";
        }
        */
    }

    void TestChemicalOdePde()
    {
        
        std::cout<<"Mass action kinetics of Schnackenberg using the AbstractChemicalOdeSystem Class"<<std::endl;
        std::cout<<"-----------------------------"<<std::endl;
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

        // r1: 2U + V -> 3U     forwardRate = 0.1
        // r2: 0 <-> U          forwardRate = 0.1 reverseRate = 0.2
        // r3: 0 -> V           forwardRate = 0.3

        double reaction_1_rate = 1.0;
        double reaction_2_forward_rate = 0.5;
        double reaction_2_reverse_rate = 2.2;
        double reaction_3_rate = 1.5;

        AbstractReaction* p_reaction_1 = new AbstractReaction(p_substrates_1, p_products_1, stoich_substrates_1, stoich_products_1,reaction_1_rate);
        AbstractReversibleReaction* p_reaction_2 = new AbstractReversibleReaction(p_substrates_2, p_products_2, stoich_substrates_2, stoich_products_2,reaction_2_forward_rate, reaction_2_reverse_rate);
        AbstractReaction* p_reaction_3 = new AbstractReaction(p_substrates_3, p_products_3, stoich_substrates_3, stoich_products_3,reaction_3_rate);

        std::vector<AbstractReaction*> p_reaction_vector;
        p_reaction_vector.push_back(p_reaction_1);
        p_reaction_vector.push_back(p_reaction_2);
        p_reaction_vector.push_back(p_reaction_3);
    
        //AbstractReactionSystem* p_reaction_system = new AbstractReactionSystem(p_system_chemistry, p_reaction_vector);

        //SchnackenbergCoupledPdeSystem<2> pde(1e-4, 1e-2, 0.5, 2.2, 1.5, 1);
    
        MassActionReaction* p_mass_action_reaction_1 = new MassActionReaction(p_substrates_1, p_products_1, stoich_substrates_1, stoich_products_1,false,false,reaction_1_rate);
        MassActionReaction* p_mass_action_reaction_2 = new MassActionReaction(p_substrates_2, p_products_2, stoich_substrates_2, stoich_products_2,false, true, reaction_2_forward_rate, reaction_2_reverse_rate);
        MassActionReaction* p_mass_action_reaction_3 = new MassActionReaction(p_substrates_3, p_products_3, stoich_substrates_3, stoich_products_3,false, false,reaction_3_rate);

        std::vector<AbstractReaction*> p_mass_action_reaction_vector;
        p_mass_action_reaction_vector.push_back(p_mass_action_reaction_1);
        p_mass_action_reaction_vector.push_back(p_mass_action_reaction_2);
        p_mass_action_reaction_vector.push_back(p_mass_action_reaction_3);

        AbstractReactionSystem* p_mass_action_reaction_system = new AbstractReactionSystem(p_system_chemistry, p_mass_action_reaction_vector);

        // form ode system
        std::cout<<"SchnackenbergCoupledPdeOdeSystem  -As AbstractChemicalODeSystem"<<std::endl;
        // system properties
        const unsigned probDim =2;
        const unsigned spaceDim=2;
        const unsigned elementDim=2;
        std::vector<double> initValues = {2.0, 0.75};
        std::vector<double> bcValues = {0.0, 0.0};

        // mesh
        HoneycombMeshGenerator generator(10, 10, 0);
        MutableMesh<elementDim,spaceDim>* p_mesh = generator.GetMesh();

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
        std::cout<<"Initial conditions"<<std::endl;
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
        std::cout<<"PETSc Vec"<<std::endl;
        Vec initial_condition = PetscTools::CreateVec(init_conds);

        
        //SchnackenbergCoupledPdeSystem<2> pde(1e-4, 1e-2, 0.5, 2.2, 1.5, 1);
        // coupled ode system
        std::cout<<"Ode loop"<<std::endl;
        std::vector<AbstractOdeSystemForCoupledPdeSystem*> odeSystem;
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++){
            // number of ode system objects must match the number of nodes, i.e the individual odes may be multi-dimensional
            odeSystem.push_back(new AbstractChemicalOdeForCoupledPdeSystem(p_mass_action_reaction_system));
        }
        std::cout<<"Number odesAtNodes: "<<p_mesh->GetNumNodes()<<std::endl;

        // pde system
        std::cout<<"Pde"<<std::endl;
        PdeSchnackenbergCoupledPdeOdeSystem<elementDim, spaceDim, probDim> pde(odeSystem[0],1e-4, 1e-2 );

        // used the explicitly defined EulerIvpOdeSolver rather than the default defined BackwardEuler method

        //EulerIvpOdeSolver euler_solver; 
        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);

        std::cout<<"Solver"<<std::endl;
        // solver
        LinearParabolicPdeSystemWithCoupledOdeSystemSolver<elementDim,spaceDim,probDim> solver(p_mesh, &pde, &bcc,odeSystem,p_solver);

        // solver properties
        double t_end = 10;
        solver.SetTimes(0, t_end);
        solver.SetTimeStep(1e-2);
        solver.SetSamplingTimeStep(1e-2);
        solver.SetOutputDirectory("TestAbstractChemicalOdeOutput_groupMeeting");
        solver.SetInitialCondition(initial_condition);
        // solve
        std::cout<<"Solve"<<std::endl;
        //solver.SolveAndWriteResultsToFile();
        solver.SolveAndWriteResultsToFile();
        std::cout<<"Clean"<<std::endl;

        // clean
        PetscTools::Destroy(initial_condition);
        
    }

    void TestReactionSystemFromFile()
    {
        /*
        std::cout<<"Reaction system from file"<<std::endl;
        std::cout<<"-----------------------------"<<std::endl;

        AbstractReaction* p_reaction = new AbstractReaction();
        std::cout<<"Before cast reaction type: "<<p_reaction -> GetReactionType()<<std::endl;
        ReactionTablet(p_reaction,"MassActionReaction");
        std::cout<<"After cast reaction type: "<<p_reaction -> GetReactionType()<<std::endl;


        std::cout<<"-----------------------------"<<std::endl;
    
        std::string reactionFilename = "/home/chaste/projects/ChemicalChaste/src/Data/SchnackenbergReactionFile.txt";
        //std::string reactionFilename = "/home/chaste/projects/ChemicalChaste/src/Data/SchnackenbergReactionFileMixed.txt";
        AbstractReactionSystemFromFile* p_file_reaction_system = new  AbstractReactionSystemFromFile(reactionFilename);

        


        std::cout<<"Mass action kinetics of Schnackenberg using the AbstractChemicalOdeSystem Class"<<std::endl;
        std::cout<<"-----------------------------"<<std::endl;
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
        // add U to reactions
        p_substrates_1.push_back(p_chemical_V);
        stoich_substrates_1.push_back(1);
        p_products_3.push_back(p_chemical_V);
        stoich_products_3.push_back(1);

        double reaction_1_rate = 0.1;
        double reaction_2_forward_rate = 0.1;
        double reaction_2_reverse_rate = 0.2;
        double reaction_3_rate = 0.3;

        AbstractReaction* p_reaction_1 = new AbstractReaction(p_substrates_1, p_products_1, stoich_substrates_1, stoich_products_1,reaction_1_rate);
        AbstractReversibleReaction* p_reaction_2 = new AbstractReversibleReaction(p_substrates_2, p_products_2, stoich_substrates_2, stoich_products_2,reaction_2_forward_rate, reaction_2_reverse_rate);
        AbstractReaction* p_reaction_3 = new AbstractReaction(p_substrates_3, p_products_3, stoich_substrates_3, stoich_products_3,reaction_3_rate);

        std::vector<AbstractReaction*> p_reaction_vector;
        p_reaction_vector.push_back(p_reaction_1);
        p_reaction_vector.push_back(p_reaction_2);
        p_reaction_vector.push_back(p_reaction_3);
    

        MassActionReaction* p_Mass_action_reaction_1 = new MassActionReaction(p_substrates_1, p_products_1, stoich_substrates_1, stoich_products_1,false,false,reaction_1_rate);
        MassActionReaction* p_Mass_action_reaction_2 = new MassActionReaction(p_substrates_2, p_products_2, stoich_substrates_2, stoich_products_2,false, true, reaction_2_forward_rate, reaction_2_reverse_rate);
        MassActionReaction* p_Mass_action_reaction_3 = new MassActionReaction(p_substrates_3, p_products_3, stoich_substrates_3, stoich_products_3,false, false,reaction_3_rate);

        std::vector<AbstractReaction*> p_mass_action_reaction_vector;
        p_mass_action_reaction_vector.push_back(p_Mass_action_reaction_1);
        p_mass_action_reaction_vector.push_back(p_Mass_action_reaction_2);
        p_mass_action_reaction_vector.push_back(p_Mass_action_reaction_3);

        AbstractReactionSystem* p_mass_action_reaction_system = new AbstractReactionSystem(p_system_chemistry, p_mass_action_reaction_vector);


        
        std::cout<<"--------------------------------------"<<std::endl;
        std::cout<<"Test the reaction systems for equality"<<std::endl;
        std::cout<<"--------------------------------------"<<std::endl;

        std::cout<<"-----------------------------"<<std::endl;
        std::cout<<"System from file"<<std::endl;
        std::cout<<"-----------------------------"<<std::endl;
        std::cout<<"Read out reaction details: "<<std::endl;
        std::cout<<"-----------------------------"<<std::endl;
        std::cout<<"Number of reactions: "<<p_file_reaction_system -> GetNumberOfReactions()<<std::endl;
        for(unsigned i=0; i<p_file_reaction_system -> GetNumberOfReactions(); i++)
        {
            std::cout<<"Reaction: "<<i<<" "<< p_file_reaction_system -> GetReactionByIndex(i) -> GetReactionType()<<std::endl;
            std::cout<<"Reaction: "<<i<<" "<< dynamic_cast<MassActionReaction*>(p_file_reaction_system -> GetReactionByIndex(i)) -> GetForwardReactionRateConstant()<<std::endl;
            std::cout<<"Reaction: "<<i<<" "<< dynamic_cast<MassActionReaction*>(p_file_reaction_system -> GetReactionByIndex(i)) -> GetReverseReactionRateConstant()<<std::endl;
        }

        std::vector<std::string> chemNames = p_file_reaction_system-> GetSystemChemistry() -> GetChemicalNames();
        std::cout<<"System chemical names:"<<std::endl;
        for(unsigned i=0; i<chemNames.size();i++)
        {
            std::cout<<chemNames[i]<<std::endl;
        }

        for(unsigned i=0; i<p_file_reaction_system-> GetNumberOfReactions(); i++ )
        {
            AbstractReaction* p_reaction = p_file_reaction_system-> GetReactionByIndex(i);

            std::cout<<"Reaction type: "<<p_reaction ->GetReactionType()<<std::endl;
            for(unsigned j=0; j<p_reaction -> GetNumberOfSubstrates(); j++)
            {
                std::cout<<"Substrate: "<<p_reaction -> GetSubstratesByIndex(j) -> GetChemicalName()<<" Stoich: "<<p_reaction ->GetStoichSubstratesByIndex(j)<<std::endl;
            }

            for(unsigned j=0; j<p_reaction -> GetNumberOfProducts(); j++)
            {
                std::cout<<"Product: "<<p_reaction -> GetProductsByIndex(j) -> GetChemicalName()<<" Stoich: "<<p_reaction ->GetStoichProductsByIndex(j)<<std::endl;
            }

        }

        std::cout<<"-----------------------------"<<std::endl;
        std::cout<<"Hard coded system"<<std::endl;
        std::cout<<"-----------------------------"<<std::endl;
        std::cout<<"Read out reaction details: "<<std::endl;
        std::cout<<"-----------------------------"<<std::endl;
        std::cout<<"Number of reactions: "<<p_mass_action_reaction_system -> GetNumberOfReactions()<<std::endl;
        for(unsigned i=0; i<p_mass_action_reaction_system -> GetNumberOfReactions(); i++)
        {
            std::cout<<"Reaction: "<<i<<" "<< p_mass_action_reaction_system -> GetReactionByIndex(i) -> GetReactionType()<<std::endl;
            std::cout<<"Reaction: "<<i<<" "<< dynamic_cast<MassActionReaction*>(p_mass_action_reaction_system -> GetReactionByIndex(i)) -> GetForwardReactionRateConstant()<<std::endl;
            std::cout<<"Reaction: "<<i<<" "<< dynamic_cast<MassActionReaction*>(p_mass_action_reaction_system -> GetReactionByIndex(i)) -> GetReverseReactionRateConstant()<<std::endl;
        }

        std::vector<std::string> chemNames_hard_coded = p_mass_action_reaction_system-> GetSystemChemistry() -> GetChemicalNames();
        std::cout<<"System chemical names:"<<std::endl;
        for(unsigned i=0; i<chemNames_hard_coded.size();i++)
        {
            std::cout<<chemNames_hard_coded[i]<<std::endl;
        }

        for(unsigned i=0; i<p_mass_action_reaction_system-> GetNumberOfReactions(); i++ )
        {
            AbstractReaction* p_reaction = p_mass_action_reaction_system-> GetReactionByIndex(i);

            std::cout<<"Reaction type: "<<p_reaction ->GetReactionType()<<std::endl;
            for(unsigned j=0; j<p_reaction -> GetNumberOfSubstrates(); j++)
            {
                std::cout<<"Substrate: "<<p_reaction -> GetSubstratesByIndex(j) -> GetChemicalName()<<" Stoich: "<<p_reaction ->GetStoichSubstratesByIndex(j)<<std::endl;
            }

            for(unsigned j=0; j<p_reaction -> GetNumberOfProducts(); j++)
            {
                std::cout<<"Product: "<<p_reaction -> GetProductsByIndex(j) -> GetChemicalName()<<" Stoich: "<<p_reaction ->GetStoichProductsByIndex(j)<<std::endl;
            }

        }



        std::cout<<"-----------------------------"<<std::endl;
        std::cout<<"Run the reaction ODE system  "<<std::endl;
        std::cout<<"-----------------------------"<<std::endl;


        AbstractChemicalOdeSystem chemicalOde(p_mass_action_reaction_system);
        EulerIvpOdeSolver euler_solver;
        std::vector<double> initial_condition = {1.0, 1.0};

        OdeSolution solutions = euler_solver.Solve(&chemicalOde, initial_condition, 0, 1, 0.1, 0.1);
        for (unsigned i=0; i<solutions.rGetTimes().size(); i++)
        {
            std::cout << solutions.rGetTimes()[i] << " " << solutions.rGetSolutions()[i][0] << " " << solutions.rGetSolutions()[i][1]<< "\n";
        }



        std::cout<<"-----------------------------"<<std::endl;
        std::cout<<"Run the file reaction system "<<std::endl;
        std::cout<<"-----------------------------"<<std::endl;
        // implicit upcast AbstractReactionSystemFromFile to AbstractReactionSystem
        AbstractChemicalOdeSystem chemical_ode_file(p_file_reaction_system);
        std::vector<double> initial_condition_file = {1.0, 1.0};
        EulerIvpOdeSolver euler_solver_file;
        OdeSolution solutions_file = euler_solver_file.Solve(&chemical_ode_file, initial_condition_file, 0, 1, 0.1, 0.1);
        for (unsigned i=0; i<solutions_file.rGetTimes().size(); i++)
        {
            std::cout << solutions_file.rGetTimes()[i] << " " << solutions_file.rGetSolutions()[i][0] << " " << solutions_file.rGetSolutions()[i][1]<< "\n";
        }
        */
    }

    void TestPdeFromFile()
    {
        /*
        std::cout<<"ConsumerProducerFromFile"<<std::endl;

        //  reaction system involving two species, A and B
        // 0 -> A   rateConstant = k1
        // B -> 0   rateConstant = k2
        // A <-> B  rateConstantForward = k3    rateConstantReverse = k_3
        // A diffuses at rate Da
        // B diffuses at rate Db

        // nodal ode of the form:   OdeConsumerProducer(k1, k2, k3, k4)
        // with nodes 22, 27, 72, 77 selected forming corners of a square domain offset from the boundary

        double Da = 1e-1;
        double Db = 5e-2;

        // system properties
        const unsigned probDim =2;
        const unsigned spaceDim=2;
        const unsigned elementDim=2;
        std::vector<double> initValues = {0, 0};
        std::vector<double> bcValues = {0, 0};

        // mesh
        HoneycombMeshGenerator generator(10, 10, 0);
        MutableMesh<elementDim,spaceDim>* p_mesh = generator.GetMesh();

        // Process Boundary Conditions
        std::cout<<"Process Boundary Conditions"<<std::endl;
        BoundaryConditionsContainer<elementDim,spaceDim,probDim> bcc;
        std::vector<bool> areNeumannBoundaryConditions(probDim, false);
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
                init_conds[probDim*i + pdeDim] =    initValues[pdeDim] ;//fabs(initValues[pdeDim] + RandomNumberGenerator::Instance()->ranf());
            }
        }
        // PETSc Vec
        std::cout<<"PETSc Vec"<<std::endl;
        Vec initial_condition = PetscTools::CreateVec(init_conds);

        // pde system
        std::cout<<"Pde"<<std::endl;
        PdeConsumerProducer<elementDim, spaceDim, probDim> pde(Da, Db);
  
        // coupled ode system
        std::cout<<"Ode loop"<<std::endl;
        std::vector<AbstractOdeSystemForCoupledPdeSystem*> odeSystem;
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++){
            // number of ode system objects must match the number of nodes, i.e the individual odes may be multi-dimensional
            if(i==22)
            {
                odeSystem.push_back(new OdeConsumerProducer(1, 0, 0, 0));

            }else if(i==77)
            {
                odeSystem.push_back(new OdeConsumerProducer(0, 1, 0, 0));  
            }else{
                odeSystem.push_back(new OdeConsumerProducer(0, 0, 0.1, 0));
            }
        }
        // 0 -> A   rateConstant = k1
        // B -> 0   rateConstant = k2
        // A <-> B  rateConstantForward = k3    rateConstantReverse = k_3
        // A diffuses at rate Da
        // B diffuses at rate Db


        std::cout<<"Number odesAtNodes: "<<p_mesh->GetNumNodes()<<std::endl;

        // used the explicitly defined EulerIvpOdeSolver rather than the default defined BackwardEuler method

        //EulerIvpOdeSolver euler_solver; 
        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);

        std::cout<<"Solver"<<std::endl;
        // solver
        LinearParabolicPdeSystemWithCoupledOdeSystemSolver<elementDim,spaceDim,probDim> solver(p_mesh, &pde, &bcc, odeSystem, p_solver);

        // solver properties
        double t_end = 10;
        solver.SetTimes(0, t_end);
        solver.SetTimeStep(1e-2);
        solver.SetSamplingTimeStep(1e-2);
        solver.SetOutputDirectory("TestConsumerProducerOutputFromFile");
        solver.SetInitialCondition(initial_condition);
        // solve
        std::cout<<"Solve"<<std::endl;
        //solver.SolveAndWriteResultsToFile();
        solver.SolveAndWriteResultsToFile();
        std::cout<<"Clean"<<std::endl;

        // clean
        PetscTools::Destroy(initial_condition);

        */
    }

    void TestCHASTESpheroidTutorial()
    {
        /*
        EXIT_IF_PARALLEL;

        HoneycombMeshGenerator generator(10, 10, 0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;

        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);

        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            SimpleOxygenBasedCellCycleModel* p_model = new SimpleOxygenBasedCellCycleModel;
            p_model->SetDimension(2);
            CellPtr p_cell(new Cell(p_state, p_model));
            p_cell->SetCellProliferativeType(p_stem_type);

            
            p_model->SetStemCellG1Duration(8.0);
            p_model->SetTransitCellG1Duration(8.0);

            
            double birth_time = - RandomNumberGenerator::Instance()->ranf() *
                                 (  p_model->GetStemCellG1Duration()
                                  + p_model->GetSG2MDuration() );
            
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        MAKE_PTR_ARGS(CellwiseSourceEllipticPde<2>, p_pde, (cell_population, -0.03));

        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (1.0));
        bool is_neumann_bc = false;

        MAKE_PTR_ARGS(EllipticGrowingDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, is_neumann_bc));
        p_pde_modifier->SetDependentVariableName("oxygen");

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.AddSimulationModifier(p_pde_modifier);

        simulator.SetOutputDirectory("SpheroidTutorial");
        simulator.SetEndTime(1.0);

        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        p_linear_force->SetCutOffLength(3);
        simulator.AddForce(p_linear_force);


        simulator.Solve();
*/


    }
    
    void TestChemicalSpheroid()
    {
        /*
        std::vector<double> initValues = {2.0, 0.75};
        std::vector<double> bcValues = {2.0, 0.75};

        // mesh
        HoneycombMeshGenerator generator(100, 100, 0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();


        SchnackenbergCoupledPdeSystem<2> pde(1e-4, 1e-2, 0.1, 0.2, 0.3, 0.1);

        BoundaryConditionsContainer<2,2,2> bcc;
        ConstBoundaryCondition<2>* p_bc_for_u = new ConstBoundaryCondition<2>(2.0);
        ConstBoundaryCondition<2>* p_bc_for_v = new ConstBoundaryCondition<2>(0.75);
        for (TetrahedralMesh<2,2>::BoundaryNodeIterator node_iter = p_mesh->GetBoundaryNodeIteratorBegin();
             node_iter != p_mesh->GetBoundaryNodeIteratorEnd();
             ++node_iter)
        {
            bcc.AddDirichletBoundaryCondition(*node_iter, p_bc_for_u, 0);
            bcc.AddDirichletBoundaryCondition(*node_iter, p_bc_for_v, 1);
        }

        LinearParabolicPdeSystemWithCoupledOdeSystemSolver<2,2,2> solver(p_mesh, &pde, &bcc);

        double t_end = 10;
        solver.SetTimes(0, t_end);
        solver.SetTimeStep(1e-2);
        solver.SetSamplingTimeStep(1e-2);
        solver.SetOutputDirectory("TestSchnackenbergSystemOnHoneycombMesh");

        std::vector<double> init_conds(2*p_mesh->GetNumNodes());
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            init_conds[2*i] = fabs(2.0 + RandomNumberGenerator::Instance()->ranf());
            init_conds[2*i + 1] = fabs(0.75 + RandomNumberGenerator::Instance()->ranf());
        }
        Vec initial_condition = PetscTools::CreateVec(init_conds);
        solver.SetInitialCondition(initial_condition);

        solver.SolveAndWriteResultsToFile();

        PetscTools::Destroy(initial_condition);

        */
    }

    void TestChemicalSpheroidPAPER()
    {
        
        std::vector<double> initValues = {2.0, 0.75};
        std::vector<double> bcValues = {2.0, 0.75};

        // mesh
        HoneycombMeshGenerator generator(3, 3, 0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();


        SchnackenbergCoupledPdeSystem<2> pde(1e-4, 1e-2, 0.5, 2.2, 1.5, 1);

        BoundaryConditionsContainer<2,2,2> bcc;
        ConstBoundaryCondition<2>* p_bc_for_u = new ConstBoundaryCondition<2>(2.0);
        ConstBoundaryCondition<2>* p_bc_for_v = new ConstBoundaryCondition<2>(0.75);
        for (TetrahedralMesh<2,2>::BoundaryNodeIterator node_iter = p_mesh->GetBoundaryNodeIteratorBegin();
             node_iter != p_mesh->GetBoundaryNodeIteratorEnd();
             ++node_iter)
        {
            bcc.AddDirichletBoundaryCondition(*node_iter, p_bc_for_u, 0);
            bcc.AddDirichletBoundaryCondition(*node_iter, p_bc_for_v, 1);
        }

        LinearParabolicPdeSystemWithCoupledOdeSystemSolver<2,2,2> solver(p_mesh, &pde, &bcc);

        double t_end = 10;
        solver.SetTimes(0, t_end);
        solver.SetTimeStep(1e-2);
        solver.SetSamplingTimeStep(1e-2);
        solver.SetOutputDirectory("TestSchnackenbergSystemOnHoneycombMesh_paper_test");

        std::vector<double> init_conds(2*p_mesh->GetNumNodes());
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            init_conds[2*i] = fabs(2.0 + RandomNumberGenerator::Instance()->ranf());
            init_conds[2*i + 1] = fabs(0.75 + RandomNumberGenerator::Instance()->ranf());
        }
        Vec initial_condition = PetscTools::CreateVec(init_conds);
        solver.SetInitialCondition(initial_condition);

        solver.SolveAndWriteResultsToFile();

        PetscTools::Destroy(initial_condition);
        
    }

    void TestSpectatorDependentReactionClass()
    {
        /*
        std::cout<<"--------------------------------------"<<std::endl;
        std::cout<<"Test Spectator dependent reaction class"<<std::endl;
        std::cout<<"--------------------------------------"<<std::endl;

        std::string reactionFilename = "/home/chaste/projects/ChemicalChaste/src/Data/SpectatorReactionFile.txt";
        AbstractReactionSystemFromFile* p_file_reaction_system = new  AbstractReactionSystemFromFile(reactionFilename);



        std::cout<<"-----------------------------"<<std::endl;
        std::cout<<"System from file"<<std::endl;
        std::cout<<"-----------------------------"<<std::endl;
        std::cout<<"Read out reaction details: "<<std::endl;
        std::cout<<"-----------------------------"<<std::endl;
        std::cout<<"Number of reactions: "<<p_file_reaction_system -> GetNumberOfReactions()<<std::endl;
        for(unsigned i=0; i<p_file_reaction_system -> GetNumberOfReactions(); i++)
        {
            std::cout<<"Reaction: "<<i<<" "<< p_file_reaction_system -> GetReactionByIndex(i) -> GetReactionType()<<std::endl;
            std::cout<<"Number of spectators: "<<i<<" "<< dynamic_cast<SpectatorDependentReaction*>(p_file_reaction_system -> GetReactionByIndex(i)) -> GetNumberOfSpectators()<<std::endl;
        }

        std::vector<std::string> chemNames = p_file_reaction_system-> GetSystemChemistry() -> GetChemicalNames();
        std::cout<<"System chemical names:"<<std::endl;
        for(unsigned i=0; i<chemNames.size();i++)
        {
            std::cout<<chemNames[i]<<std::endl;
        }

        for(unsigned i=0; i<p_file_reaction_system-> GetNumberOfReactions(); i++ )
        {
            AbstractReaction* p_reaction = p_file_reaction_system-> GetReactionByIndex(i);

            std::cout<<"Reaction type: "<<p_reaction ->GetReactionType()<<std::endl;
            for(unsigned j=0; j<p_reaction -> GetNumberOfSubstrates(); j++)
            {
                std::cout<<"Substrate: "<<p_reaction -> GetSubstratesByIndex(j) -> GetChemicalName()<<" Stoich: "<<p_reaction ->GetStoichSubstratesByIndex(j)<<std::endl;
            }

            for(unsigned j=0; j<p_reaction -> GetNumberOfProducts(); j++)
            {
                std::cout<<"Product: "<<p_reaction -> GetProductsByIndex(j) -> GetChemicalName()<<" Stoich: "<<p_reaction ->GetStoichProductsByIndex(j)<<std::endl;
            }

        }

        std::cout<<"-----------------------------"<<std::endl;
        std::cout<<"Run the file reaction system "<<std::endl;
        std::cout<<"-----------------------------"<<std::endl;
        // implicit upcast AbstractReactionSystemFromFile to AbstractReactionSystem
        AbstractChemicalOdeSystem chemical_ode_file(p_file_reaction_system);
        // read in the order Alpha, Beta, A, B, C
        // order speceis occurance in reaction system, then order of spectator species when duplicates removed
        std::vector<double> initial_condition = {1.0, 1.0,1.0,1.0,0.0};

        // initial_condition = {1.0, 1.0,1.0,1.0,0.0}, in order Alpha, Beta, A, B, C
        // C set to 0.0 makes reaction 1 have zero rate, reaction 0 occurs at coanstant rate, while reaction 2 has variable rate


        EulerIvpOdeSolver euler_solver_file;
        OdeSolution solutions_file = euler_solver_file.Solve(&chemical_ode_file, initial_condition, 0, 1, 0.1, 0.1);
        for (unsigned i=0; i<solutions_file.rGetTimes().size(); i++)
        {
            std::cout << "Time: "<<solutions_file.rGetTimes()[i] << " " << solutions_file.rGetSolutions()[i][0] << " " << solutions_file.rGetSolutions()[i][1]<< " " << solutions_file.rGetSolutions()[i][2]<< "\n";
        }

        */
    }

    void TestChemicalOdePdeConvergence()
    {
        /*
        std::cout<<"Mass action kinetics of Schnackenberg using the AbstractChemicalOdeSystem Class"<<std::endl;
        std::cout<<"-----------------------------"<<std::endl;
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

        // r1: 2U + V -> 3U     forwardRate = 0.1
        // r2: 0 <-> U          forwardRate = 0.1 reverseRate = 0.2
        // r3: 0 -> V           forwardRate = 0.3

        double reaction_1_rate = 1.0;
        double reaction_2_forward_rate = 0.5;
        double reaction_2_reverse_rate = 2.2;
        double reaction_3_rate = 1.5;

        AbstractReaction* p_reaction_1 = new AbstractReaction(p_substrates_1, p_products_1, stoich_substrates_1, stoich_products_1,reaction_1_rate);
        AbstractReversibleReaction* p_reaction_2 = new AbstractReversibleReaction(p_substrates_2, p_products_2, stoich_substrates_2, stoich_products_2,reaction_2_forward_rate, reaction_2_reverse_rate);
        AbstractReaction* p_reaction_3 = new AbstractReaction(p_substrates_3, p_products_3, stoich_substrates_3, stoich_products_3,reaction_3_rate);

        std::vector<AbstractReaction*> p_reaction_vector;
        p_reaction_vector.push_back(p_reaction_1);
        p_reaction_vector.push_back(p_reaction_2);
        p_reaction_vector.push_back(p_reaction_3);
    
        //AbstractReactionSystem* p_reaction_system = new AbstractReactionSystem(p_system_chemistry, p_reaction_vector);

        //SchnackenbergCoupledPdeSystem<2> pde(1e-4, 1e-2, 0.5, 2.2, 1.5, 1);
    
        MassActionReaction* p_mass_action_reaction_1 = new MassActionReaction(p_substrates_1, p_products_1, stoich_substrates_1, stoich_products_1,false,false,reaction_1_rate);
        MassActionReaction* p_mass_action_reaction_2 = new MassActionReaction(p_substrates_2, p_products_2, stoich_substrates_2, stoich_products_2,false, true, reaction_2_forward_rate, reaction_2_reverse_rate);
        MassActionReaction* p_mass_action_reaction_3 = new MassActionReaction(p_substrates_3, p_products_3, stoich_substrates_3, stoich_products_3,false, false,reaction_3_rate);

        std::vector<AbstractReaction*> p_mass_action_reaction_vector;
        p_mass_action_reaction_vector.push_back(p_mass_action_reaction_1);
        p_mass_action_reaction_vector.push_back(p_mass_action_reaction_2);
        p_mass_action_reaction_vector.push_back(p_mass_action_reaction_3);

        AbstractReactionSystem* p_mass_action_reaction_system = new AbstractReactionSystem(p_system_chemistry, p_mass_action_reaction_vector);

        // form ode system
        std::cout<<"SchnackenbergCoupledPdeOdeSystem  -As AbstractChemicalODeSystem"<<std::endl;
        // system properties
        const unsigned probDim =2;
        const unsigned spaceDim=2;
        const unsigned elementDim=2;
        std::vector<double> initValues = {2.0, 0.75};
        std::vector<double> bcValues = {0.0, 0.0};

        // mesh
        HoneycombMeshGenerator generator(10, 10, 0);
        MutableMesh<elementDim,spaceDim>* p_mesh = generator.GetMesh();

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
        std::cout<<"Initial conditions"<<std::endl;
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
        std::cout<<"PETSc Vec"<<std::endl;
        Vec initial_condition = PetscTools::CreateVec(init_conds);

        
        //SchnackenbergCoupledPdeSystem<2> pde(1e-4, 1e-2, 0.5, 2.2, 1.5, 1);
        // coupled ode system
        std::cout<<"Ode loop"<<std::endl;
        std::vector<AbstractOdeSystemForCoupledPdeSystem*> odeSystem;
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++){
            // number of ode system objects must match the number of nodes, i.e the individual odes may be multi-dimensional
            odeSystem.push_back(new AbstractChemicalOdeForCoupledPdeSystem(p_mass_action_reaction_system));
        }
        std::cout<<"Number odesAtNodes: "<<p_mesh->GetNumNodes()<<std::endl;

        // pde system
        std::cout<<"Pde"<<std::endl;
        PdeSchnackenbergCoupledPdeOdeSystem<elementDim, spaceDim, probDim> pde(odeSystem[0],1e-4, 1e-2 );

        // used the explicitly defined EulerIvpOdeSolver rather than the default defined BackwardEuler method

        //EulerIvpOdeSolver euler_solver; 
        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);

        double t_end = 10;
        std::cout<<"Solver 1"<<std::endl;
        // solver
        LinearParabolicPdeSystemWithCoupledOdeSystemSolver<elementDim,spaceDim,probDim> solver1(p_mesh, &pde, &bcc,odeSystem,p_solver);

        // solver properties
        
        solver1.SetTimes(0, t_end);
        solver1.SetTimeStep(1e-1);
        solver1.SetSamplingTimeStep(1e-1);
        solver1.SetOutputDirectory("TestAbstractChemicalOdeOutput_groupMeeting_1e-1");
        solver1.SetInitialCondition(initial_condition);
        // solve
        std::cout<<"Solve"<<std::endl;
        //solver.SolveAndWriteResultsToFile();
        solver1.SolveAndWriteResultsToFile();
        std::cout<<"Clean"<<std::endl;


        std::cout<<"Solver 2"<<std::endl;
        // solver
        LinearParabolicPdeSystemWithCoupledOdeSystemSolver<elementDim,spaceDim,probDim> solver2(p_mesh, &pde, &bcc,odeSystem,p_solver);

        solver2.SetTimes(0, t_end);
        solver2.SetTimeStep(1e-2);
        solver2.SetSamplingTimeStep(1e-2);
        solver2.SetOutputDirectory("TestAbstractChemicalOdeOutput_groupMeeting_1e-2");
        solver2.SetInitialCondition(initial_condition);
        // solve
        std::cout<<"Solve"<<std::endl;
        //solver.SolveAndWriteResultsToFile();
        solver2.SolveAndWriteResultsToFile();
        std::cout<<"Clean"<<std::endl;

        std::cout<<"Solver 3"<<std::endl;
        // solver
        LinearParabolicPdeSystemWithCoupledOdeSystemSolver<elementDim,spaceDim,probDim> solver3(p_mesh, &pde, &bcc,odeSystem,p_solver);

        // solver properties
        
        solver3.SetTimes(0, t_end);
        solver3.SetTimeStep(1e-3);
        solver3.SetSamplingTimeStep(1e-3);
        solver3.SetOutputDirectory("TestAbstractChemicalOdeOutput_groupMeeting_1e-3");
        solver3.SetInitialCondition(initial_condition);
        // solve
        std::cout<<"Solve"<<std::endl;
        //solver.SolveAndWriteResultsToFile();
        solver3.SolveAndWriteResultsToFile();
        std::cout<<"Clean"<<std::endl;


        std::cout<<"Solver 4"<<std::endl;
        // solver
        LinearParabolicPdeSystemWithCoupledOdeSystemSolver<elementDim,spaceDim,probDim> solver4(p_mesh, &pde, &bcc,odeSystem,p_solver);

        solver4.SetTimes(0, t_end);
        solver4.SetTimeStep(1e-4);
        solver4.SetSamplingTimeStep(1e-4);
        solver4.SetOutputDirectory("TestAbstractChemicalOdeOutput_groupMeeting_1e-4");
        solver4.SetInitialCondition(initial_condition);
        // solve
        std::cout<<"Solve"<<std::endl;
        //solver.SolveAndWriteResultsToFile();
        solver4.SolveAndWriteResultsToFile();
        std::cout<<"Clean"<<std::endl;


        std::cout<<"Solver 5"<<std::endl;
        // solver
        LinearParabolicPdeSystemWithCoupledOdeSystemSolver<elementDim,spaceDim,probDim> solver5(p_mesh, &pde, &bcc,odeSystem,p_solver);

        // solver properties
        
        solver5.SetTimes(0, t_end);
        solver5.SetTimeStep(1e-5);
        solver5.SetSamplingTimeStep(1e-5);
        solver5.SetOutputDirectory("TestAbstractChemicalOdeOutput_groupMeeting_1e-5");
        solver5.SetInitialCondition(initial_condition);
        // solve
        std::cout<<"Solve"<<std::endl;
        //solver.SolveAndWriteResultsToFile();
        solver5.SolveAndWriteResultsToFile();
        std::cout<<"Clean"<<std::endl;


        std::cout<<"Solver 6"<<std::endl;
        // solver
        LinearParabolicPdeSystemWithCoupledOdeSystemSolver<elementDim,spaceDim,probDim> solver6(p_mesh, &pde, &bcc,odeSystem,p_solver);

        solver6.SetTimes(0, t_end);
        solver6.SetTimeStep(1e-6);
        solver6.SetSamplingTimeStep(1e-6);
        solver6.SetOutputDirectory("TestAbstractChemicalOdeOutput_groupMeeting_1e-6");
        solver6.SetInitialCondition(initial_condition);
        // solve
        std::cout<<"Solve"<<std::endl;
        //solver.SolveAndWriteResultsToFile();
        solver6.SolveAndWriteResultsToFile();
        std::cout<<"Clean"<<std::endl;

        // clean
        PetscTools::Destroy(initial_condition);
        */
        
    }

    void TestChemicalOdePdeParameters()
    {
        
        std::cout<<"Mass action kinetics of Schnackenberg using the AbstractChemicalOdeSystem Class"<<std::endl;
        std::cout<<"-----------------------------"<<std::endl;
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

        // r1: 2U + V -> 3U     forwardRate = 0.1
        // r2: 0 <-> U          forwardRate = 0.1 reverseRate = 0.2
        // r3: 0 -> V           forwardRate = 0.3

        double reaction_1_rate = 0.1;
        double reaction_2_forward_rate = 0.1;
        double reaction_2_reverse_rate = 0.2;
        double reaction_3_rate = 0.3;

        AbstractReaction* p_reaction_1 = new AbstractReaction(p_substrates_1, p_products_1, stoich_substrates_1, stoich_products_1,reaction_1_rate);
        AbstractReversibleReaction* p_reaction_2 = new AbstractReversibleReaction(p_substrates_2, p_products_2, stoich_substrates_2, stoich_products_2,reaction_2_forward_rate, reaction_2_reverse_rate);
        AbstractReaction* p_reaction_3 = new AbstractReaction(p_substrates_3, p_products_3, stoich_substrates_3, stoich_products_3,reaction_3_rate);

        std::vector<AbstractReaction*> p_reaction_vector;
        p_reaction_vector.push_back(p_reaction_1);
        p_reaction_vector.push_back(p_reaction_2);
        p_reaction_vector.push_back(p_reaction_3);
    
        //AbstractReactionSystem* p_reaction_system = new AbstractReactionSystem(p_system_chemistry, p_reaction_vector);

        //SchnackenbergCoupledPdeSystem<2> pde(1e-4, 1e-2, 0.5, 2.2, 1.5, 1);
    
        MassActionReaction* p_mass_action_reaction_1 = new MassActionReaction(p_substrates_1, p_products_1, stoich_substrates_1, stoich_products_1,false,false,reaction_1_rate);
        MassActionReaction* p_mass_action_reaction_2 = new MassActionReaction(p_substrates_2, p_products_2, stoich_substrates_2, stoich_products_2,false, true, reaction_2_forward_rate, reaction_2_reverse_rate);
        MassActionReaction* p_mass_action_reaction_3 = new MassActionReaction(p_substrates_3, p_products_3, stoich_substrates_3, stoich_products_3,false, false,reaction_3_rate);

        std::vector<AbstractReaction*> p_mass_action_reaction_vector;
        p_mass_action_reaction_vector.push_back(p_mass_action_reaction_1);
        p_mass_action_reaction_vector.push_back(p_mass_action_reaction_2);
        p_mass_action_reaction_vector.push_back(p_mass_action_reaction_3);

        AbstractReactionSystem* p_mass_action_reaction_system = new AbstractReactionSystem(p_system_chemistry, p_mass_action_reaction_vector);

        // form ode system
        std::cout<<"SchnackenbergCoupledPdeOdeSystem  -As AbstractChemicalODeSystem"<<std::endl;
        // system properties
        const unsigned probDim =2;
        const unsigned spaceDim=2;
        const unsigned elementDim=2;
        std::vector<double> initValues = {2.0, 0.75};
        std::vector<double> bcValues = {0.0, 0.0};

        // mesh
        HoneycombMeshGenerator generator(10, 10, 0);
        MutableMesh<elementDim,spaceDim>* p_mesh = generator.GetMesh();

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
        std::cout<<"Initial conditions"<<std::endl;
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
        std::cout<<"PETSc Vec"<<std::endl;
        Vec initial_condition = PetscTools::CreateVec(init_conds);

        
        //SchnackenbergCoupledPdeSystem<2> pde(1e-4, 1e-2, 0.5, 2.2, 1.5, 1);
        // coupled ode system
        std::cout<<"Ode loop"<<std::endl;
        std::vector<AbstractOdeSystemForCoupledPdeSystem*> odeSystem;
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++){
            // number of ode system objects must match the number of nodes, i.e the individual odes may be multi-dimensional
            odeSystem.push_back(new AbstractChemicalOdeForCoupledPdeSystem(p_mass_action_reaction_system));
        }
        std::cout<<"Number odesAtNodes: "<<p_mesh->GetNumNodes()<<std::endl;

        // pde system
        std::cout<<"Pde"<<std::endl;
        PdeSchnackenbergCoupledPdeOdeSystem<elementDim, spaceDim, probDim> pde(odeSystem[0],1e-4, 1e-2 );

        // used the explicitly defined EulerIvpOdeSolver rather than the default defined BackwardEuler method

        //EulerIvpOdeSolver euler_solver; 
        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);

        std::cout<<"Solver"<<std::endl;
        // solver
        LinearParabolicPdeSystemWithCoupledOdeSystemSolver<elementDim,spaceDim,probDim> solver(p_mesh, &pde, &bcc,odeSystem,p_solver);

        // solver properties
        double t_end = 10;
        solver.SetTimes(0, t_end);
        solver.SetTimeStep(1e-2);
        solver.SetSamplingTimeStep(1e-2);
        solver.SetOutputDirectory("TestAbstractChemicalOdeOutput_groupMeeting_params1");
        solver.SetInitialCondition(initial_condition);
        // solve
        std::cout<<"Solve"<<std::endl;
        //solver.SolveAndWriteResultsToFile();
        solver.SolveAndWriteResultsToFile();
        std::cout<<"Clean"<<std::endl;

        // clean
        PetscTools::Destroy(initial_condition);
        
    }

    void TestChemicalOdePdeParameters2()
    {
        
        std::cout<<"Mass action kinetics of Schnackenberg using the AbstractChemicalOdeSystem Class"<<std::endl;
        std::cout<<"-----------------------------"<<std::endl;
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

        // r1: 2U + V -> 3U     forwardRate = 0.1
        // r2: 0 <-> U          forwardRate = 0.1 reverseRate = 0.2
        // r3: 0 -> V           forwardRate = 0.3

        double reaction_1_rate = 0.2;
        double reaction_2_forward_rate = 0.2;
        double reaction_2_reverse_rate = 0.4;
        double reaction_3_rate = 0.6;

        AbstractReaction* p_reaction_1 = new AbstractReaction(p_substrates_1, p_products_1, stoich_substrates_1, stoich_products_1,reaction_1_rate);
        AbstractReversibleReaction* p_reaction_2 = new AbstractReversibleReaction(p_substrates_2, p_products_2, stoich_substrates_2, stoich_products_2,reaction_2_forward_rate, reaction_2_reverse_rate);
        AbstractReaction* p_reaction_3 = new AbstractReaction(p_substrates_3, p_products_3, stoich_substrates_3, stoich_products_3,reaction_3_rate);

        std::vector<AbstractReaction*> p_reaction_vector;
        p_reaction_vector.push_back(p_reaction_1);
        p_reaction_vector.push_back(p_reaction_2);
        p_reaction_vector.push_back(p_reaction_3);
    
        //AbstractReactionSystem* p_reaction_system = new AbstractReactionSystem(p_system_chemistry, p_reaction_vector);

        //SchnackenbergCoupledPdeSystem<2> pde(1e-4, 1e-2, 0.5, 2.2, 1.5, 1);
    
        MassActionReaction* p_mass_action_reaction_1 = new MassActionReaction(p_substrates_1, p_products_1, stoich_substrates_1, stoich_products_1,false,false,reaction_1_rate);
        MassActionReaction* p_mass_action_reaction_2 = new MassActionReaction(p_substrates_2, p_products_2, stoich_substrates_2, stoich_products_2,false, true, reaction_2_forward_rate, reaction_2_reverse_rate);
        MassActionReaction* p_mass_action_reaction_3 = new MassActionReaction(p_substrates_3, p_products_3, stoich_substrates_3, stoich_products_3,false, false,reaction_3_rate);

        std::vector<AbstractReaction*> p_mass_action_reaction_vector;
        p_mass_action_reaction_vector.push_back(p_mass_action_reaction_1);
        p_mass_action_reaction_vector.push_back(p_mass_action_reaction_2);
        p_mass_action_reaction_vector.push_back(p_mass_action_reaction_3);

        AbstractReactionSystem* p_mass_action_reaction_system = new AbstractReactionSystem(p_system_chemistry, p_mass_action_reaction_vector);

        // form ode system
        std::cout<<"SchnackenbergCoupledPdeOdeSystem  -As AbstractChemicalODeSystem"<<std::endl;
        // system properties
        const unsigned probDim =2;
        const unsigned spaceDim=2;
        const unsigned elementDim=2;
        std::vector<double> initValues = {2.0, 0.75};
        std::vector<double> bcValues = {0.0, 0.0};

        // mesh
        HoneycombMeshGenerator generator(10, 10, 0);
        MutableMesh<elementDim,spaceDim>* p_mesh = generator.GetMesh();

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
        std::cout<<"Initial conditions"<<std::endl;
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
        std::cout<<"PETSc Vec"<<std::endl;
        Vec initial_condition = PetscTools::CreateVec(init_conds);

        
        //SchnackenbergCoupledPdeSystem<2> pde(1e-4, 1e-2, 0.5, 2.2, 1.5, 1);
        // coupled ode system
        std::cout<<"Ode loop"<<std::endl;
        std::vector<AbstractOdeSystemForCoupledPdeSystem*> odeSystem;
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++){
            // number of ode system objects must match the number of nodes, i.e the individual odes may be multi-dimensional
            odeSystem.push_back(new AbstractChemicalOdeForCoupledPdeSystem(p_mass_action_reaction_system));
        }
        std::cout<<"Number odesAtNodes: "<<p_mesh->GetNumNodes()<<std::endl;

        // pde system
        std::cout<<"Pde"<<std::endl;
        PdeSchnackenbergCoupledPdeOdeSystem<elementDim, spaceDim, probDim> pde(odeSystem[0],1e-4, 1e-2 );

        // used the explicitly defined EulerIvpOdeSolver rather than the default defined BackwardEuler method

        //EulerIvpOdeSolver euler_solver; 
        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);

        std::cout<<"Solver"<<std::endl;
        // solver
        LinearParabolicPdeSystemWithCoupledOdeSystemSolver<elementDim,spaceDim,probDim> solver(p_mesh, &pde, &bcc,odeSystem,p_solver);

        // solver properties
        double t_end = 10;
        solver.SetTimes(0, t_end);
        solver.SetTimeStep(1e-2);
        solver.SetSamplingTimeStep(1e-2);
        solver.SetOutputDirectory("TestAbstractChemicalOdeOutput_groupMeeting_params2");
        solver.SetInitialCondition(initial_condition);
        // solve
        std::cout<<"Solve"<<std::endl;
        //solver.SolveAndWriteResultsToFile();
        solver.SolveAndWriteResultsToFile();
        std::cout<<"Clean"<<std::endl;

        // clean
        PetscTools::Destroy(initial_condition);
        
    }


    void TestChemicalOdePdeParameters3()
    {
        
        std::cout<<"Mass action kinetics of Schnackenberg using the AbstractChemicalOdeSystem Class"<<std::endl;
        std::cout<<"-----------------------------"<<std::endl;
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

        // r1: 2U + V -> 3U     forwardRate = 0.1
        // r2: 0 <-> U          forwardRate = 0.1 reverseRate = 0.2
        // r3: 0 -> V           forwardRate = 0.3

        double reaction_1_rate = 0.3;
        double reaction_2_forward_rate = 0.3;
        double reaction_2_reverse_rate = 0.6;
        double reaction_3_rate = 0.9;

        AbstractReaction* p_reaction_1 = new AbstractReaction(p_substrates_1, p_products_1, stoich_substrates_1, stoich_products_1,reaction_1_rate);
        AbstractReversibleReaction* p_reaction_2 = new AbstractReversibleReaction(p_substrates_2, p_products_2, stoich_substrates_2, stoich_products_2,reaction_2_forward_rate, reaction_2_reverse_rate);
        AbstractReaction* p_reaction_3 = new AbstractReaction(p_substrates_3, p_products_3, stoich_substrates_3, stoich_products_3,reaction_3_rate);

        std::vector<AbstractReaction*> p_reaction_vector;
        p_reaction_vector.push_back(p_reaction_1);
        p_reaction_vector.push_back(p_reaction_2);
        p_reaction_vector.push_back(p_reaction_3);
    
        //AbstractReactionSystem* p_reaction_system = new AbstractReactionSystem(p_system_chemistry, p_reaction_vector);

        //SchnackenbergCoupledPdeSystem<2> pde(1e-4, 1e-2, 0.5, 2.2, 1.5, 1);
    
        MassActionReaction* p_mass_action_reaction_1 = new MassActionReaction(p_substrates_1, p_products_1, stoich_substrates_1, stoich_products_1,false,false,reaction_1_rate);
        MassActionReaction* p_mass_action_reaction_2 = new MassActionReaction(p_substrates_2, p_products_2, stoich_substrates_2, stoich_products_2,false, true, reaction_2_forward_rate, reaction_2_reverse_rate);
        MassActionReaction* p_mass_action_reaction_3 = new MassActionReaction(p_substrates_3, p_products_3, stoich_substrates_3, stoich_products_3,false, false,reaction_3_rate);

        std::vector<AbstractReaction*> p_mass_action_reaction_vector;
        p_mass_action_reaction_vector.push_back(p_mass_action_reaction_1);
        p_mass_action_reaction_vector.push_back(p_mass_action_reaction_2);
        p_mass_action_reaction_vector.push_back(p_mass_action_reaction_3);

        AbstractReactionSystem* p_mass_action_reaction_system = new AbstractReactionSystem(p_system_chemistry, p_mass_action_reaction_vector);

        // form ode system
        std::cout<<"SchnackenbergCoupledPdeOdeSystem  -As AbstractChemicalODeSystem"<<std::endl;
        // system properties
        const unsigned probDim =2;
        const unsigned spaceDim=2;
        const unsigned elementDim=2;
        std::vector<double> initValues = {2.0, 0.75};
        std::vector<double> bcValues = {0.0, 0.0};

        // mesh
        HoneycombMeshGenerator generator(10, 10, 0);
        MutableMesh<elementDim,spaceDim>* p_mesh = generator.GetMesh();

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
        std::cout<<"Initial conditions"<<std::endl;
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
        std::cout<<"PETSc Vec"<<std::endl;
        Vec initial_condition = PetscTools::CreateVec(init_conds);

        
        //SchnackenbergCoupledPdeSystem<2> pde(1e-4, 1e-2, 0.5, 2.2, 1.5, 1);
        // coupled ode system
        std::cout<<"Ode loop"<<std::endl;
        std::vector<AbstractOdeSystemForCoupledPdeSystem*> odeSystem;
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++){
            // number of ode system objects must match the number of nodes, i.e the individual odes may be multi-dimensional
            odeSystem.push_back(new AbstractChemicalOdeForCoupledPdeSystem(p_mass_action_reaction_system));
        }
        std::cout<<"Number odesAtNodes: "<<p_mesh->GetNumNodes()<<std::endl;

        // pde system
        std::cout<<"Pde"<<std::endl;
        PdeSchnackenbergCoupledPdeOdeSystem<elementDim, spaceDim, probDim> pde(odeSystem[0],1e-4, 1e-2 );

        // used the explicitly defined EulerIvpOdeSolver rather than the default defined BackwardEuler method

        //EulerIvpOdeSolver euler_solver; 
        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);

        std::cout<<"Solver"<<std::endl;
        // solver
        LinearParabolicPdeSystemWithCoupledOdeSystemSolver<elementDim,spaceDim,probDim> solver(p_mesh, &pde, &bcc,odeSystem,p_solver);

        // solver properties
        double t_end = 10;
        solver.SetTimes(0, t_end);
        solver.SetTimeStep(1e-2);
        solver.SetSamplingTimeStep(1e-2);
        solver.SetOutputDirectory("TestAbstractChemicalOdeOutput_groupMeeting_params3");
        solver.SetInitialCondition(initial_condition);
        // solve
        std::cout<<"Solve"<<std::endl;
        //solver.SolveAndWriteResultsToFile();
        solver.SolveAndWriteResultsToFile();
        std::cout<<"Clean"<<std::endl;

        // clean
        PetscTools::Destroy(initial_condition);
        
    }


    void TestChemicalOdePdeDiffusion4()
    {
        
        std::cout<<"Mass action kinetics of Schnackenberg using the AbstractChemicalOdeSystem Class"<<std::endl;
        std::cout<<"-----------------------------"<<std::endl;
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

        // r1: 2U + V -> 3U     forwardRate = 0.1
        // r2: 0 <-> U          forwardRate = 0.1 reverseRate = 0.2
        // r3: 0 -> V           forwardRate = 0.3

        double reaction_1_rate = 0.3;
        double reaction_2_forward_rate = 0.3;
        double reaction_2_reverse_rate = 0.6;
        double reaction_3_rate = 0.9;

        AbstractReaction* p_reaction_1 = new AbstractReaction(p_substrates_1, p_products_1, stoich_substrates_1, stoich_products_1,reaction_1_rate);
        AbstractReversibleReaction* p_reaction_2 = new AbstractReversibleReaction(p_substrates_2, p_products_2, stoich_substrates_2, stoich_products_2,reaction_2_forward_rate, reaction_2_reverse_rate);
        AbstractReaction* p_reaction_3 = new AbstractReaction(p_substrates_3, p_products_3, stoich_substrates_3, stoich_products_3,reaction_3_rate);

        std::vector<AbstractReaction*> p_reaction_vector;
        p_reaction_vector.push_back(p_reaction_1);
        p_reaction_vector.push_back(p_reaction_2);
        p_reaction_vector.push_back(p_reaction_3);
    
        //AbstractReactionSystem* p_reaction_system = new AbstractReactionSystem(p_system_chemistry, p_reaction_vector);

        //SchnackenbergCoupledPdeSystem<2> pde(1e-4, 1e-2, 0.5, 2.2, 1.5, 1);
    
        MassActionReaction* p_mass_action_reaction_1 = new MassActionReaction(p_substrates_1, p_products_1, stoich_substrates_1, stoich_products_1,false,false,reaction_1_rate);
        MassActionReaction* p_mass_action_reaction_2 = new MassActionReaction(p_substrates_2, p_products_2, stoich_substrates_2, stoich_products_2,false, true, reaction_2_forward_rate, reaction_2_reverse_rate);
        MassActionReaction* p_mass_action_reaction_3 = new MassActionReaction(p_substrates_3, p_products_3, stoich_substrates_3, stoich_products_3,false, false,reaction_3_rate);

        std::vector<AbstractReaction*> p_mass_action_reaction_vector;
        p_mass_action_reaction_vector.push_back(p_mass_action_reaction_1);
        p_mass_action_reaction_vector.push_back(p_mass_action_reaction_2);
        p_mass_action_reaction_vector.push_back(p_mass_action_reaction_3);

        AbstractReactionSystem* p_mass_action_reaction_system = new AbstractReactionSystem(p_system_chemistry, p_mass_action_reaction_vector);

        // form ode system
        std::cout<<"SchnackenbergCoupledPdeOdeSystem  -As AbstractChemicalODeSystem"<<std::endl;
        // system properties
        const unsigned probDim =2;
        const unsigned spaceDim=2;
        const unsigned elementDim=2;
        std::vector<double> initValues = {2.0, 0.75};
        std::vector<double> bcValues = {0.0, 0.0};

        // mesh
        HoneycombMeshGenerator generator(10, 10, 0);
        MutableMesh<elementDim,spaceDim>* p_mesh = generator.GetMesh();

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
        std::cout<<"Initial conditions"<<std::endl;
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
        std::cout<<"PETSc Vec"<<std::endl;
        Vec initial_condition = PetscTools::CreateVec(init_conds);

        
        //SchnackenbergCoupledPdeSystem<2> pde(1e-4, 1e-2, 0.5, 2.2, 1.5, 1);
        // coupled ode system
        std::cout<<"Ode loop"<<std::endl;
        std::vector<AbstractOdeSystemForCoupledPdeSystem*> odeSystem;
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++){
            // number of ode system objects must match the number of nodes, i.e the individual odes may be multi-dimensional
            odeSystem.push_back(new AbstractChemicalOdeForCoupledPdeSystem(p_mass_action_reaction_system));
        }
        std::cout<<"Number odesAtNodes: "<<p_mesh->GetNumNodes()<<std::endl;

        // pde system
        std::cout<<"Pde"<<std::endl;
        PdeSchnackenbergCoupledPdeOdeSystem<elementDim, spaceDim, probDim> pde(odeSystem[0],1e-3, 1e-1 );

        // used the explicitly defined EulerIvpOdeSolver rather than the default defined BackwardEuler method

        //EulerIvpOdeSolver euler_solver; 
        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);

        std::cout<<"Solver"<<std::endl;
        // solver
        LinearParabolicPdeSystemWithCoupledOdeSystemSolver<elementDim,spaceDim,probDim> solver(p_mesh, &pde, &bcc,odeSystem,p_solver);

        // solver properties
        double t_end = 10;
        solver.SetTimes(0, t_end);
        solver.SetTimeStep(1e-2);
        solver.SetSamplingTimeStep(1e-2);
        solver.SetOutputDirectory("TestAbstractChemicalOdeOutput_groupMeeting_diffusion1e-31e-1 ");
        solver.SetInitialCondition(initial_condition);
        // solve
        std::cout<<"Solve"<<std::endl;
        //solver.SolveAndWriteResultsToFile();
        solver.SolveAndWriteResultsToFile();
        std::cout<<"Clean"<<std::endl;

        // clean
        PetscTools::Destroy(initial_condition);
        
    }

    void TestChemicalOdePdeDiffusion3()
    {
        
        std::cout<<"Mass action kinetics of Schnackenberg using the AbstractChemicalOdeSystem Class"<<std::endl;
        std::cout<<"-----------------------------"<<std::endl;
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

        // r1: 2U + V -> 3U     forwardRate = 0.1
        // r2: 0 <-> U          forwardRate = 0.1 reverseRate = 0.2
        // r3: 0 -> V           forwardRate = 0.3

        double reaction_1_rate = 0.3;
        double reaction_2_forward_rate = 0.3;
        double reaction_2_reverse_rate = 0.6;
        double reaction_3_rate = 0.9;

        AbstractReaction* p_reaction_1 = new AbstractReaction(p_substrates_1, p_products_1, stoich_substrates_1, stoich_products_1,reaction_1_rate);
        AbstractReversibleReaction* p_reaction_2 = new AbstractReversibleReaction(p_substrates_2, p_products_2, stoich_substrates_2, stoich_products_2,reaction_2_forward_rate, reaction_2_reverse_rate);
        AbstractReaction* p_reaction_3 = new AbstractReaction(p_substrates_3, p_products_3, stoich_substrates_3, stoich_products_3,reaction_3_rate);

        std::vector<AbstractReaction*> p_reaction_vector;
        p_reaction_vector.push_back(p_reaction_1);
        p_reaction_vector.push_back(p_reaction_2);
        p_reaction_vector.push_back(p_reaction_3);
    
        //AbstractReactionSystem* p_reaction_system = new AbstractReactionSystem(p_system_chemistry, p_reaction_vector);

        //SchnackenbergCoupledPdeSystem<2> pde(1e-4, 1e-2, 0.5, 2.2, 1.5, 1);
    
        MassActionReaction* p_mass_action_reaction_1 = new MassActionReaction(p_substrates_1, p_products_1, stoich_substrates_1, stoich_products_1,false,false,reaction_1_rate);
        MassActionReaction* p_mass_action_reaction_2 = new MassActionReaction(p_substrates_2, p_products_2, stoich_substrates_2, stoich_products_2,false, true, reaction_2_forward_rate, reaction_2_reverse_rate);
        MassActionReaction* p_mass_action_reaction_3 = new MassActionReaction(p_substrates_3, p_products_3, stoich_substrates_3, stoich_products_3,false, false,reaction_3_rate);

        std::vector<AbstractReaction*> p_mass_action_reaction_vector;
        p_mass_action_reaction_vector.push_back(p_mass_action_reaction_1);
        p_mass_action_reaction_vector.push_back(p_mass_action_reaction_2);
        p_mass_action_reaction_vector.push_back(p_mass_action_reaction_3);

        AbstractReactionSystem* p_mass_action_reaction_system = new AbstractReactionSystem(p_system_chemistry, p_mass_action_reaction_vector);

        // form ode system
        std::cout<<"SchnackenbergCoupledPdeOdeSystem  -As AbstractChemicalODeSystem"<<std::endl;
        // system properties
        const unsigned probDim =2;
        const unsigned spaceDim=2;
        const unsigned elementDim=2;
        std::vector<double> initValues = {2.0, 0.75};
        std::vector<double> bcValues = {0.0, 0.0};

        // mesh
        HoneycombMeshGenerator generator(10, 10, 0);
        MutableMesh<elementDim,spaceDim>* p_mesh = generator.GetMesh();

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
        std::cout<<"Initial conditions"<<std::endl;
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
        std::cout<<"PETSc Vec"<<std::endl;
        Vec initial_condition = PetscTools::CreateVec(init_conds);

        
        //SchnackenbergCoupledPdeSystem<2> pde(1e-4, 1e-2, 0.5, 2.2, 1.5, 1);
        // coupled ode system
        std::cout<<"Ode loop"<<std::endl;
        std::vector<AbstractOdeSystemForCoupledPdeSystem*> odeSystem;
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++){
            // number of ode system objects must match the number of nodes, i.e the individual odes may be multi-dimensional
            odeSystem.push_back(new AbstractChemicalOdeForCoupledPdeSystem(p_mass_action_reaction_system));
        }
        std::cout<<"Number odesAtNodes: "<<p_mesh->GetNumNodes()<<std::endl;

        // pde system
        std::cout<<"Pde"<<std::endl;
        PdeSchnackenbergCoupledPdeOdeSystem<elementDim, spaceDim, probDim> pde(odeSystem[0],1e-2, 1 );

        // used the explicitly defined EulerIvpOdeSolver rather than the default defined BackwardEuler method

        //EulerIvpOdeSolver euler_solver; 
        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);

        std::cout<<"Solver"<<std::endl;
        // solver
        LinearParabolicPdeSystemWithCoupledOdeSystemSolver<elementDim,spaceDim,probDim> solver(p_mesh, &pde, &bcc,odeSystem,p_solver);

        // solver properties
        double t_end = 10;
        solver.SetTimes(0, t_end);
        solver.SetTimeStep(1e-2);
        solver.SetSamplingTimeStep(1e-2);
        solver.SetOutputDirectory("TestAbstractChemicalOdeOutput_groupMeeting_diffusion1e-21e-0 ");
        solver.SetInitialCondition(initial_condition);
        // solve
        std::cout<<"Solve"<<std::endl;
        //solver.SolveAndWriteResultsToFile();
        solver.SolveAndWriteResultsToFile();
        std::cout<<"Clean"<<std::endl;

        // clean
        PetscTools::Destroy(initial_condition);
        
    }


};


#endif