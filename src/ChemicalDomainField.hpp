#ifndef CHEMICALDOMAINFIELD_HPP
#define CHEMICALDOMAINFIELD_HPP

#include <string>
#include <vector>
#include "AbstractDomainField.hpp"
#include "AbstractDiffusiveChemistry.hpp"
#include "AbstractDiffusiveChemical.hpp"
#include "AbstractInhomogenousChemicalOdeSystemForCoupledPdeSystem.hpp"
#include "AbstractReactionSystemFromFile.hpp"


class ChemicalDomainField: public AbstractDomainField
{
private:
    using AbstractDomainField::GetDiffusionValueBasedOnPoint;
    using AbstractDomainField::ReturnDiffusionValueFromStateNameAndDomainLabel;
    using AbstractDomainField::GetFieldType;

protected:

    std::string mReactionFileRoot;

    std::vector<AbstractInhomogenousChemicalOdeSystemForCoupledPdeSystem*> mOdeSystem;

    StateVariableRegister* mpDomainRegister; // the state variable vector for the union of reacting species

    std::vector<std::string> mReactingStateVariableVector;

    AbstractDiffusiveChemistry* mpDiffusiveChemistry;

    unsigned mProbDim =0;
    const unsigned mSpaceDim=2;
    const unsigned mElementDim=2;


public:

    ChemicalDomainField(std::string reactionFileRoot="",
                        std::string domainLabelFilename="", 
                        std::string domainKeyFilename="", 
                        std::string odeLabelFilename="", 
                        std::string odeKeyFilename="",
                        std::string diffusionFilename="");

    virtual~ChemicalDomainField()
    {
    }

    // setup methods

    void FormReactionSystemAtNodes();

    void DeriveSystemProperties();

    void DeriveExtendedSystemProperties();

    // interface methods

    virtual double GetDiffusionValueBasedOnPoint(const ChastePoint<2>&, unsigned);

    virtual double ReturnDiffusionValueFromStateNameAndDomainLabel(std::string stateName, std::string domainLabel = "");

    virtual std::string GetFieldType()
    {
        return "ChemicalDomainField";
    };

    // set methods
    void SetDomainStateVariableRegister(StateVariableRegister*);

    void SetChemistry(AbstractDiffusiveChemistry*);

    void SetOdeSystem(std::vector<AbstractInhomogenousChemicalOdeSystemForCoupledPdeSystem*>);

    // get methods
    StateVariableRegister* GetDomainStateVariableRegister();

    unsigned GetProblemDimensions();

    unsigned GetSpaceDimensions();

    unsigned GetElementDimensions();

    AbstractDiffusiveChemistry* GetChemistry();

    std::vector<AbstractInhomogenousChemicalOdeSystemForCoupledPdeSystem*> GetOdeSystem();

};

ChemicalDomainField::ChemicalDomainField(   std::string reactionFileRoot,
                                            std::string domainLabelFilename, 
                                            std::string domainKeyFilename, 
                                            std::string odeLabelFilename, 
                                            std::string odeKeyFilename, 
                                            std::string diffusionFilename)
    :   AbstractDomainField(domainLabelFilename, 
                            domainKeyFilename, 
                            odeLabelFilename, 
                            odeKeyFilename, 
                            diffusionFilename),
        mReactionFileRoot(reactionFileRoot)
{

    // set up the ode system container
    std::vector<AbstractInhomogenousChemicalOdeSystemForCoupledPdeSystem*> mOdeSystem = std::vector<AbstractInhomogenousChemicalOdeSystemForCoupledPdeSystem*>();

    // map odes systems to nodes 
    FormReactionSystemAtNodes();

    // update properties of the domain field
    //DeriveSystemProperties();
    DeriveExtendedSystemProperties();
}

void ChemicalDomainField::FormReactionSystemAtNodes()
{   
    // run through each node of the mesh formed from the input file and formed from the base class
    // imbue the node with the ode system 

    // run the system update for the first node separately to the rest in order to set up the correct domain state variable register

    std::vector<std::string> nodeLabels = GetNodeLabels();
    
    // populate the node ODE, system chemistry and state variable registers for the first
    std::string reactionFilename = mReactionFileRoot + ReturnKeyValueFromNodeLabel(nodeLabels[0]);

    // pointer to read a the ODE domain specified reaction system
    AbstractReactionSystemFromFile* p_file_reaction_system = new  AbstractReactionSystemFromFile(reactionFilename);
    // convert read in reaction system to a chemical ODE system
    AbstractInhomogenousChemicalOdeSystemForCoupledPdeSystem* p_chemical_ode = new AbstractInhomogenousChemicalOdeSystemForCoupledPdeSystem(p_file_reaction_system);

    // form the state variable register for the node from the read in reaction systems chemistry
    std::vector<std::string> chemical_names = p_file_reaction_system -> GetSystemChemistry() -> GetChemicalNames();

    p_chemical_ode -> SetStateVariableRegister(new StateVariableRegister(chemical_names));

    mOdeSystem.push_back(p_chemical_ode);

    // domain state vector
    StateVariableRegister* p_domain_register  = new StateVariableRegister(chemical_names); 

    // populate the rest of the nodes
    for(unsigned node_index =1; node_index<nodeLabels.size(); node_index++)
    {

        std::string reactionFilename = mReactionFileRoot + ReturnKeyValueFromNodeLabel(nodeLabels[node_index]);

        // pointer to read a the ODE domain specified reaction system
        AbstractReactionSystemFromFile* p_file_reaction_system = new  AbstractReactionSystemFromFile(reactionFilename);

        // convert read in reaction system to a chemical ODE system
        AbstractInhomogenousChemicalOdeSystemForCoupledPdeSystem* p_chemical_ode = new AbstractInhomogenousChemicalOdeSystemForCoupledPdeSystem(p_file_reaction_system);

        // form the state variable register for the node from the read in reaction systems chemistry
        std::vector<std::string> chemical_names = p_file_reaction_system -> GetSystemChemistry() -> GetChemicalNames();

        p_chemical_ode -> SetStateVariableRegister(new StateVariableRegister(chemical_names));

        mOdeSystem.push_back(p_chemical_ode);

        p_domain_register -> AddStateVariableVector(chemical_names);

    }

    SetDomainStateVariableRegister(p_domain_register);
    
    // run through the nodes and set the complete domain state variable register so that
    // each ode has knowledge of the system size

    // edit: this was setting the p_domain_register but as the diffusion register is larger and contains
    // more species which diffuse via pde than chemically react, need to input the diffusion state register

    for(unsigned node_index =0; node_index<nodeLabels.size(); node_index++)
    {
        mOdeSystem[node_index] -> SetDomainStateVariableRegister(mStateVariableVector);
    }
   
}


void ChemicalDomainField::DeriveSystemProperties()
{
    // for when we want to track the diffusion of all speces in the system
    // this function derives the state register from the chemical reaction system
    // leads to incorrect initial conditions if the initial conditions take the form of all the diffusing 
    // species in the system

    
    // determine the system properties from the read in reaction systems
    mProbDim = mpDomainRegister -> GetNumberOfStateVariables();
    
    AbstractDiffusiveChemistry* p_diffusive_chemistry = new AbstractDiffusiveChemistry();

    // create a diffusive chemistry from the diffusion database

    for(unsigned system_number=0; system_number<mOdeSystem.size(); system_number++)
    {
 
        std::vector<AbstractChemical*> ode_chemical_vector = mOdeSystem[system_number] -> GetReactionSystem() -> GetSystemChemistry() -> rGetChemicalVector();
    
        // for each chemical in the reaction system 
        for(unsigned ode_system_chemical_index=0; ode_system_chemical_index<ode_chemical_vector.size(); ode_system_chemical_index++)
        {
            // populate their diffusive properties from the labels in the system and the diffusion database
            
            for(unsigned domain_index=0; domain_index<GetNumberOfDomains(); domain_index++)
            {
                std::string domainLabel = GetDomainLabelByIndex(domain_index);

                std::string chemicalName = ode_chemical_vector[ode_system_chemical_index] -> GetChemicalName();
                double chemicalSize = ode_chemical_vector[ode_system_chemical_index] -> GetChemicalSize();
                double chemicalMass = ode_chemical_vector[ode_system_chemical_index] -> GetChemicalMass();
                int chemicalValence = ode_chemical_vector[ode_system_chemical_index] -> GetChemicalValence();

                // diffusion database may be larger than the active chemicals in the system, the chemical vector
                bool IsInDatabase = false;

                // for a diffusive chemical object through the information given in the chemical
                // change the pointer of the original chemical to point to the diffusive chemical
                AbstractDiffusiveChemical *p_chemical = new AbstractDiffusiveChemical(chemicalName,chemicalSize,chemicalMass,chemicalValence);
             
                for(unsigned record_index=0; record_index<mDiffusionDatabase.size(); record_index++)
                {
                    if(mDiffusionDatabase[record_index][0] == chemicalName)
                    {   
                        p_chemical -> AddDiffusiveDomain(mDiffusionDatabase[record_index][1],std::stod(mDiffusionDatabase[record_index][2]));
                        IsInDatabase = true;
                    }
                    // run through the rest of the database incase another domain type is specified
                }
                // if not in the database then don't add to diffusive chemistry as cannot diffuse
                if(IsInDatabase)
                {
                    p_diffusive_chemistry -> AddChemical(p_chemical);
                }
                ode_chemical_vector[ode_system_chemical_index] = p_chemical;
            }
            
        }
    }

    SetChemistry(p_diffusive_chemistry);
}

void ChemicalDomainField::DeriveExtendedSystemProperties()
{
    // determine the system properties from the read in reaction systems
    mProbDim = GetStateVariableVector() -> GetNumberOfStateVariables();

    AbstractDiffusiveChemistry* p_diffusive_chemistry = new AbstractDiffusiveChemistry();

    // each chemical in the system deined in state variable vector, derived from diffusion database
    StateVariableRegister* p_stateVariableVector = GetStateVariableVector();

    for(unsigned chem_index=0; chem_index < p_stateVariableVector -> GetNumberOfStateVariables(); chem_index++)
    {
        // for each chemical in the diffusion database
        std::string chemicalName = p_stateVariableVector -> RetrieveStateVariableName(chem_index);

        // for each domain
        for(unsigned domain_index=0; domain_index<GetNumberOfDomains(); domain_index++)
        {
            std::string domainLabel = GetDomainLabelByIndex(domain_index);

            AbstractDiffusiveChemical *p_chemical = new AbstractDiffusiveChemical(chemicalName);
            bool IsInDatabase = false;
            // test whether the chemical and domain match, expecting at least one match
            for(unsigned record_index=0; record_index<mDiffusionDatabase.size(); record_index++)
            {
                if(mDiffusionDatabase[record_index][0] == chemicalName)
                {   
                    p_chemical -> AddDiffusiveDomain(mDiffusionDatabase[record_index][1],std::stod(mDiffusionDatabase[record_index][2]));
                    IsInDatabase = true;
                }
                // run through the rest of the database incase another domain type is specified
            }
            if(IsInDatabase)
            {
                p_diffusive_chemistry -> AddChemical(p_chemical);
            }

        SetChemistry(p_diffusive_chemistry);
        }
    }
}




double ChemicalDomainField::GetDiffusionValueBasedOnPoint(const ChastePoint<2>& chastePoint, unsigned stateIndex)
{

    // take in the point and state index then determine the state name and domain label, then determine the diffusion value
    if(stateIndex<mStateVariableVector ->GetNumberOfStateVariables())
    {
        
        // retrive state name from stateVariableRegister
        std::string stateName = mStateVariableVector -> RetrieveStateVariableName(stateIndex);

        // retrieve label from domainLabels
        std::string domainLabel = ReturnDomainLabelAtPosition(chastePoint.rGetLocation());
        std::string domainKeyName = ReturnDomainKeyFromDomainLabel(domainLabel);

        return ReturnDiffusionValueFromStateNameAndDomainLabel(stateName,domainKeyName);
    
    }
    else
    {
        std::cout<<"Error: ChemicalDomainField::GetDiffusionValueBasedOnPoint: State not in state variable"<<std::endl;
        return 0.0;
    }

}

double ChemicalDomainField::ReturnDiffusionValueFromStateNameAndDomainLabel(std::string stateName, std::string domainLabel)
{
    // look for diffusion value (from diffusiondatabase) based on state name and domain label
    if(domainLabel != "")
    {

        return mpDiffusiveChemistry -> GetDiffusivityValueByChemicalAndDomainName(stateName,domainLabel);
    }
    else
    {
        //diffusivity is non-domain specific
        return mpDiffusiveChemistry -> GetDiffusivityValueByChemicalName(stateName);
    }
}




// set methods

void ChemicalDomainField::SetDomainStateVariableRegister(StateVariableRegister* p_register)
{
    mpDomainRegister = p_register;
}

void ChemicalDomainField::SetChemistry(AbstractDiffusiveChemistry* p_chemistry)
{
    mpDiffusiveChemistry = p_chemistry;
}

void ChemicalDomainField::SetOdeSystem(std::vector<AbstractInhomogenousChemicalOdeSystemForCoupledPdeSystem*> odeSystem)
{
    mOdeSystem = odeSystem;
}

// get methods
StateVariableRegister* ChemicalDomainField::GetDomainStateVariableRegister()
{
    return mpDomainRegister;
}

unsigned ChemicalDomainField::GetProblemDimensions()
{
    return mProbDim;
}

unsigned ChemicalDomainField::GetSpaceDimensions()
{
    return mSpaceDim;
}

unsigned ChemicalDomainField::GetElementDimensions()
{
    return mElementDim;
}

AbstractDiffusiveChemistry* ChemicalDomainField::GetChemistry()
{
    return mpDiffusiveChemistry;
}

std::vector<AbstractInhomogenousChemicalOdeSystemForCoupledPdeSystem*> ChemicalDomainField::GetOdeSystem()
{
    return mOdeSystem;
}


#endif