#ifndef CHEMICALDOMAINFIELD_TEMPLATED_HPP
#define CHEMICALDOMAINFIELD_TEMPLATED_HPP

#include <string>
#include <vector>
#include "AbstractDomainField_templated.hpp"
#include "AbstractDiffusiveChemistry.hpp"
#include "AbstractDiffusiveChemical.hpp"
#include "AbstractInhomogenousChemicalOdeSystemForCoupledPdeSystem.hpp"
#include "AbstractReactionSystemFromFile.hpp"
#include "ChastePoint.hpp"

// class to handle the chemistry proeprties for a domain field additional to the AbstractDomainField class

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
class ChemicalDomainFieldTemplated : public AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>
{
private:
    using AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::GetDiffusionValueBasedOnPoint;
    using AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::ReturnDiffusionValueFromStateNameAndDomainLabel;
    using AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::GetFieldType;
    using AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetUpDomainFromFiles;

protected:

    std::string mReactionFileRoot;

    std::vector<AbstractInhomogenousChemicalOdeSystemForCoupledPdeSystem*> mOdeSystem;

    StateVariableRegister* mpDomainRegister; // the state variable vector for the union of reacting species

    std::vector<std::string> mReactingStateVariableVector;

    AbstractDiffusiveChemistry* mpDiffusiveChemistry;

public:

    ChemicalDomainFieldTemplated(
                        //TetrahedralMesh<ELEMENT_DIM,SPACE_DIM>*,
                        std::string reactionFileRoot="",
                        std::string domainLabelFilename="", 
                        std::string domainKeyFilename="", 
                        std::string odeLabelFilename="", 
                        std::string odeKeyFilename="",
                        std::string diffusionFilename="",
                        bool isHoneyCombMesh=false,
                        std::vector<double> labelOrigin = std::vector<double>(),
                        std::vector<double> cartesianCellScaleXY = std::vector<double>(),
                        std::vector<double> cartesianOdeScaleXY = std::vector<double>()
                        );

    virtual~ChemicalDomainFieldTemplated()
    {
    }

    // setup methods

    virtual void SetUpDomainFromFiles();

    void FormReactionSystemAtNodes();

    void DeriveSystemProperties();

    void DeriveExtendedSystemProperties();

    // interface methods

    virtual double GetDiffusionValueBasedOnPoint(const ChastePoint<2>&, unsigned);

    virtual double ReturnDiffusionValueFromStateNameAndDomainLabel(std::string stateName, std::string domainLabel = "");

    virtual std::string GetFieldType()
    {
        return "ChemicalDomainFieldTemplated";
    };

    // set methods
    void SetDomainStateVariableRegister(StateVariableRegister*);

    void SetChemistry(AbstractDiffusiveChemistry*);

    void SetOdeSystem(std::vector<AbstractInhomogenousChemicalOdeSystemForCoupledPdeSystem*>);

    // get methods
    StateVariableRegister* GetDomainStateVariableRegister();

    AbstractDiffusiveChemistry* GetChemistry();

    std::vector<AbstractInhomogenousChemicalOdeSystemForCoupledPdeSystem*> GetOdeSystem();

};


template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
ChemicalDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::ChemicalDomainFieldTemplated(
                                            //TetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pFeMesh,    
                                            std::string reactionFileRoot,
                                            std::string domainLabelFilename, 
                                            std::string domainKeyFilename, 
                                            std::string odeLabelFilename, 
                                            std::string odeKeyFilename, 
                                            std::string diffusionFilename,
                                            bool isHoneyCombMesh,
                                            std::vector<double> labelOrigin,
                                            std::vector<double> cartesianCellScale,
                                            std::vector<double> cartesianOdeScale)
    :   AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>(
                            //pFeMesh,
                            domainLabelFilename, 
                            domainKeyFilename, 
                            odeLabelFilename, 
                            odeKeyFilename, 
                            diffusionFilename,
                            isHoneyCombMesh,
                            labelOrigin,
                            cartesianCellScale,
                            cartesianOdeScale),
        mReactionFileRoot(reactionFileRoot)
{
    // set up the ode system container
    std::vector<AbstractInhomogenousChemicalOdeSystemForCoupledPdeSystem*> mOdeSystem = std::vector<AbstractInhomogenousChemicalOdeSystemForCoupledPdeSystem*>();
    SetUpDomainFromFiles();
}


template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
void ChemicalDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetUpDomainFromFiles()
{
    // read and set up the domain as an abstract domain, form mesh etc
    if(!this->mIsFeMeshGenerated)
    {
        AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::GenerateFeMesh();
    }
    
    // populate abstract domain with a chemistry

    // map odes systems to nodes 
    FormReactionSystemAtNodes();

    // update properties of the domain field
    //DeriveSystemProperties();
    DeriveExtendedSystemProperties();
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
void ChemicalDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::FormReactionSystemAtNodes()
{
    // run through each node of the mesh formed from the input file and formed from the base class
    // imbue the node with the ode system 

    // run the system update for the first node separately to the rest in order to set up the correct domain state variable register

    std::vector<std::string> nodeLabels = this->GetNodeLabels();

    // populate the node ODE, system chemistry and state variable registers for the first
    std::string reactionFilename = mReactionFileRoot + this->ReturnKeyValueFromNodeLabel(nodeLabels[0]);

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

        std::string reactionFilename = mReactionFileRoot + this->ReturnKeyValueFromNodeLabel(nodeLabels[node_index]);

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
        mOdeSystem[node_index] -> SetDomainStateVariableRegister(this->GetStateVariableVector());
    }
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
void ChemicalDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::DeriveSystemProperties()
{
    // for when we want to track the diffusion of all speces in the system
    // this function derives the state register from the chemical reaction system
    // leads to incorrect initial conditions if the initial conditions take the form of all the diffusing 
    // species in the system
    
    // determine the system properties from the read in reaction systems
    
    AbstractDiffusiveChemistry* p_diffusive_chemistry = new AbstractDiffusiveChemistry();

    // create a diffusive chemistry from the diffusion database

    for(unsigned system_number=0; system_number<mOdeSystem.size(); system_number++)
    {
 
        std::vector<AbstractChemical*> ode_chemical_vector = mOdeSystem[system_number] -> GetReactionSystem() -> GetSystemChemistry() -> rGetChemicalVector();
    
        // for each chemical in the reaction system 
        for(unsigned ode_system_chemical_index=0; ode_system_chemical_index<ode_chemical_vector.size(); ode_system_chemical_index++)
        {
            // populate their diffusive properties from the labels in the system and the diffusion database
            
            for(unsigned domain_index=0; domain_index<this.GetNumberOfDomains(); domain_index++)
            {
                std::string domainLabel = this.GetDomainLabelByIndex(domain_index);

                std::string chemicalName = ode_chemical_vector[ode_system_chemical_index] -> GetChemicalName();
                double chemicalSize = ode_chemical_vector[ode_system_chemical_index] -> GetChemicalSize();
                double chemicalMass = ode_chemical_vector[ode_system_chemical_index] -> GetChemicalMass();
                int chemicalValence = ode_chemical_vector[ode_system_chemical_index] -> GetChemicalValence();

                // diffusion database may be larger than the active chemicals in the system, the chemical vector
                bool IsInDatabase = false;

                // for a diffusive chemical object through the information given in the chemical
                // change the pointer of the original chemical to point to the diffusive chemical
                AbstractDiffusiveChemical *p_chemical = new AbstractDiffusiveChemical(chemicalName,chemicalSize,chemicalMass,chemicalValence);
             
                for(unsigned record_index=0; record_index<this.mDiffusionDatabase.size(); record_index++)
                {
                    if(this.mDiffusionDatabase[record_index][0] == chemicalName)
                    {   
                        p_chemical -> AddDiffusiveDomain(this.mDiffusionDatabase[record_index][1],std::stod(this.mDiffusionDatabase[record_index][2]));
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

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
void ChemicalDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::DeriveExtendedSystemProperties()
{
    // determine the system properties from the read in reaction systems
    AbstractDiffusiveChemistry* p_diffusive_chemistry = new AbstractDiffusiveChemistry();

    // each chemical in the system deined in state variable vector, derived from diffusion database
    StateVariableRegister* p_stateVariableVector = this->GetStateVariableVector();

    for(unsigned chem_index=0; chem_index < p_stateVariableVector -> GetNumberOfStateVariables(); chem_index++)
    {
        // for each chemical in the diffusion database
        std::string chemicalName = p_stateVariableVector -> RetrieveStateVariableName(chem_index);

        // for each domain
        for(unsigned domain_index=0; domain_index<this->GetNumberOfDomains(); domain_index++)
        {
            std::string domainLabel = this->GetDomainLabelByIndex(domain_index);

            AbstractDiffusiveChemical *p_chemical = new AbstractDiffusiveChemical(chemicalName);
            bool IsInDatabase = false;
            // test whether the chemical and domain match, expecting at least one match
            for(unsigned record_index=0; record_index<this->GetDiffusionDatabase().size(); record_index++)
            {
                if(this->GetDiffusionDatabase()[record_index][0] == chemicalName)
                {   
                    p_chemical -> AddDiffusiveDomain(this->GetDiffusionDatabase()[record_index][1],std::stod(this->GetDiffusionDatabase()[record_index][2]));
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


template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
double ChemicalDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::GetDiffusionValueBasedOnPoint(const ChastePoint<2>& chastePoint, unsigned stateIndex)
{

    // take in the point and state index then determine the state name and domain label, then determine the diffusion value
    if(stateIndex<this ->GetStateVariableVector() ->GetNumberOfStateVariables())
    {
        
        // retrive state name from stateVariableRegister
        std::string stateName = this ->GetStateVariableVector()-> RetrieveStateVariableName(stateIndex);

        // retrieve label from domainLabels
        std::string domainLabel = this->ReturnDomainLabelAtPosition(chastePoint.rGetLocation());
        std::string domainKeyName = this->ReturnDomainKeyFromDomainLabel(domainLabel);

        return ReturnDiffusionValueFromStateNameAndDomainLabel(stateName,domainKeyName);
    
    }
    else
    {
        std::cout<<"ChemicalDomainFieldTemplated::GetDiffusionValueBasedOnPoint: State not in state variable"<<std::endl;
        return 0.0;
    }

}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
double ChemicalDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::ReturnDiffusionValueFromStateNameAndDomainLabel(std::string stateName, std::string domainLabel)
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

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
void ChemicalDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetDomainStateVariableRegister(StateVariableRegister* p_register)
{
    mpDomainRegister = p_register;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
void ChemicalDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetChemistry(AbstractDiffusiveChemistry* p_chemistry)
{
    mpDiffusiveChemistry = p_chemistry;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
void ChemicalDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetOdeSystem(std::vector<AbstractInhomogenousChemicalOdeSystemForCoupledPdeSystem*> odeSystem)
{
    mOdeSystem = odeSystem;
}

// get methods

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
StateVariableRegister* ChemicalDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::GetDomainStateVariableRegister()
{
    return mpDomainRegister;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
AbstractDiffusiveChemistry* ChemicalDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::GetChemistry()
{
    return mpDiffusiveChemistry;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
std::vector<AbstractInhomogenousChemicalOdeSystemForCoupledPdeSystem*> ChemicalDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::GetOdeSystem()
{
    return mOdeSystem;
}


#endif