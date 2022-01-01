#ifndef ABSTRACTDOMAINFIELD_TEMPLATED_HPP
#define ABSTRACTDOMAINFIELD_TEMPLATED_HPP

#include <string>
#include <vector>
#include "ChastePoint.hpp"
#include "TetrahedralMesh.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "AbstractTetrahedralMesh.hpp"
#include "StateVariableRegister.hpp"

#include "RandomNumberGenerator.hpp"

// class containing the information required to read and define a diffusive domain system
// for diffusion pde simulations.

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
class AbstractDomainFieldTemplated
{
protected:

    bool mIsFeMeshGenerated=false;

    TetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* mpFeMesh;

    std::vector<double> mLabelOrigin;

    // file names containing the system information

    unsigned mNumberOfSlices=0;

    std::string mDomainLabelFilename;

    std::string mDomainKeyFilename;

    std::string mOdeLabelFilename;

    std::string mOdeKeyFilename;

    std::string mDiffusionFilename;


    // domain information input
    unsigned mNumberOfDomains;

    std::vector<std::vector<std::vector<std::string>>> mDomainLabels;

    std::vector<std::string> mDomainLabelVector;

    std::vector<std::vector<std::string>> mDomainKeys;

    std::vector<unsigned> mCartesianCellDimensions; // (x, y, z) of the input csv file grid

    std::vector<double> mCartesianCellScaleXY; // (Sx, Sy, Sz) scale factor of the input csv grid, default (1,1)



    // ode system information input
    unsigned mNumberOfOdeSystems;

    std::vector<std::vector<std::vector<std::string>>> mOdeLabels;

    std::vector<std::string> mOdeLabelVector;

    std::vector<std::vector<std::string>> mOdeKeys;

    std::vector<double> mCartesianOdeScaleXY; // (Sx, Sy, Sz) scale factor of the input csv grid, default (1,1)

    std::vector<unsigned> mCartesianOdeDimensions; // (x, y, z) of the input csv file grid



    // collective system information
    unsigned mNumberOfStates; // number of unique species in the diffusion database

    StateVariableRegister* mpStateVariableVector; // container of each unique species in the diffusion database

    std::vector<std::vector<std::string>> mDiffusionDatabase;

    unsigned mDiffusionDomainLabelPosition =1;



    // mapped cartesian domain 
    std::vector<unsigned> mCartesianChasteDimensions;

    std::vector<double> mCartesianChasteScaleXY; // (Sx, Sy, Sz) scale factor of the Chaste grid

    std::vector<std::vector<std::vector<std::string>>> mCartesianChasteDomain;

    std::vector<std::vector<std::vector<std::string>>> mCartesianOdeDomain;


    // mapped mesh domain
    std::vector<unsigned> mMeshDimensions;

    std::vector<double> mMeshScale;

    double mStepSize=1.0;

    HoneycombMeshGenerator* mpGenerator;

    std::vector<std::string> mLabelledNodes;

    std::vector<std::string> mDomainNodes;

    // initial and boundary conditions
    std::string mInitialConditionsFilename;

    std::string mBoundaryConditionsFilename;

    std::vector<double> mInitalNodeConditions;

    std::vector<std::string> mBoundaryConditionTypes;

    std::vector<double> mBoundaryConditionValues;


    bool mIsHoneyCombMesh;



public:

    AbstractDomainFieldTemplated(
                        //TetrahedralMesh<ELEMENT_DIM,SPACE_DIM>*,
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

    virtual ~AbstractDomainFieldTemplated()
    {
    }

    virtual void SetUpDomainFromFiles();

    virtual void GenerateFeMesh();

    virtual void LabelFeMeshNodally();

    // virtual methods involved with basic function
    virtual void SetupAndInitialiseLabelDomain();

    virtual void SetupAndInitialiseOdes();

    virtual void SetupDomainMappingDimensions();

    virtual void FormMesh();

    virtual void SetupAndInitialiseDiffusionDatabase();

    virtual std::vector<std::string> LabelNodesWithOdes();

    virtual std::vector<std::string> LabelNodesWithDomains();

    virtual void MapToChasteCartesianDomain();

    virtual void MapOdeLabelsToCartesianOdeDomain();

    virtual double GetDiffusionValueBasedOnPoint(const ChastePoint<2>& , unsigned);

    virtual double ReturnDiffusionValueFromStateNameAndDomainLabel(std::string stateName, std::string domainLabel = "");

    virtual void ParseInitialConditionsFromFile(std::string initialConditionsFilename, bool IsPerturbInitialConditions=false );

    virtual void ParseBoundaryConditionsFromFile(std::string);

    std::vector<std::vector<std::string>> ReturnStateDataWithinRegister(std::vector<std::vector<std::string>>);

    virtual std::string GetFieldType()
    {
        return "AbstractDomainFieldTemplated";
    };


    // return methods
    std::string ReturnNodeOdeLabelAtPosition(const c_vector<double,SPACE_DIM>&);

    std::string ReturnDomainLabelAtPosition(const c_vector<double,SPACE_DIM>&);

    std::string ReturnMappedDomainLabel(std::vector<double>, unsigned xindex=0, unsigned yindex=0, unsigned zindex=0);

    std::string ReturnMappedOdeLabel(std::vector<double>, unsigned xindex=0, unsigned yindex=0, unsigned zindex=0);

    std::string ReturnKeyValueFromNodeLabel(std::string);

    std::string ReturnDomainKeyFromDomainLabel(std::string);

   
    // file read methods
    std::vector<std::vector<std::string>> ReadMatrix(std::string);

    std::vector<std::string> parseMatrixLineString(std::string);

    void ReadDomainLabels();

    void ReadDomainKey();

    void ReadOdeLabels();

    void ReadOdeKey();

    void ReadDiffusionDatabase();


    // auxiliary methods
    std::vector<std::string> ReturnUnique(std::vector<std::string>);

    std::vector<std::string> ReturnUnique(std::vector<std::vector<std::string>>);

    std::vector<std::string> ReturnUnique(std::vector<std::vector<std::vector<std::string>>>);

    void printMatrix(std::vector<std::vector<std::vector<std::string>>>);

    void printMatrix(std::vector<std::vector<std::string>>);

    void printVector(std::vector<std::string>);

    bool PerturbInitialConditionTest(std::vector<std::string>);

    // display methods

    void PrintDiffusionDomain();

    void PrintMappedDiffusionDomain();

    void PrintDomainLabelKeys();

    void PrintODEDomain();

    void PrintMappedODEDomain();

    void PrintODELabelKeys();

    void PrintDiffusionDatabase();

    // set methods
    void SetNumberOfDomains(unsigned);

    void SetDomainLabels(std::vector<std::vector<std::vector<std::string>>> );

    void SetDomainLabelVector(std::vector<std::string>);

    void SetDomainKeys(std::vector<std::vector<std::string>>);    

    void SetNumberOfOdeSystems(unsigned);

    void SetOdeLabels(std::vector<std::vector<std::vector<std::string>>> );

    void SetOdeLabelVector(std::vector<std::string>);    

    void SetOdeKeys(std::vector<std::vector<std::string>>);

    void SetStateVariableVector(StateVariableRegister*);

    void SetNumberOfStates(unsigned);    
        
    void SetCartesianCellDimensions(std::vector<unsigned>);

    void SetCartesianCellScale(std::vector<double>);

    void SetCartesianOdeDimensions(std::vector<unsigned>);

    void SetCartesianOdeScale(std::vector<double>);
    
    void SetCartesianChasteScale(std::vector<double>);

    void SetCartesianChasteDomain(std::vector<std::vector<std::vector<std::string>>> );

    void SetCartesianOdeDomain(std::vector<std::vector<std::vector<std::string>>> );

    void SetCartesianChasteDimensions(std::vector<unsigned>);

    void SetMeshDimensions(std::vector<unsigned>);

    void SetDomainFeMesh(TetrahedralMesh<ELEMENT_DIM,SPACE_DIM>*&);

    void SetMeshGenerator(HoneycombMeshGenerator*);

    void SetMeshScale(std::vector<double>);

    void SetDiffusionDatabase(std::vector<std::vector<std::string>>);

    void SetNodeLabels(std::vector<std::string>);

    void SetNodeDomains(std::vector<std::string>);

    void SetInitialNodeConditions(std::vector<double>);

    void SetBoundaryConditionTypes(std::vector<std::string>);

    void SetBoundaryConditionValues(std::vector<double>);

    void SetLabelOrigin(std::vector<double>);

    // get methods
    unsigned GetNumberOfDomains();

    std::vector<std::vector<std::vector<std::string>>> GetDomainLabels();

    std::vector<std::string> GetDomainLabelVector();

    std::string GetDomainLabelByIndex(unsigned);

    std::vector<std::vector<std::string>> GetDomainKeys();

    unsigned GetNumberOfOdeSystems();

    std::vector<std::vector<std::vector<std::string>>>  GetOdeLabels();

    std::vector<std::string> GetOdeLabelVector();

    std::vector<std::vector<std::string>> GetOdeKeys();

    StateVariableRegister* GetStateVariableVector();

    unsigned GetNumberOfStates();

    std::vector<unsigned> GetCartesianCellDimensions();

    std::vector<double> GetCartesianCellScale();

    std::vector<unsigned> GetCartesianOdeDimensions();

    std::vector<double> GetCartesianOdeScale();

    std::vector<std::vector<std::vector<std::string>>>  GetCartesianChasteDomain();

    std::vector<double> GetCartesianChasteScale();

    std::vector<std::vector<std::vector<std::string>>>  GetCartesianOdeDomain();

    std::vector<unsigned> GetCartesianChasteDimensions();

    std::vector<unsigned> GetMeshDimensions();

    TetrahedralMesh<ELEMENT_DIM,SPACE_DIM>*& rGetDomainFeMesh();

    HoneycombMeshGenerator* GetMeshGenerator();

    std::vector<double> GetMeshScale();

    std::vector<std::vector<std::string>> GetDiffusionDatabase();

    std::vector<std::string> GetNodeLabels();

    std::vector<std::string> GetNodeDomains();

    std::vector<double> GetInitialNodeConditions();

    std::vector<std::string> GetBoundaryConditionTypes();

    std::vector<double> GetBoundaryConditionValues();

    std::vector<double> GetLabelOrigin();

};

//======================================================================//
//                          IMPLEMENTATION                              //  
//======================================================================//

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::AbstractDomainFieldTemplated(
                                            //TetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pFeMesh,   
                                            std::string domainLabelFilename, 
                                            std::string domainKeyFilename, 
                                            std::string odeLabelFilename, 
                                            std::string odeKeyFilename, 
                                            std::string diffusionFilename,
                                            bool isHoneyCombMesh,
                                            std::vector<double> labelOrigin,
                                            std::vector<double> cartesianCellScale,
                                            std::vector<double> cartesianOdeScale)
    :   //mpFeMesh(pFeMesh),
        mDomainLabelFilename(domainLabelFilename),
        mDomainKeyFilename(domainKeyFilename),
        mOdeLabelFilename(odeLabelFilename),
        mOdeKeyFilename(odeKeyFilename),
        mDiffusionFilename(diffusionFilename),
        mIsHoneyCombMesh(isHoneyCombMesh),
        mLabelOrigin(labelOrigin),
        mCartesianCellScaleXY(cartesianCellScale),
        mCartesianOdeScaleXY(cartesianOdeScale)
        
        
{   //std::cout<<"AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::AbstractDomainFieldTemplated"<<std::endl;

    if(mLabelOrigin.empty())
    {
        std::vector<double> origin(SPACE_DIM,0);
        mLabelOrigin = origin;
    }

    if(mCartesianCellScaleXY.empty())
    {
        // default to scale 1.0 in all directions
        std::vector<double> cellScaleXY(SPACE_DIM,1.0);
        mCartesianCellScaleXY = cellScaleXY;
    }

    if(mCartesianOdeScaleXY.empty())
    {
        // default to scale 1.0 in all directions
        std::vector<double> odeScaleXY(SPACE_DIM,1.0);
        mCartesianOdeScaleXY = odeScaleXY;
    }

    SetUpDomainFromFiles();
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
void AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetUpDomainFromFiles()
{//std::cout<<"AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetUpDomainFromFiles"<<std::endl;
    // read the domain labels and label keys form files

    SetupAndInitialiseLabelDomain();

    // scale the input file dimensions to the chaste rectilinear cartesian grid
    MapToChasteCartesianDomain();

    // translate the ode labels to ode nodes, form the ode vector
    SetupDomainMappingDimensions();

    // read in the ode information, set the number of ode systems and 
    SetupAndInitialiseOdes();

    MapOdeLabelsToCartesianOdeDomain();

    // read in the domain information, set the number of domains
    SetupAndInitialiseDiffusionDatabase();

}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
void AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::GenerateFeMesh()
{

    // generate the finite element mesh and populate the nodes with the corrtesponding ode systems
    // and solvers populated from the domain information. Function to be called whenever the domain
    // pde mesh changes in the simulation.

    // form the mesh as either a honey comb lattice or a tetrahedral mesh, dependent on mIsHoneycomb
    FormMesh();

    // label the newly formed mesh
    LabelFeMeshNodally();

    mIsFeMeshGenerated = true;

}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
void AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::LabelFeMeshNodally()
{

    // for a mesh that has been already defined, look at location of each node, compare with reference
    // chaste domain and populate node with corresponding property label

    // store the node labels as std::string type
    SetNodeDomains(LabelNodesWithDomains());

    SetNodeLabels(LabelNodesWithOdes());

}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
void AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetupAndInitialiseLabelDomain()
{
    // read in and store the domain grid labels for mapping to mesh nodes and chaste rectilinear grid
    
    std::vector<unsigned> cartesianCellDimensions;

    ReadDomainLabels();
    ReadDomainKey();

    cartesianCellDimensions.push_back(mDomainLabels[0][0].size());
    cartesianCellDimensions.push_back(mDomainLabels[0].size());
    cartesianCellDimensions.push_back(mDomainLabels.size());

    // store these values
    SetCartesianCellDimensions(cartesianCellDimensions);
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
void AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetupAndInitialiseOdes()
{
    // read in and store the ode labels from the input grid, take in the ode label keys 
    ReadOdeLabels();
    ReadOdeKey();

    std::vector<unsigned> cartesianOdeDimensions;

    cartesianOdeDimensions.push_back(mOdeLabels[0][0].size());
    cartesianOdeDimensions.push_back(mOdeLabels[0].size());
    cartesianOdeDimensions.push_back(mOdeLabels.size());

    // store these values
    SetCartesianOdeDimensions(cartesianOdeDimensions);

}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
void AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetupDomainMappingDimensions()
{
    // hard code the honey comb mesh scale factor and dimensions with a default of tetrahedral mesh
    std::vector<double> meshScale;
    if(mIsHoneyCombMesh)
    {   
        switch (SPACE_DIM)
        {
            case 2:
                meshScale.push_back(1.0);
                meshScale.push_back(0.866025);
                break;
            default:
            
                for(unsigned dim=0;dim<SPACE_DIM; dim++)
                {
                    meshScale.push_back(1.0);
                }
        }
    }
    else
    {   
        // mesh is tetrahedral
        for(unsigned dim=0;dim<SPACE_DIM; dim++)
        {
            meshScale.push_back(1.0);
        }
    }
    SetMeshScale(meshScale);
    

    std::vector<unsigned> meshDimensions(SPACE_DIM,0);
    std::vector<double> meshSegment(SPACE_DIM,0.0);

    // determine the mesh dimensions up to the maximum allowed through the scaling
    // cover the whole domain with the mesh

    for(unsigned dim=0; dim<SPACE_DIM; dim++)
    {
        unsigned count=0;
        while( meshSegment[dim] < mCartesianCellDimensions[dim]*mCartesianCellScaleXY[dim])
        {
            meshSegment[dim]  = count*meshScale[dim];
            count=count+1;
        }
        meshDimensions[dim]=count-2;
    }

    SetMeshDimensions(meshDimensions);
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
void AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::FormMesh()
{
    // use the aforementioned mesh dimension to produce a new honeycomb mesh, as mutable mesh type

    //HoneycombMeshGenerator generator(meshDimensions[0], meshDimensions[1], 0);
    if(mIsHoneyCombMesh)
    {
        HoneycombMeshGenerator* p_generator = new HoneycombMeshGenerator(mMeshDimensions[0], mMeshDimensions[1], 0);

        mpFeMesh = dynamic_cast<TetrahedralMesh<ELEMENT_DIM,SPACE_DIM>*>(p_generator -> GetMesh());

        SetMeshGenerator(p_generator);
        SetDomainFeMesh(mpFeMesh);
    }
    else
    {
        mpFeMesh = new TetrahedralMesh<ELEMENT_DIM,SPACE_DIM>();

        switch (SPACE_DIM)
        {
            case 1:
                mpFeMesh->ConstructRegularSlabMesh(mStepSize, mMeshDimensions[0]);
                break;
            case 2:
                mpFeMesh->ConstructRegularSlabMesh(mStepSize, mMeshDimensions[0], mMeshDimensions[1]);
                break;
            case 3:
                mpFeMesh->ConstructRegularSlabMesh(mStepSize, mMeshDimensions[0], mMeshDimensions[1], mMeshDimensions[2]);
                break;
            default:
                NEVER_REACHED;
        }
        SetDomainFeMesh(mpFeMesh);
    }
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
void AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetupAndInitialiseDiffusionDatabase()
{
    // to be overriden for further meshing and diffusion schemes
    ReadDiffusionDatabase();
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
std::vector<std::string> AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::LabelNodesWithDomains()
{
    // retireve the node positions and compare to domain labels in order to label the nodes
    // nodes are defined bottom-left of domain to top-right in a serialised manner

    unsigned numberOfNodes = mpFeMesh -> GetNumNodes();

    // form the serialised node label array
    std::vector<std::string> serialisedNodeDomains(numberOfNodes,"");
    
    for (typename AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>::NodeIterator iter = mpFeMesh ->GetNodeIteratorBegin();
             iter != mpFeMesh ->GetNodeIteratorEnd();
             ++iter)
    {
        // for each node, retireve index and position, label with correct ode label
        unsigned node_index = iter ->GetIndex();
        // retireve the position of the node of index node_index
        const c_vector<double,SPACE_DIM>& position = iter->rGetLocation();
        // search through the position ode labels for the correct label to give the node
        serialisedNodeDomains.at(node_index) = ReturnDomainLabelAtPosition(position);
    }

    return serialisedNodeDomains;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
std::vector<std::string> AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::LabelNodesWithOdes()
{
    // retireve the node positions and compare to domain labels in order to label the nodes
    // nodes are defined bottom-left of domain to top-right in a serialised manner
   
    unsigned numberOfNodes = mpFeMesh -> GetNumNodes();

    // form the serialised node label array
    std::vector<std::string> serialisedNodeLabels(numberOfNodes,"");
    
    for (typename AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>::NodeIterator iter = mpFeMesh ->GetNodeIteratorBegin();
             iter != mpFeMesh ->GetNodeIteratorEnd();
             ++iter)
    {
        // for each node, retireve index and position, label with correct ode label
        unsigned node_index = iter ->GetIndex();
        // retireve the position of the node of index node_index
        const c_vector<double,SPACE_DIM>& position = iter->rGetLocation();
        // search through the position ode labels for the correct label to give the node
        serialisedNodeLabels.at(node_index) = ReturnNodeOdeLabelAtPosition(position);
    }

    return serialisedNodeLabels;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
void AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::ParseInitialConditionsFromFile(std::string initialConditionsFilename, bool IsPerturbInitialConditions)
{

    // read the input intial conditions file as a matrix of strings, first column state name, 
    // 2nd column subdomain, 3rd column value.

    // if the state is not within the state register then the state doesn't evolve over time, remaining constant
    std::vector<std::vector<std::string>> fullInitialConditionsAsStrings = ReadMatrix(initialConditionsFilename);

    std::vector<std::vector<std::string>> initialConditionsAsStrings = ReturnStateDataWithinRegister(fullInitialConditionsAsStrings);

    // use state name against stateVariableRegister to remove any state value not in the system
    // replace any state in register but not in intitial conditions with value of 0.0
    
    // for the vector in order of StateVariableRegister
    unsigned numberOfStateVariables = mpStateVariableVector->GetNumberOfStateVariables();
    unsigned numberOfNodes = mpFeMesh -> GetNumNodes();
    
 
    // determine the nodal initial condition, test for perturbation
    std::vector<double> init_conds(numberOfStateVariables*numberOfNodes,0.0);
    std::vector<std::string> node_domain_labels = GetNodeDomains();
    for (unsigned node_index=0; node_index<numberOfNodes; node_index++)
    {   
        // each node set as being a vector of the read in state condition values

        // test for the domain of the node, retireve the condition for the specific domain
        for(unsigned pdeDim=0; pdeDim<numberOfStateVariables; pdeDim++)
        {   
            // determine which state the pde dimension refers to and then retrieve the condition if present
            std::string stateName = mpStateVariableVector -> RetrieveStateVariableName(pdeDim);
  
            bool IsFoundState =false;
            // run through input data, line by line
            for(unsigned inputState=0; inputState<initialConditionsAsStrings.size();inputState++)
            {
                // test whether the state variable is specified within the input data
                if(stateName==initialConditionsAsStrings[inputState][0])
                { 
                    // the state is in the input date but now test for the correct domain
                    if(initialConditionsAsStrings[inputState].size()>2)
                    {

                        // then sub domain is also specified, test for sub domain
                        if(node_domain_labels[node_index]==initialConditionsAsStrings[inputState][1]||ReturnDomainKeyFromDomainLabel(node_domain_labels[node_index])==initialConditionsAsStrings[inputState][1])
                        {  
                            // sub domain has been found, state found, fill data record
                            IsFoundState = true;
                           
                            // store the condition, serialised for nodes

                            if(PerturbInitialConditionTest(initialConditionsAsStrings[inputState])||IsPerturbInitialConditions)
                            {
                            
                                init_conds[numberOfStateVariables*node_index + pdeDim] =fabs(std::stod(initialConditionsAsStrings[inputState][2]) + RandomNumberGenerator::Instance()->ranf());
                            }
                            else
                            {
                                init_conds[numberOfStateVariables*node_index + pdeDim] =std::stod(initialConditionsAsStrings[inputState][2]);
             
                            }
                            break;
                        }
                        // subdomain specified, but not found.  State not present in particular subdomain
                    }
                    else
                    {
                        //state data found, but no sub domain is specied, assume present for all sub domains if sub domains are indeed present
                        IsFoundState = true;
                        // store the condition, serialised for nodes
                        if(IsPerturbInitialConditions)
                        {
                            init_conds[numberOfStateVariables*node_index + pdeDim] =fabs(std::stod(initialConditionsAsStrings[inputState][2]) + RandomNumberGenerator::Instance()->ranf());
                        }
                        else
                        {
                            init_conds[numberOfStateVariables*node_index + pdeDim] =std::stod(initialConditionsAsStrings[inputState][2]);
  
                        }
                        break;
                    }   
                }
                // data record isn't correct for the pde state name, if at end of data record state not found then state not present in input data
            }  
            if(!IsFoundState)
            {
                // state in system but not in the input data, default to 0.0
                // initialse the missing condition, serialised for nodes    
                if(IsPerturbInitialConditions)
                {
                    // initialise to a small perturbation, otherwise leave as intitialised 0.0
                    init_conds[numberOfStateVariables*node_index + pdeDim] =fabs(RandomNumberGenerator::Instance()->ranf());
                }
            }


            std::cout<<"Initial condition: State: "<<stateName<<" "<<init_conds[numberOfStateVariables*node_index + pdeDim]<<std::endl;
        }
        // each node pde is set
    }

    SetInitialNodeConditions(init_conds);
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
bool AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::PerturbInitialConditionTest(std::vector<std::string> inputVector)
{
    // determine whether or not the initial conditions is to be perturbed
    // perturbation on more specific basis
        
    if(inputVector.size()>3)
    {
        // then have the perturbation bool selections
        if(inputVector[3]=="true"||inputVector[3]=="True"||inputVector[3]=="TRUE")
        {
            return true;
            //IsPerturbSpecificConditions[mpStateVariableVector -> RetrieveStateVariableIndex(initialConditionsAsStrings[i][0])] = true;
        }
        else
        {
            return false;
        }
        
    }
    else
    {
        return false;
    }
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
std::string AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::ReturnKeyValueFromNodeLabel(std::string nodeLabel)
{
    bool IsFound=false;

    for(unsigned key_index=0; key_index<mOdeKeys.size(); key_index++)
    {
        //std::cout<<"key_index: "<<key_index<<std::endl;
        if(mOdeKeys[key_index][0] == nodeLabel)
        {
            return mOdeKeys[key_index][1];
            IsFound=true;
            break;
        }
    }
    if(!IsFound)
    {
        std::cout<<"Error: AbstractDomainFieldTemplated::ReturnKeyValueFromNodeLabel, node label not found"<<std::endl;
        return "Null";
    }
    return "Null";
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
std::string AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::ReturnDomainKeyFromDomainLabel(std::string domainLabel)
{
    bool IsFound=false;

    for(unsigned key_index=0; key_index<mDomainKeys.size(); key_index++)
    {
        
        if(mDomainKeys[key_index][0] == domainLabel)
        {
            return mDomainKeys[key_index][1];
            IsFound=true;
            break;
        }
    }
    if(!IsFound)
    {
        std::cout<<"Error: AbstractDomainFieldTemplated::ReturnDomainKeyFromDomainLabel, domain label not found"<<std::endl;
        return "Null";
    }
    return "Null";
}


template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
std::vector<std::vector<std::string>> AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::ReturnStateDataWithinRegister(std::vector<std::vector<std::string>> inputMatrix)
{
    // take an input matrix read from a file wherein the first column denotes a state
    // check against the StateVariableRegister for states that are in the system and remove input
    // that are not in the system
    std::string inputName;
    std::vector<std::vector<std::string>> newMatrix;

    for(unsigned i=0; i<inputMatrix.size();i++)
    {
        inputName = inputMatrix[i][0];

        if(mpStateVariableVector -> IsStateVariablePresent(inputName))
        {
            newMatrix.push_back(inputMatrix[i]);
        }
    }
    // newMatrix conatins only the data records of states in the system
    return newMatrix;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
void AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::ParseBoundaryConditionsFromFile(std::string boundaryConditionsFilename)
{
    // read the input boundary conditions file as a matrix of strings, first column state name, 
    // 2nd column boundary condition type, 3rd column value.
    std::vector<std::vector<std::string>> fullBoundaryConditionsAsStrings = ReadMatrix(boundaryConditionsFilename);
    std::vector<std::vector<std::string>> boundaryConditionsAsStrings = ReturnStateDataWithinRegister(fullBoundaryConditionsAsStrings);
    // use state name against stateVariableRegister to remove any state value not in the system
    // replace any state in register but not in boundary conditions with value of 0.0, type Derichlet

    unsigned numberOfStateVariables = GetStateVariableVector() ->GetNumberOfStateVariables();

    std::vector<std::string> boundaryConditionTypes(numberOfStateVariables,"Dirichlet");
    std::vector<double> boundaryConditionValues(numberOfStateVariables,0.0);

    for(unsigned pdeDim=0; pdeDim<numberOfStateVariables; pdeDim++)
    {   
        // determine which state the pde dimension refers to and then retrieve the condition if present
        std::string stateName = mpStateVariableVector -> RetrieveStateVariableName(pdeDim);

        for(unsigned inputState=0; inputState<boundaryConditionsAsStrings.size();inputState++)
        {
            // test whether the state variable is specified within the input data
            if(stateName==boundaryConditionsAsStrings[inputState][0])
            {
                // state found, read data
                if(boundaryConditionsAsStrings[inputState].size()>2)
                {
                    // type value is specified also
                    boundaryConditionTypes[pdeDim] = boundaryConditionsAsStrings[inputState][1];
                    boundaryConditionValues[pdeDim] = std::stod(boundaryConditionsAsStrings[inputState][2]);
                }
                else
                {
                    // type value not specified, assume Dirichlet
                    boundaryConditionValues[pdeDim] = std::stod(boundaryConditionsAsStrings[inputState][1]);
                }
                break;
            }
        }
        // if state isn't found, leave as default, 0.0 Dirichlet
        std::cout<<"Boundary condition: State: "<<stateName<<" "<<boundaryConditionValues[pdeDim]<<std::endl;

    }

    SetBoundaryConditionTypes(boundaryConditionTypes);

    SetBoundaryConditionValues(boundaryConditionValues);
}




template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
std::string AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::ReturnNodeOdeLabelAtPosition(const c_vector<double,SPACE_DIM>& position)
{
    // search through the node labels to provide the specific label at the given position

    // search through the domain labels to provide the specific label at the given position
    std::vector<unsigned> labelIndex(SPACE_DIM,0);

    // translate the simulation position from the origin to the domain grid
    double position_in_domain=0.0;

    // scale the dimension to provide the X,Y indices of the grid element containing the position under inspection
    for(unsigned dim=0; dim<SPACE_DIM; dim++)
    {
        position_in_domain = position[dim] + mLabelOrigin[dim];

        for(unsigned i=1; i<(mCartesianOdeDimensions[dim]+1); i++)
        {
            if(position_in_domain <= i*mCartesianOdeScaleXY[dim])
            {
                labelIndex[dim] = i-1;
                break;
            }
        }
    }

    // return the ode label for the grid element at the X,Y index containing the node
    return mOdeLabels[labelIndex[2]][labelIndex[1]][labelIndex[0]];
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
std::string AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::ReturnDomainLabelAtPosition(const c_vector<double,SPACE_DIM>& position)
{
    // search through the domain labels to provide the specific label at the given position
    std::vector<unsigned> labelIndex(SPACE_DIM,0);

    // translate the simulation position from the origin to the domain grid
    double position_in_domain=0.0;

    // scale the dimension to provide the X,Y indices of the grid element containing the position under inspection
    for(unsigned dim=0; dim<SPACE_DIM; dim++)
    {
        position_in_domain = position[dim] + mLabelOrigin[dim];

        for(unsigned i=1; i<(mCartesianChasteDimensions[dim]+1); i++)
        {
            if(position_in_domain <= i*mCartesianChasteScaleXY[dim])
            {
                labelIndex[dim] = i-1;
                break;
            }
        }
    }

    switch(SPACE_DIM)
    {
        case 2:
            return mCartesianChasteDomain[0][labelIndex[1]][labelIndex[0]];
            break;

        case 3:
            return mCartesianChasteDomain[labelIndex[2]][labelIndex[1]][labelIndex[0]];
            break;
        
    }
    
    // return the domain label for the grid element at the Z,Y,X index containing the position
    return "";
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
void AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::MapToChasteCartesianDomain()
{
    
    std::vector<double> cartesianChasteScaleXY;
    if(mIsHoneyCombMesh)
    {
        // is a a honeycomb mesh
        cartesianChasteScaleXY.push_back(0.25);
        cartesianChasteScaleXY.push_back(0.144338);
        cartesianChasteScaleXY.push_back(0.144338);

    }
    else
    {
        // else tetrahedral mesh
        cartesianChasteScaleXY.push_back(0.1);
        cartesianChasteScaleXY.push_back(0.1);
        cartesianChasteScaleXY.push_back(0.1);
    }
    SetCartesianChasteScale(cartesianChasteScaleXY);

    unsigned numberOfXSteps = std::round((mCartesianCellDimensions[0]*mCartesianCellScaleXY[0])/cartesianChasteScaleXY[0]); 
    unsigned numberOfYSteps = std::round((mCartesianCellDimensions[1]*mCartesianCellScaleXY[1])/cartesianChasteScaleXY[1]); 
    unsigned numberOfZSteps = std::round((mCartesianCellDimensions[2]*mCartesianCellScaleXY[2])/cartesianChasteScaleXY[2]); 
    std::vector<unsigned> cartesianChasteDimensions;
    cartesianChasteDimensions.push_back(numberOfXSteps);
    cartesianChasteDimensions.push_back(numberOfYSteps);
    cartesianChasteDimensions.push_back(numberOfZSteps);


    std::vector<std::vector<std::vector<std::string>>> cartesianChasteDomain;
    std::vector<std::vector<std::string>> tempCartesianChasteDomain;
    
    switch(SPACE_DIM)
    {
        case 2:

            // define the boundary with the edge steps included for easier comparisons of (double) spatial values
            for(unsigned yindex=0; yindex<numberOfYSteps; yindex++)
            {
                std::vector<std::string> rowVector;
                for(unsigned xindex=0; xindex<numberOfXSteps; xindex++)
                {
                    // push back the domain label of the file domain onto the mapped domain, here everything is in terms of std::string label
                    rowVector.push_back(ReturnMappedDomainLabel(mCartesianChasteScaleXY,xindex,yindex,0));

                }

                tempCartesianChasteDomain.push_back(rowVector);
            }

            cartesianChasteDomain.push_back(tempCartesianChasteDomain);
            break;

        case 3:

            // define the boundary with the edge steps included for easier comparisons of (double) spatial values
            for(unsigned zindex=0; zindex<numberOfZSteps; zindex++)
            {
                for(unsigned yindex=0; yindex<numberOfYSteps; yindex++)
                {
                    std::vector<std::string> rowVector;
                    for(unsigned xindex=0; xindex<numberOfXSteps; xindex++)
                    {
                        // push back the domain label of the file domain onto the mapped domain, here everything is in terms of std::string label
                        rowVector.push_back(ReturnMappedDomainLabel(mCartesianChasteScaleXY,xindex,yindex,zindex));

                    }

                    tempCartesianChasteDomain.push_back(rowVector);
                }

                cartesianChasteDomain.push_back(tempCartesianChasteDomain);
            }
            break;

        case 1:
            std::cout<<"Error: Abstract Domain 1D not supported"<<std::endl;
            break;
    }

    SetCartesianChasteDomain(cartesianChasteDomain);
    SetCartesianChasteDimensions(cartesianChasteDimensions);

}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
void AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::MapOdeLabelsToCartesianOdeDomain()
{
    // map the labels of the input grid file to the scaled chaste cartesian grid.
    
    std::vector<double> cartesianChasteScaleXY = GetCartesianChasteScale();

    // determine the number of steps in a mapped chaste grid based on the the number of input domain dimensions and the
    // ration of scaling parameters, round to the nearest integer value.
    unsigned numberOfXSteps = std::round((mCartesianOdeDimensions[0]*mCartesianOdeScaleXY[0])/cartesianChasteScaleXY[0]); 
    unsigned numberOfYSteps = std::round((mCartesianOdeDimensions[1]*mCartesianOdeScaleXY[1])/cartesianChasteScaleXY[1]); 
    unsigned numberOfZSteps = std::round((mCartesianOdeDimensions[2]*mCartesianOdeScaleXY[2])/cartesianChasteScaleXY[2]);

    std::vector<std::vector<std::string>> tempCartesianOdeDomain;
    std::vector<std::vector<std::vector<std::string>>> cartesianOdeDomain;

    switch(SPACE_DIM)
    {
        case 2:

            // define the boundary with the edge steps included for easier comparisons of (double) spatial values
            for(unsigned yindex=0; yindex<numberOfYSteps; yindex++)
            {
                std::vector<std::string> rowVector;
                for(unsigned xindex=0; xindex<numberOfXSteps; xindex++)
                {
                    // push back the ode label of the file domain onto the mapped domain, here everything is in terms of std::string label
                    rowVector.push_back(ReturnMappedOdeLabel(mCartesianChasteScaleXY,xindex,yindex,0));

                }

                tempCartesianOdeDomain.push_back(rowVector);
            }

            cartesianOdeDomain.push_back(tempCartesianOdeDomain);
            break;

        case 3:

            // define the boundary with the edge steps included for easier comparisons of (double) spatial values
            for(unsigned zindex=0; zindex<numberOfZSteps; zindex++)
            {
                for(unsigned yindex=0; yindex<numberOfYSteps; yindex++)
                {
                    std::vector<std::string> rowVector;
                    for(unsigned xindex=0; xindex<numberOfXSteps; xindex++)
                    {
                        // push back the ode label of the file domain onto the mapped domain, here everything is in terms of std::string label
                        rowVector.push_back(ReturnMappedOdeLabel(mCartesianChasteScaleXY,xindex,yindex,zindex));

                    }

                    tempCartesianOdeDomain.push_back(rowVector);
                }

                cartesianOdeDomain.push_back(tempCartesianOdeDomain);
            }
            break;

        case 1:
            std::cout<<"Error: Abstract Domain 1D not supported"<<std::endl;
            break;
    }

    SetCartesianOdeDomain(cartesianOdeDomain);
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
double AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::GetDiffusionValueBasedOnPoint(const ChastePoint<2>& chastePoint, unsigned stateIndex)
{

    // take in the point and state index then determine the state name and domain label, then determine the diffusion value
    if(stateIndex<mpStateVariableVector ->GetNumberOfStateVariables())
    {
        if(mDiffusionDatabase[0].size() ==1 )
        {
            // only a diffusion value vector
            return std::stod(mDiffusionDatabase[stateIndex][0]);
        }
        else
        {
            // retrive state name from stateVariableRegister
            std::string stateName = mpStateVariableVector -> RetrieveStateVariableName(stateIndex);

            // retrieve label from domainLabels
            std::string domainLabel = ReturnDomainLabelAtPosition(chastePoint.rGetLocation());
            return ReturnDiffusionValueFromStateNameAndDomainLabel(stateName,domainLabel);
        }
    }
    else
    {
        // state hasn't been found therefore return a diffusion value of 0.0
        std::cout<<"Error: AbstractDomainFieldTemplated::GetDiffusionValueBasedOnPoint: State not in state variable"<<std::endl;
        return 0.0;
    }
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
double AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::ReturnDiffusionValueFromStateNameAndDomainLabel(std::string stateName, std::string domainLabel)
{
    // look for diffusion value (from diffusiondatabase) based on state name and domain label

    unsigned numberOfDiffusiveEntries = mDiffusionDatabase.size();
    bool IsStateFound=false;
    if(mDiffusionDatabase[0].size() ==2)
    {
        // take to be only state name, diffusion value structure
        // find stateName
        for(unsigned diffusive_index=0; diffusive_index<numberOfDiffusiveEntries; diffusive_index++)
        {
            if(stateName == mDiffusionDatabase[diffusive_index][0])
            {
                IsStateFound = true;
                return std::stod(mDiffusionDatabase[diffusive_index][1]);
                break;
                
            }
        }
        // if the state isn't found in the diffusion database then state cannot diffuse, return value of 0.0
        if(!IsStateFound)
        {
            return 0.0;
        }
    }
    else
    {
        // find stateName, test for domain name
        for(unsigned diffusive_index=0; diffusive_index<numberOfDiffusiveEntries; diffusive_index++)
        {
            if(stateName == mDiffusionDatabase[diffusive_index][0])
            {
                // found state, then test for domain
                if(domainLabel == mDiffusionDatabase[diffusive_index][1])
                {
                    // state and domain found, return diffusivity
                    IsStateFound = true;
                    return std::stod(mDiffusionDatabase[diffusive_index][2]);
                    break;
                }
            }
        }
        // if the state isn't found in the diffusion database then state cannot diffuse, return value of 0.0
        if(!IsStateFound)
        {
            return 0.0;
        }
    }
    return 0.0;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
std::string AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::ReturnMappedDomainLabel(std::vector<double> scaleFactor, unsigned xindex , unsigned yindex, unsigned zindex)
{
    // compare the xindex and yindex of the new chaste domain to the label domain
    // return the label of the chaste domain, work from the top left to bottom right

    unsigned labelXIndex =0;
    unsigned labelYIndex =0;
    unsigned labelZIndex =0;

    double positionX = xindex*scaleFactor[0];
    double positionY = yindex*scaleFactor[1];
    double positionZ = zindex*scaleFactor[2];

    switch(SPACE_DIM)
    {
        case 2:

            // for the y position

            for(unsigned i=1; i<(mCartesianCellDimensions[1]+1); i++) // top to bottom
            {
                if(positionY <= i*mCartesianCellScaleXY[1])
                {
                    labelYIndex = i-1;
                    break;
                }
            }

            // for the x position

            for(unsigned i=1; i<(mCartesianCellDimensions[0]+1); i++)
            {
                if(positionX <= i*mCartesianCellScaleXY[0])
                {
                    labelXIndex = i-1;
                    break;
                }
            }

            break; 

        case 3:
            // for the z position

            for(unsigned i=1; i<(mCartesianCellDimensions[2]+1); i++) // top to bottom
            {
                if(positionZ <= i*mCartesianCellScaleXY[2])
                {
                    labelZIndex = i-1;
                    break;
                }
            }
            
            // for the y position

            for(unsigned i=1; i<(mCartesianCellDimensions[1]+1); i++) // top to bottom
            {
                if(positionY <= i*mCartesianCellScaleXY[1])
                {
                    labelYIndex = i-1;
                    break;
                }
            }

            // for the x position

            for(unsigned i=1; i<(mCartesianCellDimensions[0]+1); i++)
            {
                if(positionX <= i*mCartesianCellScaleXY[0])
                {
                    labelXIndex = i-1;
                    break;
                }
            }

            break;

        case 1:
            std::cout<<"Error: Abstract Domain 1D not supported"<<std::endl;
            break;
    }

    return mDomainLabels[labelZIndex][labelYIndex][labelXIndex];
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
std::string AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::ReturnMappedOdeLabel(std::vector<double> scaleFactor, unsigned xindex, unsigned yindex, unsigned zindex)
{
    // compare the xindex and yindex of the new chaste domain to the label ode
    // return the label of the chaste domain, work from the top left to bottom right

    unsigned labelXIndex =0;
    unsigned labelYIndex =0;
    unsigned labelZIndex =0;

    std::string label;
    double positionX = xindex*scaleFactor[0];
    double positionY = yindex*scaleFactor[1];
    double positionZ = zindex*scaleFactor[2];

    switch(SPACE_DIM)
    {
        case 2:

            // for the z position

            for(unsigned i=1; i<(mCartesianOdeDimensions[2]+1); i++) // top to bottom
            {
                if(positionZ <= i*mCartesianOdeScaleXY[2])
                {
                    labelZIndex = i-1;
                    break;
                }
            }
            
            // for the y position

            for(unsigned i=1; i<(mCartesianOdeDimensions[1]+1); i++) // top to bottom
            {
                if(positionY <= i*mCartesianOdeScaleXY[1])
                {
                    labelYIndex = i-1;
                    break;
                }
            }

            // for the x position

            for(unsigned i=1; i<(mCartesianOdeDimensions[0]+1); i++)
            {
                if(positionX <= i*mCartesianOdeScaleXY[0])
                {
                    labelXIndex = i-1;
                    break;
                }
            }

            break;

        case 3:

            // for the y position

            for(unsigned i=1; i<(mCartesianOdeDimensions[1]+1); i++) // top to bottom
            {
                if(positionY <= i*mCartesianOdeScaleXY[1])
                {
                    labelYIndex = i-1;
                    break;
                }
            }

            // for the x position

            for(unsigned i=1; i<(mCartesianOdeDimensions[0]+1); i++)
            {
                if(positionX <= i*mCartesianOdeScaleXY[0])
                {
                    labelXIndex = i-1;
                    break;
                }
            }

            break;

        case 1:
            std::cout<<"Error: Abstract Domain 1D not supported"<<std::endl;
            break;

    }

    

    return mOdeLabels[labelZIndex][labelYIndex][labelXIndex];
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
std::vector<std::vector<std::string>> AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::ReadMatrix(std::string filename)
{
    // parse a matrix file (.csv) line by line, ignore escape line,s containing file information
    // that is lines starting with '#' 
    
    std::string line;
    std::ifstream inputFile(filename);
    
    // read all data types as std::string therefore return the matrix of strings for personalised
    // methods down stream
    std::vector<std::vector<std::string>> outputMatrix = std::vector<std::vector<std::string>>();

    // check file exists and is openable
    if(inputFile.is_open()){
        // open the matrix file
        while (getline(inputFile,line)){
            // while the file still has lines not read.
            // read line left to right, top to bottom.
            if(!line.empty()){
                if(line.at(0)=='#')
                {
                    //std::cout<<"Escape line: "<<line<<std::endl;
                }
                else
                {
                    outputMatrix.push_back(parseMatrixLineString(line));
                }   
            }
        }
        inputFile.close();

        return outputMatrix;
    }else{
        std::cout<<"Error: Unable to open file: "<<filename<<std::endl;
        return outputMatrix;
    }
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
std::vector<std::string> AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::parseMatrixLineString(std::string line)
{
    // for a line string in the matrix read, parse into vector data entries based on delimiters ','
    std::vector<std::string> rowVector = std::vector<std::string>();

    // delimiter, may be modified by further methods
    std::string delim = ",";
    std::string matrixCell;

    // determine the position of the delimiter
    size_t posSnew=line.find(delim);

    bool IsEndOfLine = false;
    
    while(!IsEndOfLine)
    {
        // while not at the end of the file, sample sub strings from the posiiton of the delimiter
        if(posSnew == std::string::npos)
        {
            IsEndOfLine = true;
        }
        
        // sample substring from begining of the string to the delimiter positioon, store as data entry
        matrixCell = line.substr(0,posSnew);

        // remove the sampled entry from the string
        line = line.substr(posSnew+1,std::string::npos);

        rowVector.push_back(matrixCell);

        // update delimiter position
        posSnew=line.find(delim);
    }
    return rowVector;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
void AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::ReadDomainLabels()
{
    // read and store the information contained within the domain labels files
    std::vector<std::vector<std::vector<std::string>>> tempDomainLabels;

    // mDomainLabelFilename will eb a .csv if 2D or a directory containg .csv files if 3D
    std::string filename = mDomainLabelFilename;
    switch(SPACE_DIM)
    {
        case 2:
            tempDomainLabels.push_back(ReadMatrix(filename));
            break;

        case 3: 
            for(unsigned i=0; i<mNumberOfSlices;i++)
            {
                filename = mDomainLabelFilename + "/" + std::to_string(i) + ".csv";
                tempDomainLabels.push_back(ReadMatrix(filename));
            }
            break;

        case 1:
            std::cout<<"Error: Abstract Domain 1D not supported"<<std::endl;
            break;

    }

    SetDomainLabels(tempDomainLabels);

    SetDomainLabelVector(ReturnUnique(tempDomainLabels));
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
void AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::ReadDomainKey()
{
    // read the label keys for the domain, providing more information as to the type of the sun domain
    std::vector<std::vector<std::string>> domainKeys;

    domainKeys = ReadMatrix(mDomainKeyFilename);

    SetDomainKeys(domainKeys);
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
void AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::ReadOdeLabels()
{
    // read and store the ode labels, form a vector of the unique labels to store as a reference vector    
    std::vector<std::string> odeLabelVector;

    std::vector<std::vector<std::vector<std::string>>> tempOdeLabels;

    // mOdeLabelFilename will eb a .csv if 2D or a directory containg .csv files if 3D
    std::string filename = mOdeLabelFilename;
    switch(SPACE_DIM)
    {
        case 2:
            tempOdeLabels.push_back(ReadMatrix(filename));
            break;

        case 3: 
            for(unsigned i=0; i<mNumberOfSlices;i++)
            {
                filename = mOdeLabelFilename + "/" + std::to_string(i) + ".csv";
                tempOdeLabels.push_back(ReadMatrix(filename));
            }
            break;

        case 1:
            std::cout<<"Error: Abstract Domain 1D not supported"<<std::endl;
            break;

    }

    odeLabelVector = ReturnUnique(tempOdeLabels);

    SetOdeLabelVector(odeLabelVector);
    SetOdeLabels(tempOdeLabels);
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
void AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::ReadOdeKey()
{
    // read the ode keys from a .csv file, providing the label used within the domain and
    // the respective ode to be used
    std::vector<std::vector<std::string>> odeKeys;

    odeKeys = ReadMatrix(mOdeKeyFilename);

    //printMatrix(odeKeys);
    SetOdeKeys(odeKeys);
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
void AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::ReadDiffusionDatabase()
{
    // read and store the diffusion data as a matrix of strings 
    std::vector<std::vector<std::string>> diffusionDatabase;

    diffusionDatabase = ReadMatrix(mDiffusionFilename);

    //printMatrix(diffusionDatabase);

    // vector to store the state variables
    std::vector<std::string> stateVector;
    std::vector<std::string> domainVector;
    for(unsigned i=0; i<diffusionDatabase.size(); i++)
    {
        // the first element of the read in database is the identifying names
        // of the state variables that diffuse
        stateVector.push_back(diffusionDatabase[i][0]);
        if(diffusionDatabase[0].size()>2)
        {
            domainVector.push_back(diffusionDatabase[i][1]);
        }
    }

    SetNumberOfDomains(ReturnUnique(domainVector).size());
    // for a new state variable registers with the unique identifying names from the diffusion file
    StateVariableRegister*   p_stateVariableVector = new StateVariableRegister(ReturnUnique(stateVector));
    unsigned numberOfStates = p_stateVariableVector -> GetNumberOfStateVariables();
    

    // store the diffusion information
    SetStateVariableVector(p_stateVariableVector);

    SetNumberOfStates(numberOfStates);

    SetDiffusionDatabase(diffusionDatabase);
}


// auxillary methods

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
std::vector<std::string> AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::ReturnUnique(std::vector<std::string> candidateVector)
{
    // serach through a vector of strings and store the first occurance of each unique string
    std::vector<std::string> resultUnique = std::vector<std::string>();
    resultUnique.push_back(candidateVector.at(0));
    unsigned uniqueCount =1;
    bool IsFound = false;
    for(unsigned i=1; i<candidateVector.size(); i++)
    {
        IsFound = false;
        for(unsigned j=0; j<uniqueCount; j++)
        {

            if(candidateVector[i] == resultUnique[j])
            {
                IsFound =true;
                break;
            }
        }
        if(!IsFound)
        {
            resultUnique.push_back(candidateVector[i]);
            uniqueCount++;
        }
    }
    return resultUnique;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
std::vector<std::string> AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::ReturnUnique(std::vector<std::vector<std::string>> candidateMatrix)
{
    // serach through a matrix of strings and store the first occurance of each unique string
    std::vector<std::string> resultUnique = std::vector<std::string>();
    resultUnique.push_back(candidateMatrix[0][0]);
    unsigned uniqueCount =1;
    bool IsFound = false;

    for(unsigned j=0; j<candidateMatrix.size(); j++)
    {
        for(unsigned i=0; i<candidateMatrix[j].size(); i++)
        {
            IsFound = false;
            for(unsigned k=0; k<uniqueCount; k++)
            {
                if(candidateMatrix[j][i] == resultUnique[k])
                {
                    IsFound =true;
                    break;
                }
            }
            if(!IsFound)
            {
                resultUnique.push_back(candidateMatrix[j][i]);
                uniqueCount++;
            }
        }
    }
    return resultUnique;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
std::vector<std::string> AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::ReturnUnique(std::vector<std::vector<std::vector<std::string>>> candidateMatrix)
{
    // serach through a matrix of strings and store the first occurance of each unique string
    std::vector<std::string> resultUnique = std::vector<std::string>();
    resultUnique.push_back(candidateMatrix[0][0][0]);
    unsigned uniqueCount =1;
    bool IsFound = false;

    for(unsigned j=0; j<candidateMatrix.size(); j++)
    {
        for(unsigned i=0; i<candidateMatrix[j].size(); i++)
        {
            for(unsigned k=0; k<candidateMatrix[j][i].size();k++)
            {
                IsFound = false;
                for(unsigned l=0; l<uniqueCount; l++)
                {
                    if(candidateMatrix[j][i][k] == resultUnique[l])
                    {
                        IsFound =true;
                        break;
                    }
                }
                if(!IsFound)
                {
                    resultUnique.push_back(candidateMatrix[j][i][k]);
                    uniqueCount++;
                }
            }
            
        }
    }
    return resultUnique;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
void AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::printMatrix(std::vector<std::vector<std::vector<std::string>>> matrix)
{
    for (unsigned int k = 0; k < matrix.size(); k++)
    {
        for (unsigned int i = 0; i < matrix.size(); i++)
        {
            for (unsigned int j = 0; j < matrix[i].size(); j++)
            {
            std::cout << matrix[i][j][k]<< ' ';
            }
            std::cout<<std::endl;
        }
        std::cout<<std::endl;
    }
    return;
};

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
void AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::printMatrix(std::vector<std::vector<std::string>> matrix)
{
  for (unsigned int i = 0; i < matrix.size(); i++)
  {
    for (unsigned int j = 0; j < matrix[i].size(); j++)
    {
      std::cout << matrix[i][j]<< ' ';
    }
    std::cout<<std::endl;
  }
  std::cout<<std::endl;

    return;
};

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
void AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::printVector(std::vector<std::string> vec)
{
  for (unsigned int i = 0; i < vec.size(); i++)
  {
    std::cout << vec[i]<< ' ';
  }
  std::cout<<std::endl;

    return;
};

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
void AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::PrintDiffusionDomain()
{
    std::cout<<"Diffusion domain:"<<std::endl;
    printMatrix(mDomainLabels);
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
void AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::PrintMappedDiffusionDomain()
{
    std::cout<<"Mapped diffusion domain:"<<std::endl;
    printMatrix(mCartesianChasteDomain);
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
void AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::PrintDomainLabelKeys()
{
    std::cout<<"Domain label keys:"<<std::endl;
    printMatrix(mDomainKeys);
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
void AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::PrintODEDomain()
{
    std::cout<<"ODE domain:"<<std::endl;
    printMatrix(mOdeLabels);
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
void AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::PrintMappedODEDomain()
{
    std::cout<<"Mapped ODE domain:"<<std::endl;
    printMatrix(mCartesianOdeDomain);
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
void AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::PrintODELabelKeys()
{
    std::cout<<"ODE label keys:"<<std::endl;
    printMatrix(mOdeKeys);
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
void AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::PrintDiffusionDatabase()
{
    std::cout<<"Diffusion database:"<<std::endl;
    printMatrix(mDiffusionDatabase);
}

// set methods

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
void AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetNumberOfDomains(unsigned numberOfDomains)
{
    mNumberOfDomains = numberOfDomains;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
unsigned AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::GetNumberOfDomains()
{
    return mNumberOfDomains;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
void AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetDomainLabels(std::vector<std::vector<std::vector<std::string>>> domainLabels)
{
    mDomainLabels = domainLabels;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
void AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetDomainLabelVector(std::vector<std::string> domainLabelVector)
{
    mDomainLabelVector = domainLabelVector;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
void AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetDomainKeys(std::vector<std::vector<std::string>> domainKeys)
{
    mDomainKeys = domainKeys;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
void AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetNumberOfOdeSystems(unsigned numberSystems)
{
    mNumberOfOdeSystems = numberSystems;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
void AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetOdeLabels(std::vector<std::vector<std::vector<std::string>>>  odeLabels)
{
    mOdeLabels = odeLabels;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
void AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetOdeLabelVector(std::vector<std::string> labelVector)
{
    mOdeLabelVector = labelVector;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
void AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetOdeKeys(std::vector<std::vector<std::string>> odeKeys)
{
    mOdeKeys = odeKeys;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
void AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetStateVariableVector(StateVariableRegister* stateVector)
{
    mpStateVariableVector = stateVector;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
void AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetNumberOfStates(unsigned numberStates)
{
    mNumberOfStates = numberStates;
} 

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
void AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetCartesianCellDimensions(std::vector<unsigned> cellDims)
{
    mCartesianCellDimensions = cellDims;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
void AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetCartesianOdeDimensions(std::vector<unsigned> odeDims)
{
    mCartesianOdeDimensions = odeDims;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
void AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetCartesianCellScale(std::vector<double> cellScale)
{
    mCartesianCellScaleXY = cellScale;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
void AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetCartesianOdeScale(std::vector<double> odeScale)
{
    mCartesianOdeScaleXY = odeScale;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
void AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetCartesianChasteScale(std::vector<double> chasteScale)
{
    mCartesianChasteScaleXY = chasteScale;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
void AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetCartesianChasteDomain(std::vector<std::vector<std::vector<std::string>>> cartesianChasteDomain)
{
    mCartesianChasteDomain = cartesianChasteDomain;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
void AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetCartesianOdeDomain(std::vector<std::vector<std::vector<std::string>>> cartesianOdeDomain)
{
    mCartesianOdeDomain = cartesianOdeDomain;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
void AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetCartesianChasteDimensions(std::vector<unsigned> dimensions)
{
    mCartesianChasteDimensions = dimensions;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
void AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetMeshDimensions(std::vector<unsigned> meshDimensions)
{
    mMeshDimensions = meshDimensions;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
void AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetDomainFeMesh(TetrahedralMesh<ELEMENT_DIM,SPACE_DIM>*& p_mesh)
{
    mpFeMesh = p_mesh;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
void AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetMeshGenerator(HoneycombMeshGenerator* p_generator)
{
    mpGenerator = p_generator;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
void AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetMeshScale(std::vector<double> meshScale)
{
    mMeshScale = meshScale;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
void AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetDiffusionDatabase(std::vector<std::vector<std::string>> diffusionDatabase)
{
    mDiffusionDatabase = diffusionDatabase;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
void AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetNodeLabels(std::vector<std::string> nodeLabels)
{
    mLabelledNodes = nodeLabels;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
void AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetNodeDomains(std::vector<std::string> nodeDomains)
{
    mDomainNodes = nodeDomains;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
void AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetInitialNodeConditions(std::vector<double> initialNodeConditions)
{
    mInitalNodeConditions = initialNodeConditions;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
void AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetBoundaryConditionTypes(std::vector<std::string> boundaryConditionsTypes)
{
    mBoundaryConditionTypes = boundaryConditionsTypes;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
void AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetBoundaryConditionValues(std::vector<double> boundaryConditionValues)
{
    mBoundaryConditionValues = boundaryConditionValues;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
void AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetLabelOrigin(std::vector<double> origin)
{
    mLabelOrigin = origin;
}

// get methods

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
std::vector<std::vector<std::vector<std::string>>> AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::GetDomainLabels()
{
    return mDomainLabels;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
std::string AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::GetDomainLabelByIndex(unsigned index)
{
    if(index<mNumberOfDomains)
    {
        return mDomainLabelVector[index];
    }
    else
    {
        std::cout<<"Error: AbstractDomainFieldTemplated::GetDomainLabelByIndex: index out of bounds"<<std::endl;
        return "Error";
    }
    
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
std::vector<std::string> AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::GetDomainLabelVector()
{
    return mDomainLabelVector;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
std::vector<std::vector<std::string>> AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::GetDomainKeys()
{
    return mDomainKeys;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
unsigned AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::GetNumberOfOdeSystems()
{
    return mNumberOfOdeSystems;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
std::vector<std::vector<std::vector<std::string>>>  AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::GetOdeLabels()
{
    return mOdeLabels;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
std::vector<std::string> AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::GetOdeLabelVector()
{
    return mOdeLabelVector;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
std::vector<std::vector<std::string>> AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::GetOdeKeys()
{
    return mOdeKeys;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
StateVariableRegister* AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::GetStateVariableVector()
{
    return mpStateVariableVector;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
unsigned AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::GetNumberOfStates()
{
    return mNumberOfStates;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
std::vector<unsigned> AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::GetCartesianCellDimensions()
{
    return mCartesianCellDimensions;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
std::vector<unsigned> AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::GetCartesianOdeDimensions()
{
    return mCartesianOdeDimensions;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
std::vector<double> AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::GetCartesianCellScale()
{
    return mCartesianCellScaleXY;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
std::vector<double> AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::GetCartesianOdeScale()
{
    return mCartesianOdeScaleXY;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
std::vector<double> AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::GetCartesianChasteScale()
{
    return mCartesianChasteScaleXY;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
std::vector<std::vector<std::vector<std::string>>> AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::GetCartesianChasteDomain()
{
    return mCartesianChasteDomain;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
std::vector<std::vector<std::vector<std::string>>> AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::GetCartesianOdeDomain()
{
    return mCartesianOdeDomain;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
std::vector<unsigned> AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::GetCartesianChasteDimensions()
{
    return mCartesianChasteDimensions;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
std::vector<unsigned> AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::GetMeshDimensions()
{
    return mMeshDimensions;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
TetrahedralMesh<ELEMENT_DIM,SPACE_DIM>*& AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::rGetDomainFeMesh()
{
    return mpFeMesh;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
HoneycombMeshGenerator* AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::GetMeshGenerator()
{
    return mpGenerator;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
std::vector<double> AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::GetMeshScale()
{
    return mMeshScale;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
std::vector<std::vector<std::string>> AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::GetDiffusionDatabase()
{
    return mDiffusionDatabase;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
std::vector<std::string> AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::GetNodeLabels()
{
    return mLabelledNodes;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
std::vector<std::string> AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::GetNodeDomains()
{
    return mDomainNodes;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
std::vector<double> AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::GetInitialNodeConditions()
{
    return mInitalNodeConditions;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
std::vector<std::string> AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::GetBoundaryConditionTypes()
{
    return mBoundaryConditionTypes;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
std::vector<double> AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::GetBoundaryConditionValues()
{
    return mBoundaryConditionValues;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
std::vector<double> AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::GetLabelOrigin()
{
    return mLabelOrigin;
}


#endif