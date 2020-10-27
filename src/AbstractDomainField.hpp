#ifndef ABSTRACTDOMAINFIELD_HPP
#define ABSTRACTDOMAINFIELD_HPP

#include <string>
#include <vector>
#include "ChastePoint.hpp"
#include "TetrahedralMesh.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "AbstractTetrahedralMesh.hpp"
#include "StateVariableRegister.hpp"

/*
    Class to hold the information about the spatial domain, mesh, domain labels, ode types etc
*/
class AbstractDomainField
{
protected:

    // file names containing the system information
    std::string mDomainLabelFilename;

    std::string mDomainKeyFilename;

    std::string mOdeLabelFilename;

    std::string mOdeKeyFilename;

    std::string mDiffusionFilename;


    // domain information input
    unsigned mNumberOfDomains;

    std::vector<std::vector<std::string>> mDomainLabels;

    std::vector<std::string> mDomainLabelVector;

    std::vector<std::vector<std::string>> mDomainKeys;

    std::vector<unsigned> mCartesianCellDimensions; // (x, y) of the input csv file grid

    std::vector<double> mCartesianCellScaleXY; // (Sx, Sy) scale factor of the input csv grid, default (1,1)



    // ode system information input
    unsigned mNumberOfOdeSystems;

    std::vector<std::vector<std::string>> mOdeLabels;

    std::vector<std::string> mOdeLabelVector;

    std::vector<std::vector<std::string>> mOdeKeys;

    std::vector<double> mCartesianOdeScaleXY; // (Sx, Sy) scale factor of the input csv grid, default (1,1)

    std::vector<unsigned> mCartesianOdeDimensions; // (x, y) of the input csv file grid



    // collective system information
    unsigned mNumberOfStates; // numebr of unique species in the diffusion database

    StateVariableRegister* mStateVariableVector; // container of each unique species in the diffusion database

    std::vector<std::vector<std::string>> mDiffusionDatabase;

    unsigned mDiffusionDomainLabelPosition =1;



    // mapped cartesian domain 
    std::vector<unsigned> mCartesianChasteDimensions;

    std::vector<double> mCartesianChasteScaleXY; // (Sx, Sy) scale factor of the Chaste grid

    std::vector<std::vector<std::string>> mCartesianChasteDomain;

    std::vector<std::vector<std::string>> mCartesianOdeDomain;


    // mapped mesh domain
    std::vector<unsigned> mMeshDimensions;

    std::vector<double> mMeshScale;

    MutableMesh<2,2>* mpMesh;

    HoneycombMeshGenerator* mpGenerator;

    std::vector<std::string> mLabelledNodes;

    std::vector<std::string> mDomainNodes;

    // initial and boundary conditions
    std::string mInitialConditionsFilename;

    std::string mBoundaryConditionsFilename;

    std::vector<double> mInitalNodeConditions;

    std::vector<std::string> mBoundaryConditionTypes;

    std::vector<double> mBoundaryConditionValues;




    // system properties
    unsigned mProbDim =0; // to be overridden for the class object
    const unsigned mSpaceDim=2;
    const unsigned mElementDim=2;



public:

    AbstractDomainField(std::string domainLabelFilename="", 
                        std::string domainKeyFilename="", 
                        std::string odeLabelFilename="", 
                        std::string odeKeyFilename="", 
                        std::string diffusionFilename="");

    virtual ~AbstractDomainField()
    {
    }

    // virtual methods involved with basic function
    virtual void SetupAndInitialiseLabelDomain();

    virtual void SetupAndInitialiseOdes();

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
        return "AbstractDomainField";
    };

    // display methods

    void PrintDiffusionDomain();

    void PrintMappedDiffusionDomain();

    void PrintDomainLabelKeys();

    void PrintODEDomain();

    void PrintMappedODEDomain();

    void PrintODELabelKeys();

    void PrintDiffusionDatabase();

    // return methods
    std::string ReturnNodeOdeLabelAtPosition(const c_vector<double,2>&);

    std::string ReturnDomainLabelAtPosition(const c_vector<double,2>&);

    std::string ReturnMappedDomainLabel(unsigned, unsigned, std::vector<double>);

    std::string ReturnMappedOdeLabel(unsigned, unsigned,std::vector<double>);

   
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

    void printMatrix(std::vector<std::vector<std::string>>);

    void printVector(std::vector<std::string>);

    bool PerturbInitialConditionTest(std::vector<std::string>);

    // set methods
    void SetNumberOfDomains(unsigned);

    void SetDomainLabels(std::vector<std::vector<std::string>>);

    void SetDomainLabelVector(std::vector<std::string>);

    void SetDomainKeys(std::vector<std::vector<std::string>>);    

    void SetNumberOfOdeSystems(unsigned);

    void SetOdeLabels(std::vector<std::vector<std::string>>);

    void SetOdeLabelVector(std::vector<std::string>);    

    void SetOdeKeys(std::vector<std::vector<std::string>>);

    void SetStateVariableVector(StateVariableRegister*);

    void SetNumberOfStates(unsigned);    
        
    void SetCartesianCellDimensions(std::vector<unsigned>);

    void SetCartesianCellScale(std::vector<double>);

    void SetCartesianOdeDimensions(std::vector<unsigned>);

    void SetCartesianOdeScale(std::vector<double>);
    
    void SetCartesianChasteScale(std::vector<double>);

    void SetCartesianChasteDomain(std::vector<std::vector<std::string>>);

    void SetCartesianOdeDomain(std::vector<std::vector<std::string>>);

    void SetCartesianChasteDimensions(std::vector<unsigned>);

    void SetMeshDimensions(std::vector<unsigned>);

    void SetDomainMesh(MutableMesh<2,2>*);

    void SetMeshGenerator(HoneycombMeshGenerator*);

    void SetMeshScale(std::vector<double>);

    void SetDiffusionDatabase(std::vector<std::vector<std::string>>);

    void SetNodeLabels(std::vector<std::string>);

    void SetNodeDomains(std::vector<std::string>);

    void SetInitialNodeConditions(std::vector<double>);

    void SetBoundaryConditionTypes(std::vector<std::string>);

    void SetBoundaryConditionValues(std::vector<double>);


    // get methods
    unsigned GetNumberOfDomains();

    std::vector<std::vector<std::string>> GetDomainLabels();

    std::vector<std::string> GetDomainLabelVector();

    std::string GetDomainLabelByIndex(unsigned);

    std::vector<std::vector<std::string>> GetDomainKeys();

    unsigned GetNumberOfOdeSystems();

    std::vector<std::vector<std::string>> GetOdeLabels();

    std::vector<std::string> GetOdeLabelVector();

    std::vector<std::vector<std::string>> GetOdeKeys();

    StateVariableRegister* GetStateVariableVector();

    unsigned GetNumberOfStates();

    std::vector<unsigned> GetCartesianCellDimensions();

    std::vector<double> GetCartesianCellScale();

    std::vector<unsigned> GetCartesianOdeDimensions();

    std::vector<double> GetCartesianOdeScale();

    std::vector<std::vector<std::string>> GetCartesianChasteDomain();

    std::vector<double> GetCartesianChasteScale();

    std::vector<std::vector<std::string>> GetCartesianOdeDomain();

    std::vector<unsigned> GetCartesianChasteDimensions();

    std::vector<unsigned> GetMeshDimensions();

    MutableMesh<2,2>* GetDomainMesh();

    HoneycombMeshGenerator* GetMeshGenerator();

    std::vector<double> GetMeshScale();

    std::vector<std::vector<std::string>> GetDiffusionDatabase();

    std::vector<std::string> GetNodeLabels();

    std::vector<std::string> GetNodeDomains();

    std::vector<double> GetInitialNodeConditions();

    std::vector<std::string> GetBoundaryConditionTypes();

    std::vector<double> GetBoundaryConditionValues();



    std::string ReturnKeyValueFromNodeLabel(std::string);

    std::string ReturnDomainKeyFromDomainLabel(std::string);

};

#endif