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

    // Fe mesh scaling
    bool mIsScaleBy  =true;

    double mTargetWidth;

    double mTargetHeight;

    double mTargetDepth;

    double mScaleWidth;

    double mScaleHeight;

    double mScaleDepth;

    std::string mInitialConditionsFilename;

    std::string mBoundaryConditionsFilename;

    boost::shared_ptr<BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM> > mpBoundaryConditionsContainer; 

    std::vector<bool> mPiecewiseBCs;

public:

    ChemicalDomainFieldTemplated(
                        //TetrahedralMesh<ELEMENT_DIM,SPACE_DIM>*,
                        std::string reactionFileRoot="",
                        std::string domainLabelFilename="", 
                        std::string domainKeyFilename="", 
                        std::string odeLabelFilename="", 
                        std::string odeKeyFilename="",
                        std::string diffusionFilename="",
                        std::string initialConditionsFilename = "",
                        std::string boundaryConditionsFilename = "",
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

    virtual boost::shared_ptr<BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>> ProcessBoundaryConditions();

    virtual void ParseBoundaryConditionsFromFile(std::string);

    void FormReactionSystemAtNodes();

    void DeriveSystemProperties();

    void DeriveExtendedSystemProperties();

    void FeMeshScaling( double targetWidth = 0.0,
                        double targetHeight= 0.0,
                        double targetDepth= 0.0,
                        double scaleWidth  = 1.0,
                        double scaleHeight = 1.0,
                        double scaleDepth = 1.0
    );

    void ScaleFeMesh();

    // interface methods

    virtual double GetDiffusionValueBasedOnPoint(const ChastePoint<2>&, unsigned);

    virtual double ReturnDiffusionValueFromStateNameAndDomainLabel(std::string stateName, std::string domainLabel = "");

    virtual std::string GetFieldType()
    {
        return "ChemicalDomainFieldTemplated";
    };

    boost::shared_ptr<BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM> > ReturnSharedPtrBoundaryConditionsContainer();

    // set methods
    void SetDomainStateVariableRegister(StateVariableRegister*);

    void SetChemistry(AbstractDiffusiveChemistry*);

    void SetOdeSystem(std::vector<AbstractInhomogenousChemicalOdeSystemForCoupledPdeSystem*>);

    void SetPiecewiseBCs(std::vector<bool>);

    void SetBoundaryConditionsContainer(BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM> bcc);

    // get methods
    StateVariableRegister* GetDomainStateVariableRegister();

    AbstractDiffusiveChemistry* GetChemistry();

    std::vector<AbstractInhomogenousChemicalOdeSystemForCoupledPdeSystem*> GetOdeSystem();

    std::vector<bool> GetPiecewiseBCs();

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
                                            std::string initialConditionsFilename,
                                            std::string boundaryConditionsFilename,
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
        mReactionFileRoot(reactionFileRoot),
        mInitialConditionsFilename(initialConditionsFilename),
        mBoundaryConditionsFilename(boundaryConditionsFilename),
        mpBoundaryConditionsContainer(boost::shared_ptr<BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>>())
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
    
    if(mInitialConditionsFilename != "")
    {
        // process initial conditions
        this -> ParseInitialConditionsFromFile(mInitialConditionsFilename);
    }
    

    if(mBoundaryConditionsFilename != "")
    {
        // process boundary conditions
        this -> ParseBoundaryConditionsFromFile(mBoundaryConditionsFilename);

        mpBoundaryConditionsContainer = ProcessBoundaryConditions();
    }

    // populate abstract domain with a chemistry
    // map odes systems to nodes 
    FormReactionSystemAtNodes();
    // update properties of the domain field
    //DeriveSystemProperties();
    DeriveExtendedSystemProperties();
}


template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
void ChemicalDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::ParseBoundaryConditionsFromFile(std::string boundaryConditionsFilename)
{
    // read the input boundary conditions file as a matrix of strings, first column state name, 
    // 2nd column boundary condition type, 3rd column value.
    std::vector<std::vector<std::string>> fullBoundaryConditionsAsStrings = this->ReadMatrix(boundaryConditionsFilename);
    std::vector<std::vector<std::string>> boundaryConditionsAsStrings = this->ReturnStateDataWithinRegister(fullBoundaryConditionsAsStrings);
    // use state name against stateVariableRegister to remove any state value not in the system
    // replace any state in register but not in boundary conditions with value of 0.0, type Derichlet

    unsigned numberOfStateVariables = this -> GetStateVariableVector() ->GetNumberOfStateVariables();

    std::vector<std::string> boundaryConditionTypes;
    std::vector<bool> piecewiseBCs(numberOfStateVariables,false);
    std::vector<double> boundaryConditionValues;
    std::string BCtype ="";
    for(unsigned pdeDim=0; pdeDim<numberOfStateVariables; pdeDim++)
    {   
        // determine which state the pde dimension refers to and then retrieve the condition if present
        std::string stateName = this -> mpStateVariableVector -> RetrieveStateVariableName(pdeDim);

        for(unsigned inputState=0; inputState<boundaryConditionsAsStrings.size();inputState++)
        {
            std::cout<<"input: "<<boundaryConditionsAsStrings[inputState][0]<<std::endl;
            // test whether the state variable is specified within the input data
            if(stateName==boundaryConditionsAsStrings[inputState][0])
            {
                // state found, read data
                if(boundaryConditionsAsStrings[inputState].size()>2)
                {
                    unsigned numberOfVariables = static_cast<unsigned>((boundaryConditionsAsStrings[inputState].size()-1)/2);
                    std::cout<<"numberOfVariables: "<<numberOfVariables<<std::endl;
                    if(numberOfVariables>1)
                    {
                        piecewiseBCs[pdeDim]= true;
                    }
                    for(unsigned i=0; i<numberOfVariables;i++)
                    {
                        // type value is specified also
                        BCtype = boundaryConditionsAsStrings[inputState][1+2*i];
                        std::cout<<"BC type: "<<BCtype<<" BCValue: "<<boundaryConditionsAsStrings[inputState][2+2*i]<<std::endl;
                        // standardise the input
                        if((BCtype=="Dirichlet"||BCtype=="dirichlet"||BCtype=="D"||BCtype=="d"))
                        {
                            BCtype="Dirichlet";
                        }
                        else if((BCtype=="Neumann"||BCtype=="neumann"||BCtype=="N"||BCtype=="n"))
                        {
                            BCtype="Neumann";
                        }
                        else
                        {
                            std::cout<<"Error: ChemicalDomainFieldForCellCoupling::ParseBoundaryConditionsFromFile BC value missing - "<<BCtype<<std::endl;
                        }
                        boundaryConditionTypes.push_back(BCtype);
                        boundaryConditionValues.push_back(std::stod(boundaryConditionsAsStrings[inputState][2+2*i]));
                    }
                }
                else
                {
                    // type value not specified, assume Dirichlet
                    boundaryConditionValues.push_back(std::stod(boundaryConditionsAsStrings[inputState][1]));
                    boundaryConditionTypes.push_back("Dirichlet");
                    piecewiseBCs[pdeDim] = false;

                }
                break;
            }
        }
        // if state isn't found, leave as default, 0.0 Dirichlet
    }

    this -> SetBoundaryConditionTypes(boundaryConditionTypes);

    this -> SetBoundaryConditionValues(boundaryConditionValues);

    SetPiecewiseBCs(piecewiseBCs);
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
boost::shared_ptr<BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>> ChemicalDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::ProcessBoundaryConditions()
{

    // retrieive the boundary conditions information derived from the files inputs
    std::vector<std::string> boundaryConditionTypes = this -> GetBoundaryConditionTypes();
    std::vector<double> boundaryConditionValues = this -> GetBoundaryConditionValues();
    std::vector<bool> piecewiseBCs = this -> GetPiecewiseBCs();



    // create containeers to store the boundary conditions, assume each boundary condition is constant over the simulation
    boost::shared_ptr<BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>> p_bcc(new BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>());
    std::vector<ConstBoundaryCondition<SPACE_DIM>*> vectorConstBCs;

    for(unsigned pdeDim=0; pdeDim<boundaryConditionValues.size(); pdeDim++){
        vectorConstBCs.push_back(new ConstBoundaryCondition<SPACE_DIM>(boundaryConditionValues[pdeDim]));
    }



    for(unsigned pdeDim=0; pdeDim<PROBLEM_DIM; pdeDim++)
    {
        if(piecewiseBCs[pdeDim])
        {
            double maxNodeX=0.0;
            double maxNodeY=0.0;
            double maxNodeZ=0.0;
            double maxElementNodeX=0.0;
            double maxElementNodeY=0.0;
            double maxElementNodeZ=0.0;

            double prevNodeX=0.0;
            double prevNodeY=0.0;

            double distanceBetweenNodesX=10.0;
            double distanceBetweenNodesY=10.0;

            switch (SPACE_DIM)
            {
                case 2:

               

                    // run through Boundary nodes
                    for (typename TetrahedralMesh<ELEMENT_DIM,SPACE_DIM>::BoundaryNodeIterator node_iter = this->mpFeMesh->GetBoundaryNodeIteratorBegin();
                    node_iter != this->mpFeMesh->GetBoundaryNodeIteratorEnd();
                    ++node_iter)
                    {
                        double x = (*node_iter)->GetPoint()[0];
                        double y = (*node_iter)->GetPoint()[1];
                        if(x>maxNodeX)
                        {
                            maxNodeX = x;
                        }
                        if(y>maxNodeY)
                        {
                            maxNodeY = y;
                        }

                        if(fabs(prevNodeX-x)<distanceBetweenNodesX)
                        {
                            distanceBetweenNodesX=fabs(prevNodeX-x);
                        }

                        if(fabs(prevNodeY-y)<distanceBetweenNodesY)
                        {
                            distanceBetweenNodesY=fabs(prevNodeY-y);
                        }

                        prevNodeX=x;
                        prevNodeY=y;
                    }

                    // run through Neumann BCs
                    for (typename TetrahedralMesh<ELEMENT_DIM,SPACE_DIM>::BoundaryElementIterator boundary_iter = this->mpFeMesh->GetBoundaryElementIteratorBegin();
                    boundary_iter != this->mpFeMesh->GetBoundaryElementIteratorEnd();
                    boundary_iter++)
                    {
                        unsigned node_index = (*boundary_iter)->GetNodeGlobalIndex(0);
                        double x = this->mpFeMesh->GetNode(node_index)->GetPoint()[0];
                        double y = this->mpFeMesh->GetNode(node_index)->GetPoint()[1];
                        if(x>maxElementNodeX)
                        {
                            maxElementNodeX = x;
                        }
                        if(y>maxElementNodeY)
                        {
                            maxElementNodeY = y;
                        }


                    }
                    // run through Dirichlet BCs
                    for (typename TetrahedralMesh<ELEMENT_DIM,SPACE_DIM>::BoundaryNodeIterator node_iter = this->mpFeMesh->GetBoundaryNodeIteratorBegin();
                    node_iter != this->mpFeMesh->GetBoundaryNodeIteratorEnd();
                    ++node_iter)
                    {

                        double x = (*node_iter)->GetPoint()[0];
                        double y = (*node_iter)->GetPoint()[1];
                        if (x==0)
                        {
                            if(boundaryConditionTypes[pdeDim+0]=="Dirichlet")
                            {
                                p_bcc->AddDirichletBoundaryCondition(*node_iter, vectorConstBCs[pdeDim+0], pdeDim);
                            }
                        }
                        else if(y==0)
                        {
                            if(boundaryConditionTypes[pdeDim+1]=="Dirichlet")
                            {
                                p_bcc->AddDirichletBoundaryCondition(*node_iter, vectorConstBCs[pdeDim+1], pdeDim);
                            }
                        }
                        //else if(fabs(x-this -> mMeshDimensions[0])<1e-6)
                        else if(fabs(x-maxNodeX-distanceBetweenNodesX)<1e-6)
                        {
                            if(boundaryConditionTypes[pdeDim+2]=="Dirichlet")
                            {
                                p_bcc->AddDirichletBoundaryCondition(*node_iter, vectorConstBCs[pdeDim+2], pdeDim);
                            }
                        }
                        //else if(fabs(y-this -> mMeshDimensions[1])<1e-6)
                        else if(fabs(y-maxNodeY-distanceBetweenNodesY)<1e-6)
                        {
                            if(boundaryConditionTypes[pdeDim+3]=="Dirichlet")
                            {
                                p_bcc->AddDirichletBoundaryCondition(*node_iter, vectorConstBCs[pdeDim+3], pdeDim);
                            }
                        }
                        else
                        {
                            std::cout<<"Error: ChemicalDomainFieldForCellCoupling::ProcessBoundaryConditions() Dirichlet node not on boundary"<<std::endl;
                        }
                
                        node_iter++;
                    } 
                    // run through Neumann BCs
                    for (typename TetrahedralMesh<ELEMENT_DIM,SPACE_DIM>::BoundaryElementIterator boundary_iter = this->mpFeMesh->GetBoundaryElementIteratorBegin();
                    boundary_iter != this->mpFeMesh->GetBoundaryElementIteratorEnd();
                    boundary_iter++)
                    {
                        unsigned node_index = (*boundary_iter)->GetNodeGlobalIndex(0);
                        double x = this->mpFeMesh->GetNode(node_index)->GetPoint()[0];
                        double y = this->mpFeMesh->GetNode(node_index)->GetPoint()[1];
    
                        if (x==0)
                        {
                            if(boundaryConditionTypes[pdeDim+0]=="Neumann")
                            {
                                p_bcc->AddNeumannBoundaryCondition(*boundary_iter, vectorConstBCs[pdeDim+0], pdeDim);
                            }
                        }
                        else if(y==0)
                        {
                            if(boundaryConditionTypes[pdeDim+1]=="Neumann")
                            {
                                p_bcc->AddNeumannBoundaryCondition(*boundary_iter, vectorConstBCs[pdeDim+1], pdeDim);
                            }
                        }
                        //else if(fabs(x-this -> mMeshDimensions[0])<1e-6)
                        else if(fabs(x-maxElementNodeX)<1e-6)
                        {
                            if(boundaryConditionTypes[pdeDim+2]=="Neumann")
                            {
                                p_bcc->AddNeumannBoundaryCondition(*boundary_iter, vectorConstBCs[pdeDim+2], pdeDim);
                            }
                        }
                        //else if(fabs(y-this -> mMeshDimensions[1])<1e-6)
                        else if(fabs(y-maxElementNodeY)<1e-6)
                        {
                            if(boundaryConditionTypes[pdeDim+3]=="Neumann")
                            {
                                p_bcc->AddNeumannBoundaryCondition(*boundary_iter, vectorConstBCs[pdeDim+3], pdeDim);
                            }
                        }
                        else
                        {
                            std::cout<<"Error: ChemicalDomainFieldForCellCoupling::ProcessBoundaryConditions() Neumann node not on boundary"<<std::endl;
                        }
                
                        boundary_iter++;
                    } 

                    break;

                case 3:

                    

                    // run through Boundary nodes
                    for (typename TetrahedralMesh<ELEMENT_DIM,SPACE_DIM>::BoundaryNodeIterator node_iter = this->mpFeMesh->GetBoundaryNodeIteratorBegin();
                    node_iter != this->mpFeMesh->GetBoundaryNodeIteratorEnd();
                    ++node_iter)
                    {
                        double x = (*node_iter)->GetPoint()[0];
                        double y = (*node_iter)->GetPoint()[1];
                        double z = (*node_iter)->GetPoint()[2];
                        if(x>maxNodeX)
                        {
                            maxNodeX = x;
                        }
                        if(y>maxNodeY)
                        {
                            maxNodeY = y;
                        }
                        if(z>maxNodeZ)
                        {
                            maxNodeZ = z;
                        }
                    }

                    // run through Neumann BCs
                    for (typename TetrahedralMesh<ELEMENT_DIM,SPACE_DIM>::BoundaryElementIterator boundary_iter = this->mpFeMesh->GetBoundaryElementIteratorBegin();
                    boundary_iter != this->mpFeMesh->GetBoundaryElementIteratorEnd();
                    boundary_iter++)
                    {
                        unsigned node_index = (*boundary_iter)->GetNodeGlobalIndex(0);
                        double x = this->mpFeMesh->GetNode(node_index)->GetPoint()[0];
                        double y = this->mpFeMesh->GetNode(node_index)->GetPoint()[1];
                        double z = this->mpFeMesh->GetNode(node_index)->GetPoint()[2];

                        if(x>maxElementNodeX)
                        {
                            maxElementNodeX = x;
                        }
                        if(y>maxElementNodeY)
                        {
                            maxElementNodeY = y;
                        }
                        if(z>maxElementNodeZ)
                        {
                            maxElementNodeZ = z;
                        }

                    }



                    // run through Dirichlet BCs
                    for (typename TetrahedralMesh<ELEMENT_DIM,SPACE_DIM>::BoundaryNodeIterator node_iter = this->mpFeMesh->GetBoundaryNodeIteratorBegin();
                    node_iter != this->mpFeMesh->GetBoundaryNodeIteratorEnd();
                    ++node_iter)
                    {

                        double x = (*node_iter)->GetPoint()[0];
                        double y = (*node_iter)->GetPoint()[1];
                        double z = (*node_iter)->GetPoint()[2];
                        if ((x==0) && (z==0))
                        {
                            if(boundaryConditionTypes[pdeDim+0]=="Dirichlet")
                            {
                                p_bcc->AddDirichletBoundaryCondition(*node_iter, vectorConstBCs[pdeDim+0], pdeDim);
                            }
                        }
                        else if((y==0) && (z==0))
                        {
                            if(boundaryConditionTypes[pdeDim+1]=="Dirichlet")
                            {
                                p_bcc->AddDirichletBoundaryCondition(*node_iter, vectorConstBCs[pdeDim+1], pdeDim);
                            }
                        }
                        //else if((fabs(x-this -> mMeshDimensions[0])<1e-6) && (z==0))
                        else if((fabs(x-maxNodeX)<1e-6) && (z==0))
                        {
                            if(boundaryConditionTypes[pdeDim+2]=="Dirichlet")
                            {
                                p_bcc->AddDirichletBoundaryCondition(*node_iter, vectorConstBCs[pdeDim+2], pdeDim);
                            }
                        }
                        //else if((fabs(y-this -> mMeshDimensions[1])<1e-6) && (z==0))
                        else if((fabs(y-maxNodeY)<1e-6) && (z==0))
                        {
                            if(boundaryConditionTypes[pdeDim+3]=="Dirichlet")
                            {
                                p_bcc->AddDirichletBoundaryCondition(*node_iter, vectorConstBCs[pdeDim+3], pdeDim);
                            }
                        }

                        else if ((x==0) && (y==0))
                        {
                            if(boundaryConditionTypes[pdeDim+4]=="Dirichlet")
                            {
                                p_bcc->AddDirichletBoundaryCondition(*node_iter, vectorConstBCs[pdeDim+4], pdeDim);
                            }
                        }
                        //else if((fabs(x-this -> mMeshDimensions[0])<1e-6) && (y==0))
                        else if((fabs(x-maxNodeX)<1e-6) && (y==0))
                        {
                            if(boundaryConditionTypes[pdeDim+5]=="Dirichlet")
                            {
                                p_bcc->AddDirichletBoundaryCondition(*node_iter, vectorConstBCs[pdeDim+5], pdeDim);
                            }
                        }
                        //else if((fabs(x-this -> mMeshDimensions[0])<1e-6) && (fabs(y-this -> mMeshDimensions[1])<1e-6))
                        else if((fabs(x-maxNodeX)<1e-6) && (fabs(y-maxNodeY)<1e-6))
                        {
                            if(boundaryConditionTypes[pdeDim+6]=="Dirichlet")
                            {
                                p_bcc->AddDirichletBoundaryCondition(*node_iter, vectorConstBCs[pdeDim+6], pdeDim);
                            }
                        }
                        //else if((x==0) && (fabs(y-this -> mMeshDimensions[1])<1e-6))
                        else if((x==0) && (fabs(y-maxNodeY)<1e-6))
                        {
                            if(boundaryConditionTypes[pdeDim+7]=="Dirichlet")
                            {
                                p_bcc->AddDirichletBoundaryCondition(*node_iter, vectorConstBCs[pdeDim+7], pdeDim);
                            }
                        }

                        //else if ((x==0) && (fabs(z-this -> mMeshDimensions[2])<1e-6))
                        else if ((x==0) && (fabs(z-maxNodeZ)<1e-6))
                        {
                            if(boundaryConditionTypes[pdeDim+8]=="Dirichlet")
                            {
                                p_bcc->AddDirichletBoundaryCondition(*node_iter, vectorConstBCs[pdeDim+8], pdeDim);
                            }
                        }
                        //else if((y==0) && (fabs(z-this -> mMeshDimensions[2])<1e-6))
                        else if((y==0) && (fabs(z-maxNodeZ)<1e-6))
                        {
                            if(boundaryConditionTypes[pdeDim+9]=="Dirichlet")
                            {
                                p_bcc->AddDirichletBoundaryCondition(*node_iter, vectorConstBCs[pdeDim+9], pdeDim);
                            }
                        }
                        //else if((fabs(x-this -> mMeshDimensions[0])<1e-6) && (fabs(z-this -> mMeshDimensions[2])<1e-6))
                        else if((fabs(x-maxNodeX)<1e-6) && (fabs(z-maxNodeZ)<1e-6))
                        {
                            if(boundaryConditionTypes[pdeDim+10]=="Dirichlet")
                            {
                                p_bcc->AddDirichletBoundaryCondition(*node_iter, vectorConstBCs[pdeDim+10], pdeDim);
                            }
                        }
                        //else if((fabs(y-this -> mMeshDimensions[1])<1e-6) && (fabs(z-this -> mMeshDimensions[2])<1e-6))
                        else if((fabs(y-maxNodeY)<1e-6) && (fabs(z-maxNodeZ)<1e-6))
                        {
                            if(boundaryConditionTypes[pdeDim+11]=="Dirichlet")
                            {
                                p_bcc->AddDirichletBoundaryCondition(*node_iter, vectorConstBCs[pdeDim+11], pdeDim);
                            }
                        }


                        else
                        {
                            std::cout<<"Error: ChemicalDomainFieldForCellCoupling::ProcessBoundaryConditions() Dirichlet node not on boundary"<<std::endl;
                        }
                
                        node_iter++;
                    } 

                    // run through Neumann BCs
                    for (typename TetrahedralMesh<ELEMENT_DIM,SPACE_DIM>::BoundaryElementIterator boundary_iter = this->mpFeMesh->GetBoundaryElementIteratorBegin();
                    boundary_iter != this->mpFeMesh->GetBoundaryElementIteratorEnd();
                    boundary_iter++)
                    {
                        unsigned node_index = (*boundary_iter)->GetNodeGlobalIndex(0);
                        double x = this->mpFeMesh->GetNode(node_index)->GetPoint()[0];
                        double y = this->mpFeMesh->GetNode(node_index)->GetPoint()[1];
                        double z = this->mpFeMesh->GetNode(node_index)->GetPoint()[2];
                
                        if ((x==0) && (z==0))
                        {
                            if(boundaryConditionTypes[pdeDim+0]=="Neumann")
                            {
                                p_bcc->AddNeumannBoundaryCondition(*boundary_iter, vectorConstBCs[pdeDim+0], pdeDim);
                            }
                        }
                        else if((y==0) && (z==0))
                        {
                            if(boundaryConditionTypes[pdeDim+1]=="Neumann")
                            {
                                p_bcc->AddNeumannBoundaryCondition(*boundary_iter, vectorConstBCs[pdeDim+1], pdeDim);
                            }
                        }
                        //else if((fabs(x-this -> mMeshDimensions[0])<1e-6) && (z==0))
                        else if((fabs(x-maxElementNodeX)<1e-6) && (z==0))
                        {
                            if(boundaryConditionTypes[pdeDim+2]=="Neumann")
                            {
                                p_bcc->AddNeumannBoundaryCondition(*boundary_iter, vectorConstBCs[pdeDim+2], pdeDim);
                            }
                        }
                        //else if((fabs(y-this -> mMeshDimensions[1])<1e-6) && (z==0))
                        else if((fabs(y-maxElementNodeY)<1e-6) && (z==0))
                        {
                            if(boundaryConditionTypes[pdeDim+3]=="Neumann")
                            {
                                p_bcc->AddNeumannBoundaryCondition(*boundary_iter, vectorConstBCs[pdeDim+3], pdeDim);
                            }
                        }

                        else if ((x==0) && (y==0))
                        {
                            if(boundaryConditionTypes[pdeDim+4]=="Neumann")
                            {
                                p_bcc->AddNeumannBoundaryCondition(*boundary_iter, vectorConstBCs[pdeDim+4], pdeDim);
                            }
                        }
                        //else if((fabs(x-this -> mMeshDimensions[0])<1e-6) && (y==0))
                        else if((fabs(x-maxElementNodeX)<1e-6) && (y==0))
                        {
                            if(boundaryConditionTypes[pdeDim+5]=="Neumann")
                            {
                                p_bcc->AddNeumannBoundaryCondition(*boundary_iter, vectorConstBCs[pdeDim+5], pdeDim);
                            }
                        }
                        //else if((fabs(x-this -> mMeshDimensions[0])<1e-6) && (fabs(y-this -> mMeshDimensions[1])<1e-6))
                        else if((fabs(x-maxElementNodeX)<1e-6) && (fabs(y-maxElementNodeY)<1e-6))
                        {
                            if(boundaryConditionTypes[pdeDim+6]=="Neumann")
                            {
                                p_bcc->AddNeumannBoundaryCondition(*boundary_iter, vectorConstBCs[pdeDim+6], pdeDim);
                            }
                        }
                        //else if((x==0) && (fabs(y-this -> mMeshDimensions[1])<1e-6))
                        else if((x==0) && (fabs(y-maxElementNodeY)<1e-6))
                        {
                            if(boundaryConditionTypes[pdeDim+7]=="Neumann")
                            {
                                p_bcc->AddNeumannBoundaryCondition(*boundary_iter, vectorConstBCs[pdeDim+7], pdeDim);
                            }
                        }

                        //else if ((x==0) && (fabs(z-this -> mMeshDimensions[2])<1e-6))
                        else if ((x==0) && (fabs(z-maxElementNodeZ)<1e-6))
                        {
                            if(boundaryConditionTypes[pdeDim+8]=="Neumann")
                            {
                                p_bcc->AddNeumannBoundaryCondition(*boundary_iter, vectorConstBCs[pdeDim+8], pdeDim);
                            }
                        }
                        //else if((y==0) && (fabs(z-this -> mMeshDimensions[2])<1e-6))
                        else if((y==0) && (fabs(z-maxElementNodeZ)<1e-6))
                        {
                            if(boundaryConditionTypes[pdeDim+9]=="Neumann")
                            {
                                p_bcc->AddNeumannBoundaryCondition(*boundary_iter, vectorConstBCs[pdeDim+9], pdeDim);
                            }
                        }
                        //else if((fabs(x-this -> mMeshDimensions[0])<1e-6) && (fabs(z-this -> mMeshDimensions[2])<1e-6))
                        else if((fabs(x-maxElementNodeX)<1e-6) && (fabs(z-maxElementNodeZ)<1e-6))
                        {
                            if(boundaryConditionTypes[pdeDim+10]=="Neumann")
                            {
                                p_bcc->AddNeumannBoundaryCondition(*boundary_iter, vectorConstBCs[pdeDim+10], pdeDim);
                            }
                        }
                        //else if((fabs(y-this -> mMeshDimensions[1])<1e-6) && (fabs(z-this -> mMeshDimensions[2])<1e-6))
                        else if((fabs(y-maxElementNodeY)<1e-6) && (fabs(z-maxElementNodeZ)<1e-6))
                        {
                            if(boundaryConditionTypes[pdeDim+11]=="Neumann")
                            {
                                p_bcc->AddNeumannBoundaryCondition(*boundary_iter, vectorConstBCs[pdeDim+11], pdeDim);
                            }
                        }


                        else
                        {
                            std::cout<<"Error: ChemicalDomainFieldForCellCoupling::ProcessBoundaryConditions() Neumann node not on boundary"<<std::endl;
                        }
                
                        boundary_iter++;
                    } 

                    break;

                case 1:

                    break;

                default:
                    // error
                    std::cout<<"Error: ChemicalDomainFieldForCellCoupling::ProcessBoundaryConditions() default boundary condition dimension reached"<<std::endl;

            }


        }else
        {
            
            // the whole boundary perimeter is of the same BC type
            if(boundaryConditionTypes[pdeDim]=="Dirichlet"||boundaryConditionTypes[pdeDim]=="dirichlet"||boundaryConditionTypes[pdeDim]=="D"||boundaryConditionTypes[pdeDim]=="d")
            {
                // standardise
                boundaryConditionTypes[pdeDim] = "Dirichlet";

                // run though each bondary node and set boundary condition
                for (typename TetrahedralMesh<ELEMENT_DIM,SPACE_DIM>::BoundaryNodeIterator node_iter = this->mpFeMesh->GetBoundaryNodeIteratorBegin();
                node_iter != this->mpFeMesh->GetBoundaryNodeIteratorEnd();
                ++node_iter)
                {
                    p_bcc->AddDirichletBoundaryCondition(*node_iter, vectorConstBCs[pdeDim], pdeDim);
                }
            }else if(boundaryConditionTypes[pdeDim]=="Neumann"||boundaryConditionTypes[pdeDim]=="neumann"||boundaryConditionTypes[pdeDim]=="N"||boundaryConditionTypes[pdeDim]=="n")
            {
                // standardise
                boundaryConditionTypes[pdeDim] = "Neumann";

                // run though each bondary node and set boundary condition
                for (typename TetrahedralMesh<ELEMENT_DIM,SPACE_DIM>::BoundaryElementIterator boundary_iter = this->mpFeMesh->GetBoundaryElementIteratorBegin();
                boundary_iter != this->mpFeMesh->GetBoundaryElementIteratorEnd();
                boundary_iter++)
                {
                    p_bcc->AddNeumannBoundaryCondition(*boundary_iter, vectorConstBCs[pdeDim], pdeDim);
                }
            }

        }
    }
    // update boundary condition types
    this -> SetBoundaryConditionTypes(boundaryConditionTypes);

    return p_bcc;
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

    SetDomainStateVariableRegister(this->GetStateVariableVector());//p_domain_register);
    
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
void ChemicalDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::FeMeshScaling( double targetWidth,
                        double scaleWidth,
                        double targetHeight,
                        double scaleHeight,
                        double targetDepth,
                        double scaleDepth)
{

    this->mTargetWidth = targetWidth;
    this->mTargetHeight = targetHeight;
    this->mTargetDepth = targetDepth;
    this->mScaleWidth = scaleWidth;
    this->mScaleHeight = scaleHeight;
    this->mScaleHeight = scaleDepth;

    double tol =1e-10;
    if(std::abs(targetWidth)>tol)
    {
        mIsScaleBy=false;
    }
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
void ChemicalDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::ScaleFeMesh()
{

    switch(SPACE_DIM)
    {
        case 2:

            std::cout<<"2D: pre-scaled width: "<<this->mpFeMesh->GetWidth(0)<<" pre-scaled height: "<<this->mpFeMesh->GetWidth(1)<<std::endl;
            /*
            if(mIsScaleBy)
            {
                std::cout<<"mIsScaleBy"<<std::endl;
                // apply constant factor to the scale in x and y
                this->mpFeMesh->Scale(mScaleWidth,mScaleHeight);
            }
            else
            {
                std::cout<<"scale to"<<std::endl;
                // scale to a given target width/height
                this->mpFeMesh->Scale(mTargetWidth/this->mpFeMesh->GetWidth(0),mTargetHeight/this->mpFeMesh ->GetWidth(1));
            }
            */

            std::cout<<"2D: Scaled width: "<<this->mpFeMesh->GetWidth(0)<<" Scaled height: "<<this->mpFeMesh->GetWidth(1)<<std::endl;

            break;

        case 3:

            std::cout<<"3D: pre-scaled width: "<<this->mpFeMesh->GetWidth(0)<<" pre-scaled height: "<<this->mpFeMesh->GetWidth(1)<<" pre-scaled depth: "<<this->mpFeMesh->GetWidth(2)<<std::endl;

            if(mIsScaleBy)
            {
                std::cout<<"mIsScaleBy"<<std::endl;
                // apply constant factor to the scale in x and y
                this->mpFeMesh->Scale(mScaleWidth,mScaleHeight,mScaleDepth);
            }
            else
            {
                std::cout<<"scale to"<<std::endl;
                // scale to a given target width/height
                this->mpFeMesh->Scale(mTargetWidth/this->mpFeMesh->GetWidth(0),mTargetHeight/this->mpFeMesh ->GetWidth(1),mTargetDepth/this->mpFeMesh ->GetWidth(2));
            }


            std::cout<<"3D: Scaled width: "<<this->mpFeMesh->GetWidth(0)<<" Scaled height: "<<this->mpFeMesh->GetWidth(1)<<" Scaled depth: "<<this->mpFeMesh->GetWidth(2)<<std::endl;

            break;


        case 1:

            std::cout<<"Error ChemicalDomainField: 1D domain not supported"<<std::endl;
            break;

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
        std::cout<<"Error: ChemicalDomainFieldTemplated::GetDiffusionValueBasedOnPoint: State not in state variable"<<std::endl;
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

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
void ChemicalDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetPiecewiseBCs(std::vector<bool> piecewiseBCs)
{
    mPiecewiseBCs = piecewiseBCs;
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

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
std::vector<bool> ChemicalDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::GetPiecewiseBCs()
{
    return mPiecewiseBCs;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
boost::shared_ptr<BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM> > ChemicalDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::ReturnSharedPtrBoundaryConditionsContainer()
{
    return mpBoundaryConditionsContainer;
}
#endif