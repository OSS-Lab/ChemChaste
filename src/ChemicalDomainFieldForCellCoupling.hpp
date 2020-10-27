#ifndef CHEMICALDOMAINFIELDFORCELLCOUPLING_HPP
#define CHEMICALDOMAINFIELDFORCELLCOUPLING_HPP

#include <string>
#include <vector>
#include "PetscSetupAndFinalize.hpp"
#include "InhomogenousParabolicPdeForCoupledOdeSystem_templated.hpp"

#include "ChemicalDomainField_templated.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "TetrahedralMesh.hpp"
#include "EulerIvpOdeSolver.hpp"

#include "BoundaryConditionsContainer_extended.hpp"
#include "ConstBoundaryCondition.hpp"


// class to handle the formation of a chemcially active domain ready for coupling to a cell 
// simulation

// need to template over <ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM> for PDE and BCC

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
class ChemicalDomainFieldForCellCoupling: public ChemicalDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>
{
private:

    using ChemicalDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::GetFieldType;
    using ChemicalDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetUpDomainFromFiles;

protected:

    std::string mCellLabelFilename;

    std::string mCellKeyFilename;

    std::string mInitialConditionsFilename;

    std::string mBoundaryConditionsFilename;

    boost::shared_ptr<BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM> > mpBoundaryConditionsContainer; 

    boost::shared_ptr<InhomogenousParabolicPdeForCoupledOdeSystemTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>> mpPdeSystem; 

    std::vector<boost::shared_ptr<AbstractIvpOdeSolver>> mOdeSolverSystem;

    std::vector<AbstractInhomogenousOdeSystemForCoupledPdeSystem*> mNodalOdeSystem;

    std::vector<double> mMeshDomainLower; // for locating the origin of the cell layer for cell based simulations

    std::vector<double> mMeshDomainUpper;

    std::vector<double> mMeshCentre;


public:


    ChemicalDomainFieldForCellCoupling( std::string reactionFileRoot="",
                                        std::string cellLabelFilename="",
                                        std::string cellKeyFilename="",
                                        std::string domainLabelFilename="", 
                                        std::string domainKeyFilename="", 
                                        std::string odeLabelFilename="", 
                                        std::string odeKeyFilename="",
                                        std::string diffusionFilename="",
                                        std::string initialConditionsFilename="",
                                        std::string boundaryConditionsFilename="",
                                        bool isHoneyCombMesh=false,
                                        std::vector<double> labelOrigin = std::vector<double>(),
                                        std::vector<double> cartesianCellScaleXY = std::vector<double>(),
                                        std::vector<double> cartesianOdeScaleXY = std::vector<double>(),
                                        std::vector<double> meshDomainLower = std::vector<double>(),
                                        std::vector<double> meshDomainUpper = std::vector<double>(),
                                        std::vector<double> meshCentre = std::vector<double>()
                                        );

    virtual ~ChemicalDomainFieldForCellCoupling()
    {
    }

    virtual void SetUpDomainFromFiles();

    virtual std::string GetFieldType()
    {
        return "ChemicalDomainFieldForCellCoupling";
    };

    virtual boost::shared_ptr<BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>> ProcessBoundaryConditions();


    virtual void SetUpPdeSystem();

    virtual void SetUpNodalOdeSystems();

    void CalculateFeMeshCentreOffset();

    boost::shared_ptr<InhomogenousParabolicPdeForCoupledOdeSystemTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>> ReturnSharedPtrPdeSystem();
  
    boost::shared_ptr<BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM> > ReturnSharedPtrBoundaryConditionsContainer();

    // set methods

    void SetBoundaryConditionsContainer(BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM> bcc);

    void SetPdeSystem(boost::shared_ptr<InhomogenousParabolicPdeForCoupledOdeSystemTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>> pdeSystem);

    void SetNodalOdeSystems(std::vector<AbstractInhomogenousOdeSystemForCoupledPdeSystem*> nodalOdeSystem);

    void SetNodalOdeSolvers(std::vector<boost::shared_ptr<AbstractIvpOdeSolver>> odeSolverSystem);

    void SetMeshDomainLower(std::vector<double> lower);

    void SetMeshDomainUpper(std::vector<double> upper);

    void SetMeshCentre(std::vector<double> centre);

    // get methods

    std::vector<AbstractInhomogenousOdeSystemForCoupledPdeSystem*> GetNodalOdeSystems();

    std::vector<boost::shared_ptr<AbstractIvpOdeSolver>> GetNodalOdeSolvers();

    std::vector<double> GetMeshDomainLower();

    std::vector<double> GetMeshDomainUpper();

    std::vector<double> GetMeshDomainCentre();
    

};


template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
ChemicalDomainFieldForCellCoupling<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::ChemicalDomainFieldForCellCoupling( 
                                        std::string reactionFileRoot,
                                        std::string cellLabelFilename,
                                        std::string cellKeyFilename,
                                        std::string domainLabelFilename, 
                                        std::string domainKeyFilename, 
                                        std::string odeLabelFilename, 
                                        std::string odeKeyFilename,
                                        std::string diffusionFilename,
                                        std::string initialConditionsFilename,
                                        std::string boundaryConditionsFilename,
                                        bool isHoneyCombMesh,
                                        std::vector<double> labelOrigin,
                                        std::vector<double> cartesianCellScaleXY,
                                        std::vector<double> cartesianOdeScaleXY,
                                        std::vector<double> meshDomainLower,
                                        std::vector<double> meshDomainUpper,
                                        std::vector<double> meshCentre
                                        )

    :   ChemicalDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>(
                            reactionFileRoot,
                            domainLabelFilename, 
                            domainKeyFilename, 
                            odeLabelFilename, 
                            odeKeyFilename, 
                            diffusionFilename,
                            isHoneyCombMesh,
                            labelOrigin,
                            cartesianCellScaleXY,
                            cartesianOdeScaleXY),
        mCellLabelFilename(cellLabelFilename),
        mCellKeyFilename(cellKeyFilename),
        mInitialConditionsFilename(initialConditionsFilename),
        mBoundaryConditionsFilename(boundaryConditionsFilename),
        mpBoundaryConditionsContainer(boost::shared_ptr<BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>>()),
        mpPdeSystem(boost::shared_ptr<InhomogenousParabolicPdeForCoupledOdeSystemTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>> ()),
        mMeshDomainLower(meshDomainLower),
        mMeshDomainUpper(meshDomainUpper),
        mMeshCentre(meshCentre)
    
    {
        if(mMeshDomainLower.empty())
        {
            std::vector<double> lower(SPACE_DIM,0.0);
            mMeshDomainLower = lower;
        }
        if(mMeshDomainUpper.empty())
        {
            std::vector<double> upper(SPACE_DIM,0.0);
            mMeshDomainUpper = upper;
        }
        if(mMeshCentre.empty())
        {
            std::vector<double> centre(SPACE_DIM,0.0);
            mMeshCentre = centre;
        }
        SetUpDomainFromFiles();

    }

    template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
    void ChemicalDomainFieldForCellCoupling<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetUpDomainFromFiles()
    {
        // set up a domain with a chemistry, form mesh and reaction nodes etc

        // process initial conditions

        this -> ParseInitialConditionsFromFile(mInitialConditionsFilename);
        // process boundary conditions
        this -> ParseBoundaryConditionsFromFile(mBoundaryConditionsFilename);

        mpBoundaryConditionsContainer = ProcessBoundaryConditions();

        // convert pde and nodal ode systems for interface with cell based population simulation 
        SetUpPdeSystem();

        SetUpNodalOdeSystems();

        CalculateFeMeshCentreOffset();
    }

    template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
    boost::shared_ptr<BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>> ChemicalDomainFieldForCellCoupling<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::ProcessBoundaryConditions()
    {

        // retireive the boundary conditions information derived from the files inputs
        std::vector<std::string> boundaryConditionTypes = this -> GetBoundaryConditionTypes();
        std::vector<double> boundaryConditionValues = this -> GetBoundaryConditionValues();

        // create containeers to store the boundary conditions, assume each boundary condition is constant over the simulation
        boost::shared_ptr<BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>> p_bcc(new BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>());
        std::vector<ConstBoundaryCondition<SPACE_DIM>*> vectorConstBCs;
 
        for(unsigned pdeDim=0; pdeDim<PROBLEM_DIM; pdeDim++){
            vectorConstBCs.push_back(new ConstBoundaryCondition<SPACE_DIM>(boundaryConditionValues[pdeDim]));
        }

        for(unsigned pdeDim=0; pdeDim<PROBLEM_DIM; pdeDim++)
        {
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
        // update boundary condition types
        this -> SetBoundaryConditionTypes(boundaryConditionTypes);

        return p_bcc;
    }

    
    
    template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
    void ChemicalDomainFieldForCellCoupling<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetUpPdeSystem()
    { 
        boost::shared_ptr<InhomogenousParabolicPdeForCoupledOdeSystemTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>> pde( new InhomogenousParabolicPdeForCoupledOdeSystemTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>(this));
        SetPdeSystem(pde);
    }


    template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
    void ChemicalDomainFieldForCellCoupling<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetUpNodalOdeSystems()
    {
        // used the explicitly defined EulerIvpOdeSolver rather than the default defined BackwardEuler method
        std::vector<boost::shared_ptr<AbstractIvpOdeSolver>> odeSolverSystem;
        std::vector<AbstractInhomogenousOdeSystemForCoupledPdeSystem*> odeSystem;


        for(unsigned i=0; i<this->mpFeMesh->GetNumNodes(); i++)
        {
            odeSystem.push_back(this->GetOdeSystem()[i]);
            boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);
            odeSolverSystem.push_back(p_solver);
        }


        SetNodalOdeSystems(odeSystem);
        SetNodalOdeSolvers(odeSolverSystem);
    }


    template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
    void ChemicalDomainFieldForCellCoupling<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::CalculateFeMeshCentreOffset()
    {
        std::vector<double> mesh_lengths(SPACE_DIM,0.0);

        for(unsigned dim=0; dim<SPACE_DIM; dim++)
        {
            mesh_lengths[dim] = this-> mMeshDimensions[dim] * this-> mMeshScale[dim];
            mMeshCentre[dim] = 0.5*mesh_lengths[dim]; 
        }
    }

    template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
    boost::shared_ptr<InhomogenousParabolicPdeForCoupledOdeSystemTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>> ChemicalDomainFieldForCellCoupling<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::ReturnSharedPtrPdeSystem()
    {
        return mpPdeSystem;
    }
  
    template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
    boost::shared_ptr<BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM> > ChemicalDomainFieldForCellCoupling<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::ReturnSharedPtrBoundaryConditionsContainer()
    {
        return mpBoundaryConditionsContainer;
    }

    // set methods

    template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
    void ChemicalDomainFieldForCellCoupling<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetPdeSystem(boost::shared_ptr<InhomogenousParabolicPdeForCoupledOdeSystemTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>> p_pdeSystem)
    {
        mpPdeSystem = p_pdeSystem;
    }

    template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
    void ChemicalDomainFieldForCellCoupling<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetNodalOdeSystems(std::vector<AbstractInhomogenousOdeSystemForCoupledPdeSystem*> nodalOdeSystem)
    {
        mNodalOdeSystem = nodalOdeSystem;
    }

    template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
    void ChemicalDomainFieldForCellCoupling<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetNodalOdeSolvers(std::vector<boost::shared_ptr<AbstractIvpOdeSolver>> odeSolverSystem)
    {
        mOdeSolverSystem = odeSolverSystem;
    }

    template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
    void ChemicalDomainFieldForCellCoupling<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetMeshDomainLower(std::vector<double> lower)
    {
        mMeshDomainLower = lower;
    }

    template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
    void ChemicalDomainFieldForCellCoupling<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetMeshDomainUpper(std::vector<double> upper)
    {
        mMeshDomainUpper = upper;
    }

    template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
    void ChemicalDomainFieldForCellCoupling<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetMeshCentre(std::vector<double> meshCentre)
    {
        mMeshCentre = meshCentre;
    }

    // get methods

    template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
    std::vector<AbstractInhomogenousOdeSystemForCoupledPdeSystem*> ChemicalDomainFieldForCellCoupling<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::GetNodalOdeSystems()
    {
        return mNodalOdeSystem;
    }

    template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
    std::vector<boost::shared_ptr<AbstractIvpOdeSolver>> ChemicalDomainFieldForCellCoupling<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::GetNodalOdeSolvers()
    {
        return mOdeSolverSystem;
    }

    template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
    std::vector<double> ChemicalDomainFieldForCellCoupling<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::GetMeshDomainLower()
    {
        return mMeshDomainLower;
    }

    template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
    std::vector<double> ChemicalDomainFieldForCellCoupling<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::GetMeshDomainUpper()
    {
        return mMeshDomainUpper;
    }

    template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
    std::vector<double> ChemicalDomainFieldForCellCoupling<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::GetMeshDomainCentre()
    {
        return mMeshCentre;
    }

#endif