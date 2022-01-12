#ifndef CHEMICALDOMAINFIELDFORCELLCOUPLING_HPP
#define CHEMICALDOMAINFIELDFORCELLCOUPLING_HPP

#include <string>
#include <vector>
//#include "PetscSetupAndFinalize.hpp"
#include "InhomogenousParabolicPdeForCoupledOdeSystem_templated.hpp"
#include "BoundaryConditionsContainer_extended.hpp"
#include "ConstBoundaryCondition.hpp"

#include "ChemicalDomainField_templated.hpp"
#include "TetrahedralMesh.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "AbstractTetrahedralMesh.hpp"
#include "EulerIvpOdeSolver.hpp"




// class to handle the formation of a chemcially active domain ready for coupling to a cell 
// simulation

// need to template over <ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM> for PDE and BCC

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
class ChemicalDomainFieldForCellCoupling: public ChemicalDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>
{
private:

    using ChemicalDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::GetFieldType;
    using ChemicalDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetUpDomainFromFiles;
    //using ChemicalDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::ParseBoundaryConditionsFromFile;

protected:

    bool mIsCellLayerFileBased=false;

    bool mIsCellMeshGenerated=false;

    TetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* mpCellMesh;

    std::string mCellLabelFilename;

    unsigned mNumberOfCellSlices=0;

    std::string mCellKeyFilename;

    boost::shared_ptr<InhomogenousParabolicPdeForCoupledOdeSystemTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>> mpPdeSystem; 

    std::vector<boost::shared_ptr<AbstractIvpOdeSolver>> mOdeSolverSystem;

    std::vector<AbstractInhomogenousOdeSystemForCoupledPdeSystem*> mNodalOdeSystem;

    std::vector<double> mMeshDomainLower; // for locating the origin of the cell layer for cell based simulations

    std::vector<double> mMeshDomainUpper;

    std::vector<double> mMeshCentre;

    // cell information input
    unsigned mNumberOfCellTypes;

    std::vector<std::vector<std::vector<std::string>>> mCellLabels;

    std::vector<std::string> mCellLabelVector;

    std::vector<std::vector<std::string>> mCellKeys;

    std::vector<std::string> mCellKeyVector;

    std::vector<unsigned> mCartesianCellLayerDimensions; // (x, y) of the input csv file grid

    std::vector<double> mCartesianCellLayerScaleXY; // (Sx, Sy) scale factor of the input csv grid, default (1,1)

    // mapped cartesian domain 
    std::vector<unsigned> mCartesianChasteCellDimensions;

    std::vector<double> mCartesianChasteCellScaleXY; // (Sx, Sy) scale factor of the Chaste grid

    std::vector<std::vector<std::string>> mCartesianChasteCellLayer;

    // mapped mesh domain
    std::vector<unsigned> mCellMeshDimensions;

    std::vector<double> mCellMeshScale;

    double mCellStepSize=1.0;

    HoneycombMeshGenerator* mpCellMeshGenerator;

    std::vector<std::string> mCellNodes;

    bool mIsCellHoneyCombMesh;

    std::vector<double> mCellLabelOrigin;

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
                                        std::vector<double> meshCentre = std::vector<double>(),
                                        bool isCellHoneyCombMesh=false,
                                        std::vector<double> cellLabelOrigin = std::vector<double>(),
                                        std::vector<double> cartesianCellLayerScaleXY = std::vector<double>()
                                        );

    virtual ~ChemicalDomainFieldForCellCoupling()
    {
    }

    virtual void SetUpCellLayer();

    virtual void SetUpDomainFromFiles();

    virtual std::string GetFieldType()
    {
        return "ChemicalDomainFieldForCellCoupling";
    };

    virtual void SetUpPdeSystem();

    virtual void SetUpNodalOdeSystems();

    void CalculateFeMeshCentreOffset();

    void SetupAndInitialiseLabelCellLayer();

    void SetupCellLayerMappingDimensions();

    void FormCellMesh();

    void LabelCellMeshNodally();

    std::vector<std::string> LabelNodesWithCells();

    void SetNodeCells(std::vector<std::string>);

    void ReadCellLabels();

    void ReadCellKey();

    void MapToChasteCartesianCellLayer();

    std::string ReturnCellKeyFromCellLabel(std::string);

    std::string ReturnMappedCellLabel(unsigned, unsigned, unsigned, std::vector<double>);

    std::string ReturnCellLabelAtPosition(const c_vector<double,SPACE_DIM>&);

    unsigned ReturnUnsignedIDFromCellKeyString(std::string);

    void SetCellLabels(std::vector<std::vector<std::vector<std::string>>>);

    void SetCellLabelVector(std::vector<std::string>);

    void SetCellKeys(std::vector<std::vector<std::string>>);

    void SetCellKeyVector(std::vector<std::string>);

    void SetCartesianCellLayerDimensions(std::vector<unsigned>);

    void SetCartesianChasteCellScale(std::vector<double>);

    void SetCartesianChasteCellLayer(std::vector<std::vector<std::string>>);

    std::vector<std::vector<std::vector<std::string>>> GetCellLabels();

    std::string GetCellLabelByIndex(unsigned);

    std::vector<std::string> GetCellLabelVector();

    std::vector<std::vector<std::string>> GetCellKeys();

    std::vector<std::string> GetCellKeyVector();

    void SetCartesianChasteCellDimensions(std::vector<unsigned>);

    void SetCellMesh(TetrahedralMesh<ELEMENT_DIM,SPACE_DIM>*&);

    void SetCellMeshGenerator(HoneycombMeshGenerator*);

    void SetCellMeshScale(std::vector<double>);

    void SetCellMeshDimensions(std::vector<unsigned>);

    unsigned GetNumberOfCellTypes();

    std::vector<std::string> GetNodeCells();


    boost::shared_ptr<InhomogenousParabolicPdeForCoupledOdeSystemTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>> ReturnSharedPtrPdeSystem();


    void PrintCellIDWithKey();

    void PrintCellIDWithCellLabel();

    // set methods

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

    TetrahedralMesh<ELEMENT_DIM,SPACE_DIM>*& rGetCellMesh();

    HoneycombMeshGenerator* GetCellMeshGenerator();
    

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
                                        std::vector<double> meshCentre,
                                        bool isCellHoneyCombMesh,
                                        std::vector<double> cellLabelOrigin,
                                        std::vector<double> cartesianCellLayerScaleXY
                                        )

    :   ChemicalDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>(
                            reactionFileRoot,
                            domainLabelFilename, 
                            domainKeyFilename, 
                            odeLabelFilename, 
                            odeKeyFilename, 
                            diffusionFilename,
                            initialConditionsFilename,
                            boundaryConditionsFilename,
                            isHoneyCombMesh,
                            labelOrigin,
                            cartesianCellScaleXY,
                            cartesianOdeScaleXY),
        mCellLabelFilename(cellLabelFilename),
        mCellKeyFilename(cellKeyFilename),
        mpPdeSystem(boost::shared_ptr<InhomogenousParabolicPdeForCoupledOdeSystemTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>> ()),
        mMeshDomainLower(meshDomainLower),
        mMeshDomainUpper(meshDomainUpper),
        mMeshCentre(meshCentre),
        mIsCellHoneyCombMesh(isCellHoneyCombMesh),
        mCellLabelOrigin(cellLabelOrigin),
        mCartesianCellLayerScaleXY(cartesianCellLayerScaleXY)
    
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

        if(mCartesianCellLayerScaleXY.empty())
        {
            // default to scale 1.0 in all directions
            std::vector<double> cellLayerScaleXY(SPACE_DIM,1.0);
            mCartesianCellLayerScaleXY = cellLayerScaleXY;
        }

        if(mCellLabelOrigin.empty())
        {
            std::vector<double> origin(SPACE_DIM,0);
            mCellLabelOrigin = origin;
        }

        if(cellLabelFilename != cellKeyFilename)
        {
            mIsCellLayerFileBased = true;
            SetUpCellLayer();
        }

        SetUpDomainFromFiles();

        PrintCellIDWithKey();

        PrintCellIDWithCellLabel();
    }

    template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
    void ChemicalDomainFieldForCellCoupling<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetUpCellLayer()
    {
        SetupAndInitialiseLabelCellLayer();

        // scale the input file dimensions to the chaste rectilinear cartesian grid
        MapToChasteCartesianCellLayer();

        // translate the ode labels to ode nodes, form the ode vector
        SetupCellLayerMappingDimensions();

        FormCellMesh();
        // label the newly formed mesh
        LabelCellMeshNodally();

        PrintCellIDWithCellLabel();

        mIsCellMeshGenerated = true;
    }


    template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
    void ChemicalDomainFieldForCellCoupling<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetUpDomainFromFiles()
    {
        // set up a domain with a chemistry, form mesh and reaction nodes etc


        // convert pde and nodal ode systems for interface with cell based population simulation 
        SetUpPdeSystem();

        SetUpNodalOdeSystems();

        CalculateFeMeshCentreOffset();
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
    void ChemicalDomainFieldForCellCoupling<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetupAndInitialiseLabelCellLayer()
    {
        // read in and store the cell grid labels for mapping to mesh nodes and chaste rectilinear grid
        ReadCellLabels();
        ReadCellKey();
        
        std::vector<unsigned> cartesianCellLayerDimensions;

        cartesianCellLayerDimensions.push_back(mCellLabels[0][0].size());
        cartesianCellLayerDimensions.push_back(mCellLabels[0].size());
        cartesianCellLayerDimensions.push_back(mCellLabels.size());

        // store these values
        SetCartesianCellLayerDimensions(cartesianCellLayerDimensions);
    }

    template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
    void ChemicalDomainFieldForCellCoupling<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetupCellLayerMappingDimensions()
    {
        // hard code the honey comb mesh scale factor and dimensions with a default of tetrahedral mesh
        std::vector<double> meshScale;
        if(mIsCellHoneyCombMesh)
        {   
            switch (SPACE_DIM)
            {
                case 2:
                    
                    meshScale.push_back(1.0);
                    meshScale.push_back(0.866025);
                    meshScale.push_back(1.0);
                    break;
                default:
                    // not sure as to honeycomb dimensions
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
        SetCellMeshScale(meshScale);
        

        std::vector<unsigned> meshDimensions(SPACE_DIM,0);
        std::vector<double> meshSegment(SPACE_DIM,0.0);
        double meshElementX = 0.0;
        double meshElementY = 0.0;
        double meshElementZ = 0.0;
        
        // determine the mesh dimensions up to the maximum allowed through the scaling
        // cover the whole cell set with the mesh

        for(unsigned dim=0; dim<SPACE_DIM; dim++)
        {
            unsigned count=0;
            while( meshSegment[dim] < mCartesianCellLayerDimensions[dim]*mCartesianCellLayerScaleXY[dim])
            {
                meshSegment[dim]  = count*meshScale[dim];
                count=count+1;
            }
            meshDimensions[dim]=count-2;
        }

        SetCellMeshDimensions(meshDimensions);
    }

    template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
    void ChemicalDomainFieldForCellCoupling<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::FormCellMesh()
    {
        // use the aforementioned mesh dimension to produce a new honeycomb mesh, as mutable mesh type
        std::cout<<"ChemicalDomainFieldForCellCoupling<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::FormCellMesh()"<<std::endl;
        switch (SPACE_DIM)
        {
            case 2:
                if(mIsCellHoneyCombMesh)
                {
                    HoneycombMeshGenerator* p_generator = new HoneycombMeshGenerator(mCellMeshDimensions[0], mCellMeshDimensions[1], 0,mCellStepSize);

                    mpCellMesh = dynamic_cast<TetrahedralMesh<ELEMENT_DIM,SPACE_DIM>*>(p_generator -> GetMesh());

                    SetCellMeshGenerator(p_generator);
        
                    SetCellMesh(mpCellMesh);
                }
                else
                {
                    mpCellMesh = new TetrahedralMesh<ELEMENT_DIM,SPACE_DIM>();

                    mpCellMesh->ConstructRegularSlabMesh(mCellStepSize, mCellMeshDimensions[0], mCellMeshDimensions[1], mCellMeshDimensions[2]);

                    SetCellMesh(mpCellMesh);
                }

                break;

            case 3:

                mpCellMesh = new TetrahedralMesh<ELEMENT_DIM,SPACE_DIM>();
                mpCellMesh->ConstructRegularSlabMesh(mCellStepSize, mCellMeshDimensions[0], mCellMeshDimensions[1]);
                SetCellMesh(mpCellMesh);
                break;

            case 1:

                mpCellMesh = new TetrahedralMesh<ELEMENT_DIM,SPACE_DIM>();
                mpCellMesh->ConstructRegularSlabMesh(mCellStepSize, mCellMeshDimensions[0]);
                SetCellMesh(mpCellMesh);
                break;

        }
        
    }

    template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
    std::vector<std::string> ChemicalDomainFieldForCellCoupling<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::LabelNodesWithCells()
    {
        // retireve the node positions and compare to domain labels in order to label the nodes
        // nodes are defined bottom-left of domain to top-right in a serialised manner

        unsigned numberOfNodes = mpCellMesh -> GetNumNodes();

        // form the serialised node label array
        std::vector<std::string> serialisedNodeCells(numberOfNodes,"");
        
        for (typename AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>::NodeIterator iter = mpCellMesh ->GetNodeIteratorBegin();
                iter != mpCellMesh ->GetNodeIteratorEnd();
                ++iter)
        {
            // for each node, retireve index and position, label with correct ode label
            unsigned node_index = iter ->GetIndex();
            // retireve the position of the node of index node_index
            const c_vector<double,SPACE_DIM>& position = iter->rGetLocation();
            // search through the position ode labels for the correct label to give the node

            serialisedNodeCells.at(node_index) = ReturnCellLabelAtPosition(position);
        }

        return serialisedNodeCells;
    }


    template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
    std::string ChemicalDomainFieldForCellCoupling<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::ReturnCellLabelAtPosition(const c_vector<double,SPACE_DIM>& position)
    {
        // search through the cell labels to provide the specific label at the given position
        std::vector<unsigned> labelIndex(SPACE_DIM,0);


        double position_in_domain=0.0;

        for(unsigned dim=0; dim<SPACE_DIM; dim++)
        {
            position_in_domain = position[dim] + mCellLabelOrigin[dim];

            for(unsigned i=0; i<(mCartesianChasteCellDimensions[dim]); i++)
            {
                if(position_in_domain <= i*mCartesianChasteCellScaleXY[dim])
                {
                    labelIndex[dim] = i;
                    break;
                }
            }
        }
        
        // return the domain label for the grid element at the Y,X index containing the position
        return mCartesianChasteCellLayer[labelIndex[1]][labelIndex[0]];
    }


    template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
    unsigned ChemicalDomainFieldForCellCoupling<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::ReturnUnsignedIDFromCellKeyString(std::string keyString)
    {
        std::vector<std::string> keyVector = GetCellKeyVector();

        for(unsigned i=0; i<keyVector.size();i++)
        {
            if(keyVector[i]==keyString)
            {
                return i;
            }
        }
        std::cout<<"Error: unsigned ChemicalDomainFieldForCellCoupling<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::ReturnUnsignedIDFromCellKeyString(std::string) index out of bounds"<<std::endl;
        return keyVector.size()+1;
    }


    template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
    void ChemicalDomainFieldForCellCoupling<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::ReadCellLabels()
    {
        // read and store the information contained within the cell labels topology file(s)

        std::vector<std::vector<std::vector<std::string>>> tempCellLabels;

        // mCellLabelFilename will be a .csv if 2D or a directory containg .csv files if 3D
        std::string filename = mCellLabelFilename;
        switch(SPACE_DIM)
        {
            case 2:
                tempCellLabels.push_back(AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::ReadMatrix(filename));
                break;

            case 3:
                for(unsigned i=0; i<mNumberOfCellSlices;i++)
                {
                    filename = mCellLabelFilename + "/" + std::to_string(i) + ".csv";
                    tempCellLabels.push_back(AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::ReadMatrix(filename));
                }
                break;

            case 1:
                std::cout<<"Error: ChemicalDomainFieldForCellCoupling 1D domain not supported"<<std::endl;
                break;



        }


        SetCellLabels(tempCellLabels);

        SetCellLabelVector(AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::ReturnUnique(tempCellLabels));
    }

    template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
    void ChemicalDomainFieldForCellCoupling<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::ReadCellKey()
    {
        // read the label keys for the cell layer, providing more information as to the type of the cell
        std::vector<std::vector<std::string>> cellKeys;

        cellKeys = AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::ReadMatrix(mCellKeyFilename);

        // set the cell key vector for identifying the cell/strain type
        std::vector<std::string> cellKeysUnique = AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::ReturnUnique(cellKeys);
        std::vector<std::string> cellKeyNames;

        for(unsigned i=0; i<cellKeysUnique.size(); i++)
        {
            // take every other value to grab the names of the cell types
            if(i%2 ==1)
            {
                cellKeyNames.push_back(cellKeysUnique[i]);
            }

        }

        SetCellKeyVector(cellKeyNames);

        SetCellKeys(cellKeys);
    }

    template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
    void ChemicalDomainFieldForCellCoupling<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::MapToChasteCartesianCellLayer()
    {
        
        std::vector<double> cartesianChasteCellScaleXY;
        if(mIsCellHoneyCombMesh)
        {
            // is a a honeycomb mesh
            //cartesianChasteCellScaleXY.push_back(0.25);
            //cartesianChasteCellScaleXY.push_back(0.144338);
            cartesianChasteCellScaleXY.push_back(1);
            cartesianChasteCellScaleXY.push_back(1);
        }
        else
        {
            // else tetrahedral mesh
            //cartesianChasteCellScaleXY.push_back(0.166666);
            //cartesianChasteCellScaleXY.push_back(0.166666);
            cartesianChasteCellScaleXY.push_back(1);
            cartesianChasteCellScaleXY.push_back(1);
        }
        SetCartesianChasteCellScale(cartesianChasteCellScaleXY);

        unsigned numberOfXSteps = std::round((mCartesianCellLayerDimensions[0]*mCartesianCellLayerScaleXY[0])/cartesianChasteCellScaleXY[0]); 
        unsigned numberOfYSteps = std::round((mCartesianCellLayerDimensions[1]*mCartesianCellLayerScaleXY[1])/cartesianChasteCellScaleXY[1]); 
    
        std::vector<std::vector<std::string>> cartesianChasteCellLayer;
        std::vector<unsigned> cartesianChasteCellDimensions;
        cartesianChasteCellDimensions.push_back(numberOfXSteps);
        cartesianChasteCellDimensions.push_back(numberOfYSteps);

        // define the boundary with the edge steps included for easier comparisons of (double) spatial values
        for(unsigned yindex=0; yindex<numberOfYSteps; yindex++)
        {
            std::vector<std::string> rowVector;
            for(unsigned xindex=0; xindex<numberOfXSteps; xindex++)
            {
                // push back the domain label of the file domain onto the mapped domain, here everything is in terms of std::string label
                rowVector.push_back(ReturnMappedCellLabel(xindex,yindex,0,cartesianChasteCellScaleXY));

            }

            cartesianChasteCellLayer.push_back(rowVector);
        }

        AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::printMatrix(mCellLabels);
        AbstractDomainFieldTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::printMatrix(cartesianChasteCellLayer);

        SetCartesianChasteCellLayer(cartesianChasteCellLayer);
        SetCartesianChasteCellDimensions(cartesianChasteCellDimensions);

    }

    template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
    std::string ChemicalDomainFieldForCellCoupling<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::ReturnCellKeyFromCellLabel(std::string cellLabel)
    {
        bool IsFound=false;

        for(unsigned key_index=0; key_index<mCellKeys.size(); key_index++)
        {
            if(mCellKeys[key_index][0] == cellLabel)
            {
                return mCellKeys[key_index][1];
                IsFound=true;
                break;
            }
        }
        if(!IsFound)
        {
            std::cout<<"Error: ChemicalDomainFieldForCellCoupling::ReturnCellKeyFromCellLabel, cell label not found"<<std::endl;
            return "Null";
        }
        return "Null";
    }



    template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
    std::string ChemicalDomainFieldForCellCoupling<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::ReturnMappedCellLabel(unsigned xindex, unsigned yindex, unsigned zindex, std::vector<double> scaleFactor)
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

                for(unsigned i=1; i<(mCartesianCellLayerDimensions[1]+1); i++) // top to bottom
                {
                    if(positionY <= i*mCartesianCellLayerScaleXY[1])
                    {
                        labelYIndex = i-1;
                        break;
                    }
                }

                // for the x position

                for(unsigned i=1; i<(mCartesianCellLayerDimensions[0]+1); i++)
                {
                    if(positionX <= i*mCartesianCellLayerScaleXY[0])
                    {
                        labelXIndex = i-1;
                        break;
                    }
                }

                break; 

            case 3:
                // for the z position

                for(unsigned i=1; i<(mCartesianCellLayerDimensions[2]+1); i++) // top to bottom
                {
                    if(positionZ <= i*mCartesianCellLayerScaleXY[2])
                    {
                        labelZIndex = i-1;
                        break;
                    }
                }
                
                // for the y position

                for(unsigned i=1; i<(mCartesianCellLayerDimensions[1]+1); i++) // top to bottom
                {
                    if(positionY <= i*mCartesianCellLayerScaleXY[1])
                    {
                        labelYIndex = i-1;
                        break;
                    }
                }

                // for the x position

                for(unsigned i=1; i<(mCartesianCellLayerDimensions[0]+1); i++)
                {
                    if(positionX <= i*mCartesianCellLayerScaleXY[0])
                    {
                        labelXIndex = i-1;
                        break;
                    }
                }

                break;

            case 1:
                std::cout<<"Error: ChemicalDomainFieldForCellCoupling A 1D cell layer not supported"<<std::endl;
                break;
        }


        return mCellLabels[labelZIndex][labelYIndex][labelXIndex];
    }


    template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
    void ChemicalDomainFieldForCellCoupling<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::LabelCellMeshNodally()
    {
        SetNodeCells(LabelNodesWithCells());
    }


    template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
    void ChemicalDomainFieldForCellCoupling<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::PrintCellIDWithKey()
    {
        std::vector<std::string> keyVector = GetCellKeyVector();
        std::cout<<"CellID,CellKey"<<std::endl;
        for(unsigned i=0; i<keyVector.size(); i++)
        {
            std::cout<<ReturnUnsignedIDFromCellKeyString(keyVector[i])<<","<<keyVector[i]<<std::endl;
        }
    }

    template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
    void ChemicalDomainFieldForCellCoupling<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::PrintCellIDWithCellLabel()
    {
        std::vector<std::string> labelVector = GetCellLabelVector();
        std::cout<<"CellID,CellKey"<<std::endl;
        for(unsigned i=0; i<labelVector.size(); i++)
        {
            std::cout<<ReturnUnsignedIDFromCellKeyString(ReturnCellKeyFromCellLabel(labelVector[i]))<<","<<labelVector[i]<<std::endl;
        }
    }

    template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
    void ChemicalDomainFieldForCellCoupling<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetNodeCells(std::vector<std::string> nodeCells)
    {
        mCellNodes = nodeCells;
    }


    template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
    void ChemicalDomainFieldForCellCoupling<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetCellLabels(std::vector<std::vector<std::vector<std::string>>> cellLabels)
    {
        mCellLabels = cellLabels;
    }

    template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
    void ChemicalDomainFieldForCellCoupling<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetCellLabelVector(std::vector<std::string> cellLabelVector)
    {
        mCellLabelVector = cellLabelVector;
    }

    template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
    void ChemicalDomainFieldForCellCoupling<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetCellKeys(std::vector<std::vector<std::string>> cellKeys)
    {
        mCellKeys = cellKeys;
    }

    template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
    void ChemicalDomainFieldForCellCoupling<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetCellKeyVector(std::vector<std::string> cellKeyVector)
    {
        mCellKeyVector = cellKeyVector;
    }

    template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
    void ChemicalDomainFieldForCellCoupling<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetCartesianCellLayerDimensions(std::vector<unsigned> cellLayerDims)
    {
        mCartesianCellLayerDimensions = cellLayerDims;
    }

    template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
    void ChemicalDomainFieldForCellCoupling<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetCartesianChasteCellScale(std::vector<double> chasteCellScale)
    {
        mCartesianChasteCellScaleXY = chasteCellScale;
    }

    template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
    void ChemicalDomainFieldForCellCoupling<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetCartesianChasteCellLayer(std::vector<std::vector<std::string>> cartesianChasteCellLayer)
    {
        mCartesianChasteCellLayer = cartesianChasteCellLayer;
    }

    template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
    std::vector<std::vector<std::vector<std::string>>> ChemicalDomainFieldForCellCoupling<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::GetCellLabels()
    {
        return mCellLabels;
    }

    template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
    std::string ChemicalDomainFieldForCellCoupling<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::GetCellLabelByIndex(unsigned index)
    {
        if(index<mCellNodes.size())
        {
            return mCellNodes[index];
        }
        else
        {
            std::cout<<"Error: ChemicalDomainFieldForCellCoupling::GetCellLabelByIndex: index out of bounds"<<std::endl;
            return "Error";
        }
        
    }

    template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
    std::vector<std::string> ChemicalDomainFieldForCellCoupling<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::GetCellLabelVector()
    {
        return mCellLabelVector;
    }

    template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
    std::vector<std::vector<std::string>> ChemicalDomainFieldForCellCoupling<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::GetCellKeys()
    {
        return mCellKeys;
    }

    template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
    std::vector<std::string> ChemicalDomainFieldForCellCoupling<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::GetCellKeyVector()
    {
        return mCellKeyVector;
    }

    template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
    void ChemicalDomainFieldForCellCoupling<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetCartesianChasteCellDimensions(std::vector<unsigned> dimensions)
    {
        mCartesianChasteCellDimensions = dimensions;
    }


    template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
    void ChemicalDomainFieldForCellCoupling<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetCellMesh(TetrahedralMesh<ELEMENT_DIM,SPACE_DIM>*& p_mesh)
    {
        mpCellMesh = p_mesh;
    }

    template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
    void ChemicalDomainFieldForCellCoupling<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetCellMeshGenerator(HoneycombMeshGenerator* p_generator)
    {
        mpCellMeshGenerator = p_generator;
    }


    template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
    void ChemicalDomainFieldForCellCoupling<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetCellMeshScale(std::vector<double> cellMeshScale)
    {
        mCellMeshScale = cellMeshScale;
    }

    template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
    void ChemicalDomainFieldForCellCoupling<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetCellMeshDimensions(std::vector<unsigned> meshDimensions)
    {
        mCellMeshDimensions = meshDimensions;
    }

    template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
    unsigned ChemicalDomainFieldForCellCoupling<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::GetNumberOfCellTypes()
    {
        return mNumberOfCellTypes;
    }

    template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
    TetrahedralMesh<ELEMENT_DIM,SPACE_DIM>*& ChemicalDomainFieldForCellCoupling<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::rGetCellMesh()
    {
        return mpCellMesh;
    }

    template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
    HoneycombMeshGenerator* ChemicalDomainFieldForCellCoupling<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::GetCellMeshGenerator()
    {
        return mpCellMeshGenerator;
    }

    template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
    std::vector<std::string> ChemicalDomainFieldForCellCoupling<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::GetNodeCells()
    {
        return mCellNodes;
    }


    template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
    boost::shared_ptr<InhomogenousParabolicPdeForCoupledOdeSystemTemplated<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>> ChemicalDomainFieldForCellCoupling<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::ReturnSharedPtrPdeSystem()
    {
        return mpPdeSystem;
    }
  
    //template<unsigned ELEMENT_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
    //boost::shared_ptr<BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM> > ChemicalDomainFieldForCellCoupling<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::ReturnSharedPtrBoundaryConditionsContainer()
    //{
     //   return mpBoundaryConditionsContainer;
    //}

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