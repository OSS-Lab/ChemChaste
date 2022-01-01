#ifndef CHEMICALCELLFROMFILE_HPP
#define CHEMICALCELLFROMFILE_HPP

#include "ChemicalCell.hpp"
#include "SimpleChemicalThresholdCellCycleFromFile.hpp"
#include "AbstractSrnModel.hpp"
#include "ChemicalSRNFromFile.hpp"
#include "AbstractCellCycleModel.hpp"
#include "NoCellCycleModel.hpp"
#include "InitialCellConditionsFromFile.hpp"
#include "TransportCellPropertyFromFile.hpp"
#include "MembraneCellPropertyFromFile.hpp"
#include "CellAnalyticsPropertyFromCellID.hpp"


class ChemicalCellFromFile
{
protected:

    CellPtr mpCell;

    CellPropertyCollection mPropertyCollection;

    boost::shared_ptr<ChemicalCellProperty> mp_cell_chemical;

    boost::shared_ptr<MembraneCellProperty> mp_cell_membrane;

    boost::shared_ptr<TransportCellProperty> mp_cell_transport;

    boost::shared_ptr<CellAnalyticsProperty> mp_cell_analytics;

    std::string mCellCycleFilename;

    bool mIsCellCycleSet = false;

    std::string mSrnFilename;

    bool mIsSRNSet = false;

    std::string mInitialConditionsFilename;

    bool mIsInitConditionsSet = false;

    std::string mTransportPropertyFilename;

    bool mIsTransportPropertySet = false;

    std::string mMembranePropertyFilename;

    bool mIsMembranePropertySet = false;

    unsigned mCellID;

    bool mIsCellIDSet = false;



    StateVariableRegister* mpFullChemicalStateRegister; 

    SimpleChemicalThresholdCellCycleModel*  mpSimpleChemicalThresholdCellCycleModel;

    ChemicalSrnModel*   mpChemicalSrnModel;

public:

    ChemicalCellFromFile(   std::string cellCycleFilename="", 
                            std::string srnFilename="",
                            std::string initialConditionFilename="",
                            std::string transportPropertyFilename="",
                            std::string membranePropertyFilename="",
                            unsigned cellID =0,
                            bool isCellIDSet = false);

    virtual ~ChemicalCellFromFile()
    {
    };

    void SetUpSRNandCellCycle();

    void SetUpCellObject();

    void SetUpCellProperties();

    void SetUpCellInitialConditions(CellPtr, std::vector<std::string>, std::vector<double>);


    // get methods

    CellPtr GetCellPtr();

    CellPropertyCollection GetCellPropertyCollection();
    
    boost::shared_ptr<ChemicalCellProperty> GetChemicalCellProperty();

    boost::shared_ptr<MembraneCellProperty> GetMembraneCellProperty();

    boost::shared_ptr<TransportCellProperty> GetTransportCellProperty();
    
    boost::shared_ptr<CellAnalyticsProperty> GetCellAnalyticsProperty();

    ChemicalSrnModel* GetChemicalSrnModel();

    SimpleChemicalThresholdCellCycleModel* GetChemicalCellCycleModel();

    std::string GetCellCycleFilename();

    bool GetIsCellCycleSet();

    std::string GetSrnFilename();

    bool GetIsSrnSet();

    std::string GetInitialConditionsFilename();

    bool GetInitConditionsSet();

    std::string GetTransportPropertyFilename();

    bool GetIsTransportPropertySet();

    std::string GetMembranePropertyFilename();

    bool GetIsMembranePropertySet();

    unsigned GetCellID();

    bool GetIsCellIDSet();



    StateVariableRegister* GetFullChemicalStateRegister();

    std::vector<std::string> GetFullChemicalNamesVector();


    // Set methods

    void SetCellPtr(CellPtr);

    void SetCellPropertyCollection(CellPropertyCollection);

    void SetChemicalCellProperty(boost::shared_ptr<ChemicalCellProperty>);

    void SetMembraneCellProperty(boost::shared_ptr<MembraneCellProperty>);

    void SetTransportCellProperty(boost::shared_ptr<TransportCellProperty>);

    void SetCellAnalyticsProperty(boost::shared_ptr<CellAnalyticsProperty>);

    void SetChemicalSrnModel(ChemicalSrnModel*);

    void SetChemicalCellCycleModel(SimpleChemicalThresholdCellCycleModel*);

    void SetCellCycleFilename(std::string);

    void SetSrnFilename(std::string);

    void SetInitialConditionsFilename(std::string);

    void SetTransportPropertyFilename(std::string);

    void SetMembranePropertyFilename(std::string);


    void SetFullChemicalStateRegister(StateVariableRegister*);


};

ChemicalCellFromFile::ChemicalCellFromFile(
                        std::string cellCycleFilename,
                        std::string srnFilename,
                        std::string initialConditionFilename,
                        std::string transportPropertyFilename,
                        std::string membranePropertyFilename,
                        unsigned cellID,
                        bool isCellIDSet)
        :   mCellCycleFilename(cellCycleFilename),
            mSrnFilename(srnFilename),
            mInitialConditionsFilename(initialConditionFilename),
            mTransportPropertyFilename(transportPropertyFilename),
            mMembranePropertyFilename(membranePropertyFilename),
            mCellID(cellID),
            mIsCellIDSet(isCellIDSet)
{

    if(mCellCycleFilename != "")
    {
        mIsCellCycleSet = true;
    }

    if(mSrnFilename != "")
    {
        mIsSRNSet = true;
    }

    if(mInitialConditionsFilename != "")
    {
        mIsInitConditionsSet = true;
        boost::shared_ptr<ChemicalCellProperty> p_cell_chemical(new ChemicalCellProperty());

        SetChemicalCellProperty(p_cell_chemical);
    }

    if(mTransportPropertyFilename != "")
    {
        mIsTransportPropertySet = true;
        boost::shared_ptr<TransportCellProperty> p_cell_transport(new TransportCellProperty());

        SetTransportCellProperty(p_cell_transport);
    }

    if(mMembranePropertyFilename != "")
    {
        mIsMembranePropertySet = true;
        boost::shared_ptr<MembraneCellProperty> p_cell_membrane(new MembraneCellProperty());

        SetMembraneCellProperty(p_cell_membrane);
    }

    
    if(mIsCellIDSet)
    {
        boost::shared_ptr<CellAnalyticsProperty> p_cell_analytics(new CellAnalyticsProperty());

        SetCellAnalyticsProperty(p_cell_analytics);
    }


    SetUpSRNandCellCycle();

    SetUpCellObject();

    SetUpCellProperties();

    // update the initial conditons of the cell

    // superset of all internal cell chemical species, if not defined in intial conditions then set concentration to 0

    // chemical species from initial conditions
    
    std::vector<std::string> cellChemicalNames = mp_cell_chemical -> GetStateVariableRegister() -> GetStateVariableRegisterVector();
    std::vector<double> cellConcentrationVector = mp_cell_chemical -> GetCellConcentrationVector();

    StateVariableRegister* pFullStateVariableRegister = new StateVariableRegister(cellChemicalNames);


    // chemical species from srn
    if(mIsSRNSet)
    {
        std::vector<std::string> srn_chemical_names = mpChemicalSrnModel -> GetCellChemistry() -> GetChemicalNames();
        pFullStateVariableRegister -> AddStateVariableVector(srn_chemical_names);
    }  

    // chemical species from cell cycle
    if(mIsCellCycleSet)
    {
        std::vector<std::string> cell_cycle_chemical_names = mpSimpleChemicalThresholdCellCycleModel -> GetThresholdChemistry() -> GetChemicalNames();
        pFullStateVariableRegister -> AddStateVariableVector(cell_cycle_chemical_names);
    }    

    // chemical species from membrane property
    if(mIsMembranePropertySet)
    {
        std::vector<std::string> cell_membrane_chemical_names = mp_cell_membrane -> GetCellStateVariableRegister() -> GetStateVariableRegisterVector();
        pFullStateVariableRegister -> AddStateVariableVector(cell_membrane_chemical_names);
    }    

    // chemical species from transport property
    if(mIsTransportPropertySet)
    {
        std::vector<std::string> cell_transport_chemical_names = mp_cell_transport -> GetCellStateVariableRegister() -> GetStateVariableRegisterVector();
        pFullStateVariableRegister -> AddStateVariableVector(cell_transport_chemical_names);
    } 


    // ensure all species are accounted for in initial conditions

    std::vector<std::string> fullChemicalNames = pFullStateVariableRegister -> GetStateVariableRegisterVector();
    if(cellChemicalNames.size() != fullChemicalNames.size())
    {
        for(unsigned i=cellChemicalNames.size(); i<fullChemicalNames.size(); i++)
        {
            cellConcentrationVector.push_back(0.0);
        }
    }
    
    SetFullChemicalStateRegister(pFullStateVariableRegister);

    SetUpCellInitialConditions(mpCell, cellChemicalNames, cellConcentrationVector);
}



void ChemicalCellFromFile::SetUpCellProperties()
{
    //std::cout<<"ChemicalCellFromFile::SetUpCellProperties - start"<<std::endl;
    if(mIsInitConditionsSet)
    {

        InitialCellConditionsFromFile* p_init_conditions_from_file = new InitialCellConditionsFromFile(mInitialConditionsFilename);
        
        StateVariableRegister* p_state_register = new StateVariableRegister(p_init_conditions_from_file -> GetChemicalNamesVector());

        mp_cell_chemical -> InitialiseCell(p_state_register,p_init_conditions_from_file -> GetConcentrationVector());

        SetUpCellInitialConditions(mpCell, p_state_register -> GetStateVariableRegisterVector(), p_init_conditions_from_file -> GetConcentrationVector());
  
    }

    if(mIsTransportPropertySet)
    {
        TransportCellPropertyFromFile* p_transport_property_from_file = new TransportCellPropertyFromFile(mTransportPropertyFilename);

        p_transport_property_from_file -> SetUpTransportProperty(mpCell);

        mp_cell_transport = p_transport_property_from_file -> GetTransportProperty();
    }

    if(mIsMembranePropertySet)
    {
        MembraneCellPropertyFromFile* p_membrane_property_from_file = new MembraneCellPropertyFromFile(mMembranePropertyFilename);

        p_membrane_property_from_file -> SetUpMembraneProperty(mpCell);

        mp_cell_membrane = p_membrane_property_from_file -> GetMembraneProperty();

        mp_cell_membrane -> SetMembraneThickness(5.0);
    }

    if(mIsCellIDSet)
    {
        CellAnalyticsPropertyFromCellID* p_cell_analytics_property_from_cellID = new CellAnalyticsPropertyFromCellID(mCellID);

        p_cell_analytics_property_from_cellID -> SetUpCellAnalyticsProperty(mpCell);

        mp_cell_analytics = p_cell_analytics_property_from_cellID -> GetCellAnalyticsProperty();
    }



    //std::cout<<"ChemicalCellFromFile::SetUpCellProperties - end"<<std::endl;
}

void ChemicalCellFromFile::SetUpSRNandCellCycle()
{
    if(mIsSRNSet)
    {
 
        ChemicalSRNFromFile* p_srn_reaction_system_from_file = new ChemicalSRNFromFile(mSrnFilename);
        SetChemicalSrnModel(p_srn_reaction_system_from_file -> GetChemicalSrnModel());


        AbstractChemistry* this_cell_chemistry = mpChemicalSrnModel->GetCellChemistry();
        unsigned numberOfChemicals = this_cell_chemistry->GetNumberChemicals();
  
    }


    if(mIsCellCycleSet)
    {
        SimpleChemicalThresholdCellCycleFromFile* pCellCycle = new SimpleChemicalThresholdCellCycleFromFile(mCellCycleFilename);

        pCellCycle -> SetUp();
        SetChemicalCellCycleModel(pCellCycle);
    }

}

void ChemicalCellFromFile::SetUpCellObject()
{
    std::cout<<"ChemicalCellFromFile::SetUpCellObject - start"<<std::endl;
    // form cell
    CellPropertyCollection collection;
    MAKE_PTR(WildTypeCellMutationState, p_state);
    if(mIsInitConditionsSet)
    {
        collection.AddProperty(mp_cell_chemical);
    }

    if(mIsTransportPropertySet)
    {
        collection.AddProperty(mp_cell_transport);
    }
    
    if(mIsMembranePropertySet)
    {
        collection.AddProperty(mp_cell_membrane);
    }

    if(mIsCellIDSet)
    {
        collection.AddProperty(mp_cell_analytics);
    }

    // check if srn and cell cycle are set, if not provide default models
    AbstractSrnModel* pSrnModel=nullptr;

    if(mIsSRNSet)
    {
        pSrnModel = GetChemicalSrnModel();
    }
    
    AbstractCellCycleModel* pCellCycleModel = new NoCellCycleModel();

    if(mIsCellCycleSet)
    {
        pCellCycleModel = GetChemicalCellCycleModel();
    }


    // call cell constructor
    CellPtr p_cell(new ChemicalCell(p_state, pCellCycleModel, pSrnModel, false, collection));

    // at oresent cell state and proliferation type is default
    MAKE_PTR(StemCellProliferativeType, p_stem_type);
    p_cell->SetCellProliferativeType(p_stem_type);
    
    SetCellPtr(p_cell);

    ChemicalSrnModel *p_chemical_srn = static_cast<ChemicalSrnModel*>(pSrnModel);
    std::cout<<"Distribute ########################"<<std::endl;
    p_chemical_srn->GetReactionSystem()->SetCell(p_cell); // should work or else type cast to ChemicalSrnModel*
    p_chemical_srn->GetReactionSystem()->DistributeCellPtr();
    std::cout<<"ChemicalCellFromFile::SetUpCellObject - end"<<std::endl;
}


void ChemicalCellFromFile::SetUpCellInitialConditions(CellPtr p_cell, std::vector<std::string> speciesNames, std::vector<double> initValue)
{
    for(unsigned i=0; i<speciesNames.size();i++)
    {
        p_cell->GetCellData()->SetItem(speciesNames[i],initValue[i]);
    }
}

// get methods

CellPtr ChemicalCellFromFile::GetCellPtr()
{
    return mpCell;
}

CellPropertyCollection ChemicalCellFromFile::GetCellPropertyCollection()
{
    return mPropertyCollection;
}

boost::shared_ptr<ChemicalCellProperty> ChemicalCellFromFile::GetChemicalCellProperty()
{
    return mp_cell_chemical;
}

boost::shared_ptr<MembraneCellProperty> ChemicalCellFromFile::GetMembraneCellProperty()
{
    return mp_cell_membrane;
}

boost::shared_ptr<TransportCellProperty> ChemicalCellFromFile::GetTransportCellProperty()
{
    return mp_cell_transport;
}

boost::shared_ptr<CellAnalyticsProperty> ChemicalCellFromFile::GetCellAnalyticsProperty()
{
    return mp_cell_analytics;
}


ChemicalSrnModel* ChemicalCellFromFile::GetChemicalSrnModel()
{
    return mpChemicalSrnModel;
}

SimpleChemicalThresholdCellCycleModel* ChemicalCellFromFile::GetChemicalCellCycleModel()
{
    return mpSimpleChemicalThresholdCellCycleModel;
}

std::string ChemicalCellFromFile::GetCellCycleFilename()
{
    return mCellCycleFilename;
}

bool ChemicalCellFromFile::GetIsCellCycleSet()
{
    return mIsCellCycleSet;
}

std::string ChemicalCellFromFile::GetSrnFilename()
{
    return mSrnFilename;
}

bool ChemicalCellFromFile::GetIsSrnSet()
{
    return mIsSRNSet;
}

std::string ChemicalCellFromFile::GetInitialConditionsFilename()
{
    return mInitialConditionsFilename;
}

bool ChemicalCellFromFile::GetInitConditionsSet()
{
    return mIsInitConditionsSet;
}

std::string ChemicalCellFromFile::GetTransportPropertyFilename()
{
    return mTransportPropertyFilename;
}

bool ChemicalCellFromFile::GetIsTransportPropertySet()
{
    return mIsTransportPropertySet;
}

std::string ChemicalCellFromFile::GetMembranePropertyFilename()
{
    return mMembranePropertyFilename;
}

bool ChemicalCellFromFile::GetIsMembranePropertySet()
{
    return mIsMembranePropertySet;
}

unsigned ChemicalCellFromFile::GetCellID()
{
    return mCellID;
}

bool ChemicalCellFromFile::GetIsCellIDSet()
{
    return mIsCellIDSet;
}



StateVariableRegister* ChemicalCellFromFile::GetFullChemicalStateRegister()
{
    return mpFullChemicalStateRegister;
}

std::vector<std::string> ChemicalCellFromFile::GetFullChemicalNamesVector()
{
    return mpFullChemicalStateRegister -> GetStateVariableRegisterVector();
}


// set methods

void ChemicalCellFromFile::SetCellPtr(CellPtr p_cell)
{
    mpCell = p_cell;
}

void ChemicalCellFromFile::SetCellPropertyCollection(CellPropertyCollection propertyCollection)
{
    mPropertyCollection = propertyCollection;
}

void ChemicalCellFromFile::SetChemicalCellProperty(boost::shared_ptr<ChemicalCellProperty> p_chemical)
{
    mp_cell_chemical = p_chemical;
}

void ChemicalCellFromFile::SetMembraneCellProperty(boost::shared_ptr<MembraneCellProperty> p_membrane)
{
    mp_cell_membrane = p_membrane;
}

void ChemicalCellFromFile::SetTransportCellProperty(boost::shared_ptr<TransportCellProperty> p_transport)
{
    mp_cell_transport = p_transport;
}

void ChemicalCellFromFile::SetCellAnalyticsProperty(boost::shared_ptr<CellAnalyticsProperty> p_cellAnalytics)
{
    mp_cell_analytics = p_cellAnalytics;
}

void ChemicalCellFromFile::SetChemicalSrnModel(ChemicalSrnModel* p_chemicalSrn)
{
    mpChemicalSrnModel = p_chemicalSrn;
}

void ChemicalCellFromFile::SetChemicalCellCycleModel(SimpleChemicalThresholdCellCycleModel* p_chemicalCellCycleModel)
{
    mpSimpleChemicalThresholdCellCycleModel = p_chemicalCellCycleModel;
}




void ChemicalCellFromFile::SetCellCycleFilename(std::string filename)
{
    mCellCycleFilename = filename;
    mIsCellCycleSet = true;
}

void ChemicalCellFromFile::SetSrnFilename(std::string filename)
{
    mSrnFilename = filename;
    mIsSRNSet = true;
}

void ChemicalCellFromFile::SetInitialConditionsFilename(std::string filename)
{
    mInitialConditionsFilename = filename;
    mIsInitConditionsSet = true;
}

void ChemicalCellFromFile::SetTransportPropertyFilename(std::string filename)
{
    mTransportPropertyFilename = filename;
    mIsTransportPropertySet = true;
}

void ChemicalCellFromFile::SetMembranePropertyFilename(std::string filename)
{
    mMembranePropertyFilename = filename;
    mIsMembranePropertySet = true;
}

void ChemicalCellFromFile::SetFullChemicalStateRegister(StateVariableRegister* p_register)
{
    mpFullChemicalStateRegister = p_register;
}


#endif