#ifndef COMPLEXCELLFROMFILE_HPP
#define COMPLEXCELLFROMFILE_HPP

#include "ComplexCell.hpp"
#include "SimpleChemicalThresholdCellCycleFromFile.hpp"
#include "AbstractSrnModel.hpp"
#include "ChemicalSRNFromFile.hpp"
#include "AbstractCellCycleModel.hpp"
#include "NoCellCycleModel.hpp"
#include "InitialCellConditionsFromFile.hpp"
#include "TransportCellPropertyFromFile.hpp"
#include "MembraneCellPropertyFromFile.hpp"
#include "EnvironmentCellPropertyFromFile.hpp"
#include "CellAnalyticsPropertyFromCellID.hpp"


class ComplexCellFromFile
{
protected:

    CellPtr mpCell;

    CellPropertyCollection mPropertyCollection;

    boost::shared_ptr<ChemicalCellProperty> mp_cell_chemical;

    boost::shared_ptr<MembraneCellProperty> mp_cell_membrane;

    boost::shared_ptr<TransportCellProperty> mp_cell_transport;

    boost::shared_ptr<CellAnalyticsProperty> mp_cell_analytics;

    boost::shared_ptr<EnvironmentCellProperty> mp_cell_environment;

    std::string mCellFileRoot;

    std::string mCellCycleFilename;

    bool mIsCellCycleSet = false;

    std::string mDivisionRulesFilename;

    bool mIsDivisionRulesSet = false;

    std::string mSrnFilename;

    bool mIsSRNSet = false;

    std::string mInitialConditionsFilename;

    bool mIsInitConditionsSet = false;

    std::string mTransportPropertyFilename;

    bool mIsTransportPropertySet = false;

    std::string mMembranePropertyFilename;

    bool mIsMembranePropertySet = false;

    std::string mEnvironmentPropertyFilename;

    bool mIsEnvironmentPropertySet = false;

    unsigned mCellID;

    bool mIsCellIDSet = false;

    std::string mCellTypeName="";

    StateVariableRegister* mpFullChemicalStateRegister; 

    SimpleChemicalThresholdCellCycleModel*  mpSimpleChemicalThresholdCellCycleModel;

    ChemicalSrnModel*   mpChemicalSrnModel;

    double mSplitRatio = 0.5;

public:

    ComplexCellFromFile(    std::string cellFileRoot="",
                            std::string cellCycleFilename="",
                            std::string divisionRulesFilename="", 
                            std::string srnFilename="",
                            std::string initialConditionFilename="",
                            std::string transportPropertyFilename="",
                            std::string membranePropertyFilename="",
                            std::string environmentPropertyFilename="",
                            unsigned cellID =0,
                            bool isCellIDSet = false);

    virtual ~ComplexCellFromFile()
    {
    };

    virtual void SetUp();

    virtual void BuildCellPropertyCollection(CellPropertyCollection& collection);

    virtual void SetUpSRNandCellCycle();

    void SetUpCellObject();

    void SetUpCellProperties();

    void SetUpCellInitialConditions(CellPtr, std::vector<std::string>, std::vector<double>);

    virtual void SetUpCellDivisionRules(CellPtr);

    std::vector<std::string> parseMatrixLineString(std::string);

    std::vector<std::vector<std::string>> ReadMatrix(std::string);

    // get methods
    CellPtr GetCellPtr();

    CellPropertyCollection GetCellPropertyCollection();
    
    boost::shared_ptr<ChemicalCellProperty> GetChemicalCellProperty();

    boost::shared_ptr<MembraneCellProperty> GetMembraneCellProperty();

    boost::shared_ptr<EnvironmentCellProperty> GetEnvironmentCellProperty();

    boost::shared_ptr<TransportCellProperty> GetTransportCellProperty();
    
    boost::shared_ptr<CellAnalyticsProperty> GetCellAnalyticsProperty();

    ChemicalSrnModel* GetChemicalSrnModel();

    SimpleChemicalThresholdCellCycleModel* GetChemicalCellCycleModel();

    std::string GetCellFileRoot();

    std::string GetCellCycleFilename();

    bool GetIsCellCycleSet();

    std::string GetDivisionRulesFilename();

    bool GetIsDivisionRulesSet();

    std::string GetSrnFilename();

    bool GetIsSrnSet();

    std::string GetInitialConditionsFilename();

    bool GetInitConditionsSet();

    std::string GetTransportPropertyFilename();

    bool GetIsTransportPropertySet();

    std::string GetMembranePropertyFilename();

    bool GetIsMembranePropertySet();

    std::string GetEnvironmentPropertyFilename();

    bool GetIsEnvironmentPropertySet();

    unsigned GetCellID();

    bool GetIsCellIDSet();

    StateVariableRegister* GetFullChemicalStateRegister();

    std::vector<std::string> GetFullChemicalNamesVector();


    // Set methods
    void SetCellPtr(CellPtr);

    void SetCellPropertyCollection(CellPropertyCollection);

    void SetChemicalCellProperty(boost::shared_ptr<ChemicalCellProperty>);

    void SetMembraneCellProperty(boost::shared_ptr<MembraneCellProperty>);

    void SetEnvironmentCellProperty(boost::shared_ptr<EnvironmentCellProperty>);

    void SetTransportCellProperty(boost::shared_ptr<TransportCellProperty>);

    void SetCellAnalyticsProperty(boost::shared_ptr<CellAnalyticsProperty>);

    void SetChemicalSrnModel(ChemicalSrnModel*);

    void SetChemicalCellCycleModel(SimpleChemicalThresholdCellCycleModel*);

    void SetCellFileRoot(std::string);

    void SetCellCycleFilename(std::string);

    void SetDivisionRulesFilename(std::string);

    void SetSrnFilename(std::string);

    void SetInitialConditionsFilename(std::string);

    void SetTransportPropertyFilename(std::string);

    void SetMembranePropertyFilename(std::string);

    void SetEnvironmentPropertyFilename(std::string);

    void SetFullChemicalStateRegister(StateVariableRegister*);

    void SetCellTypeName(std::string);

};

ComplexCellFromFile::ComplexCellFromFile(
                        std::string cellFileRoot,
                        std::string cellCycleFilename,
                        std::string divisionRulesFilename,
                        std::string srnFilename,
                        std::string initialConditionFilename,
                        std::string transportPropertyFilename,
                        std::string membranePropertyFilename,
                        std::string environmentPropertyFilename,
                        unsigned cellID,
                        bool isCellIDSet)
        :   mCellFileRoot(cellFileRoot),
            mCellCycleFilename(cellCycleFilename),
            mDivisionRulesFilename(divisionRulesFilename),
            mSrnFilename(srnFilename),
            mInitialConditionsFilename(initialConditionFilename),
            mTransportPropertyFilename(transportPropertyFilename),
            mMembranePropertyFilename(membranePropertyFilename),
            mEnvironmentPropertyFilename(environmentPropertyFilename),
            mCellID(cellID),
            mIsCellIDSet(isCellIDSet)
{
    SetFullChemicalStateRegister(new StateVariableRegister(std::vector<std::string>()));
}

void ComplexCellFromFile::SetUp()
{

    if(mCellCycleFilename != "")
    {
        mIsCellCycleSet = true;
        mCellCycleFilename = mCellFileRoot+mCellCycleFilename;
    }
    //std::cout<<"division filename: "<<mDivisionRulesFilename<<std::endl;
    if(mDivisionRulesFilename != "")
    {
        mIsDivisionRulesSet = true;
        mDivisionRulesFilename = mCellFileRoot+mDivisionRulesFilename;
    }


    if(mSrnFilename != "")
    {
        mIsSRNSet = true;
        mSrnFilename = mCellFileRoot+mSrnFilename;
    }

    if(mInitialConditionsFilename != "")
    {
        
        mIsInitConditionsSet = true;
        mInitialConditionsFilename = mCellFileRoot+mInitialConditionsFilename;
        boost::shared_ptr<ChemicalCellProperty> p_cell_chemical(new ChemicalCellProperty());
        std::cout<<"mInitialConditionsFilename: "<<mInitialConditionsFilename<<std::endl;
        SetChemicalCellProperty(p_cell_chemical);
    }

    if(mTransportPropertyFilename != "")
    {
        
  //      std::cout<<"################test if there is and set transprot"<<std::endl;
        mIsTransportPropertySet = true;
        mTransportPropertyFilename = mCellFileRoot+mTransportPropertyFilename;
        boost::shared_ptr<TransportCellProperty> p_cell_transport(new TransportCellProperty());
        std::cout<<"mTransportPropertyFilename: "<<mTransportPropertyFilename<<std::endl;
        SetTransportCellProperty(p_cell_transport);
    }

    if(mMembranePropertyFilename != "")
    {
   //     std::cout<<"################test if there is and set membrane"<<std::endl;
        mIsMembranePropertySet = true;
        mMembranePropertyFilename = mCellFileRoot+mMembranePropertyFilename;
        boost::shared_ptr<MembraneCellProperty> p_cell_membrane(new MembraneCellProperty());
        std::cout<<"mMembranePropertyFilename: "<<mMembranePropertyFilename<<std::endl;
        SetMembraneCellProperty(p_cell_membrane);
    }

    if(mEnvironmentPropertyFilename != "")
    {
        mIsEnvironmentPropertySet = true;
        mEnvironmentPropertyFilename = mCellFileRoot+mEnvironmentPropertyFilename;
        boost::shared_ptr<EnvironmentCellProperty> p_cell_environment(new EnvironmentCellProperty());
        std::cout<<"mEnvironmentPropertyFilename: "<<mEnvironmentPropertyFilename<<std::endl;
        SetEnvironmentCellProperty(p_cell_environment);
    }

    
    if(mIsCellIDSet)
    {
        boost::shared_ptr<CellAnalyticsProperty> p_cell_analytics(new CellAnalyticsProperty());

        SetCellAnalyticsProperty(p_cell_analytics);
    }

    // all good to here

    SetUpSRNandCellCycle();

    SetUpCellObject();

    SetUpCellProperties();

    // update the initial conditons of the cell

    // superset of all internal cell chemical species, if not defined in intial conditions then set concentration to 0

    // chemical species from initial conditions
    
    std::vector<std::string> cellChemicalNames = mp_cell_chemical -> GetStateVariableRegister() -> GetStateVariableRegisterVector();
    std::vector<double> cellConcentrationVector = mp_cell_chemical -> GetCellConcentrationVector();

    StateVariableRegister* pFullStateVariableRegister = GetFullChemicalStateRegister();
    
    pFullStateVariableRegister -> AddStateVariableVector(cellChemicalNames);


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


void ComplexCellFromFile::SetUpCellProperties()
{
    if(mIsInitConditionsSet)
    {

        InitialCellConditionsFromFile* p_init_conditions_from_file = new InitialCellConditionsFromFile(mInitialConditionsFilename);
        
        StateVariableRegister* p_state_register = new StateVariableRegister(p_init_conditions_from_file -> GetChemicalNamesVector());

        mp_cell_chemical -> InitialiseCell(p_state_register,p_init_conditions_from_file -> GetConcentrationVector());

        SetUpCellInitialConditions(mpCell, p_state_register -> GetStateVariableRegisterVector(), p_init_conditions_from_file -> GetConcentrationVector());
    }

    if(mIsTransportPropertySet)
    {
      //  std::cout<<"here transport property is set"<<std::endl;
        TransportCellPropertyFromFile* p_transport_property_from_file = new TransportCellPropertyFromFile(mTransportPropertyFilename);

        p_transport_property_from_file -> SetUpTransportProperty(mpCell);

        mp_cell_transport = p_transport_property_from_file -> GetTransportProperty();
    }

    if(mIsMembranePropertySet)
    {
     //   std::cout<<"here membrane property is set"<<std::endl;
        MembraneCellPropertyFromFile* p_membrane_property_from_file = new MembraneCellPropertyFromFile(mMembranePropertyFilename);

        p_membrane_property_from_file -> SetUpMembraneProperty(mpCell);

        mp_cell_membrane = p_membrane_property_from_file -> GetMembraneProperty();

        mp_cell_membrane -> SetMembraneThickness(5.0);
    }

    if(mIsEnvironmentPropertySet)
    {
        //std::cout<<"here environment property is set"<<std::endl;
        EnvironmentCellPropertyFromFile* p_environment_property_from_file = new EnvironmentCellPropertyFromFile(mEnvironmentPropertyFilename);

        p_environment_property_from_file -> SetUpEnvironmentProperty(mpCell);

        mp_cell_environment = p_environment_property_from_file -> GetEnvironmentProperty();
        mp_cell_environment->SetCellPtr(mpCell);

    }

    if(mIsCellIDSet)
    {
        //std::cout<<"CellID: "<<mCellID<<std::endl;
        CellAnalyticsPropertyFromCellID* p_cell_analytics_property_from_cellID = new CellAnalyticsPropertyFromCellID(mCellID,mCellTypeName);

        p_cell_analytics_property_from_cellID -> SetUpCellAnalyticsProperty(mpCell);

        mp_cell_analytics = p_cell_analytics_property_from_cellID -> GetCellAnalyticsProperty();
    }

}

void ComplexCellFromFile::SetUpSRNandCellCycle()
{

    if(mIsSRNSet)
    {
        if(mSrnFilename != "")
        {
            std::cout<<"mSrnFileName: "<<mSrnFilename<<std::endl;
            ChemicalSRNFromFile* p_srn_reaction_system_from_file = new ChemicalSRNFromFile(mSrnFilename);
            SetChemicalSrnModel(p_srn_reaction_system_from_file -> GetChemicalSrnModel());
        }
    }


    if(mIsCellCycleSet && mCellCycleFilename != "")
    {
        SimpleChemicalThresholdCellCycleFromFile* pCellCycle = new SimpleChemicalThresholdCellCycleFromFile(mCellCycleFilename);

        pCellCycle -> SetUp();
        SetChemicalCellCycleModel(pCellCycle);
    }

}

void ComplexCellFromFile::SetUpCellObject()
{

    // form cell
    CellPropertyCollection collection;
    MAKE_PTR(WildTypeCellMutationState, p_state);
    

    BuildCellPropertyCollection(collection);

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
    CellPtr p_cell(new ComplexCell(p_state, pCellCycleModel, pSrnModel, false, collection));
    
    // at present cell state and proliferation type is default
    MAKE_PTR(StemCellProliferativeType, p_stem_type);
    p_cell->SetCellProliferativeType(p_stem_type);



    ChemicalSrnModel *p_chemical_srn = static_cast<ChemicalSrnModel*>(pSrnModel);
    
    p_chemical_srn->GetReactionSystem()->SetCell(p_cell); // should work or else type cast to ChemicalSrnModel*
    p_chemical_srn->GetReactionSystem()->DistributeCellPtr();


    if(mIsDivisionRulesSet)
    {
        SetUpCellDivisionRules(p_cell);
    }

    
    SetCellPtr(p_cell);

}


void ComplexCellFromFile::BuildCellPropertyCollection(CellPropertyCollection& collection)
{

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

    if(mIsEnvironmentPropertySet)
    {
        collection.AddProperty(mp_cell_environment);
    }

    if(mIsCellIDSet)
    {
        collection.AddProperty(mp_cell_analytics);
    }
}


void ComplexCellFromFile::SetUpCellInitialConditions(CellPtr p_cell, std::vector<std::string> speciesNames, std::vector<double> initValue)
{
    for(unsigned i=0; i<speciesNames.size();i++)
    {
        p_cell->GetCellData()->SetItem(speciesNames[i],initValue[i]);
    }
}

void ComplexCellFromFile::SetUpCellDivisionRules(CellPtr p_cell)
{
    // function to read the division rules from a file and supply them to the cell model

    std::vector<std::string> chemicalNames;
    std::vector<std::string> chemicalDivisionRules;

    std::vector<std::vector<std::string>> divisionDatabase;

    divisionDatabase = ReadMatrix(mDivisionRulesFilename);

    unsigned numberOfChemicals = divisionDatabase.size();
    //unsigned numberOfitems = divisionDatabase[0].size();

    for(unsigned row=0; row<numberOfChemicals;row++)
    {
        //for(unsigned item=0; item<numberOfitems; item++)
        //{
        chemicalNames.push_back(divisionDatabase[row][0]);
        chemicalDivisionRules.push_back(divisionDatabase[row][1]);
        //}
    }
    //ComplexCell* p_cell_new = dynamic_cast<ComplexCell*>(p_cell);
    boost::shared_ptr<ComplexCell> p_cell_new = boost::static_pointer_cast<ComplexCell>(p_cell);
    p_cell_new->SetChemicalNames(chemicalNames);

    p_cell_new->SetChemicalDivsionRules(chemicalDivisionRules);
}

std::vector<std::vector<std::string>> ComplexCellFromFile::ReadMatrix(std::string filename)
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


std::vector<std::string> ComplexCellFromFile::parseMatrixLineString(std::string line)
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


// get methods

CellPtr ComplexCellFromFile::GetCellPtr()
{
    return mpCell;
}

CellPropertyCollection ComplexCellFromFile::GetCellPropertyCollection()
{
    return mPropertyCollection;
}

boost::shared_ptr<ChemicalCellProperty> ComplexCellFromFile::GetChemicalCellProperty()
{
    return mp_cell_chemical;
}

boost::shared_ptr<MembraneCellProperty> ComplexCellFromFile::GetMembraneCellProperty()
{
    return mp_cell_membrane;
}

boost::shared_ptr<EnvironmentCellProperty> ComplexCellFromFile::GetEnvironmentCellProperty()
{
    return mp_cell_environment;
}

boost::shared_ptr<TransportCellProperty> ComplexCellFromFile::GetTransportCellProperty()
{
    return mp_cell_transport;
}

boost::shared_ptr<CellAnalyticsProperty> ComplexCellFromFile::GetCellAnalyticsProperty()
{
    return mp_cell_analytics;
}


ChemicalSrnModel* ComplexCellFromFile::GetChemicalSrnModel()
{
    return mpChemicalSrnModel;
}

SimpleChemicalThresholdCellCycleModel* ComplexCellFromFile::GetChemicalCellCycleModel()
{
    return mpSimpleChemicalThresholdCellCycleModel;
}

std::string ComplexCellFromFile::GetCellFileRoot()
{
    return mCellFileRoot;
}

std::string ComplexCellFromFile::GetCellCycleFilename()
{
    return mCellCycleFilename;
}

bool ComplexCellFromFile::GetIsCellCycleSet()
{
    return mIsCellCycleSet;
}

std::string ComplexCellFromFile::GetSrnFilename()
{
    return mSrnFilename;
}

bool ComplexCellFromFile::GetIsSrnSet()
{
    return mIsSRNSet;
}

std::string ComplexCellFromFile::GetInitialConditionsFilename()
{
    return mInitialConditionsFilename;
}

bool ComplexCellFromFile::GetInitConditionsSet()
{
    return mIsInitConditionsSet;
}

std::string ComplexCellFromFile::GetTransportPropertyFilename()
{
    return mTransportPropertyFilename;
}

bool ComplexCellFromFile::GetIsTransportPropertySet()
{
    return mIsTransportPropertySet;
}

std::string ComplexCellFromFile::GetMembranePropertyFilename()
{
    return mMembranePropertyFilename;
}

bool ComplexCellFromFile::GetIsMembranePropertySet()
{
    return mIsMembranePropertySet;
}

std::string ComplexCellFromFile::GetEnvironmentPropertyFilename()
{
    return mEnvironmentPropertyFilename;
}

bool ComplexCellFromFile::GetIsEnvironmentPropertySet()
{
    return mIsEnvironmentPropertySet;
}

unsigned ComplexCellFromFile::GetCellID()
{
    return mCellID;
}

bool ComplexCellFromFile::GetIsCellIDSet()
{
    return mIsCellIDSet;
}



StateVariableRegister* ComplexCellFromFile::GetFullChemicalStateRegister()
{
    return mpFullChemicalStateRegister;
}

std::vector<std::string> ComplexCellFromFile::GetFullChemicalNamesVector()
{
    return mpFullChemicalStateRegister -> GetStateVariableRegisterVector();
}


// set methods

void ComplexCellFromFile::SetCellPtr(CellPtr p_cell)
{
    mpCell = p_cell;
}

void ComplexCellFromFile::SetCellPropertyCollection(CellPropertyCollection propertyCollection)
{
    mPropertyCollection = propertyCollection;
}

void ComplexCellFromFile::SetChemicalCellProperty(boost::shared_ptr<ChemicalCellProperty> p_chemical)
{
    mp_cell_chemical = p_chemical;
}

void ComplexCellFromFile::SetMembraneCellProperty(boost::shared_ptr<MembraneCellProperty> p_membrane)
{
    mp_cell_membrane = p_membrane;
}

void ComplexCellFromFile::SetEnvironmentCellProperty(boost::shared_ptr<EnvironmentCellProperty> p_environment)
{
    mp_cell_environment = p_environment;
}

void ComplexCellFromFile::SetTransportCellProperty(boost::shared_ptr<TransportCellProperty> p_transport)
{
    mp_cell_transport = p_transport;
}

void ComplexCellFromFile::SetCellAnalyticsProperty(boost::shared_ptr<CellAnalyticsProperty> p_cellAnalytics)
{
    mp_cell_analytics = p_cellAnalytics;
}

void ComplexCellFromFile::SetChemicalSrnModel(ChemicalSrnModel* p_chemicalSrn)
{
    mpChemicalSrnModel = p_chemicalSrn;
}

void ComplexCellFromFile::SetChemicalCellCycleModel(SimpleChemicalThresholdCellCycleModel* p_chemicalCellCycleModel)
{
    mpSimpleChemicalThresholdCellCycleModel = p_chemicalCellCycleModel;
}


void ComplexCellFromFile::SetCellFileRoot(std::string fileRoot)
{
    mCellFileRoot = fileRoot;
}

void ComplexCellFromFile::SetCellCycleFilename(std::string filename)
{
    mCellCycleFilename = filename;
    mIsCellCycleSet = true;
}

void ComplexCellFromFile::SetSrnFilename(std::string filename)
{
    mSrnFilename = filename;
    mIsSRNSet = true;
}

void ComplexCellFromFile::SetInitialConditionsFilename(std::string filename)
{
    mInitialConditionsFilename = filename;
    mIsInitConditionsSet = true;
}

void ComplexCellFromFile::SetTransportPropertyFilename(std::string filename)
{
    mTransportPropertyFilename = filename;
    mIsTransportPropertySet = true;
}

void ComplexCellFromFile::SetMembranePropertyFilename(std::string filename)
{
    mMembranePropertyFilename = filename;
    mIsMembranePropertySet = true;
}

void ComplexCellFromFile::SetEnvironmentPropertyFilename(std::string filename)
{
    mEnvironmentPropertyFilename = filename;
    mIsEnvironmentPropertySet = true;
}

void ComplexCellFromFile::SetFullChemicalStateRegister(StateVariableRegister* p_register)
{
    mpFullChemicalStateRegister = p_register;
}

void ComplexCellFromFile::SetCellTypeName(std::string typeName)
{
    mCellTypeName = typeName;
}



#endif