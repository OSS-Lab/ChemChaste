
/*#include "ComplexCell.hpp"
#include "Cell_virtual.hpp"
#include <vector>
#include <string>
#include <boost/utility.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/enable_shared_from_this.hpp>
#include "CellId.hpp"
#include "ChasteSerialization.hpp"
#include <boost/serialization/shared_ptr.hpp>

#include "AbstractCellMutationState.hpp"
#include "AbstractCellProliferativeType.hpp"
#include "CellData.hpp"
#include "CellVecData.hpp"
#include "ChemicalCell.hpp"

#include "AbstractCellCycleModel.hpp"
#include "AbstractSrnModel.hpp"
#include "CellPropertyCollection.hpp"
#include "ChemicalSrnModel.hpp"
#include "SimpleChemicalThresholdCellCycleModel.hpp"

#include "AbstractChemistry.hpp"
#include "MembraneCellProperty.hpp"
#include "TransportCellProperty.hpp"
#include "CellAnalyticsProperty.hpp"
#include "EnvironmentCellProperty.hpp"
//class StateSwitchingCellProperty;
//#include "StateSwitchingCellProperty.hpp"


ComplexCell::ComplexCell(
        boost::shared_ptr<AbstractCellProperty> pMutationState,
        AbstractCellCycleModel* pCellCycleModel,
        AbstractSrnModel* pSrnModel,
        bool archiving,
        CellPropertyCollection cellPropertyCollection)
    : ChemicalCell(pMutationState,pCellCycleModel,pSrnModel,archiving,cellPropertyCollection),
        mSplitRatio(0.5)

         {

         };

CellPtr ComplexCell::Divide()
{
    std::cout<<"ComplexCell::Divide()#######################################"<<std::endl;
    // Check we're allowed to divide
    assert(!IsDead());
    assert(mCanDivide);
    mCanDivide = false;

    // Reset properties of parent cell
    mpCellCycleModel->ResetForDivision();
    mpSrnModel->ResetForDivision();

    // Create copy of cell property collection to modify for daughter cell
    CellPropertyCollection daughter_property_collection = mCellPropertyCollection;

    // Remove the CellId from the daughter cell, as a new one will be assigned in the constructor
    daughter_property_collection.RemoveProperty<CellId>();

    // copy cell data
    std::cout<<"here0"<<std::endl;
    // Copy all cell data (note we create a new object not just copying the pointer)
    assert(daughter_property_collection.HasPropertyType<CellData>());
    // Get the existing copy of the cell data and remove it from the daughter cell
    boost::shared_ptr<CellData> p_cell_data = GetCellData();
    daughter_property_collection.RemoveProperty(p_cell_data);

    boost::shared_ptr<CellData> p_daughter_cell_data(new CellData(*p_cell_data));
    // Create a new cell data object using the copy constructor and add this to the daughter cell
    //MAKE_PTR_ARGS(CellData, p_daughter_cell_data, (*p_cell_data));
    daughter_property_collection.AddProperty(p_daughter_cell_data);

    // Copy all cell Vec data (note we create a new object not just copying the pointer)
    if (daughter_property_collection.HasPropertyType<CellVecData>())
    {
        // Get the existing copy of the cell data and remove it from the daughter cell
        boost::shared_ptr<CellVecData> p_cell_vec_data = GetCellVecData();
        daughter_property_collection.RemoveProperty(p_cell_vec_data);
        // Create a new cell data object using the copy constructor and add this to the daughter cell
        boost::shared_ptr<CellVecData> p_daughter_cell_vec_data(new CellVecData(*p_cell_vec_data));
        //MAKE_PTR_ARGS(CellVecData, p_daughter_cell_vec_data, (*p_cell_vec_data));
        daughter_property_collection.AddProperty(p_daughter_cell_vec_data);
    }
std::cout<<"here1"<<std::endl;

    // record the cell chemistry for splitting cell data
    AbstractChemistry* cellChemistry = new AbstractChemistry();

    if(static_cast<ChemicalSrnModel*>(mpSrnModel)->SRNType()=="Chemical")
    {
        ChemicalSrnModel* p_srn_model = static_cast<ChemicalSrnModel*>(mpSrnModel);
        cellChemistry -> AddChemistry(p_srn_model->GetCellChemistry()); // from SRN
    }
std::cout<<"here2"<<std::endl;
    if(static_cast<SimpleChemicalThresholdCellCycleModel*>(mpCellCycleModel)->CellCycleType()=="Chemical")
    {
        SimpleChemicalThresholdCellCycleModel* p_cc_model = static_cast<SimpleChemicalThresholdCellCycleModel*>(mpCellCycleModel);
        cellChemistry -> AddChemistry(p_cc_model->GetThresholdChemistry()); // from cell cycle model
    }
std::cout<<"here3"<<std::endl;
    // transport property
    if(mCellPropertyCollection.HasProperty<TransportCellProperty>())
    {
        boost::shared_ptr<TransportCellProperty> transport_cell_property = boost::static_pointer_cast<TransportCellProperty>(mCellPropertyCollection.GetPropertiesType<TransportCellProperty>().GetProperty());
        AbstractChemistry* transportChemistry = transport_cell_property ->GetTransportReactionSystem()->GetCellChemistry();
        // add the cell chemistry due to transport
        cellChemistry -> AddChemistry(transportChemistry);
    }
std::cout<<"here4"<<std::endl;
    // membrane property
    if(mCellPropertyCollection.HasProperty<MembraneCellProperty>())
    {
        boost::shared_ptr<MembraneCellProperty> membrane_cell_property = boost::static_pointer_cast<MembraneCellProperty>(mCellPropertyCollection.GetPropertiesType<MembraneCellProperty>().GetProperty());
        AbstractChemistry* membraneChemistry = membrane_cell_property ->GetMembraneReactionSystem()->GetCellChemistry();
        // add the cell chemistry due to membrane
        cellChemistry -> AddChemistry(membraneChemistry);
    }

std::cout<<"here5"<<std::endl;
    // based on the parent cell data determine the split ratio
    DetermineSplitRatio();

    // share the two cellDatas between the two cells, use cellChemistry
    // split chemical cell data
    double parent_species_concentration=0.0;
    double new_parent_species_concentration=0.0;
    double daughter_species_concentration=0.0;
    unsigned numberOfChemicals = cellChemistry -> GetNumberChemicals();
std::cout<<"here6"<<std::endl;
    for(unsigned i=0; i<numberOfChemicals; i++)
    {
        
        parent_species_concentration = this->GetCellData()->GetItem(cellChemistry -> GetChemicalNamesByIndex(i));
        

        if(IsChemicalShared(cellChemistry -> GetChemicalNamesByIndex(i)))
        {
            std::cout<<"shared"<<std::endl;
            // chemical concentration is shared upon division
            new_parent_species_concentration = SplitParentCellData(parent_species_concentration);

            daughter_species_concentration = parent_species_concentration - new_parent_species_concentration;
        }
        else
        {
            std::cout<<"duplicate"<<std::endl;
            // chemical concentration is duplicated upon division
            new_parent_species_concentration = parent_species_concentration;
            daughter_species_concentration = parent_species_concentration;
        }
        
        // check for positive concentrations
        if(daughter_species_concentration<0.0)
        {
            // must be positive concentration
            daughter_species_concentration =0.0;
        }
        
        this->GetCellData()->SetItem(cellChemistry -> GetChemicalNamesByIndex(i), new_parent_species_concentration);
        p_daughter_cell_data->SetItem(cellChemistry -> GetChemicalNamesByIndex(i), daughter_species_concentration);
    }
std::cout<<"here7"<<std::endl;

    // run through cell properties and create new objects for them 
    
    //daughter_property_collection = mCellPropertyCollection;


    // transport property
    if(mCellPropertyCollection.HasProperty<TransportCellProperty>())
    {
        // parent property
        boost::shared_ptr<TransportCellProperty> transport_cell_property = boost::static_pointer_cast<TransportCellProperty>(mCellPropertyCollection.GetPropertiesType<TransportCellProperty>().GetProperty());
        AbstractChemistry* transportChemistry = transport_cell_property ->GetTransportReactionSystem()->GetCellChemistry();

        // create new transport cell property
        daughter_property_collection.RemoveProperty(transport_cell_property);
        // use copy construtor
        boost::shared_ptr<TransportCellProperty> p_daughter_transport_property(new TransportCellProperty(*transport_cell_property));
        daughter_property_collection.AddProperty(p_daughter_transport_property);

        // split the properties betwene the two cells
        transport_cell_property->PreparePostDivisionParent(mSplitRatio);
        p_daughter_transport_property->PreparePostDivisionDaughter(*transport_cell_property, mSplitRatio);

    }

    // membrane property
    if(mCellPropertyCollection.HasProperty<MembraneCellProperty>())
    {
        boost::shared_ptr<MembraneCellProperty> membrane_cell_property = boost::static_pointer_cast<MembraneCellProperty>(mCellPropertyCollection.GetPropertiesType<MembraneCellProperty>().GetProperty());


        // create new membrane cell property
        daughter_property_collection.RemoveProperty(membrane_cell_property);
        boost::shared_ptr<MembraneCellProperty> p_daughter_membrane_property(new MembraneCellProperty(*membrane_cell_property));
        daughter_property_collection.AddProperty(p_daughter_membrane_property);

        // split the properties betwene the two cells
        membrane_cell_property->PreparePostDivisionParent(mSplitRatio);
        p_daughter_membrane_property->PreparePostDivisionDaughter(*membrane_cell_property, mSplitRatio);
    }

    // cellAnalytics property
    if(mCellPropertyCollection.HasProperty<CellAnalyticsProperty>())
    {
        // parent property
        boost::shared_ptr<CellAnalyticsProperty> cellAnalytics_cell_property = boost::static_pointer_cast<CellAnalyticsProperty>(mCellPropertyCollection.GetPropertiesType<CellAnalyticsProperty>().GetProperty());

        // create new transport cell property
        daughter_property_collection.RemoveProperty(cellAnalytics_cell_property);
        // use copy construtor
        boost::shared_ptr<CellAnalyticsProperty> p_daughter_cellAnalytics_property(new CellAnalyticsProperty(*cellAnalytics_cell_property));
        daughter_property_collection.AddProperty(p_daughter_cellAnalytics_property);

        // split the properties betwene the two cells
        cellAnalytics_cell_property->PreparePostDivisionParent(mSplitRatio);
        p_daughter_cellAnalytics_property->PreparePostDivisionDaughter(*cellAnalytics_cell_property, mSplitRatio);

    }
std::cout<<"here8"<<std::endl;
    // environment property
    if(mCellPropertyCollection.HasProperty<EnvironmentCellProperty>())
    {
        // parent property
        boost::shared_ptr<EnvironmentCellProperty> environment_cell_property = boost::static_pointer_cast<EnvironmentCellProperty>(mCellPropertyCollection.GetPropertiesType<EnvironmentCellProperty>().GetProperty());

        // create new transport cell property
        daughter_property_collection.RemoveProperty(environment_cell_property);
        // use copy construtor
        boost::shared_ptr<EnvironmentCellProperty> p_daughter_environments_property(new EnvironmentCellProperty(*environment_cell_property));
        daughter_property_collection.AddProperty(p_daughter_environments_property);

        // split the properties betwene the two cells
        environment_cell_property->PreparePostDivisionParent(mSplitRatio);
        p_daughter_environments_property->PreparePostDivisionDaughter(*environment_cell_property, mSplitRatio);

    }

    // state switching property
    
    if(mCellPropertyCollection.HasProperty<StateSwitchingCellProperty>())
    {
        // parent property
        boost::shared_ptr<StateSwitchingCellProperty> state_cell_property = boost::static_pointer_cast<StateSwitchingCellProperty>(mCellPropertyCollection.GetPropertiesType<StateSwitchingCellProperty>().GetProperty());

        // create new transport cell property
        daughter_property_collection.RemoveProperty(state_cell_property);
        // use copy construtor
        boost::shared_ptr<StateSwitchingCellProperty> p_daughter_state_property(new StateSwitchingCellProperty(*state_cell_property));
        daughter_property_collection.AddProperty(state_cell_property);

        // split the properties betwene the two cells
        state_cell_property->PreparePostDivisionParent(mSplitRatio);
        p_daughter_state_property->PreparePostDivisionDaughter(*state_cell_property, mSplitRatio);

    }

std::cout<<"here9"<<std::endl;
    // create new chemical cell
    // Create daughter cell with modified cell property collection
    CellPtr p_new_cell(new ComplexCell(GetMutationState(), mpCellCycleModel->CreateCellCycleModel(), mpSrnModel->CreateSrnModel(), false, daughter_property_collection));
    // Initialise properties of daughter cell
std::cout<<"here10"<<std::endl;
    p_new_cell->GetCellCycleModel()->InitialiseDaughterCell();
//    if(static_cast<SimpleChemicalThresholdCellCycleModel*>(mpCellCycleModel)->CellCycleType()=="Chemical")
//   {

   // p_new_cell->GetCellCycleModel()->SetUp(static_cast<SimpleChemicalThresholdCellCycleModel*>(mpCellCycleModel)->GetThresholdChemistry());
//    }

    p_new_cell->GetSrnModel()->InitialiseDaughterCell();



    // sort out the new setup of the daugheter properties

    CellPropertyCollection new_daughter_property_collection = p_new_cell-> rGetCellPropertyCollection();
    
    
    if(new_daughter_property_collection.HasProperty<EnvironmentCellProperty>())
    {
        boost::shared_ptr<EnvironmentCellProperty> new_environment_cell_property = boost::static_pointer_cast<EnvironmentCellProperty>(new_daughter_property_collection.GetPropertiesType<EnvironmentCellProperty>().GetProperty());
        new_environment_cell_property->SetCellPtr(p_new_cell);
    }
    
    if(new_daughter_property_collection.HasProperty<StateSwitchingCellProperty>())
    {
        boost::shared_ptr<StateSwitchingCellProperty> new_state_cell_property = boost::static_pointer_cast<StateSwitchingCellProperty>(new_daughter_property_collection.GetPropertiesType<StateSwitchingCellProperty>().GetProperty());
        new_state_cell_property->SetCellPtr(p_new_cell);
    }
    
    if(new_daughter_property_collection.HasProperty<TransportCellProperty>())
    {
        boost::shared_ptr<TransportCellProperty> new_transport_cell_property = boost::static_pointer_cast<TransportCellProperty>(new_daughter_property_collection.GetPropertiesType<TransportCellProperty>().GetProperty());
        new_transport_cell_property->SetCellPtr(p_new_cell);
    }
    if(new_daughter_property_collection.HasProperty<MembraneCellProperty>())
    {
        boost::shared_ptr<MembraneCellProperty> new_membrane_cell_property = boost::static_pointer_cast<MembraneCellProperty>(new_daughter_property_collection.GetPropertiesType<MembraneCellProperty>().GetProperty());
        new_membrane_cell_property->SetCellPtr(p_new_cell);
    }
    if(new_daughter_property_collection.HasProperty<CellAnalyticsProperty>())
    {
        boost::shared_ptr<CellAnalyticsProperty> new_analytics_cell_property = boost::static_pointer_cast<CellAnalyticsProperty>(new_daughter_property_collection.GetPropertiesType<CellAnalyticsProperty>().GetProperty());
        new_analytics_cell_property->SetCellPtr(p_new_cell);
    }




std::cout<<"here11"<<std::endl;
    // Set the daughter cell to inherit the apoptosis time of the parent cell
    p_new_cell->SetApoptosisTime(mApoptosisTime);

    boost::static_pointer_cast<ComplexCell>(p_new_cell)->SetChemicalNames(mChemicalNames);

    boost::static_pointer_cast<ComplexCell>(p_new_cell)->SetChemicalDivsionRules(mChemicalDivsionRules);

std::cout<<"here12"<<std::endl;
    std::vector<double> prime_cell_threshold_species_concentrations(static_cast<SimpleChemicalThresholdCellCycleModel*>(this->GetCellCycleModel())->GetNumberThresholdSpecies(),0.0);
    std::vector<double> daughter_cell_threshold_species_concentrations(static_cast<SimpleChemicalThresholdCellCycleModel*>(p_new_cell->GetCellCycleModel())->GetNumberThresholdSpecies(),0.0);

    AbstractChemistry* thresholdChemistry = static_cast<SimpleChemicalThresholdCellCycleModel*>(this->GetCellCycleModel())->GetThresholdChemistry();
    std::string chemicalName;

std::cout<<"here13"<<std::endl;
    for(unsigned i=0; i<numberOfChemicals; i++)
    {
        chemicalName = cellChemistry -> GetChemicalNamesByIndex(i);

        if(!thresholdChemistry->CheckChemical(new AbstractChemical(chemicalName)))
        {
            prime_cell_threshold_species_concentrations[thresholdChemistry->GetChemicalIndexByName(chemicalName)] = p_cell_data->GetItem(chemicalName);
            //std::cout<<"parent concentration: "<<chemicalName<<" : "<<p_cell_data->GetItem(chemicalName)<<std::endl;
        }
    }


    for(unsigned i=0; i<numberOfChemicals; i++)
    {
        chemicalName = cellChemistry -> GetChemicalNamesByIndex(i);
        //std::cout<<chemicalName<<": "<<p_daughter_cell_data->GetItem(chemicalName)<<std::endl;
        if(!thresholdChemistry->CheckChemical(new AbstractChemical(chemicalName)))
        {
            daughter_cell_threshold_species_concentrations[thresholdChemistry->GetChemicalIndexByName(chemicalName)] = p_daughter_cell_data->GetItem(chemicalName);
            //std::cout<<"daughter concentration: "<<chemicalName<<" : "<<p_daughter_cell_data->GetItem(chemicalName)<<std::endl;
        }
    }
std::cout<<"here14"<<std::endl;
    // update chemical cell cycles for each of the daughter cells
     
  //  static_cast<SimpleChemicalThresholdCellCycleModel*>(p_new_cell->GetCellCycleModel())->SetUp(static_cast<SimpleChemicalThresholdCellCycleModel*>(mpCellCycleModel)->GetThresholdChemistry());
    static_cast<SimpleChemicalThresholdCellCycleModel*>(this->GetCellCycleModel())->SetSpeciesConcentrations(prime_cell_threshold_species_concentrations);
    static_cast<SimpleChemicalThresholdCellCycleModel*>(p_new_cell->GetCellCycleModel())->SetSpeciesConcentrations(daughter_cell_threshold_species_concentrations);
std::cout<<"here15"<<std::endl;
    // update Ode from cell data for each of the daughter cells
    static_cast<ChemicalSrnModel*>(this->GetSrnModel())->UpdateOdeStatesFromCellData();
    static_cast<ChemicalSrnModel*>(p_new_cell->GetSrnModel())->UpdateOdeStatesFromCellData();
std::cout<<"here16"<<std::endl;
    ChemicalSrnModel *p_chemical_srn = static_cast<ChemicalSrnModel*>(p_new_cell->GetSrnModel());
    p_chemical_srn->GetReactionSystem()->SetCell(p_new_cell); 
    p_chemical_srn->GetReactionSystem()->DistributeCellPtr();
    std::cout<<"here17"<<std::endl;
    return p_new_cell;
}


bool ComplexCell::ReadyToDivide()
{
    std::cout<<"ComplexCell::ReadyToDivide()"<<std::endl;
    bool readyToDivide = Cell::ReadyToDivide();
    std::cout<<"ReadtToDivide? "<<readyToDivide<<std::endl;
    return readyToDivide;
}




void ComplexCell::DetermineSplitRatio()
{
    //mSplitRatio =0.5;
}


double ComplexCell::SplitParentCellData(double current_parent_value)
{
    return current_parent_value*GetSplitRatio();
}

bool ComplexCell::IsChemicalShared(std::string chemical_name)
{
    // determine whether the chemical concetration is shared between parnt and daughter after division
    bool isChemicalShared=false;

    for(unsigned i=0; i<mChemicalNames.size(); i++)
    {
        std::cout<<mChemicalNames[i]<<std::endl;
        // for each cellular chemical
        if(mChemicalNames[i]==chemical_name)
        {
            for(unsigned j=0; j<mShareKey.size(); j++)
            {
                // test if the chemical is to be shared
                if(mChemicalDivsionRules[i]==mShareKey[j])
                {
                    isChemicalShared=true;
                }
            }
        }
    }

    return isChemicalShared; 
}

double ComplexCell::GetSplitRatio()
{
    return mSplitRatio;
}

void ComplexCell::SetSplitRatio(double split_ratio)
{
    std::cout<<"ComplexCell::SetSplitRatio(double split_ratio) - start"<<std::endl;
    mSplitRatio = split_ratio;
    std::cout<<"ComplexCell::SetSplitRatio(double split_ratio) - end"<<std::endl;
}

std::vector<std::string> ComplexCell::GetShareKey()
{
    return mShareKey;
}

void ComplexCell::SetShareKey(std::vector<std::string> shareKey)
{
    mShareKey=shareKey;
}

std::vector<std::string> ComplexCell::GetChemicalNames()
{
    return mChemicalNames;
}

void ComplexCell::SetChemicalNames(std::vector<std::string> chemicalNames)
{
    mChemicalNames=chemicalNames;
}

std::vector<std::string> ComplexCell::GetChemicalDivsionRules()
{
    return mChemicalDivsionRules;
}

void ComplexCell::SetChemicalDivsionRules(std::vector<std::string> chemicalRules)
{
    mChemicalDivsionRules=chemicalRules;
}



#endif
*/