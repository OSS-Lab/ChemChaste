#ifndef CHEMICALCELL_HPP_
#define CHEMICALCELL_HPP_

#include <boost/utility.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/enable_shared_from_this.hpp>

#include "ChasteSerialization.hpp"
#include <boost/serialization/shared_ptr.hpp>

#include "AbstractCellMutationState.hpp"
#include "AbstractCellProliferativeType.hpp"
#include "CellData.hpp"
#include "CellVecData.hpp"
#include "Cell.hpp"

#include "AbstractCellCycleModel.hpp"
#include "AbstractSrnModel.hpp"
#include "CellPropertyCollection.hpp"

#include "AbstractChemistry.hpp"
#include "MembraneCellProperty.hpp"
#include "TransportCellProperty.hpp"

class AbstractCellCycleModel; // Circular definition (cells need to know about cycle models and vice-versa).
class AbstractSrnModel; // Circular definition (cells need to know about subcellular reaction network models and vice-versa).
class Cell;

// version of cell.hpp in which the cell division properties may be overridden

class ChemicalCell : public Cell
{
protected:
    using Cell::Divide;

    double mSplitRatio =0.5; // proportion of parent cell volume retained, rest goes to daughter 

public:

    ChemicalCell(boost::shared_ptr<AbstractCellProperty> pMutationState,
         AbstractCellCycleModel* pCellCycleModel,
         AbstractSrnModel* pSrnModel=nullptr,
         bool archiving=false,
         CellPropertyCollection cellPropertyCollection=CellPropertyCollection());

    virtual ~ChemicalCell()
    {
    };

    virtual CellPtr Divide();

    virtual void DetermineSplitRation();

    virtual double SplitParentCellData(double);

    double GetSplitRation();

    void SetSplitRation(double);

};

ChemicalCell::ChemicalCell(
        boost::shared_ptr<AbstractCellProperty> pMutationState,
        AbstractCellCycleModel* pCellCycleModel,
        AbstractSrnModel* pSrnModel,
        bool archiving,
        CellPropertyCollection cellPropertyCollection)

    : Cell(pMutationState,pCellCycleModel,pSrnModel,archiving,cellPropertyCollection)
         {

         };

CellPtr ChemicalCell::Divide()
{

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

    // Copy all cell data (note we create a new object not just copying the pointer)
    assert(daughter_property_collection.HasPropertyType<CellData>());
    // Get the existing copy of the cell data and remove it from the daughter cell
    boost::shared_ptr<CellData> p_cell_data = GetCellData();
    daughter_property_collection.RemoveProperty(p_cell_data);
    // Create a new cell data object using the copy constructor and add this to the daughter cell
    MAKE_PTR_ARGS(CellData, p_daughter_cell_data, (*p_cell_data));
    daughter_property_collection.AddProperty(p_daughter_cell_data);
    // Copy all cell Vec data (note we create a new object not just copying the pointer)
    if (daughter_property_collection.HasPropertyType<CellVecData>())
    {
        // Get the existing copy of the cell data and remove it from the daughter cell
        boost::shared_ptr<CellVecData> p_cell_vec_data = GetCellVecData();
        daughter_property_collection.RemoveProperty(p_cell_vec_data);
        // Create a new cell data object using the copy constructor and add this to the daughter cell
        MAKE_PTR_ARGS(CellVecData, p_daughter_cell_vec_data, (*p_cell_vec_data));
        daughter_property_collection.AddProperty(p_daughter_cell_vec_data);
    }

    // record the cell chemistry for splitting cell data
    AbstractChemistry* cellChemistry = new AbstractChemistry();

    if(static_cast<ChemicalSrnModel*>(mpSrnModel)->SRNType()=="Chemical")
    {
        ChemicalSrnModel* p_srn_model = static_cast<ChemicalSrnModel*>(mpSrnModel);
        cellChemistry -> AddChemistry(p_srn_model->GetCellChemistry()); // from SRN
    }
    //std::cout<<cellChemistry -> GetNumberChemicals()<<std::endl;
    if(static_cast<SimpleChemicalThresholdCellCycleModel*>(mpCellCycleModel)->CellCycleType()=="Chemical")
    {
        SimpleChemicalThresholdCellCycleModel* p_cc_model = static_cast<SimpleChemicalThresholdCellCycleModel*>(mpCellCycleModel);
        cellChemistry -> AddChemistry(p_cc_model->GetThresholdChemistry()); // from cell cycle model
    }
    
    //std::cout<<cellChemistry -> GetNumberChemicals()<<std::endl;
    // transport property
    if(mCellPropertyCollection.HasProperty<TransportCellProperty>())
    {
        boost::shared_ptr<TransportCellProperty> transport_cell_property = boost::static_pointer_cast<TransportCellProperty>(mCellPropertyCollection.GetPropertiesType<TransportCellProperty>().GetProperty());
        AbstractChemistry* transportChemistry = transport_cell_property ->GetTransportReactionSystem()->GetCellChemistry();
        // add the cell chemistry due to transport
        cellChemistry -> AddChemistry(transportChemistry);
    }
    //std::cout<<cellChemistry -> GetNumberChemicals()<<std::endl;
    // membrane property
    if(mCellPropertyCollection.HasProperty<MembraneCellProperty>())
    {
        boost::shared_ptr<MembraneCellProperty> membrane_cell_property = boost::static_pointer_cast<MembraneCellProperty>(mCellPropertyCollection.GetPropertiesType<MembraneCellProperty>().GetProperty());
        AbstractChemistry* membraneChemistry = membrane_cell_property ->GetMembraneReactionSystem()->GetCellChemistry();
        // add the cell chemistry due to membrane
        cellChemistry -> AddChemistry(membraneChemistry);
    }


    //std::cout<<"Cell divide"<<std::endl;
    //std::cout<<"Parent cell data"<<std::endl;

    unsigned numberOfChemicals = cellChemistry -> GetNumberChemicals();
    std::string chemicalName;

    for(unsigned i=0; i<numberOfChemicals; i++)
    {
        chemicalName = cellChemistry -> GetChemicalNamesByIndex(i);
        //std::cout<<chemicalName<<": "<<this->GetCellData()->GetItem(chemicalName)<<std::endl;
    }


    // based on the parent cell data detemrine the split ratio
    DetermineSplitRation();

    // share the two cellDatas between the two cells, use cellChemistry
    // halve chemical cell data
    double parent_species_concentration=0.0;
    double new_parent_species_concentration=0.0;
    double daughter_species_concentration=0.0;
    for(unsigned i=0; i<numberOfChemicals; i++)
    {
        parent_species_concentration = this->GetCellData()->GetItem(cellChemistry -> GetChemicalNamesByIndex(i));
        
        new_parent_species_concentration = SplitParentCellData(parent_species_concentration);

        daughter_species_concentration = parent_species_concentration - new_parent_species_concentration;

        if(daughter_species_concentration<0.0)
        {
            // must be positive concentration
            daughter_species_concentration =0.0;
        }
        
        this->GetCellData()->SetItem(cellChemistry -> GetChemicalNamesByIndex(i), new_parent_species_concentration);
        p_daughter_cell_data->SetItem(cellChemistry -> GetChemicalNamesByIndex(i), daughter_species_concentration);
    }


    // run through cell properties and create new objects for them 
    // transport property
    daughter_property_collection = mCellPropertyCollection;
    if(mCellPropertyCollection.HasProperty<TransportCellProperty>())
    {
        boost::shared_ptr<TransportCellProperty> transport_cell_property = boost::static_pointer_cast<TransportCellProperty>(mCellPropertyCollection.GetPropertiesType<TransportCellProperty>().GetProperty());
        AbstractChemistry* transportChemistry = transport_cell_property ->GetTransportReactionSystem()->GetCellChemistry();

        // create new transport cell property
        daughter_property_collection.RemoveProperty(transport_cell_property);
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
        AbstractChemistry* membraneChemistry = membrane_cell_property ->GetMembraneReactionSystem()->GetCellChemistry();

        // create new membrane cell property
        daughter_property_collection.RemoveProperty(membrane_cell_property);
        boost::shared_ptr<MembraneCellProperty> p_daughter_membrane_property(new MembraneCellProperty(*membrane_cell_property));
        daughter_property_collection.AddProperty(p_daughter_membrane_property);

        // split the properties betwene the two cells
        membrane_cell_property->PreparePostDivisionParent(mSplitRatio);
        p_daughter_membrane_property->PreparePostDivisionDaughter(*membrane_cell_property, mSplitRatio);
    }


    // create new cell
    // Create daughter cell with modified cell property collection
    CellPtr p_new_cell(new Cell(GetMutationState(), mpCellCycleModel->CreateCellCycleModel(), mpSrnModel->CreateSrnModel(), false, daughter_property_collection));
    // Initialise properties of daughter cell

    p_new_cell->GetCellCycleModel()->InitialiseDaughterCell();

    p_new_cell->GetSrnModel()->InitialiseDaughterCell();

    // Set the daughter cell to inherit the apoptosis time of the parent cell
    p_new_cell->SetApoptosisTime(mApoptosisTime);




    std::vector<double> prime_cell_threshold_species_concentrations(static_cast<SimpleChemicalThresholdCellCycleModel*>(this->GetCellCycleModel())->GetNumberThresholdSpecies(),0.0);
    std::vector<double> daughter_cell_threshold_species_concentrations(static_cast<SimpleChemicalThresholdCellCycleModel*>(p_new_cell->GetCellCycleModel())->GetNumberThresholdSpecies(),0.0);

    AbstractChemistry* thresholdChemistry = static_cast<SimpleChemicalThresholdCellCycleModel*>(this->GetCellCycleModel())->GetThresholdChemistry();

    for(unsigned i=0; i<numberOfChemicals; i++)
    {
        chemicalName = cellChemistry -> GetChemicalNamesByIndex(i);
        //std::cout<<chemicalName<<": "<<p_cell_data->GetItem(chemicalName)<<std::endl;

        if(thresholdChemistry->CheckChemical(new AbstractChemical(chemicalName)))
        {
            prime_cell_threshold_species_concentrations[thresholdChemistry->GetChemicalIndexByName(chemicalName)] = p_cell_data->GetItem(chemicalName);
        }
    }


    for(unsigned i=0; i<numberOfChemicals; i++)
    {
        chemicalName = cellChemistry -> GetChemicalNamesByIndex(i);
        //std::cout<<chemicalName<<": "<<p_daughter_cell_data->GetItem(chemicalName)<<std::endl;
        if(thresholdChemistry->CheckChemical(new AbstractChemical(chemicalName)))
        {
            daughter_cell_threshold_species_concentrations[thresholdChemistry->GetChemicalIndexByName(chemicalName)] = p_daughter_cell_data->GetItem(chemicalName);
        }
    }

    // update chemical cell cycles for each of the daughter cells
    static_cast<SimpleChemicalThresholdCellCycleModel*>(this->GetCellCycleModel())->SetSpeciesConcentrations(prime_cell_threshold_species_concentrations);
    static_cast<SimpleChemicalThresholdCellCycleModel*>(p_new_cell->GetCellCycleModel())->SetSpeciesConcentrations(daughter_cell_threshold_species_concentrations);

    // update Ode from cell data for each of the daughter cells
    static_cast<ChemicalSrnModel*>(this->GetSrnModel())->UpdateOdeStatesFromCellData();
    static_cast<ChemicalSrnModel*>(p_new_cell->GetSrnModel())->UpdateOdeStatesFromCellData();



    return p_new_cell;
}

void ChemicalCell::DetermineSplitRation()
{
    mSplitRatio =0.5;
}


double ChemicalCell::SplitParentCellData(double current_parent_value)
{
    return current_parent_value*GetSplitRation();
}

double ChemicalCell::GetSplitRation()
{
    return mSplitRatio;
}

void ChemicalCell::SetSplitRation(double split_ratio)
{
    mSplitRatio = split_ratio;
}

#endif