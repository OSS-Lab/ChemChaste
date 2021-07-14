
#include "CellProliferativeTypesWriter_ChemChaste.hpp"

#include "AbstractCellPopulation.hpp"
#include "ApoptoticCellProperty.hpp"
#include "CellLabel.hpp"
#include "WildTypeCellMutationState.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
CellProliferativeTypesWriter<ELEMENT_DIM, SPACE_DIM>::CellProliferativeTypesWriter()
    : AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>("results.vizcelltypes")
{
    this->mVtkCellDataName = "Cell types";
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double CellProliferativeTypesWriter<ELEMENT_DIM, SPACE_DIM>::GetCellDataForVtkOutput(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    double colour = pCell->GetCellProliferativeType()->GetColour();

    // Set colour dependent on cell mutation state
    if (!pCell->GetMutationState()->IsType<WildTypeCellMutationState>())
    {
        colour = pCell->GetMutationState()->GetColour();
    }
    if (pCell->rGetCellPropertyCollection().HasProperty<CellLabel>())
    {
        CellPropertyCollection collection = pCell->rGetCellPropertyCollection().GetProperties<CellLabel>();
        boost::shared_ptr<CellLabel> p_label = boost::static_pointer_cast<CellLabel>(collection.GetProperty());
        colour = p_label->GetColour();
    }
    if (pCell->rGetCellPropertyCollection().HasProperty<ApoptoticCellProperty>() || pCell->HasApoptosisBegun())
    {
        ///\todo: replace this hard-coded 6 with the ApoptoticCellProperty member mColour?
        colour = 6.0;
    }

    return colour;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellProliferativeTypesWriter<ELEMENT_DIM, SPACE_DIM>::VisitCell(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
std::cout<<"CellProliferativeTypesWriter::VisitCell"<<std::endl;
    unsigned colour = pCell->GetCellProliferativeType()->GetColour();
std::cout<<"CPTW 0"<<std::endl;
    // Set colour dependent on cell mutation state
    if (!pCell->GetMutationState()->IsType<WildTypeCellMutationState>())
    {
std::cout<<"CPTW 0 in"<<std::endl;
        colour = pCell->GetMutationState()->GetColour();
    }
std::cout<<"CPTW 1"<<std::endl;
    if (pCell->rGetCellPropertyCollection().HasProperty<CellLabel>())
    {
std::cout<<"CPTW 1 in"<<std::endl;
        CellPropertyCollection collection = pCell->rGetCellPropertyCollection().GetProperties<CellLabel>();
        boost::shared_ptr<CellLabel> p_label = boost::static_pointer_cast<CellLabel>(collection.GetProperty());
        colour = p_label->GetColour();
    }
std::cout<<"CPTW 2"<<std::endl;
    if (pCell->rGetCellPropertyCollection().HasProperty<ApoptoticCellProperty>() || pCell->HasApoptosisBegun())
    {
std::cout<<"CPTW 2 in"<<std::endl;
        ///\todo: replace this hard-coded 6 with the ApoptoticCellProperty member mColour? (#2512)
        colour = 6;
    }
std::cout<<"CPTW 3"<<std::endl;
    *this->mpOutStream << colour << " ";
std::cout<<"CPTW 4"<<std::endl;
}

// Explicit instantiation
template class CellProliferativeTypesWriter<1,1>;
template class CellProliferativeTypesWriter<1,2>;
template class CellProliferativeTypesWriter<2,2>;
template class CellProliferativeTypesWriter<1,3>;
template class CellProliferativeTypesWriter<2,3>;
template class CellProliferativeTypesWriter<3,3>;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(CellProliferativeTypesWriter)
