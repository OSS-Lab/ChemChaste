#ifndef CELLNEIGHBOURGETISORDWRITER_HPP_
#define CELLNEIGHBOURGETISORDWRITER_HPP_

#include "AbstractCellWriter.hpp"


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class CellNeighbourGetisOrdWriter : public AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>
{

public:

    CellNeighbourGetisOrdWriter();

    double GetCellDataForVtkOutput(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation);

    virtual void VisitCell(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation);

    double CalculateNeighbourGetisOrdIndex(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation);

};

#endif