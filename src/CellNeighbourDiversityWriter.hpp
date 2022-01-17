#ifndef CELLNEIGHBOURDIVERSITYWRITER_HPP_
#define CELLNEIGHBOURDIVERSITYWRITER_HPP_

#include "AbstractCellWriter.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class CellNeighbourDiversityWriter : public AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>
{

public:

    CellNeighbourDiversityWriter();

    double GetCellDataForVtkOutput(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation);

    virtual void VisitCell(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation);

    double CalculateNeighbourDiversityIndex(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation);

};

#endif
