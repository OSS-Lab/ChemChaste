#ifndef CELLNEIGHBOURDIVERSITYWRITER_HPP_
#define CELLNEIGHBOURDIVERSITYWRITER_HPP_

#include "AbstractCellWriter.hpp"

/**
 * A class written using the visitor pattern for writing cell ages to file.
 *
 * The output file is called cellages.dat by default. If VTK is switched on,
 * then the writer also specifies the VTK output for each cell, which is stored
 * in the VTK cell data "Ages" by default.
 */
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