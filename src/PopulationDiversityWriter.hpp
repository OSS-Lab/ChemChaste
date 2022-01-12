#ifndef POPULATIONDIVERSITYWRITER_HPP_
#define POPULATIONDIVERSITYWRITER_HPP_

#include "AbstractCellPopulationCountWriter.hpp"
#include <boost/serialization/base_object.hpp>
#include <string>
#include <vector>

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class PopulationDiversityWriter : public AbstractCellPopulationCountWriter<ELEMENT_DIM, SPACE_DIM>
{
private:

    std::vector<std::string> mCellTypesString;

    std::vector<unsigned> mCellTypesUnsigned;

    std::vector<unsigned> mNumberOfCellsOfType;


public:

    PopulationDiversityWriter();

    double CalculateDiversityIndex(unsigned);

    virtual void WriteHeader(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation);

    virtual void Visit(NodeBasedCellPopulation<SPACE_DIM>* pCellPopulation);

    virtual void Visit(MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation);

    virtual void Visit(CaBasedCellPopulation<SPACE_DIM>* pCellPopulation);

    virtual void Visit(PottsBasedCellPopulation<SPACE_DIM>* pCellPopulation);

    virtual void Visit(VertexBasedCellPopulation<SPACE_DIM>* pCellPopulation);

    void SetNumberOfCellsOfType(std::vector<unsigned>);

    void SetCellTypesString(std::vector<std::string>);

    std::vector<unsigned> GetNumberOfCellsOfType();

};

#endif