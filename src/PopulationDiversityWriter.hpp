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

    std::vector<double> mProbabilityOfType;

    unsigned mNumberOfCells=0;



public:

    PopulationDiversityWriter();

    double CalculateDiversityIndex(unsigned);

    double CalculateShannonIndex();

    double CalculateGiniSimpsonIndex();

    double CalculateLeeOyburnIndex(NodeBasedCellPopulation<SPACE_DIM>* pCellPopulation);

    double CalculateMoranIndex(NodeBasedCellPopulation<SPACE_DIM>* pCellPopulation);

    double CalculateGearyIndex(NodeBasedCellPopulation<SPACE_DIM>* pCellPopulation);

    virtual void WriteHeader(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation);

    virtual void Visit(NodeBasedCellPopulation<SPACE_DIM>* pCellPopulation);

    virtual void Visit(MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation);

    virtual void Visit(CaBasedCellPopulation<SPACE_DIM>* pCellPopulation);

    virtual void Visit(PottsBasedCellPopulation<SPACE_DIM>* pCellPopulation);

    virtual void Visit(VertexBasedCellPopulation<SPACE_DIM>* pCellPopulation);

    void SetNumberOfCellsOfType(std::vector<unsigned>);

    void SetProbabilityOfType(std::vector<double>);

    void SetCellTypesString(std::vector<std::string>);

    std::vector<unsigned> GetNumberOfCellsOfType();

    std::vector<double> GetProbabilityOfType();

};

#endif