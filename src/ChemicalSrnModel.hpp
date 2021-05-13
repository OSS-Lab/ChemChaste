#ifndef CHEMICALSRNMODEL_HPP_
#define CHEMICALSRNMODEL_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <string>

#include "AbstractChemicalOdeSystem.hpp"

#include "AbstractCellCycleModelOdeSolver.hpp"
#include "CellCycleModelOdeSolver.hpp"
#include "AbstractOdeSrnModel.hpp"

// class to define the subcellular reaction network (SRN) belonging to a cell. Solves an Ode system to model the 
// chemical reactions occuring within the cell. The reaction sutilise the cell chemical concentrations which
// are updated from the cellData class, updated in the ChemicalTrackingModifier class.

class ChemicalSrnModel : public AbstractOdeSrnModel
{
protected:

    AbstractReactionSystem* mpReactionSystem;

    AbstractChemistry* mpCellChemistry;

public:

    ChemicalSrnModel(AbstractReactionSystem* pReactionSystem = new AbstractReactionSystem(),boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver = boost::shared_ptr<AbstractCellCycleModelOdeSolver>());

    virtual ~ChemicalSrnModel()
    {}

    virtual std::string SRNType()
    {
        return "Chemical";
    }

    ChemicalSrnModel(const ChemicalSrnModel& rModel);

    AbstractSrnModel* CreateSrnModel();

    void Initialise(); 

    void SimulateToCurrentTime();

    void SetReactionSystem(AbstractReactionSystem*);

    AbstractReactionSystem* GetReactionSystem();

    void SetCellChemistry(AbstractChemistry*);

    AbstractChemistry* GetCellChemistry();

    void UpdateOdeStatesFromCellData();

    void UpdateOdeParameters();

    double GetStateValueByIndex(unsigned);
    
    double GetStateValueByName(std::string);

    void OutputSrnModelParameters(out_stream& rParamsFile);
    
};

#endif