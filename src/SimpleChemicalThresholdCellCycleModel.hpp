#ifndef SIMPLECHEMICALTHRESHOLDCELLCYCLEMODEL_HPP
#define SIMPLECHEMICALTHRESHOLDCELLCYCLEMODEL_HPP

#include "AbstractSimplePhaseBasedCellCycleModel.hpp"
#include "AbstractChemistry.hpp"

// simple cell cycle model wherein the concentrations of species within the cell are compared to 
// maximum and minimum thresholds. If a single species is above the threshols then cell division
// is initilaed while if below cell death is triggered

class SimpleChemicalThresholdCellCycleModel : public AbstractSimplePhaseBasedCellCycleModel
{
protected:
    
    double mCurrentStarvationDuration;
    double mCurrentStarvationOnsetTime;
    double mCriticalStarvationDuration;

    // for tracking properties
    AbstractChemistry* mpThresholdChemistry;

    std::vector<double> mSpeciesConcentrations;

    unsigned mNumberThresholdSpecies;

    // for control properties
    std::vector<double> mMaxThresholdConcentrationValues;

    std::vector<bool> mIsMaximumThresholdSet;

    std::vector<double> mMinThresholdConcentrationValues;

    std::vector<bool> mIsMinimumThresholdSet;

public:

    SimpleChemicalThresholdCellCycleModel();

    SimpleChemicalThresholdCellCycleModel(const SimpleChemicalThresholdCellCycleModel& rModel);
    
    virtual ~SimpleChemicalThresholdCellCycleModel()
    {}

    virtual std::string CellCycleType()
    {
        return "Chemical";
    }

    void SetUp(AbstractChemistry*);

    void Initialise();

    void InitialiseDaughterCell();

    double GetCurrentStarvationDuration() const;

    double GetCurrentStarvationOnsetTime() const;

    void UpdateCellCyclePhase();

    bool ReadyToDivide();

    AbstractCellCycleModel* CreateCellCycleModel();

    void UpdateStarvationDuration();

    double CellDeathProbability(std::vector<double>);

    bool ConcentrationAboveMaxThreshold();

    bool ConcentrationBelowMinThreshold();

    void PrepareForDivision();

    void SetMaximumSpeciesThreshold(std::vector<double>);

    void SetMinimumSpeciesThreshold(std::vector<double>);

    void SetMaximumThresholdByIndex(double, unsigned);

    void SetMinimumThresholdByIndex(double, unsigned);

    void SetSpeciesConcentrations(std::vector<double>);

    void SetSpeciesConcentrationsByIndex(double, unsigned);

    void SetNumberThresholdSpecies(unsigned);

    void SetMaximumThresholdCheck(std::vector<bool>);

    void SetMinimumThresholdCheck(std::vector<bool>);

    void SetMaximumThresholdCheckByIndex(bool, unsigned);

    void SetMinimumThresholdCheckByIndex(bool, unsigned);

    void SetThresholdChemistry(AbstractChemistry*);

    AbstractChemistry* GetThresholdChemistry();

    std::vector<double> GetMaximumSpeciesThreshold();

    std::vector<double> GetMinimumSpeciesThreshold();

    double GetMaximumThresholdByIndex(unsigned);

    double GetMinimumThresholdByIndex(unsigned);

    std::vector<double> GetSpeciesConcentrations();

    double GetSpeciesConcentrationsByIndex(unsigned);

    unsigned GetNumberThresholdSpecies();

    std::vector<bool> GetMaximumThresholdCheck();

    std::vector<bool> GetMinimumThresholdCheck();

    bool GetMaximumThresholdCheckByIndex(unsigned);

    bool GetMinimumThresholdCheckByIndex(unsigned);

    double GetCriticalStarvationDuration() const;

    void SetCriticalStarvationDuration(double criticalStarvationDuration);

    void SetCurrentStarvationOnsetTime(double currentStarvationOnsetTime);
    
};
#endif