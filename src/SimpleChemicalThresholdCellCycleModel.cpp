#include "SimpleChemicalThresholdCellCycleModel.hpp"
#include "RandomNumberGenerator.hpp"
#include "ApoptoticCellProperty.hpp"
#include "CellData.hpp"

SimpleChemicalThresholdCellCycleModel::SimpleChemicalThresholdCellCycleModel()
{ 
}

SimpleChemicalThresholdCellCycleModel::SimpleChemicalThresholdCellCycleModel(const SimpleChemicalThresholdCellCycleModel& rModel)
:   AbstractSimplePhaseBasedCellCycleModel(rModel),
    mCurrentStarvationDuration(rModel.mCurrentStarvationDuration),
    mCurrentStarvationOnsetTime(rModel.mCurrentStarvationOnsetTime),
    mCriticalStarvationDuration(rModel.mCriticalStarvationDuration),
    mpThresholdChemistry(rModel.mpThresholdChemistry),
    mNumberThresholdSpecies(rModel.mNumberThresholdSpecies),
    mMaxThresholdConcentrationValues(rModel.mMaxThresholdConcentrationValues),
    mIsMaximumThresholdSet(rModel.mIsMaximumThresholdSet),
    mMinThresholdConcentrationValues(rModel.mMinThresholdConcentrationValues),
    mIsMinimumThresholdSet(rModel.mIsMinimumThresholdSet)
{
    std::cout<<"Call copy constructor ---       SimpleChemicalThresholdCellCycleModel"<<std::endl;
    std::vector<double> mSpeciesConcentrations(mNumberThresholdSpecies,0.0);
    SetSpeciesConcentrations(mSpeciesConcentrations);
    mCurrentStarvationOnsetTime = SimulationTime::Instance()->GetTime();
    mIsSetUp = true;
}

void SimpleChemicalThresholdCellCycleModel::SetUp(AbstractChemistry* thresholdChemistry)
{
    std::cout<<"---------------------------------------------------setup ccm-------------"<<std::endl;
    SetThresholdChemistry(thresholdChemistry);
    SetNumberThresholdSpecies(thresholdChemistry->GetNumberChemicals());

    std::vector<double> maxThresholdConcentrationValues(mNumberThresholdSpecies,0.0);
    std::vector<bool> isMaximumThresholdSet(mNumberThresholdSpecies,false);
    std::vector<double> minThresholdConcentrationValues(mNumberThresholdSpecies,0.0);
    std::vector<bool> isMinimumThresholdSet(mNumberThresholdSpecies,false);

    std::vector<double> mSpeciesConcentrations(mNumberThresholdSpecies,0.0);

    SetMaximumSpeciesThreshold(maxThresholdConcentrationValues);
    SetMaximumThresholdCheck(isMaximumThresholdSet);
    SetMinimumSpeciesThreshold(minThresholdConcentrationValues);
    SetMinimumThresholdCheck(isMinimumThresholdSet);

    SetSpeciesConcentrations(mSpeciesConcentrations);
    mIsSetUp = true;
    mCurrentStarvationOnsetTime = SimulationTime::Instance()->GetTime();
}

void SimpleChemicalThresholdCellCycleModel::Initialise()
{
    AbstractSimplePhaseBasedCellCycleModel::Initialise();
}

void SimpleChemicalThresholdCellCycleModel::InitialiseDaughterCell()
{
    AbstractSimplePhaseBasedCellCycleModel::InitialiseDaughterCell();
}

double SimpleChemicalThresholdCellCycleModel::GetCurrentStarvationDuration() const
{
    return mCurrentStarvationDuration;
}

double SimpleChemicalThresholdCellCycleModel::GetCurrentStarvationOnsetTime() const
{
    return mCurrentStarvationOnsetTime;
}

void SimpleChemicalThresholdCellCycleModel::UpdateCellCyclePhase()
{
  //  std::cout<<"SimpleChemicalThresholdCellCycleModel::UpdateCellCyclePhase() - start"<<std::endl;
    // mG1Duration is set when the cell-cycle model is given a cell
    //bool cell_is_apoptotic = mpCell->HasCellProperty<ApoptoticCellProperty>();


    bool cell_is_apoptotic = mpCell->rGetCellPropertyCollection().HasProperty<ApoptoticCellProperty>();
    if (!cell_is_apoptotic)
    {

        UpdateStarvationDuration();

        AbstractSimplePhaseBasedCellCycleModel::UpdateCellCyclePhase();

        if (mCurrentCellCyclePhase == G_ONE_PHASE)
        {
            // Update G1 duration based on oxygen concentration
            double dt = SimulationTime::Instance()->GetTimeStep();
            mG1Duration += dt; // modified from simpleOxygen... as no mQuiescentConcentration   
        }
    }
  //  std::cout<<"SimpleChemicalThresholdCellCycleModel::UpdateCellCyclePhase() - end"<<std::endl;
}

bool SimpleChemicalThresholdCellCycleModel::ReadyToDivide()
{
  //  std::cout<<"SimpleChemicalThresholdCellCycleModel::ReadyToDivide() - start"<<std::endl;
    assert(mpCell != nullptr);
    
    if (!mReadyToDivide)
    {
 //       std::cout<<"UpdateCellCyclePhase"<<std::endl;
        UpdateCellCyclePhase();
  //      std::cout<<"Check concenrations above thresh"<<std::endl;
        if ((mCurrentCellCyclePhase != G_ZERO_PHASE) &&
            (ConcentrationAboveMaxThreshold()))
        {
            mReadyToDivide = true;
            //PrepareForDivision();
            
        }
    }
//    std::cout<<"SimpleChemicalThresholdCellCycleModel::ReadyToDivide() - end "<<mReadyToDivide<<std::endl;
    return mReadyToDivide;
}

AbstractCellCycleModel* SimpleChemicalThresholdCellCycleModel::CreateCellCycleModel()
{
    return new SimpleChemicalThresholdCellCycleModel(*this);
}

void SimpleChemicalThresholdCellCycleModel::UpdateStarvationDuration()
{  
 //   std::cout<<"SimpleChemicalThresholdCellCycleModel::UpdateStarvationDuration() - start"<<std::endl;
    //assert(!(mpCell->HasCellProperty<ApoptoticCellProperty>()));

    assert(!(mpCell->rGetCellPropertyCollection().HasProperty<ApoptoticCellProperty>()));
    assert(!mpCell->HasApoptosisBegun());
//std::cout<<"get chem names"<<std::endl;
    std::vector<std::string> chemicalNames = mpThresholdChemistry->GetChemicalNames();
//std::cout<<"get concentrations"<<std::endl;
    // get concentration vector
   // if(!mIsSetUp)
    //{
     //   boost::shared_ptr<AbstractCellProperty> p_apoptotic_property =
      //      mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<ApoptoticCellProperty>();
      //      mpCell->AddCellProperty(p_apoptotic_property);
        
    //}
    //else{

        for(unsigned species=0; species<mNumberThresholdSpecies; species++)
        {

          //  std::cout<<"can access cell data?"<<std::endl;
            boost::shared_ptr<CellData> cell_data = boost::static_pointer_cast<CellData>(mpCell->rGetCellPropertyCollection().GetPropertiesType<CellData>().GetProperty());
        //    std::cout<<"number of cell data items: "<<cell_data->GetNumItems()<<std::endl;
      //      std::cout<<"look for species:"<<std::endl;
    //        std::cout<<cell_data->GetItem(mpThresholdChemistry->GetChemicalNamesByIndex(species))<<std::endl;



  //          std::cout<<"species: "<<species<<" : "<<mpThresholdChemistry->GetChemicalNamesByIndex(species)<<std::endl;
  //          std::cout<<"Concentration: "<<mpCell->GetCellData()->GetItem(mpThresholdChemistry->GetChemicalNamesByIndex(species))<<std::endl;
            SetSpeciesConcentrationsByIndex(mpCell->GetCellData()->GetItem(mpThresholdChemistry->GetChemicalNamesByIndex(species)), species);
        }

        // check each species for starvation
        bool isStarving=false;
        for(unsigned species=0; species<mNumberThresholdSpecies; species++)
        {
   // std::cout<<"here0"<<std::endl;
            if(GetMinimumThresholdCheckByIndex(species) && GetSpeciesConcentrationsByIndex(species)<GetMinimumThresholdByIndex(species))
            {
     //           std::cout<<"here1"<<std::endl;
                mCurrentStarvationDuration = (SimulationTime::Instance()->GetTime() - mCurrentStarvationOnsetTime);
   // std::cout<<"here2"<<std::endl;
                if (mCurrentStarvationDuration > mCriticalStarvationDuration && RandomNumberGenerator::Instance()->ranf() < CellDeathProbability({GetSpeciesConcentrationsByIndex(species),GetMinimumThresholdByIndex(species)}))
                {
   //             std::cout<<"here3"<<std::endl;
                    boost::shared_ptr<AbstractCellProperty> p_apoptotic_property =
                    mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<ApoptoticCellProperty>();
                    mpCell->AddCellProperty(p_apoptotic_property);
                    isStarving=true;

                    break;
                }
            }
        }
 //       std::cout<<"if is not starving"<<std::endl;
        if(!isStarving)
        {
            // Reset the cell's Starvation duration and update the time at which the onset of Starvation occurs
            mCurrentStarvationDuration = 0.0;
            mCurrentStarvationOnsetTime = SimulationTime::Instance()->GetTime();
        }
    
 //   std::cout<<"SimpleChemicalThresholdCellCycleModel::UpdateStarvationDuration() - end"<<std::endl;
}

double SimpleChemicalThresholdCellCycleModel::CellDeathProbability(std::vector<double> parameters)
{

    //double species_concentration = parameters[0];
    //double starvation_concentration = parameters[1];
    double prob_of_death = 1.0;// - 0.5*(species_concentration/starvation_concentration);

    return prob_of_death;
}

bool SimpleChemicalThresholdCellCycleModel::ConcentrationAboveMaxThreshold()
{
    // search through all the threhsold species for whether nay are above the maximum threshold
//    std::cout<<"SimpleChemicalThresholdCellCycleModel::ConcentrationAboveMaxThreshold() - START"<<std::endl;
//    std::cout<<"number threshold species: "<<mNumberThresholdSpecies<<std::endl;
//    if(!mIsSetUp)
//    { return false;
//    }
    
    for(unsigned species=0; species<mNumberThresholdSpecies; species++)
    {
//        std::cout<<"Species: "<<species<<std::endl;
//        std::cout<< GetMaximumThresholdCheckByIndex(species) <<std::endl;
//        std::cout<<GetSpeciesConcentrationsByIndex(species) <<std::endl;
//        std::cout<< GetMaximumThresholdByIndex(species)<<std::endl;
        if(GetMaximumThresholdCheckByIndex(species) && GetSpeciesConcentrationsByIndex(species)>GetMaximumThresholdByIndex(species))
        {
 //           std::cout<<"SimpleChemicalThresholdCellCycleModel::ConcentrationAboveMaxThreshold() - end true"<<std::endl;
            return true;
        }
    }
//    std::cout<<"SimpleChemicalThresholdCellCycleModel::ConcentrationAboveMaxThreshold() - end false"<<std::endl;
    return false;
}

bool SimpleChemicalThresholdCellCycleModel::ConcentrationBelowMinThreshold()
{
    // search through all the threhsold species for whether nay are above the maximum threshold
    for(unsigned species=0; species<mNumberThresholdSpecies; species++)
    {
        if(GetMinimumThresholdCheckByIndex(species) && GetSpeciesConcentrationsByIndex(species)<GetMinimumThresholdByIndex(species))
        {
            return true;
        }
    }
    return false;
}

void SimpleChemicalThresholdCellCycleModel::PrepareForDivision()
{
    
}

void SimpleChemicalThresholdCellCycleModel::SetMaximumSpeciesThreshold(std::vector<double> maxThreshold)
{
    mMaxThresholdConcentrationValues = maxThreshold;
}

void SimpleChemicalThresholdCellCycleModel::SetMinimumSpeciesThreshold(std::vector<double> minThreshold)
{
    mMinThresholdConcentrationValues = minThreshold;
}

void SimpleChemicalThresholdCellCycleModel::SetMaximumThresholdByIndex(double threshold, unsigned index)
{
    if(index<mNumberThresholdSpecies)
    {
        mMaxThresholdConcentrationValues[index] = threshold;
    }
    else
    {
        std::cout<<"Error: SimpleChemicalThresholdCellCycleModel::SetMaximumThresholdByIndex, index out of bounds"<<std::endl;
    }
}

void SimpleChemicalThresholdCellCycleModel::SetMinimumThresholdByIndex(double threshold, unsigned index)
{
    if(index<mNumberThresholdSpecies)
    {
        mMinThresholdConcentrationValues[index] = threshold;
    }
    else
    {
        std::cout<<"Error: SimpleChemicalThresholdCellCycleModel::SetMinimumThresholdByIndex, index out of bounds"<<std::endl;
    }
}

void SimpleChemicalThresholdCellCycleModel::SetSpeciesConcentrations(std::vector<double> concentrations)
{
    mSpeciesConcentrations = concentrations;
}

void SimpleChemicalThresholdCellCycleModel::SetSpeciesConcentrationsByIndex(double concentration, unsigned index)
{
 //   std::cout<<"SimpleChemicalThresholdCellCycleModel::SetSpeciesConcentrationsByIndex"<<std::endl;
//    std::cout<<concentration<<std::endl;
//    std::cout<<index<<std::endl;
//    std::cout<<mNumberThresholdSpecies<<std::endl;
   
    if(index<mNumberThresholdSpecies)
    {
        mSpeciesConcentrations[index] = concentration;
    }
    else
    {
        std::cout<<"Error: SimpleChemicalThresholdCellCycleModel::SetSpeciesConcentrationsByIndex, index out of bounds"<<std::endl;
    }
}

void SimpleChemicalThresholdCellCycleModel::SetNumberThresholdSpecies(unsigned speciesNumber)
{
    mNumberThresholdSpecies = speciesNumber;
}

void SimpleChemicalThresholdCellCycleModel::SetMaximumThresholdCheck(std::vector<bool> thresholdCheck)
{
    mIsMaximumThresholdSet = thresholdCheck;
}

void SimpleChemicalThresholdCellCycleModel::SetMinimumThresholdCheck(std::vector<bool> thresholdCheck)
{
    mIsMinimumThresholdSet = thresholdCheck;
}

void SimpleChemicalThresholdCellCycleModel::SetMaximumThresholdCheckByIndex(bool thresholdCheck, unsigned index)
{
    if(index<mNumberThresholdSpecies)
    {
        mIsMaximumThresholdSet[index] = thresholdCheck;
    }
    else
    {
        std::cout<<"Error: SimpleChemicalThresholdCellCycleModel::SetMaximumThresholdCheckByIndex, index out of bounds"<<std::endl;
    }
}

void SimpleChemicalThresholdCellCycleModel::SetMinimumThresholdCheckByIndex(bool thresholdCheck, unsigned index)
{
    if(index<mNumberThresholdSpecies)
    {
        mIsMinimumThresholdSet[index] = thresholdCheck;
    }
    else
    {
        std::cout<<"Error: SimpleChemicalThresholdCellCycleModel::SetMinimumThresholdCheckByIndex, index out of bounds"<<std::endl;
    }
}

void SimpleChemicalThresholdCellCycleModel::SetThresholdChemistry(AbstractChemistry* p_chemistry)
{
    mpThresholdChemistry = p_chemistry;
}

AbstractChemistry* SimpleChemicalThresholdCellCycleModel::GetThresholdChemistry()
{
    return mpThresholdChemistry;
}

std::vector<double> SimpleChemicalThresholdCellCycleModel::GetMaximumSpeciesThreshold()
{
    return mMaxThresholdConcentrationValues;
}

std::vector<double> SimpleChemicalThresholdCellCycleModel::GetMinimumSpeciesThreshold()
{
    return mMinThresholdConcentrationValues;
}

double SimpleChemicalThresholdCellCycleModel::GetMaximumThresholdByIndex(unsigned index)
{
    if(index<mNumberThresholdSpecies)
    {
        return mMaxThresholdConcentrationValues[index];
    }
    std::cout<<"Error: SimpleChemicalThresholdCellCycleModel::GetMaximumThresholdByIndex, index out of bounds"<<std::endl;
    return 0.0;
}

double SimpleChemicalThresholdCellCycleModel::GetMinimumThresholdByIndex(unsigned index)
{
    if(index<mNumberThresholdSpecies)
    {
        return mMinThresholdConcentrationValues[index];
    }
    std::cout<<"Error: SimpleChemicalThresholdCellCycleModel::GetMinimumThresholdByIndex, index out of bounds"<<std::endl;
    return 0.0;
}

std::vector<double> SimpleChemicalThresholdCellCycleModel::GetSpeciesConcentrations()
{
    return mSpeciesConcentrations;
}

double SimpleChemicalThresholdCellCycleModel::GetSpeciesConcentrationsByIndex(unsigned index)
{
    if(index<mNumberThresholdSpecies)
    {
        return mSpeciesConcentrations[index];
    }
    std::cout<<"Error: SimpleChemicalThresholdCellCycleModel::GetSpeciesConcentrationsByIndex, index out of bounds"<<std::endl;
    return 0.0;
}

unsigned SimpleChemicalThresholdCellCycleModel::GetNumberThresholdSpecies()
{
    return mNumberThresholdSpecies;
}

std::vector<bool> SimpleChemicalThresholdCellCycleModel::GetMaximumThresholdCheck()
{
    return mIsMaximumThresholdSet;
}

std::vector<bool> SimpleChemicalThresholdCellCycleModel::GetMinimumThresholdCheck()
{
    return mIsMinimumThresholdSet;
}

bool SimpleChemicalThresholdCellCycleModel::GetMaximumThresholdCheckByIndex(unsigned index)
{
    if(index<mNumberThresholdSpecies)
    {
        return mIsMaximumThresholdSet[index];
    }
    std::cout<<"Error: SimpleChemicalThresholdCellCycleModel::GetMaximumThresholdCheckByIndex, index out of bounds"<<std::endl;
    return false;
}

bool SimpleChemicalThresholdCellCycleModel::GetMinimumThresholdCheckByIndex(unsigned index)
{
    if(index<mNumberThresholdSpecies)
    {
        return mIsMinimumThresholdSet[index];
    }
    std::cout<<"Error: SimpleChemicalThresholdCellCycleModel::GetMinimumThresholdCheckByIndex, index out of bounds"<<std::endl;
    return false;
}

double SimpleChemicalThresholdCellCycleModel::GetCriticalStarvationDuration() const
{
    return mCriticalStarvationDuration;
}

void SimpleChemicalThresholdCellCycleModel::SetCriticalStarvationDuration(double criticalStarvationDuration)
{
    assert(criticalStarvationDuration >= 0.0);
    mCriticalStarvationDuration = criticalStarvationDuration;
}

void SimpleChemicalThresholdCellCycleModel::SetCurrentStarvationOnsetTime(double currentStarvationOnsetTime)
{
    assert(currentStarvationOnsetTime >= 0.0);
    mCurrentStarvationOnsetTime = currentStarvationOnsetTime;
}
