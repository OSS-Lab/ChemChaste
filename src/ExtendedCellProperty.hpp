#ifndef EXTENDEDCELLPROPERTY_HPP
#define EXTENDEDCELLPROPERTY_HPP

//general includes
#include <vector>
#include <string>
#include <math.h>
#include <boost/shared_ptr.hpp>

// chaste includes
#include "ChemicalCellProperty.hpp"

// base case assume spherical geometry

template<unsigned SPACE_DIM>
class ExtendedCellProperty : public ChemicalCellProperty
{
private:

    // inhertied:
    //StateVariableRegister* mpStateVariableRegister;
    //std::vector<double> mConcentrationVector;

    using ChemicalCellProperty::UpdateCellConcentrationVector;
    using ChemicalCellProperty::InitialiseCell;
protected:

    // base case; one value equal to the radius
    std::vector<double> mCellDimensions = {1.0}; 
    
    // vector of minimal dimensions before dimension is assumed to be null.
    std::vector<double> mMinimalDimensions = {1e-5};

    double mCellVolume =0.0;

    double mVolumePerFineMeshElement=0.0;

    bool mIsCellVolumeDefined = false;

    bool mIsMeshElementVolumeDefined = false;

    unsigned mTotalNumberMeshVoxels=0;

    unsigned mNumberMeshVoxelsInCell=0;

    unsigned mNumberMeshVoxelsOnBoundary=0;

    // is the cell permitted to grow and/or change in dimensions
    bool mIsCellDynamic = false;

    std::vector<double> mChangeInCellDimensions ={0.0};

    std::vector<double> mMeshDomainScale = std::vector<double>();

    // option switchs for interpolating into the bulk FeMesh
    bool mIncludeOdeInterpolationInCell=false;

    bool mAverageInternalCellStates=true;

    // member variables for the cell boundary voxels

    std::vector<std::vector<double>> mVectorOfExternalBoundaryStateVariables;

    std::vector<std::vector<double>> mVectorOfInternalBoundaryStateVariables;

    std::vector<ChastePoint<SPACE_DIM>> mVectorOfBoundaryLocations;

    // double in [0.0,1.0] relating the proportion of cell data avaliable for transport processes
    double mTransportConcentrationAvaliability = 1.0; 

    // resultant cell boundary concentration vector; total contribution to the cell concentration for all the boundary voxels
    std::vector<double> mNextTimestepConcentrationVector;


public:

    ExtendedCellProperty();

    virtual ~ExtendedCellProperty();


    // virtual methods

    virtual void SetUpExtendedCell(StateVariableRegister*, std::vector<double>, std::vector<double>);

    virtual void SetUpExtendedCell(StateVariableRegister*, std::vector<double>, double radius =0.0);

    virtual void UpdateExtendedCell();

    virtual void InitialiseCell(StateVariableRegister*, std::vector<double>);

    virtual void CellGrowth();

    virtual bool IsPointInCell(const ChastePoint<SPACE_DIM>&, ChastePoint<SPACE_DIM>&);

    virtual bool IsPointOnCellBoundary(const ChastePoint<SPACE_DIM>&, ChastePoint<SPACE_DIM>&);

    virtual bool CheckCellBoundary(const ChastePoint<SPACE_DIM>&, ChastePoint<SPACE_DIM>&);

    virtual void CalculateVolumeOfCell();

    virtual void CalculateNumberOfCellMeshElements(ChastePoint<SPACE_DIM>);

    virtual void UpdateCellConcentrationVector(std::vector<double>&);

    virtual double TransportConcentrationAvaliabilityThisState(std::string state_name ="", unsigned state_index =0);


    // concrete methods

    void StepThroughDirectionalXBoundaryVoxels(ChastePoint<SPACE_DIM>&,ChastePoint<SPACE_DIM>&, unsigned);

    void StepThroughDirectionalXInternalVoxels(ChastePoint<SPACE_DIM>&,ChastePoint<SPACE_DIM>&, unsigned);

    void StepThroughXVoxels(ChastePoint<SPACE_DIM>&,ChastePoint<SPACE_DIM>&);

    void StepThroughDirectionalXYBoundaryVoxels(ChastePoint<SPACE_DIM>&,ChastePoint<SPACE_DIM>&, unsigned);

    void StepThroughDirectionalXYInternalVoxels(ChastePoint<SPACE_DIM>&,ChastePoint<SPACE_DIM>&, unsigned);

    void StepThroughXYVoxels(ChastePoint<SPACE_DIM>&,ChastePoint<SPACE_DIM>&);

    void StepThroughDirectionalXYZBoundaryVoxels(ChastePoint<SPACE_DIM>&,ChastePoint<SPACE_DIM>&, unsigned);

    void StepThroughDirectionalXYZInternalVoxels(ChastePoint<SPACE_DIM>&,ChastePoint<SPACE_DIM>&, unsigned);

    void StepThroughXYZVoxels(ChastePoint<SPACE_DIM>&,ChastePoint<SPACE_DIM>&);

    double RetrieveBoundarySourceByStateName(std::string, ChastePoint<SPACE_DIM>);

    double RetrieveInternalCellSourceByStateName(std::string);

    double RetrieveInternalCellSourceByStateName(std::string, ChastePoint<SPACE_DIM>);

    double RetrieveBoundaryCellSourceByStateNameAndLocation(std::string, ChastePoint<SPACE_DIM>);

    double RetrieveBoundaryInternalCellSourceByStateNameAndLocation(std::string, ChastePoint<SPACE_DIM>);

    std::vector<double>& RetrieveBoundaryCellSourceByLocation(ChastePoint<SPACE_DIM>);

    std::vector<double>& RetrieveInternalBoundaryCellSourceByLocation(ChastePoint<SPACE_DIM>);

    bool CheckChastePointsForEquality(ChastePoint<SPACE_DIM>, ChastePoint<SPACE_DIM>);

    void ResetNextTimestepConcentrationVector();

    void ResetNextTimestepConcentrationVector(unsigned);

    void ResetVectorOfBoundaryStateVariables();

    void ResetVectorOfBoundaryLocations();

    void ResetVectorOfInternalBoundaryStateVariables();

    void AddToVectorOfBoundaryStateVariables(std::vector<double>);

    void AddToVectorOfInternalBoundaryStateVariables(std::vector<double>);

    void AddToVectorOfBoundaryLocations(ChastePoint<SPACE_DIM>);

    void AppendInternalCellBoundaryConcentrations(std::vector<double>&, unsigned);

    void ReplaceBoundaryStateVariables(unsigned, std::vector<double>&);

    void RecordLocationAndStateVariable(ChastePoint<SPACE_DIM>, std::vector<double>);

    // Set methods

    void SetCellRadius(double);

    void SetCellDimensions(std::vector<double>);

    void SetMinimalCellDimensions(std::vector<double>);

    void SetCellVolume(double);

    void SetIsCellDynamic(bool);

    void SetMeshDomainScale(std::vector<double>);

    void SetTotalNumberMeshVoxels(unsigned);

    void SetNumberMeshVoxelsInCell(unsigned);

    void SetNumberMeshVoxelsOnBoundary(unsigned);

    void SetChangeInCellDimensions(std::vector<double>);

    void SetIncludeOdeInterpolationInCell(bool);

    void SetAverageInternalCellStates(bool);

    void SetVectorOfBoundaryStateVariables(std::vector<std::vector<double>>);

    void SetVectorOfBoundaryStateVariablesByIndex(unsigned,std::vector<double>);

    void SetVectorOfInternalBoundaryStateVariables(std::vector<std::vector<double>>);

    void SetVectorOfInternalBoundaryStateVariablesByIndex(unsigned,std::vector<double>);

    void SetVectorOfBoundaryLocations(std::vector<ChastePoint<SPACE_DIM>>);

    void SetVectorOfBoundaryLocationsByIndex(unsigned,ChastePoint<SPACE_DIM>);

    void SetNextTimestepConcentrationVector(std::vector<double>);

    // Get methods

    double GetCellRadius();

    std::vector<double> GetCellDimensions();

    std::vector<double> GetMinimalCellDimensions();

    double GetCellVolume();

    bool GetIsCellVolumeDefined();

    bool GetIsMeshElementVolumeDefined();

    bool GetIsCellDynamic();

    std::vector<double> GetMeshDomainScale();

    unsigned GetTotalNumberMeshVoxels();

    unsigned GetNumberMeshVoxelsInCell();

    unsigned GetNumberMeshVoxelsOnBoundary();
    
    std::vector<double> GetChangeInCellDimensions();

    bool GetIncludeOdeInterpolationInCell();

    bool GetAverageInternalCellStates();

    std::vector<std::vector<double>> GetVectorOfBoundaryStateVariables();

    std::vector<double> GetVectorOfBoundaryStateVariablesByLocationIndex(unsigned);

    std::vector<std::vector<double>> GetVectorOfInternalBoundaryStateVariables();

    std::vector<double> GetVectorOfInternalBoundaryStateVariablesByLocationIndex(unsigned);

    std::vector<ChastePoint<SPACE_DIM>> GetVectorOfBoundaryLocations();

    ChastePoint<SPACE_DIM> GetVectorOfBoundaryLocationsByIndex(unsigned);

    std::vector<double> GetNextTimestepConcentrationVector();

};

template<unsigned SPACE_DIM>
ExtendedCellProperty<SPACE_DIM>::ExtendedCellProperty()
: ChemicalCellProperty()
{

}

template<unsigned SPACE_DIM>
ExtendedCellProperty<SPACE_DIM>::~ExtendedCellProperty()
{
}


// virtual methods

template<unsigned SPACE_DIM>
void ExtendedCellProperty<SPACE_DIM>::SetUpExtendedCell(StateVariableRegister* p_StateVariableRegister, std::vector<double> concentrationVector, std::vector<double> cellDimensions)
{
    InitialiseCell(p_StateVariableRegister, concentrationVector);
    SetCellDimensions(cellDimensions);

    std::vector<double> minimalDimensions(SPACE_DIM,mMinimalDimensions[0]);
    SetMinimalCellDimensions(minimalDimensions);

    // instantansiate the cell space
    ResetNextTimestepConcentrationVector(concentrationVector.size());

    std::vector<ChastePoint<SPACE_DIM>> dummyPoint;
    std::vector<std::vector<double>> vectorDummy;
    SetVectorOfBoundaryStateVariables(vectorDummy);
    SetVectorOfBoundaryLocations(dummyPoint);
    CalculateVolumeOfCell();
    CalculateNumberOfCellMeshElements(ChastePoint<SPACE_DIM>());
}

template<unsigned SPACE_DIM>
void ExtendedCellProperty<SPACE_DIM>::SetUpExtendedCell(StateVariableRegister* p_StateVariableRegister, std::vector<double> concentrationVector, double radius)
{
    InitialiseCell(p_StateVariableRegister, concentrationVector);
    std::vector<double> cellDimensions(1,radius);
    SetCellDimensions(cellDimensions);

    // instantansiate the cell space
    ResetNextTimestepConcentrationVector(concentrationVector.size());

    std::vector<ChastePoint<SPACE_DIM>> dummyPoint;
    std::vector<std::vector<double>> vectorDummy;
    SetVectorOfBoundaryStateVariables(vectorDummy);
    SetVectorOfBoundaryLocations(dummyPoint);
    CalculateVolumeOfCell();
    CalculateNumberOfCellMeshElements(ChastePoint<SPACE_DIM>());

}

template<unsigned SPACE_DIM>
void ExtendedCellProperty<SPACE_DIM>::UpdateExtendedCell()
{
    if(mIsCellDynamic)
    {
        CellGrowth();
    }
}

template<unsigned SPACE_DIM>
void ExtendedCellProperty<SPACE_DIM>::InitialiseCell(StateVariableRegister* p_StateVariableRegister, std::vector<double> concentrationVector)
{
    ChemicalCellProperty::InitialiseCell(p_StateVariableRegister, concentrationVector);
}

template<unsigned SPACE_DIM>
void ExtendedCellProperty<SPACE_DIM>::CellGrowth()
{
    // apply the change in cell dimensions to the cell dimension vector
    for(unsigned i=0; i<mCellDimensions.size();i++)
    {
        mCellDimensions[i] += mChangeInCellDimensions[i];
    }
    
    CalculateVolumeOfCell();
    CalculateNumberOfCellMeshElements(ChastePoint<SPACE_DIM>());
}

template<unsigned SPACE_DIM>
bool ExtendedCellProperty<SPACE_DIM>::IsPointInCell(const ChastePoint<SPACE_DIM>& cellCentrePoint, ChastePoint<SPACE_DIM>& test_point)
{
    // assume the cell is SPACE_DIM spherical with radius mCellDimensions[0]

    std::vector<double> difference;
    double square_difference_distance=0.0;
    for(unsigned dim=0; dim<SPACE_DIM; dim++)
    {
        square_difference_distance += pow((cellCentrePoint.rGetLocation()[dim] - test_point.rGetLocation()[dim]),2.0);
    }

    if(sqrt(square_difference_distance)<=mCellDimensions[0])
    {
        return true;
    }

    return false;
}

template<unsigned SPACE_DIM>
bool ExtendedCellProperty<SPACE_DIM>::IsPointOnCellBoundary(const ChastePoint<SPACE_DIM>& cellCentrePoint, ChastePoint<SPACE_DIM>& test_point)
{
    // assume the cell is SPACE_DIM spherical with radius mCellDimensions[0]
    // if the test_point is within the cell but a secondary test point either side is not 
    // within the cell we can say the test_point is a boundary point. Use the mMeshDomainScale
    // to determine the locations of the points either side. Stepsize in x mMeshDomainScale[0] 
    // stepsize in y mMeshDomainScale[1], stepsize in z mMeshDomainScale[2] .

    if(IsPointInCell(cellCentrePoint,test_point))
    {
        return CheckCellBoundary(cellCentrePoint,test_point);
    }
    // test_point is not within the cell so cannot possibly be a boundary point
    return false;
}

template<unsigned SPACE_DIM>
bool ExtendedCellProperty<SPACE_DIM>::CheckCellBoundary(const ChastePoint<SPACE_DIM>& cellCentrePoint, ChastePoint<SPACE_DIM>& test_point)
{
    // the number of test cases depends on SPACE_DIM
    // SPACE_DIM=1 the single x dimension; 2 cases positive and negative x
    // SPACE_DIM=2 4 cases; positive and negative x takes precedence but also positive and negative y
    // SPACE_DIM=3 6 cases; positive and negative x, positive and negative y, positive and negative z

    ChastePoint<SPACE_DIM> secondary_test_point = test_point;
    for(unsigned dim=0; dim<SPACE_DIM; dim++)
    {
        secondary_test_point.SetCoordinate(dim, test_point[dim]+ mCellDimensions[dim]);
        if(!IsPointInCell(cellCentrePoint, secondary_test_point))
        {
            return true;
            break;
        }
        secondary_test_point.SetCoordinate(dim, test_point[dim] - mCellDimensions[dim]);
        if(!IsPointInCell(cellCentrePoint, secondary_test_point))
        {
            return true;
            break;
        }
        // reset the secondary test point
        secondary_test_point = test_point;
    }

    return false;
}

template<unsigned SPACE_DIM>
void ExtendedCellProperty<SPACE_DIM>::CalculateVolumeOfCell()
{
    if(mCellDimensions[0]<mMinimalDimensions[0])
    {
        // i.e the cell radius is below the minimal dimension threshold suggesting
        // radius is true zero or at least negligible. Set the volume ot be zero.
        mCellVolume = 0.0;
        mIsCellVolumeDefined =false;
    }
    else
    {
        // have a finite radius so calculate the volume in the spatial dimensions
        switch(SPACE_DIM)
        {
            case 2:
                mCellVolume = M_PI*pow(mCellDimensions[0],2.0);
                mIsCellVolumeDefined =true;
                break;

            case 3:
                mCellVolume = (4/3)*M_PI*pow(mCellDimensions[0],3.0);
                mIsCellVolumeDefined =true;
                break;
            default:
                mCellVolume = 0.0;
                mIsCellVolumeDefined =false;
        }
    }
}

template<unsigned SPACE_DIM>
void ExtendedCellProperty<SPACE_DIM>::CalculateNumberOfCellMeshElements(ChastePoint<SPACE_DIM> cellCentrePoint)
{

    // over ride previous calculations
    mIsMeshElementVolumeDefined = false;

    mTotalNumberMeshVoxels=0;

    mNumberMeshVoxelsInCell=0;

    mNumberMeshVoxelsOnBoundary=0;

    if(mCellDimensions[0]<mMinimalDimensions[0])
    {
        // i.e the cell radius is below the minimal dimension threshold suggesting
        // radius is true zero or at least negligible.

        // cell must exist at a point which is considered a boundary for transport phenomena
        mTotalNumberMeshVoxels =1;
        mNumberMeshVoxelsOnBoundary =1;
        mIsMeshElementVolumeDefined =true;
    }
    else
    {
        // have a finite radius so calculate the volume in terms of mesh voxels in the spatial dimensions
        if(CheckCellBoundary(cellCentrePoint,cellCentrePoint))
        {
            // if the cell centre is a boundary point then the central voxel is the only voxel within the cell 
            mTotalNumberMeshVoxels =1;
            mNumberMeshVoxelsOnBoundary =1;
            mIsMeshElementVolumeDefined =true;
        }
        else
        {
            
            switch(SPACE_DIM)
            {
                case 1:
                    // line segement; some outer boundary voxels enclose some cell voxels
                    // step along the line from the point start_point. Update the cell and boundary voxels
                    
                    StepThroughXVoxels(cellCentrePoint,cellCentrePoint);

                    break;

                case 2:
                    // 2D plane segment; a 2D shell of boundary voxels surrounding some internal cell voxels
                    // step along the y direction and run through the whole x line for a given y
                    
                    StepThroughXYVoxels(cellCentrePoint,cellCentrePoint);

                    break;

                case 3:
                    // 3D volume segment; a 3D shell of boundary voxels surrounding some internal cell voxels
                    // step along the z direction and run through the whole xy plane for a given z
                    
                    StepThroughXYZVoxels(cellCentrePoint,cellCentrePoint);

                    break;
                default:
                    // cell must exist at a point which is considered a boundary for transport phenomena
                    mTotalNumberMeshVoxels =1;
                    mNumberMeshVoxelsOnBoundary =1;
            }
            mIsMeshElementVolumeDefined =true;
            mTotalNumberMeshVoxels = mNumberMeshVoxelsOnBoundary + mNumberMeshVoxelsInCell;
        }
        
    }

}

template<unsigned SPACE_DIM>
void ExtendedCellProperty<SPACE_DIM>::UpdateCellConcentrationVector(std::vector<double>& cell_concentration_vector_from_cell_data)
{
    // input a concentration vector of the size and ordering of the internal cell concentration vector
    std::vector<double> concentrationVector(cell_concentration_vector_from_cell_data.size(),0.0);

    double constant_state_avaliability = TransportConcentrationAvaliabilityThisState();

    for(unsigned i=0; i<cell_concentration_vector_from_cell_data.size(); i++)
    {
        concentrationVector[i] =  constant_state_avaliability*cell_concentration_vector_from_cell_data[i];
        cell_concentration_vector_from_cell_data[i] -= concentrationVector[i];
    }

    mConcentrationVector = concentrationVector;
}

template<unsigned SPACE_DIM>
double ExtendedCellProperty<SPACE_DIM>::TransportConcentrationAvaliabilityThisState(std::string, unsigned)
{
    // return double in [0.0,1.0] representing the proportion of the cell data that may be avaliable for cell boundary
    // transport.  May be dependent of state properties such as intra cell diffusivities, charges etc.
    // default return mTransportConcentrationAvaliability;
    return mTransportConcentrationAvaliability;
}






// concrete methods

template<unsigned SPACE_DIM>
void ExtendedCellProperty<SPACE_DIM>::StepThroughDirectionalXBoundaryVoxels(ChastePoint<SPACE_DIM>& cellCentrePoint,ChastePoint<SPACE_DIM>& start_point, unsigned polarity)
{
    // start point is on the boundary, step along the boundary until the voxel no longer belongs in the cell
    
    // check for contigual boundary points in both positive nad negative x directions
    // positive x direction
    bool test_point_is_on_boundary = true;
    ChastePoint<SPACE_DIM> test_point = start_point;
    test_point.SetCoordinate(0, start_point[0]);
    unsigned x_step_counter=1;
    double deltaX = mMeshDomainScale[0]; // step size of the FeMesh
    if(polarity>0.0)
    {
        while(test_point_is_on_boundary)
        {
            // while we are on the boundary need to check for contigual voxels also being on the boundary
            // will occur when x plane aligns with boundary tangent. Check that next voxel in x is in cell
            test_point.SetCoordinate(0, start_point[0] + x_step_counter*deltaX);
            test_point_is_on_boundary = IsPointInCell(cellCentrePoint,test_point);
            mNumberMeshVoxelsOnBoundary += (unsigned)test_point_is_on_boundary;
            x_step_counter +=1;
        }
    }
    else
    {
        while(test_point_is_on_boundary)
        {
            // while we are on the boundary need to check for contigual voxels also being on the boundary
            // will occur when x plane aligns with boundary tangent. Check that next voxel in x is in cell
            test_point.SetCoordinate(0, start_point[0] - x_step_counter*deltaX);
            test_point_is_on_boundary = IsPointInCell(cellCentrePoint,test_point);
            mNumberMeshVoxelsOnBoundary += (unsigned)test_point_is_on_boundary;
            x_step_counter +=1;
        }
    }
}

template<unsigned SPACE_DIM>
void ExtendedCellProperty<SPACE_DIM>::StepThroughDirectionalXInternalVoxels(ChastePoint<SPACE_DIM>& cellCentrePoint,ChastePoint<SPACE_DIM>& start_point, unsigned polarity)
{
    // starting point is in the cell
    bool test_point_is_on_boundary = false;
    ChastePoint<SPACE_DIM> test_point = start_point;
    test_point.SetCoordinate(0, start_point[0]);
    unsigned x_step_counter=1;
    double deltaX = mMeshDomainScale[0]; // step size in x of FeMesh
    if(polarity>0.0)
    {
        // positive direction
        while(!test_point_is_on_boundary)
        {
            test_point.SetCoordinate(0, start_point[0] + x_step_counter*deltaX);
            test_point_is_on_boundary = CheckCellBoundary(cellCentrePoint,test_point);
            mNumberMeshVoxelsInCell += (unsigned)!test_point_is_on_boundary;
            x_step_counter +=1;
        }

    }
    else
    {
        // negative direction
        while(!test_point_is_on_boundary)
        {
            test_point.SetCoordinate(0, start_point[0] - x_step_counter*deltaX);
            test_point_is_on_boundary = CheckCellBoundary(cellCentrePoint,test_point);
            mNumberMeshVoxelsInCell += (unsigned)!test_point_is_on_boundary;
            x_step_counter +=1;
        }
    }
    
}

template<unsigned SPACE_DIM>
void ExtendedCellProperty<SPACE_DIM>::StepThroughXVoxels(ChastePoint<SPACE_DIM>& cellCentrePoint,ChastePoint<SPACE_DIM>& start_point)
{
    // perturb the centre point across x dimension until hit the boundary

    // assume the start point is a voxel of the cell

    // initialise procedure variables
    ChastePoint<SPACE_DIM> test_point = start_point;
    bool test_point_is_on_boundary = false;
    unsigned x_step_counter=1; // x_step_counter = 0 refers to the original start point

    // check if the start point is on the boundary
    if(CheckCellBoundary(cellCentrePoint,start_point))
    {
        mNumberMeshVoxelsOnBoundary +=1;
        test_point.SetCoordinate(0, start_point[0]);
        StepThroughDirectionalXBoundaryVoxels(cellCentrePoint,test_point,1);

        // negative x direction
         // re-center
        test_point.SetCoordinate(0, start_point[0]);
        StepThroughDirectionalXBoundaryVoxels(cellCentrePoint,test_point,-1);
    }
    else
    {   // starting within the cell body, begin stepping procedure
        mNumberMeshVoxelsInCell +=1; // for the start point

        // positive x direction
        test_point.SetCoordinate(0, start_point[0]);
        StepThroughDirectionalXInternalVoxels(cellCentrePoint,test_point,1);
        // ends when we hit a boundary cell
        mNumberMeshVoxelsOnBoundary +=1;
        // continue the stepping but along the boundary, don't recenter the test_point
        StepThroughDirectionalXBoundaryVoxels(cellCentrePoint,test_point,1);

        
        // negative x direction
        test_point.SetCoordinate(0, start_point[0]); // re set starting test_point
        StepThroughDirectionalXInternalVoxels(cellCentrePoint,test_point,-1);
        // ends when we hit a boundary cell
        mNumberMeshVoxelsOnBoundary +=1;
        // continue the stepping but along the boundary, don't recenter the test_point
        StepThroughDirectionalXBoundaryVoxels(cellCentrePoint,test_point,-1);

    }
}

template<unsigned SPACE_DIM>
void ExtendedCellProperty<SPACE_DIM>::StepThroughDirectionalXYBoundaryVoxels(ChastePoint<SPACE_DIM>& cellCentrePoint,ChastePoint<SPACE_DIM>& start_point, unsigned polarity)
{
    // start point is on the boundary
    // step along the boundary in the y direction and step through the x direction 
    // until the voxel no longer belongs in the cell
    
    // check for contigual boundary points in both positive nad negative x directions
    // positive x direction
    bool test_point_is_on_boundary = true;
    ChastePoint<SPACE_DIM> test_point = start_point;
    test_point.SetCoordinate(1, start_point[1]);
    unsigned y_step_counter=1;
    double deltaY = mMeshDomainScale[1]; // step size in y of FeMesh
    if(polarity>0.0)
    {
        while(test_point_is_on_boundary)
        {
            // while we are on the boundary need to check for contigual voxels also being on the boundary
            // will occur when XY plane aligns with boundary tangent. Check that next voxel in Y is in cell
            test_point.SetCoordinate(1, start_point[1] + y_step_counter*deltaY);
            test_point_is_on_boundary = IsPointInCell(cellCentrePoint,test_point);
            
            
            if(test_point_is_on_boundary)
            {
                // next voxel step in Y is in boundary, run through the X slice from this test_point
                StepThroughXVoxels(cellCentrePoint,start_point);
            }
            // if next point is not on the boundary then it isn't in the cell and subsequent voxels will also not
            // be in the cell with respect to x

            y_step_counter +=1;
        }
    }
    else
    {
        while(test_point_is_on_boundary)
        {
            // while we are on the boundary need to check for contigual voxels also being on the boundary
            // will occur when XY plane aligns with boundary tangent. Check that next voxel in Y is in cell
            test_point.SetCoordinate(1, start_point[1] - y_step_counter*deltaY);
            test_point_is_on_boundary = IsPointInCell(cellCentrePoint,test_point);
            
            
            if(test_point_is_on_boundary)
            {
                // next voxel step in Y is in boundary, run through the X slice from this test_point
                StepThroughXVoxels(cellCentrePoint,start_point);
            }
            // if next point is not on the boundary then it isn't in the cell and subsequent voxels will also not
            // be in the cell with respect to x
            
            y_step_counter +=1;
        }
    }
}

template<unsigned SPACE_DIM>
void ExtendedCellProperty<SPACE_DIM>::StepThroughDirectionalXYInternalVoxels(ChastePoint<SPACE_DIM>& cellCentrePoint,ChastePoint<SPACE_DIM>& start_point, unsigned polarity)
{
    // starting point is in the cell
    bool test_point_is_on_boundary = false;
    ChastePoint<SPACE_DIM> test_point = start_point;
    test_point.SetCoordinate(1, start_point[1]);
    unsigned y_step_counter=1;
    double deltaY = mMeshDomainScale[1]; // step size in y of FeMesh
    if(polarity>0.0)
    {
        // positive direction
        while(!test_point_is_on_boundary)
        {
            test_point.SetCoordinate(1, start_point[1] + y_step_counter*deltaY);
            test_point_is_on_boundary = CheckCellBoundary(cellCentrePoint,test_point);
            
            // if the next y step voxel is not on the boundary run through the x slice
            if(!test_point_is_on_boundary)
            {
                // next voxel step in Y is in boundary, run through the X slice from this test_point
                StepThroughXVoxels(cellCentrePoint,test_point);
            }
            // if the next y step voxel is on the boundary end
            y_step_counter +=1;
            
        }
    }
    else
    {
        // negative direction
        while(!test_point_is_on_boundary)
        {
            test_point.SetCoordinate(1, start_point[1] - y_step_counter*deltaY);
            test_point_is_on_boundary = CheckCellBoundary(cellCentrePoint,test_point);
            
            // if the next y step voxel is not on the boundary run through the x slice
            if(!test_point_is_on_boundary)
            {
                // next voxel step in Y is in boundary, run through the X slice from this test_point
                StepThroughXVoxels(cellCentrePoint,test_point);
            }
            // if the next y step voxel is on the boundary end
            y_step_counter +=1;
        }
    }
    
}

template<unsigned SPACE_DIM>
void ExtendedCellProperty<SPACE_DIM>::StepThroughXYVoxels(ChastePoint<SPACE_DIM>& cellCentrePoint,ChastePoint<SPACE_DIM>& start_point)
{
    // consider the 2D case, (X,Y). For each step in Y we need to step through the corresponding x 

    // assume the start point is a voxel of the cell

    // initialise procedure variables
    ChastePoint<SPACE_DIM> test_point = start_point;

    // check if the start point is on the boundary
    if(CheckCellBoundary(cellCentrePoint,start_point))
    {
        mNumberMeshVoxelsOnBoundary +=1;
        // positive y direction
        test_point.SetCoordinate(1, start_point[1]);
        StepThroughDirectionalXYBoundaryVoxels(cellCentrePoint,test_point,1);

        // negative y direction
        // re-center
        test_point.SetCoordinate(1, start_point[1]);
        StepThroughDirectionalXYBoundaryVoxels(cellCentrePoint,test_point,-1);
    }
    else
    {   // starting within the cell body, begin stepping procedure
        mNumberMeshVoxelsInCell +=1; // for the start point

        // positive x direction
        test_point.SetCoordinate(1, start_point[1]);
        StepThroughDirectionalXYInternalVoxels(cellCentrePoint,test_point,1);
        // ends when we hit a boundary cell
        mNumberMeshVoxelsOnBoundary +=1;
        // continue the stepping but along the boundary, don't recenter the test_point
        StepThroughDirectionalXYBoundaryVoxels(cellCentrePoint,test_point,1);

        
        // negative x direction
        test_point.SetCoordinate(1, start_point[1]); // re set starting test_point
        StepThroughDirectionalXYInternalVoxels(cellCentrePoint,test_point,-1);
        // ends when we hit a boundary cell
        mNumberMeshVoxelsOnBoundary +=1;
        // continue the stepping but along the boundary, don't recenter the test_point
        StepThroughDirectionalXYBoundaryVoxels(cellCentrePoint,test_point,-1);

    }

}

template<unsigned SPACE_DIM>
void ExtendedCellProperty<SPACE_DIM>::StepThroughDirectionalXYZBoundaryVoxels(ChastePoint<SPACE_DIM>& cellCentrePoint,ChastePoint<SPACE_DIM>& start_point, unsigned polarity)
{
    // start point is on the boundary
    // step along the boundary in the z direction and step through the x and Y directions recursively
    // until the voxel no longer belongs in the cell
    
    // check for contigual boundary points in both positive and negative y directions
    // positive y direction
    bool test_point_is_on_boundary = true;
    ChastePoint<SPACE_DIM> test_point = start_point;
    test_point.SetCoordinate(2, start_point[2]);
    unsigned z_step_counter=1;
    double deltaZ = mMeshDomainScale[2]; // step size in z of FeMesh
    if(polarity>0.0)
    {
        while(test_point_is_on_boundary)
        {
            // while we are on the boundary need to check for contigual voxels also being on the boundary
            // will occur when XY plane aligns with boundary tangent. Check that next voxel in Y is in cell
            test_point.SetCoordinate(2, start_point[2] + z_step_counter*deltaZ);
            test_point_is_on_boundary = IsPointInCell(cellCentrePoint,test_point);
            
            
            if(test_point_is_on_boundary)
            {
                // next voxel step in Z is in boundary, run through the XY slice from this test_point
                StepThroughXYVoxels(cellCentrePoint,start_point);
            }
            // if next point is not on the boundary then it isn't in the cell and subsequent voxels will also not
            // be in the cell with respect to Y

            z_step_counter +=1;
        }
    }
    else
    {
        while(test_point_is_on_boundary)
        {
            // while we are on the boundary need to check for contigual voxels also being on the boundary
            // will occur when XYZ shell aligns with boundary tangent. Check that next voxel in Z is in cell
            test_point.SetCoordinate(2, start_point[2] - z_step_counter*deltaZ);
            test_point_is_on_boundary = IsPointInCell(cellCentrePoint,test_point);
            
            
            if(test_point_is_on_boundary)
            {
                // next voxel step in Z is in boundary, run through the Y slice from this test_point
                StepThroughXYVoxels(cellCentrePoint,start_point);
            }
            // if next point is not on the boundary then it isn't in the cell and subsequent voxels will also not
            // be in the cell with respect to y
            
            z_step_counter +=1;
        }
    }
}

template<unsigned SPACE_DIM>
void ExtendedCellProperty<SPACE_DIM>::StepThroughDirectionalXYZInternalVoxels(ChastePoint<SPACE_DIM>& cellCentrePoint,ChastePoint<SPACE_DIM>& start_point, unsigned polarity)
{
    // starting point is in the cell
    bool test_point_is_on_boundary = false;
    ChastePoint<SPACE_DIM> test_point = start_point;
    test_point.SetCoordinate(2, start_point[2]);
    unsigned z_step_counter=1;
    double deltaZ = mMeshDomainScale[1]; // step size in z of FeMesh
    if(polarity>0.0)
    {
        // positive direction
        while(!test_point_is_on_boundary)
        {
            test_point.SetCoordinate(2, start_point[2] + z_step_counter*deltaZ);
            test_point_is_on_boundary = CheckCellBoundary(cellCentrePoint,test_point);
            
            // if the next y step voxel is not on the boundary run through the x slice
            if(!test_point_is_on_boundary)
            {
                // next voxel step in Z is in boundary, run through the Y slice from this test_point
                StepThroughXYVoxels(cellCentrePoint,test_point);
            }
            // if the next z step voxel is on the boundary end
            z_step_counter +=1;
            
        }
    }
    else
    {
        // negative direction
        while(!test_point_is_on_boundary)
        {
            test_point.SetCoordinate(2, start_point[2] - z_step_counter*deltaZ);
            test_point_is_on_boundary = CheckCellBoundary(cellCentrePoint,test_point);
            
            // if the next z step voxel is not on the boundary run through the y slice
            if(!test_point_is_on_boundary)
            {
                // next voxel step in Z is in boundary, run through the Y slice from this test_point
                StepThroughXYVoxels(cellCentrePoint,test_point);
            }
            // if the next z step voxel is on the boundary end
            z_step_counter +=1;
        }
    }
    
}

template<unsigned SPACE_DIM>
void ExtendedCellProperty<SPACE_DIM>::StepThroughXYZVoxels(ChastePoint<SPACE_DIM>& cellCentrePoint,ChastePoint<SPACE_DIM>& start_point)
{
    // consider the 3D case, (X,Y,Z). For each step in Z we need to step through the corresponding Y and X slices

    // assume the start point is a voxel of the cell

    // initialise procedure variables
    ChastePoint<SPACE_DIM> test_point = start_point;

    // check if the start point is on the boundary
    if(CheckCellBoundary(cellCentrePoint,start_point))
    {
        mNumberMeshVoxelsOnBoundary +=1;
        // positive y direction
        test_point.SetCoordinate(2, start_point[2]);
        StepThroughDirectionalXYZBoundaryVoxels(cellCentrePoint,test_point,1);

        // negative y direction
        // re-center
        test_point.SetCoordinate(2, start_point[2]);
        StepThroughDirectionalXYZBoundaryVoxels(cellCentrePoint,test_point,-1);
    }
    else
    {   // starting within the cell body, begin stepping procedure
        mNumberMeshVoxelsInCell +=1; // for the start point

        // positive x direction
        test_point.SetCoordinate(2, start_point[2]);
        StepThroughDirectionalXYZInternalVoxels(cellCentrePoint,test_point,1);
        // ends when we hit a boundary cell
        mNumberMeshVoxelsOnBoundary +=1;
        // continue the stepping but along the boundary, don't recenter the test_point
        StepThroughDirectionalXYZBoundaryVoxels(cellCentrePoint,test_point,1);

        
        // negative x direction
        test_point.SetCoordinate(2, start_point[2]); // re set starting test_point
        StepThroughDirectionalXYZInternalVoxels(cellCentrePoint,test_point,-1);
        // ends when we hit a boundary cell
        mNumberMeshVoxelsOnBoundary +=1;
        // continue the stepping but along the boundary, don't recenter the test_point
        StepThroughDirectionalXYZBoundaryVoxels(cellCentrePoint,test_point,-1);

    }
}

template<unsigned SPACE_DIM>
double ExtendedCellProperty<SPACE_DIM>::RetrieveBoundarySourceByStateName(std::string stateName, ChastePoint<SPACE_DIM> rX)
{
    // if cell is on boundary add the result of the transport Ode system
    unsigned index = mpStateVariableRegister -> RetrieveStateVariableIndex(stateName);

    // retrieve the extended cell property
    //boost::shared_ptr<ExtendedCellProperty<SPACE_DIM>> extended_cell_property = boost::static_pointer_cast<ExtendedCellProperty<SPACE_DIM>>(mThis_cellPtr->rGetCellPropertyCollection().GetPropertiesType<ExtendedCellProperty<SPACE_DIM>>().GetProperty());

    return RetrieveBoundaryCellSourceByStateNameAndLocation(stateName, rX);
}

template<unsigned SPACE_DIM>
double ExtendedCellProperty<SPACE_DIM>::RetrieveInternalCellSourceByStateName(std::string stateName)
{
    // the value to give the variable stateName when the point is in the cell

    unsigned index = mpStateVariableRegister -> RetrieveStateVariableIndex(stateName);

    return mConcentrationVector[index];
}


template<unsigned SPACE_DIM>
double ExtendedCellProperty<SPACE_DIM>::RetrieveBoundaryCellSourceByStateNameAndLocation(std::string stateName, ChastePoint<SPACE_DIM> rX_test)
{
    // run through the registered boundary locations to find the one that is congruent with the test point

    for(unsigned index=0; index<mVectorOfBoundaryLocations.size();index++)
    {
        if(CheckChastePointsForEquality(rX_test,mVectorOfBoundaryLocations[index]))
        {
            // return the corresponding state variable
            return mVectorOfExternalBoundaryStateVariables[index][mpStateVariableRegister -> RetrieveStateVariableIndex(stateName)];
        }
    }
    return 0.0;
}


template<unsigned SPACE_DIM>
double ExtendedCellProperty<SPACE_DIM>::RetrieveBoundaryInternalCellSourceByStateNameAndLocation(std::string stateName, ChastePoint<SPACE_DIM> rX_test)
{
    // run through the registered boundary locations to find the one that is congruent with the test point

    for(unsigned index=0; index<mVectorOfBoundaryLocations.size;index++)
    {
        if(CheckChastePointsForEquality(rX_test,mVectorOfBoundaryLocations[index]))
        {
            // return the corresponding state variable
            return mVectorOfInternalBoundaryStateVariables[index][mpStateVariableRegister -> RetrieveStateVariableIndex(stateName)];
        }
    }
    return 0.0;
}


template<unsigned SPACE_DIM>
std::vector<double>& ExtendedCellProperty<SPACE_DIM>::RetrieveBoundaryCellSourceByLocation(ChastePoint<SPACE_DIM> rX_test)
{
    // run through the registered boundary locations to find the one that is congruent with the test point

    for(unsigned index=0; index<mVectorOfBoundaryLocations.size;index++)
    {
        if(CheckChastePointsForEquality(rX_test,mVectorOfBoundaryLocations[index]))
        {
            // return the corresponding state variable vector
            return mVectorOfExternalBoundaryStateVariables[index];
        }
    }
    std::vector<double> empty_return(mVectorOfExternalBoundaryStateVariables[0].size(),0.0);

    return empty_return;
}

template<unsigned SPACE_DIM>
std::vector<double>& ExtendedCellProperty<SPACE_DIM>::RetrieveInternalBoundaryCellSourceByLocation(ChastePoint<SPACE_DIM> rX_test)
{
    // run through the registered boundary locations to find the one that is congruent with the test point

    for(unsigned index=0; index<mVectorOfBoundaryLocations.size;index++)
    {
        if(CheckChastePointsForEquality(rX_test,mVectorOfBoundaryLocations[index]))
        {
            // return the corresponding state variable vector
            return mVectorOfInternalBoundaryStateVariables[index];
        }
    }
    std::vector<double> empty_return(mVectorOfInternalBoundaryStateVariables[0].size(),0.0);

    return empty_return;
}


template<unsigned SPACE_DIM>
bool ExtendedCellProperty<SPACE_DIM>::CheckChastePointsForEquality(ChastePoint<SPACE_DIM> rX1,ChastePoint<SPACE_DIM> rX2)
{
    // function to check whether two chaste points are congruent
    // use half the scale dimension as an error

    bool equality=true;

    for(unsigned i=0; i<SPACE_DIM; i++)
    {
        if(rX1[i] < rX2[i] + 0.5*mMeshDomainScale[i] && rX1[i] >= rX2[i] - 0.5*mMeshDomainScale[i])
        {
            // both points are congruent within 0.5 of mesh scale in the ith dimension
            equality=true;
        }
        else
        {
            // as soon as a dimension where the points are disaparate is found end the function return false
            equality=false;
            break;
        }
        
    }
    return equality;
}


template<unsigned SPACE_DIM>
void ExtendedCellProperty<SPACE_DIM>::ResetNextTimestepConcentrationVector()
{   
    std::vector<double> resetVector(mNextTimestepConcentrationVector.size(),0.0);
    mNextTimestepConcentrationVector = resetVector;
}

template<unsigned SPACE_DIM>
void ExtendedCellProperty<SPACE_DIM>::ResetNextTimestepConcentrationVector(unsigned sizeOfConcentrationVector)
{   
    std::vector<double> resetVector(sizeOfConcentrationVector,0.0);
    mNextTimestepConcentrationVector = resetVector;
}

template<unsigned SPACE_DIM>
void ExtendedCellProperty<SPACE_DIM>::ResetVectorOfBoundaryStateVariables()
{
    mVectorOfExternalBoundaryStateVariables.clear();
}

template<unsigned SPACE_DIM>
void ExtendedCellProperty<SPACE_DIM>::ResetVectorOfBoundaryLocations()
{
    mVectorOfBoundaryLocations.clear();
}

template<unsigned SPACE_DIM>
void ExtendedCellProperty<SPACE_DIM>::ResetVectorOfInternalBoundaryStateVariables()
{
    mVectorOfInternalBoundaryStateVariables.clear();
}

template<unsigned SPACE_DIM>
void ExtendedCellProperty<SPACE_DIM>::AddToVectorOfBoundaryStateVariables(std::vector<double> boundaryStateVariables)
{
    mVectorOfExternalBoundaryStateVariables.push_back(boundaryStateVariables);
}

template<unsigned SPACE_DIM>
void ExtendedCellProperty<SPACE_DIM>::AddToVectorOfInternalBoundaryStateVariables(std::vector<double> boundaryStateVariables)
{
    mVectorOfInternalBoundaryStateVariables.push_back(boundaryStateVariables);
}

template<unsigned SPACE_DIM>
void ExtendedCellProperty<SPACE_DIM>::AddToVectorOfBoundaryLocations(ChastePoint<SPACE_DIM> point)
{
    mVectorOfBoundaryLocations.push_back(point);
}

template<unsigned SPACE_DIM>
void ExtendedCellProperty<SPACE_DIM>::AppendInternalCellBoundaryConcentrations(std::vector<double>& rY, unsigned boundary_index)
{
    

    if(mVectorOfBoundaryLocations.empty())
    {
        if(mNumberMeshVoxelsOnBoundary>0)
        {
            // scale the cell concentration by the number of boundary voxels shared
            for(unsigned i=0; i<mConcentrationVector.size(); i++)
            {
                rY.push_back(mConcentrationVector[i]/mNumberMeshVoxelsOnBoundary);
            }
        }
        else
        {
            // insert without scaling
            rY.insert(rY.end(), mConcentrationVector.begin(), mConcentrationVector.end());
        }
    }
    else
    {
        // use the location index to determine the internal cell concentration
        std::vector<double> internalStates = GetVectorOfInternalBoundaryStateVariablesByLocationIndex(boundary_index);

        rY.insert(rY.end(), internalStates.begin(), internalStates.end());

    }
    
}


template<unsigned SPACE_DIM>
void ExtendedCellProperty<SPACE_DIM>::ReplaceBoundaryStateVariables(unsigned boundary_index, std::vector<double>& rY)
{
    unsigned numberOfStates = mVectorOfExternalBoundaryStateVariables[boundary_index].size();
    unsigned numberOfCellStates = rY.size() - numberOfStates;
    std::vector<double> boundaryExternalStateVariables(numberOfStates,0.0);
    
    for(unsigned i=0; i<numberOfStates;i++)
    {
        boundaryExternalStateVariables[i] = rY[i];
    }

    SetVectorOfBoundaryStateVariablesByIndex(boundary_index, boundaryExternalStateVariables);

    std::vector<double> boundaryInternalStateVariables(numberOfCellStates,0.0);

    for(unsigned i=0; i<numberOfCellStates;i++)
    {
        boundaryInternalStateVariables[i] = rY[i + numberOfStates];
        mNextTimestepConcentrationVector[i] += rY[i + numberOfStates];
    }

    SetVectorOfInternalBoundaryStateVariablesByIndex(boundary_index, boundaryInternalStateVariables);
}

template<unsigned SPACE_DIM>
void ExtendedCellProperty<SPACE_DIM>::RecordLocationAndStateVariable(ChastePoint<SPACE_DIM> rX, std::vector<double> rU)
{
    AddToVectorOfBoundaryStateVariables(rU);
    AddToVectorOfBoundaryLocations(rX);
}



// set methods

template<unsigned SPACE_DIM>
void ExtendedCellProperty<SPACE_DIM>::SetCellRadius(double radius)
{   
    mCellDimensions[0] = radius;
}

template<unsigned SPACE_DIM>
void ExtendedCellProperty<SPACE_DIM>::SetCellDimensions(std::vector<double> cellDimensions)
{
    mCellDimensions = cellDimensions;
}

template<unsigned SPACE_DIM>
void ExtendedCellProperty<SPACE_DIM>::SetMinimalCellDimensions(std::vector<double> minimalDimensions)
{
    mMinimalDimensions = minimalDimensions;
}

template<unsigned SPACE_DIM>
void ExtendedCellProperty<SPACE_DIM>::SetCellVolume(double volume)
{
    mCellVolume = volume;
}

template<unsigned SPACE_DIM>
void ExtendedCellProperty<SPACE_DIM>::SetIsCellDynamic(bool isCellDynamic)
{
    mIsCellDynamic = isCellDynamic;
}

template<unsigned SPACE_DIM>
void ExtendedCellProperty<SPACE_DIM>::SetMeshDomainScale(std::vector<double> meshDomainScale)
{
    mMeshDomainScale = meshDomainScale;

    switch(SPACE_DIM)
    {
        case 2:
            mVolumePerFineMeshElement = meshDomainScale[0]*meshDomainScale[1];
            mIsMeshElementVolumeDefined = true;
            break;
        case 3:
            mVolumePerFineMeshElement = meshDomainScale[0]*meshDomainScale[1]*meshDomainScale[2];
            mIsMeshElementVolumeDefined = true;
            break;
        default:
            mVolumePerFineMeshElement = meshDomainScale[0];
            mIsMeshElementVolumeDefined = true;
    }
}

template<unsigned SPACE_DIM>
void ExtendedCellProperty<SPACE_DIM>::SetTotalNumberMeshVoxels(unsigned totalNumberMeshVoxels)
{
    mTotalNumberMeshVoxels = totalNumberMeshVoxels;
}

template<unsigned SPACE_DIM>
void ExtendedCellProperty<SPACE_DIM>::SetNumberMeshVoxelsInCell(unsigned numberInCell)
{
    mNumberMeshVoxelsInCell = numberInCell;
}

template<unsigned SPACE_DIM>
void ExtendedCellProperty<SPACE_DIM>::SetNumberMeshVoxelsOnBoundary(unsigned numberOnBoundary)
{
    mNumberMeshVoxelsOnBoundary = numberOnBoundary;
}

template<unsigned SPACE_DIM>
void ExtendedCellProperty<SPACE_DIM>::SetChangeInCellDimensions(std::vector<double> changeCellDimensions)
{
    mChangeInCellDimensions = changeCellDimensions;
}

template<unsigned SPACE_DIM>
void ExtendedCellProperty<SPACE_DIM>::SetIncludeOdeInterpolationInCell(bool includeOdeInterpolationInCell)
{
    mIncludeOdeInterpolationInCell = includeOdeInterpolationInCell;
}

template<unsigned SPACE_DIM>
void ExtendedCellProperty<SPACE_DIM>::SetAverageInternalCellStates(bool averageInternalCellStates)
{
    mAverageInternalCellStates = averageInternalCellStates;
}


template<unsigned SPACE_DIM>
void ExtendedCellProperty<SPACE_DIM>::SetVectorOfBoundaryStateVariables(std::vector<std::vector<double>> boundaryStateVariables)
{
    mVectorOfExternalBoundaryStateVariables = boundaryStateVariables;
}

template<unsigned SPACE_DIM>
void ExtendedCellProperty<SPACE_DIM>::SetVectorOfBoundaryStateVariablesByIndex(unsigned index, std::vector<double> boundaryStateVariables)
{
    if(index<mVectorOfExternalBoundaryStateVariables.size())
    {
        mVectorOfExternalBoundaryStateVariables[index] = boundaryStateVariables;
    }
    else
    {
        std::cout<<"Error: ExtendedCellProperty<SPACE_DIM>::SetVectorOfBoundaryStateVariablesByIndex index out of bounds"<<std::endl;
    }
    
}

template<unsigned SPACE_DIM>
void ExtendedCellProperty<SPACE_DIM>::SetVectorOfInternalBoundaryStateVariables(std::vector<std::vector<double>> boundaryStateVariables)
{
    mVectorOfInternalBoundaryStateVariables = boundaryStateVariables;
}

template<unsigned SPACE_DIM>
void ExtendedCellProperty<SPACE_DIM>::SetVectorOfInternalBoundaryStateVariablesByIndex(unsigned index, std::vector<double> boundaryStateVariables)
{
    if(index<mVectorOfInternalBoundaryStateVariables.size())
    {
        mVectorOfInternalBoundaryStateVariables[index] = boundaryStateVariables;
    }
    else
    {
        std::cout<<"Error: ExtendedCellProperty<SPACE_DIM>::SetVectorOfInternalBoundaryStateVariablesByIndex index out of bounds"<<std::endl;
    }
    
}

template<unsigned SPACE_DIM>
void ExtendedCellProperty<SPACE_DIM>::SetVectorOfBoundaryLocations(std::vector<ChastePoint<SPACE_DIM>> vecPoint)
{
    mVectorOfBoundaryLocations = vecPoint;
}

template<unsigned SPACE_DIM>
void ExtendedCellProperty<SPACE_DIM>::SetVectorOfBoundaryLocationsByIndex(unsigned index, ChastePoint<SPACE_DIM> point)
{
    if(index<mVectorOfBoundaryLocations.size())
    {
        mVectorOfBoundaryLocations[index] = point;
    }
    else
    {
        std::cout<<"Error: ExtendedCellProperty<SPACE_DIM>::SetVectorOfBoundaryLocationsByIndex index out of bounds"<<std::endl;
    }
}

template<unsigned SPACE_DIM>
void ExtendedCellProperty<SPACE_DIM>::SetNextTimestepConcentrationVector(std::vector<double> concentrations)
{
    mNextTimestepConcentrationVector = concentrations;
}

// Get methods

template<unsigned SPACE_DIM>
double ExtendedCellProperty<SPACE_DIM>::GetCellRadius()
{
    return mCellDimensions[0];
}

template<unsigned SPACE_DIM>
std::vector<double> ExtendedCellProperty<SPACE_DIM>::GetCellDimensions()
{
    return mCellDimensions;
}

template<unsigned SPACE_DIM>
std::vector<double> ExtendedCellProperty<SPACE_DIM>::GetMinimalCellDimensions()
{
    return mMinimalDimensions;
}

template<unsigned SPACE_DIM>
double ExtendedCellProperty<SPACE_DIM>::GetCellVolume()
{
    return mCellVolume;
}

template<unsigned SPACE_DIM>
bool ExtendedCellProperty<SPACE_DIM>::GetIsCellVolumeDefined()
{
    return mIsCellVolumeDefined;
}

template<unsigned SPACE_DIM>
bool ExtendedCellProperty<SPACE_DIM>::GetIsMeshElementVolumeDefined()
{
    return mIsMeshElementVolumeDefined;
}

template<unsigned SPACE_DIM>
bool ExtendedCellProperty<SPACE_DIM>::GetIsCellDynamic()
{
    return mIsCellDynamic;
}

template<unsigned SPACE_DIM>
std::vector<double> ExtendedCellProperty<SPACE_DIM>::GetMeshDomainScale()
{
    return mMeshDomainScale;
}

template<unsigned SPACE_DIM>
unsigned ExtendedCellProperty<SPACE_DIM>::GetTotalNumberMeshVoxels()
{
    return mTotalNumberMeshVoxels;
}

template<unsigned SPACE_DIM>
unsigned ExtendedCellProperty<SPACE_DIM>::GetNumberMeshVoxelsInCell()
{
    return mNumberMeshVoxelsInCell;
}

template<unsigned SPACE_DIM>
unsigned ExtendedCellProperty<SPACE_DIM>::GetNumberMeshVoxelsOnBoundary()
{
    return mNumberMeshVoxelsOnBoundary;
}

template<unsigned SPACE_DIM>
std::vector<double> ExtendedCellProperty<SPACE_DIM>::GetChangeInCellDimensions()
{
    return mChangeInCellDimensions;
}

template<unsigned SPACE_DIM>
bool ExtendedCellProperty<SPACE_DIM>::GetIncludeOdeInterpolationInCell()
{
    return mIncludeOdeInterpolationInCell;
}

template<unsigned SPACE_DIM>
bool ExtendedCellProperty<SPACE_DIM>::GetAverageInternalCellStates()
{
    return mAverageInternalCellStates;
}

template<unsigned SPACE_DIM>
std::vector<std::vector<double>> ExtendedCellProperty<SPACE_DIM>::GetVectorOfBoundaryStateVariables()
{
    return mVectorOfExternalBoundaryStateVariables;
}

template<unsigned SPACE_DIM>
std::vector<double> ExtendedCellProperty<SPACE_DIM>::GetVectorOfBoundaryStateVariablesByLocationIndex(unsigned index)
{
    if(index < mVectorOfExternalBoundaryStateVariables.size())
    {
        return mVectorOfExternalBoundaryStateVariables[index];
    }
    
    return std::vector<double>();
}

template<unsigned SPACE_DIM>
std::vector<std::vector<double>> ExtendedCellProperty<SPACE_DIM>::GetVectorOfInternalBoundaryStateVariables()
{
    return mVectorOfInternalBoundaryStateVariables;
}

template<unsigned SPACE_DIM>
std::vector<double> ExtendedCellProperty<SPACE_DIM>::GetVectorOfInternalBoundaryStateVariablesByLocationIndex(unsigned index)
{
    if(index < mVectorOfInternalBoundaryStateVariables.size())
    {
        return mVectorOfInternalBoundaryStateVariables[index];
    }
    
    return std::vector<double>();
}

template<unsigned SPACE_DIM>
std::vector<ChastePoint<SPACE_DIM>> ExtendedCellProperty<SPACE_DIM>::GetVectorOfBoundaryLocations()
{
    return mVectorOfBoundaryLocations;
}

template<unsigned SPACE_DIM>
ChastePoint<SPACE_DIM> ExtendedCellProperty<SPACE_DIM>::GetVectorOfBoundaryLocationsByIndex(unsigned index)
{
    if(index < mVectorOfBoundaryLocations.size())
    {
        return mVectorOfBoundaryLocations[index];
    }
    
    return ChastePoint<SPACE_DIM>();
}


template<unsigned SPACE_DIM>
std::vector<double> ExtendedCellProperty<SPACE_DIM>::GetNextTimestepConcentrationVector()
{
    return mNextTimestepConcentrationVector;
}

#endif