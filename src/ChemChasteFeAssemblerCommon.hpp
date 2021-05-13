
#ifndef ABSTRACTFEASSEMBLERCOMMON_HPP_
#define ABSTRACTFEASSEMBLERCOMMON_HPP_

#include "AbstractFeAssemblerInterface.hpp"
#include "ReplicatableVector.hpp"
#include "DistributedVector.hpp"
#include "HeartEventHandler.hpp"
#include "LinearBasisFunction.hpp"
#include "PetscTools.hpp"
#include "AbstractTetrahedralMesh.hpp"

/**
 * Enumeration for defining how much interpolation (onto quadrature points) is
 * required by the concrete class.
 *
 * CARDIAC: only interpolates the first component of the unknown (ie the voltage)
 * NORMAL: interpolates the position X and all components of the unknown u
 * NONLINEAR: interpolates X, u and grad(u). Also computes the gradient of the
 *   basis functions when assembling vectors.
 */
typedef enum InterpolationLevel_
{
    CARDIAC = 0,
    NORMAL,
    NONLINEAR
} InterpolationLevel;

/**
 *   A base class for AbstractFeVolumeIntegralAssembler (the main abstract assembler class), AbstractSurfaceFeObjectAssembler, and
 *   AbstractCableFeObjectAssembler.
 *
 *   The base class of this, AbstractFeAssemblerInterface, defines the interface for these assembler classes. This class
 *   just defines a few pde-folder-specific (ie not continuum-mechanics-related) extra methods.
 *
 *   See AbstractFeVolumeIntegralAssembler documentation for info on these assembler classes.
 */
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM, bool CAN_ASSEMBLE_VECTOR, bool CAN_ASSEMBLE_MATRIX, InterpolationLevel INTERPOLATION_LEVEL>
class AbstractFeAssemblerCommon : public AbstractFeAssemblerInterface<CAN_ASSEMBLE_VECTOR,CAN_ASSEMBLE_MATRIX>
{
protected:
    /**
     * If the matrix or vector will be dependent on a current solution, say,
     * this is where that information is put.
     */
    ReplicatableVector mCurrentSolutionOrGuessReplicated;

    /**
     * @return an entry from the solution vector.
     *
     * @param nodeIndex node index
     * @param indexOfUnknown index of unknown
     */
    virtual double GetCurrentSolutionOrGuessValue(unsigned nodeIndex, unsigned indexOfUnknown)
    {
        return mCurrentSolutionOrGuessReplicated[ PROBLEM_DIM*nodeIndex + indexOfUnknown];
    }

    /**
     * The concrete subclass can overload this and IncrementInterpolatedQuantities()
     * if there are some quantities which need to be computed at each Gauss point.
     * They are called in AssembleOnElement().
     */
    virtual void ResetInterpolatedQuantities()
    {}

    /**
     * The concrete subclass can overload this and ResetInterpolatedQuantities()
     * if there are some quantities which need to be computed at each Gauss point.
     * They are called in AssembleOnElement().
     *
     * Note that this method is called over assembly of elements, surface elements and cables.
     *
     * @param phiI
     * @param pNode pointer to a node
     */
    virtual void IncrementInterpolatedQuantities(double phiI, const Node<SPACE_DIM>* pNode)
    {}

    /**
     * The concrete subclass can overload this and ResetInterpolatedQuantities()
     * if there are some gradient dependent quantities which need to be computed at each Gauss point.
     * A matrix of all the basis function gradients at the quad point is passed for efficiency reasons.
     * To access the gradient vector use of the current basis function use rGradPhi(:, phi_index);
     * They are called in AssembleOnElement().
     *
     * Note that this method is ONLY called during assembly of elements. NOT during assembly of surface elements or cables.
     *
     * Further, it is ONLY called in the cases where rGradPhi has been computed.  Currently these cases are
     *     - When mAssembleMatrix is set
     *  or - When INTERPOLATION_LEVEL==NONLINEAR
     * If there are use-cases, then allow other cases where interpolated gradients are needed in righthand-side assembly
     * (see #2075)
     *
     * @param rGradPhi A matrix containing the gradient of all the basis functions at this Gauss point.
     * @param phiIndex The index of the current basis function in the rGradPhi matrix.
     * @param pNode pointer to the node associated with the current basis function
     */
    virtual void IncrementInterpolatedGradientQuantities(const c_matrix<double, SPACE_DIM, ELEMENT_DIM+1>& rGradPhi, unsigned phiIndex, const Node<SPACE_DIM>* pNode)
    {}
public:

    /**
     * Constructor.
     */
    AbstractFeAssemblerCommon();

    /**
     * Set a current solution vector that will be used in AssembleOnElement and can passed
     * up to ComputeMatrixTerm() or ComputeVectorTerm().
     *
     * @param currentSolution Current solution vector.
     */
    void SetCurrentSolution(Vec currentSolution);

    /**
     * Destructor.
     */
    virtual ~AbstractFeAssemblerCommon()
    {
    }

    ChastePoint<SPACE_DIM> mX;

    c_vector<double,PROBLEM_DIM> mU;

    ChastePoint<SPACE_DIM> GetPosition();

    c_vector<double,PROBLEM_DIM> GetStateVariable();


    void SetPosition(ChastePoint<SPACE_DIM>);

    void SetStateVariable(c_vector<double,PROBLEM_DIM>);
};

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM, bool CAN_ASSEMBLE_VECTOR, bool CAN_ASSEMBLE_MATRIX, InterpolationLevel INTERPOLATION_LEVEL>
AbstractFeAssemblerCommon<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM, CAN_ASSEMBLE_VECTOR, CAN_ASSEMBLE_MATRIX, INTERPOLATION_LEVEL>::AbstractFeAssemblerCommon()
    : AbstractFeAssemblerInterface<CAN_ASSEMBLE_VECTOR,CAN_ASSEMBLE_MATRIX>()
{
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM, bool CAN_ASSEMBLE_VECTOR, bool CAN_ASSEMBLE_MATRIX, InterpolationLevel INTERPOLATION_LEVEL>
void AbstractFeAssemblerCommon<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM, CAN_ASSEMBLE_VECTOR, CAN_ASSEMBLE_MATRIX, INTERPOLATION_LEVEL>::SetCurrentSolution(Vec currentSolution)
{
    assert(currentSolution != nullptr);

    // Replicate the current solution and store so can be used in AssembleOnElement
    HeartEventHandler::BeginEvent(HeartEventHandler::COMMUNICATION);
    mCurrentSolutionOrGuessReplicated.ReplicatePetscVector(currentSolution);
    HeartEventHandler::EndEvent(HeartEventHandler::COMMUNICATION);

    // The AssembleOnElement type methods will determine if a current solution or
    // current guess exists by looking at the size of the replicated vector, so
    // check the size is zero if there isn't a current solution.
    assert(mCurrentSolutionOrGuessReplicated.GetSize() > 0);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
ChastePoint<SPACE_DIM> AbstractFeAssemblerCommon<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::GetPosition()
{
    return mX;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
c_vector<double,PROBLEM_DIM> AbstractFeAssemblerCommon<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::GetStateVariable()
{
    return mU;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractFeAssemblerCommon<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetPosition(ChastePoint<SPACE_DIM> X)
{
    mX = X;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractFeAssemblerCommon<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetStateVariable(c_vector<double,PROBLEM_DIM> U)
{
    mU = U;
}


#endif 
