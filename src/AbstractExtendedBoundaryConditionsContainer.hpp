#ifndef ABSTRACTEXTENDEDBOUNDARYCONDITIONSCONTAINER_HPP_
#define ABSTRACTEXTENDEDBOUNDARYCONDITIONSCONTAINER_HPP_

#include <map>
#include <set>
#include <boost/utility.hpp>
#include "AbstractBoundaryCondition.hpp"
#include "Node.hpp"
#include "PetscTools.hpp"


//#include "AbstractExtendedBoundaryConditionsContainerImplementation.hpp"


/**
 * Helper struct storing an operator for computing whether one node
 * has a lower index than another.
 */

// redefinition of struct in AbstrctBoundaryConditionsContainer as expected
template<unsigned SPACE_DIM>
struct LessThanNode
{
    /**
     * Less-then node index comparison operator.
     * @return true if n1<n2
     * @param n1 pointer to a node
     * @param n2 pointer to a node
     */
    bool operator()(const Node<SPACE_DIM> * const &n1, const Node<SPACE_DIM> * const &n2) const
    {
        return (n1->GetIndex() < n2->GetIndex() );
    }
};




/**
 * Abstract boundary conditions container.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
class AbstractExtendedBoundaryConditionsContainer : boost::noncopyable
{
protected:

    /** To save typing */
    typedef typename std::map< const Node<SPACE_DIM> *, const AbstractBoundaryCondition<SPACE_DIM>*, LessThanNode<SPACE_DIM> >
        DirichletMapType;
    DirichletMapType* mpDirichletMap[PROBLEM_DIM]; /**< List (map) of Dirichlet boundary conditions */

    /** To save typing */
    typedef typename std::map< const Node<SPACE_DIM> *, const AbstractBoundaryCondition<SPACE_DIM>*, LessThanNode<SPACE_DIM> >::const_iterator
        DirichletIteratorType;
    DirichletIteratorType mDirichIterator; /**< Internal iterator over Dirichlet boundary conditions */

    /** @return true if there are any Dirichlet BCs anywhere on the mesh*/
    bool mHasDirichletBCs;

    /** @return true if we calculated mHasDirichletBCs. */
    bool mCheckedAndCommunicatedIfDirichletBcs;

    /** @return true if we need to delete BCs in destructor. */
    bool mDeleteConditions;

    /**
     * Delete the list of Dirichlet boundary conditions.
     *
     * @note This should stay as a protected method to avoid it being called with default arguments and causing seg faults
     *  (requires careful bookkeeping when calling this method).
     * @param alreadyDeletedConditions  This is a set of BCs that have already been deleted that we should avoid trying
     *  to delete inside this method. (defaults to empty = delete everything)
     */
    void DeleteDirichletBoundaryConditions(std::set<const AbstractBoundaryCondition<SPACE_DIM>*> alreadyDeletedConditions
                                            = std::set<const AbstractBoundaryCondition<SPACE_DIM>*>());

public:

    /**
     * Constructor allocates memory for the Dirichlet boundary conditions lists.
     *
     * @param deleteConditions whether to delete BCs in destructor (defaults to true)
     */
    AbstractExtendedBoundaryConditionsContainer(bool deleteConditions=true);

    /**
     * Destructor.
     */
    ~AbstractExtendedBoundaryConditionsContainer();

    /**
     * @return whether any Dirichlet conditions are defined (for ANY of the unknowns, on ANY of the processes).
     * Must be called collectively. The first time this is called, the result is communicated to all processes
     * and then cached locally (the bool mHasDirichletBCs). If this needs recalculating
     * mCheckedAndCommunicatedIfDirichletBcs must be reset to zero.
     */
    bool HasDirichletBoundaryConditions();

    /**
     * @return value of Dirichlet boundary condition at specified node.
     *
     * This is unlikely to be needed by the user, the methods ApplyDirichletToLinearProblem or
     * ApplyDirichletToNonlinearProblem can be called instead to apply all Dirichlet boundary conditions
     * at the same time.
     *
     * @param pBoundaryNode pointer to a boundary node
     * @param indexOfUnknown index of the unknown for which to obtain the value of the boundary condition (defaults to 0)
     */
    double GetDirichletBCValue(const Node<SPACE_DIM>* pBoundaryNode, unsigned indexOfUnknown = 0);

    /**
     * @return true if there is a Dirichlet boundary condition defined on the given node.
     *
     * @param pNode pointer to a node
     * @param indexOfUnknown index of the unknown for which to obtain the value of the boundary condition (defaults to 0)
     */
    bool HasDirichletBoundaryCondition(const Node<SPACE_DIM>* pNode, unsigned indexOfUnknown = 0);

    /**
     * When Dirichlet boundary conditions are likely to be added on one or more processes then we should call this
     * method collectively in order to ensure that all processes do a collective communication on the next call
     * to HasDirichletBoundaryConditions()
     */
    void ResetDirichletCommunication()
    {
        mCheckedAndCommunicatedIfDirichletBcs = false;
    }
};


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
AbstractExtendedBoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::AbstractExtendedBoundaryConditionsContainer(bool deleteConditions)
    : mHasDirichletBCs(false),
      mCheckedAndCommunicatedIfDirichletBcs(false),
      mDeleteConditions(deleteConditions)
{
    for (unsigned index_of_unknown=0; index_of_unknown<PROBLEM_DIM; index_of_unknown++)
    {
        mpDirichletMap[index_of_unknown] = new std::map< const Node<SPACE_DIM> *, const AbstractBoundaryCondition<SPACE_DIM>*, LessThanNode<SPACE_DIM> >;
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
AbstractExtendedBoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::~AbstractExtendedBoundaryConditionsContainer()
{
    DeleteDirichletBoundaryConditions();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
bool AbstractExtendedBoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::HasDirichletBoundaryConditions()
{
    if (!mCheckedAndCommunicatedIfDirichletBcs)
    {
        bool i_have_dirichlet=false;
        for (unsigned i=0; i<PROBLEM_DIM; i++)
        {
            if (!mpDirichletMap[i]->empty())
            {
                i_have_dirichlet=true;
                break;
            }
        }
        mHasDirichletBCs = PetscTools::ReplicateBool(i_have_dirichlet);
        mCheckedAndCommunicatedIfDirichletBcs = true;
    }
    return mHasDirichletBCs;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractExtendedBoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::DeleteDirichletBoundaryConditions(std::set<const AbstractBoundaryCondition<SPACE_DIM>*> alreadyDeletedConditions)
{
    for (unsigned i=0; i<PROBLEM_DIM; i++)
    {
        if (mpDirichletMap[i])
        {
            mDirichIterator = mpDirichletMap[i]->begin();
            while (mDirichIterator != mpDirichletMap[i]->end() )
            {
                if (alreadyDeletedConditions.count(mDirichIterator->second) == 0)
                {
                    alreadyDeletedConditions.insert(mDirichIterator->second);
                    if (mDeleteConditions)
                    {
                        delete mDirichIterator->second;
                    }
                }
                mDirichIterator++;
            }

            delete(mpDirichletMap[i]);
            mpDirichletMap[i] = nullptr;
        }
    }

    // Recommunicate that Dirichlet BCs have changed (next time we ask)
    ResetDirichletCommunication();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
double AbstractExtendedBoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::GetDirichletBCValue(const Node<SPACE_DIM>* pBoundaryNode, unsigned indexOfUnknown)
{
    assert(indexOfUnknown < PROBLEM_DIM);
    //assert(pBoundaryNode->IsBoundaryNode());

    mDirichIterator = mpDirichletMap[indexOfUnknown]->find(pBoundaryNode);
    assert(mDirichIterator != mpDirichletMap[indexOfUnknown]->end());

    return mDirichIterator->second->GetValue(pBoundaryNode->GetPoint());
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
bool AbstractExtendedBoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::HasDirichletBoundaryCondition(const Node<SPACE_DIM>* pNode, unsigned indexOfUnknown)
{
    assert(indexOfUnknown < PROBLEM_DIM);

    this->mDirichIterator = this->mpDirichletMap[indexOfUnknown]->find(pNode);

    return (this->mDirichIterator != this->mpDirichletMap[indexOfUnknown]->end());
}



#endif