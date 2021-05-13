#ifndef MEMBRANEODESYSTEMINFORMATION_HPP_
#define MEMBRANEODESYSTEMINFORMATION_HPP_

#include <boost/shared_ptr.hpp>
#include "AbstractOdeSystemInformation.hpp"
#include "AbstractMembraneReactionSystem.hpp"


template<class ODE_SYSTEM>
class MembraneOdeSystemInformation : public AbstractOdeSystemInformation
{
private:

    static boost::shared_ptr<MembraneOdeSystemInformation<ODE_SYSTEM> > mpInstance;

    AbstractMembraneReactionSystem* mp_reaction_system;

public:

    /**
     * Default constructor; calls Initialise.
     *
     * Designed to be used as follows by ODE system classes in their constructors:
     *   mpSystemInfo.reset(new CellwiseOdeSystemInformation<CLASS>);
     */
    
    static boost::shared_ptr<MembraneOdeSystemInformation<ODE_SYSTEM> > Instance();

    MembraneOdeSystemInformation(AbstractMembraneReactionSystem*);

protected:

    MembraneOdeSystemInformation(const MembraneOdeSystemInformation<ODE_SYSTEM>&);

    MembraneOdeSystemInformation& operator= (const MembraneOdeSystemInformation<ODE_SYSTEM>&);

    void Initialise();

    AbstractMembraneReactionSystem* GetReactionSystem();

    void SetReactionSystem(AbstractMembraneReactionSystem*);

};

template<class ODE_SYSTEM>
MembraneOdeSystemInformation<ODE_SYSTEM>::MembraneOdeSystemInformation(AbstractMembraneReactionSystem* p_reaction_system) : mp_reaction_system(p_reaction_system)
{
    assert(mpInstance == nullptr);
    
    MembraneOdeSystemInformation<ODE_SYSTEM>::Initialise();
}

template<class ODE_SYSTEM>
void MembraneOdeSystemInformation<ODE_SYSTEM>::Initialise()
{
}


template<class ODE_SYSTEM>
boost::shared_ptr<MembraneOdeSystemInformation<ODE_SYSTEM> > MembraneOdeSystemInformation<ODE_SYSTEM>::Instance()
{
    if (!mpInstance)
    {
        mpInstance.reset(new MembraneOdeSystemInformation<ODE_SYSTEM>);
        
        mpInstance->Initialise();
    }
    return mpInstance;
}

template<class ODE_SYSTEM>
AbstractMembraneReactionSystem* MembraneOdeSystemInformation<ODE_SYSTEM>::GetReactionSystem()
{
    return mp_reaction_system;
}

template<class ODE_SYSTEM>
void MembraneOdeSystemInformation<ODE_SYSTEM>::SetReactionSystem(AbstractMembraneReactionSystem* p_reaction_system)
{
    mp_reaction_system  = p_reaction_system;
}


/**
 * Definition of the instance static member.
 */
template<class ODE_SYSTEM>
boost::shared_ptr<MembraneOdeSystemInformation<ODE_SYSTEM> > MembraneOdeSystemInformation<ODE_SYSTEM>::mpInstance;

#endif