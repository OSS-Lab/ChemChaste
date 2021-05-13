#ifndef TRANSPORTODESYSTEMINFORMATION_HPP_
#define TRANSPORTODESYSTEMINFORMATION_HPP_

#include <boost/shared_ptr.hpp>
#include "AbstractOdeSystemInformation.hpp"
#include "AbstractTransportReactionSystem.hpp"


template<class ODE_SYSTEM>
class TransportOdeSystemInformation : public AbstractOdeSystemInformation
{
private:

    static boost::shared_ptr<TransportOdeSystemInformation<ODE_SYSTEM> > mpInstance;

    AbstractTransportReactionSystem* mp_reaction_system;

public:

    /**
     * Default constructor; calls Initialise.
     *
     * Designed to be used as follows by ODE system classes in their constructors:
     *   mpSystemInfo.reset(new CellwiseOdeSystemInformation<CLASS>);
     */
    
    static boost::shared_ptr<TransportOdeSystemInformation<ODE_SYSTEM> > Instance();

    TransportOdeSystemInformation(AbstractTransportReactionSystem*);

protected:

    TransportOdeSystemInformation(const TransportOdeSystemInformation<ODE_SYSTEM>&);

    TransportOdeSystemInformation& operator= (const TransportOdeSystemInformation<ODE_SYSTEM>&);

    void Initialise();

    AbstractTransportReactionSystem* GetReactionSystem();

    void SetReactionSystem(AbstractTransportReactionSystem*);

};

template<class ODE_SYSTEM>
TransportOdeSystemInformation<ODE_SYSTEM>::TransportOdeSystemInformation(AbstractTransportReactionSystem* p_reaction_system) : mp_reaction_system(p_reaction_system)
{
    assert(mpInstance == nullptr);
    
    TransportOdeSystemInformation<ODE_SYSTEM>::Initialise();
}

template<class ODE_SYSTEM>
void TransportOdeSystemInformation<ODE_SYSTEM>::Initialise()
{
}


template<class ODE_SYSTEM>
boost::shared_ptr<TransportOdeSystemInformation<ODE_SYSTEM> > TransportOdeSystemInformation<ODE_SYSTEM>::Instance()
{
    if (!mpInstance)
    {
        mpInstance.reset(new TransportOdeSystemInformation<ODE_SYSTEM>);
        
        mpInstance->Initialise();
    }
    return mpInstance;
}

template<class ODE_SYSTEM>
AbstractTransportReactionSystem* TransportOdeSystemInformation<ODE_SYSTEM>::GetReactionSystem()
{
    return mp_reaction_system;
}

template<class ODE_SYSTEM>
void TransportOdeSystemInformation<ODE_SYSTEM>::SetReactionSystem(AbstractTransportReactionSystem* p_reaction_system)
{
    mp_reaction_system  = p_reaction_system;
}


/**
 * Definition of the instance static member.
 */
template<class ODE_SYSTEM>
boost::shared_ptr<TransportOdeSystemInformation<ODE_SYSTEM> > TransportOdeSystemInformation<ODE_SYSTEM>::mpInstance;

#endif