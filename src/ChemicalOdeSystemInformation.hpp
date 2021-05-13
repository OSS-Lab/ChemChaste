#ifndef CHEMICALODESYSTEMINFORMATION_HPP_
#define CHEMICALODESYSTEMINFORMATION_HPP_

#include <boost/shared_ptr.hpp>
#include "AbstractOdeSystemInformation.hpp"
#include "AbstractReactionSystem.hpp"


template<class ODE_SYSTEM>
class ChemicalOdeSystemInformation : public AbstractOdeSystemInformation
{
private:

    static boost::shared_ptr<ChemicalOdeSystemInformation<ODE_SYSTEM> > mpInstance;

    AbstractReactionSystem* mp_reaction_system;

public:

    /**
     * Default constructor; calls Initialise.
     *
     * Designed to be used as follows by ODE system classes in their constructors:
     *   mpSystemInfo.reset(new CellwiseOdeSystemInformation<CLASS>);
     */
    
    static boost::shared_ptr<ChemicalOdeSystemInformation<ODE_SYSTEM> > Instance();

    ChemicalOdeSystemInformation(AbstractReactionSystem*);

protected:

    ChemicalOdeSystemInformation(const ChemicalOdeSystemInformation<ODE_SYSTEM>&);

    ChemicalOdeSystemInformation& operator= (const ChemicalOdeSystemInformation<ODE_SYSTEM>&);

    void Initialise();

    AbstractReactionSystem* GetReactionSystem();

    void SetReactionSystem(AbstractReactionSystem*);

};

template<class ODE_SYSTEM>
ChemicalOdeSystemInformation<ODE_SYSTEM>::ChemicalOdeSystemInformation(AbstractReactionSystem* p_reaction_system) : mp_reaction_system(p_reaction_system)
{
    assert(mpInstance == nullptr);
    
    ChemicalOdeSystemInformation<ODE_SYSTEM>::Initialise();
}

template<class ODE_SYSTEM>
void ChemicalOdeSystemInformation<ODE_SYSTEM>::Initialise()
{
}


template<class ODE_SYSTEM>
boost::shared_ptr<ChemicalOdeSystemInformation<ODE_SYSTEM> > ChemicalOdeSystemInformation<ODE_SYSTEM>::Instance()
{
    if (!mpInstance)
    {
        mpInstance.reset(new ChemicalOdeSystemInformation<ODE_SYSTEM>);
        
        mpInstance->Initialise();
    }
    return mpInstance;
}

template<class ODE_SYSTEM>
AbstractReactionSystem* ChemicalOdeSystemInformation<ODE_SYSTEM>::GetReactionSystem()
{
    return mp_reaction_system;
}

template<class ODE_SYSTEM>
void ChemicalOdeSystemInformation<ODE_SYSTEM>::SetReactionSystem(AbstractReactionSystem* p_reaction_system)
{
    mp_reaction_system  = p_reaction_system;
}


/**
 * Definition of the instance static member.
 */
template<class ODE_SYSTEM>
boost::shared_ptr<ChemicalOdeSystemInformation<ODE_SYSTEM> > ChemicalOdeSystemInformation<ODE_SYSTEM>::mpInstance;

#endif