// chemical ode includes
#include "AbstractReactionSystem.hpp"
#include "MassActionReaction.hpp"
#include "AbstractChemicalOdeSystem.hpp"
#include "AbstractChemicalOdeForCoupledPdeSystem.hpp"


// transport includes
#include "AbstractTransportReactionSystem.hpp"
#include "AbstractTransportReaction.hpp"
#include "AbstractTransportOutReaction.hpp"
#include "MassActionTransportReaction.hpp"
//#include "AbstractTransportReactionSystemFromFile.hpp"
#include "AbstractReversibleTransportReaction.hpp"
#include "AbstractTransportOdeSystem.hpp"

// membrane includes
#include "AbstractMembraneReactionSystem.hpp"
#include "AbstractMembraneReaction.hpp"
#include "MassActionCoupledMembraneReaction.hpp"
//#include "AbstractMembraneReactionSystemFromFile.hpp"
#include "AbstractReversibleMembraneReaction.hpp"
#include "AbstractMembraneOdeSystem.hpp"

// SRN includes
#include "SchnackenbergSrnModel.hpp"
#include "ChemicalSrnModel.hpp"
#include "NullSrnModel.hpp"

// cell cycle includes
#include "SimpleChemicalThresholdCellCycleModel.hpp"

// tracking includes
#include "ChemicalTrackingModifier.hpp"

#include "ChemicalCell.hpp"
#include "ComplexCell.hpp"

#include "AbstractBoxDomainPdeSystemModifier.hpp"
#include "AbstractPdeSystemModifier.hpp"

// added classes
//#include "AbstractDomainField_templated.hpp"
//#include "ParabolicBoxDomainPdeSystemModifier.hpp"
//#include "ChemicalDomainField_templated.hpp"
#include "ParabolicBoxDomainPdeSystemModifier.hpp"
#include "ChemicalDomainFieldForCellCoupling.hpp"

#include "ChemicalDomainField_templated.hpp"
//#include "ChemicalDomainFieldForCellCoupling.hpp"

#include "InhomogenousCoupledPdeOdeSolver_templated.hpp"
#include "InhomogenousParabolicPdeForCoupledOdeSystem_templated.hpp"
#include "InhomogenousParabolicPdeForCoupledOdeSystemInhibitedDiffusion_templated.hpp"


// cell property includes
#include "TransportCellProperty.hpp"
#include "ChemicalCellProperty.hpp"
#include "ExtendedCellProperty.hpp"
#include "MembraneCellProperty.hpp"


// chemical includes
#include "AbstractChemical.hpp"
#include "AbstractDiffusiveChemical.hpp"
#include "AbstractChemistry.hpp"
#include "AbstractDiffusiveChemistry.hpp" 

// reaction-diffusion includes
#include "AbstractReaction.hpp"
#include "AbstractReversibleReaction.hpp"
#include "AbstractReactionSystem.hpp"
#include "MassActionReaction.hpp"
#include "AbstractReactionSystemFromFile.hpp"
#include "ReactionTypeDatabase.hpp"


// custom pdeOde includes
#include "PdeSchnackenbergCoupledPdeOdeSystem.hpp"
#include "OdeSchnackenbergCoupledPdeOdeSystem.hpp"
#include "PdeConsumerProducer.hpp"
#include "OdeConsumerProducer.hpp"


#include "AbstractChemicalOdeSystem.hpp"
#include "AbstractChemicalOdeForCoupledPdeSystem.hpp"


#include "SchnackenbergCoupledPdeSystem.hpp"

// writers
#include "CellProliferativeTypesWriter_ChemChaste.hpp"

// cell boundary conditions
//#include "FixedCellBoundaryConditions.hpp"