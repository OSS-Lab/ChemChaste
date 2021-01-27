// general headers
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <tuple>
#include <cmath>

// chaste headers
#include "AbstractCellBasedTestSuite.hpp"
#include "RandomNumberGenerator.hpp"

#include "HoneycombMeshGenerator.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "ConstBoundaryCondition.hpp"
#include "OffLatticeSimulation.hpp"

#include "ApoptoticCellProperty.hpp"

#include "HoneycombVertexMeshGenerator.hpp"
#include "ReplicatableVector.hpp"

#include "WildTypeCellMutationState.hpp"
#include "StemCellProliferativeType.hpp"
#include "UniformCellCycleModel.hpp"
#include "NoCellCycleModel.hpp"

#include "EulerIvpOdeSolver.hpp"

#include "ConstBoundaryCondition.hpp"
#include "OutputFileHandler.hpp"

#include "TrianglesMeshReader.hpp"

// cell simulation
#include "NodeBasedCellPopulation.hpp"

// cell simulation writers
#include "VoronoiDataWriter.hpp"
#include "CellDataItemWriter.hpp"
#include "CellLocationIndexWriter.hpp"
#include "CellIdWriter.hpp"
#include "CellAgesWriter.hpp"
#include "CellAncestorWriter.hpp"
#include "CellAppliedForceWriter.hpp"
#include "CellAnalyticsWriter.hpp"


// chemChaste headers
#include "BoundaryConditionsContainer_extended.hpp"


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




#include "InhomogenousParabolicPdeForCoupledOdeSystem_templated.hpp"

#include "InhomogenousCoupledPdeOdeSolver_templated.hpp"

#include "InhomogenousParabolicPdeForCoupledOdeSystemInhibitedDiffusion_templated.hpp"


// cell property includes
#include "TransportCellProperty.hpp"
#include "ChemicalCellProperty.hpp"
#include "ExtendedCellProperty.hpp"
#include "MembraneCellProperty.hpp"
#include "CellAnalyticsProperty.hpp"

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

// Fisher equation
#include "InhomogenousFisherPde.hpp"
#include "InhomogenousFisherDiffusiveInhibitionPde.hpp"


#include "AbstractChemicalOdeSystem.hpp"
#include "AbstractChemicalOdeForCoupledPdeSystem.hpp"


#include "ChemicalCellFromFile.hpp"
#include "ComplexCellFromFile.hpp"


#include "AbstractPdeSystemModifier.hpp"


#include "AbstractBoxDomainPdeSystemModifier.hpp"



#include "ParabolicBoxDomainPdeSystemModifier.hpp"

#include "ChemicalDomainFieldForCellCoupling.hpp"
#include "ChemicalDomainField_templated.hpp"

