//------------------------------------------------------------------------------
// Construct a flaw distribution for a set of nodes according to a Weibull 
// (power-law) distribution.
//------------------------------------------------------------------------------
#ifndef __Spheral_weibullFlawDistribution__
#define __Spheral_weibullFlawDistribution__

#include <vector>

// Foward declarations.
namespace Spheral {
  template<typename Dimension> class FluidNodeList;
  template<typename Dimension, typename DataType> class Field;
}

namespace Spheral {

//------------------------------------------------------------------------------
// Implements the Benz-Asphaug algorithm, starting with a minimum based
// on the volume of the simulation.
//------------------------------------------------------------------------------
template<typename Dimension>
Field<Dimension, std::vector<double> >
weibullFlawDistributionBenzAsphaug(double volume,
                                   const double volumeStretchFactor,
                                   const unsigned seed,
                                   const double kWeibull,
                                   const double mWeibull,
                                   const FluidNodeList<Dimension>& nodeList,
                                   const int minFlawsPerNode,
                                   const int minTotalFlaws);

//------------------------------------------------------------------------------
// Implements the Owen algorithm, stochastically seeding flaws with a maximum
// value per node chosen based on the volume of the node.
//------------------------------------------------------------------------------
template<typename Dimension>
Field<Dimension, std::vector<double> >
weibullFlawDistributionOwen(const unsigned seed,
                            const double kWeibull,
                            const double mWeibull,
                            const FluidNodeList<Dimension>& nodeList,
                            const int minFlawsPerNode,
                            const double volumeMultiplier = 1.0);

}

#endif
