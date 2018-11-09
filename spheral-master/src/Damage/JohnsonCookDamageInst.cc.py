text = """
//------------------------------------------------------------------------------
// Instantiations.
//------------------------------------------------------------------------------
#include "JohnsonCookDamage.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class JohnsonCookDamage<Dim< %(ndim)s > >;
}
"""
