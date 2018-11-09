text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "IsothermalEquationOfState.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class IsothermalEquationOfState< Dim< %(ndim)s > >;
}
"""
