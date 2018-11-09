text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "LinearAcceleration.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class LinearAcceleration< Dim< %(ndim)s > >;
}
"""
