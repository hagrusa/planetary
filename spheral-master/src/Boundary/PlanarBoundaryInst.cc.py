text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "PlanarBoundary.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class PlanarBoundary< Dim< %(ndim)s > >;
}
"""
