text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "nodeOrdering.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template FieldList<Dim< %(ndim)s >, int> nodeOrdering<Dim< %(ndim)s > >(const FieldList<Dim< %(ndim)s >, KeyTraits::Key>&);
}
"""
