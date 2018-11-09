text = """
//---------------------------------Spheral++----------------------------------//
// PiGaussianKernel -- The gaussian interpolation kernel.
//
// Created by JMO, Wed Dec  1 14:38:51 PST 1999
//----------------------------------------------------------------------------//
#include "Kernel.hh"
#include "PiGaussianKernel.hh"

#include <math.h>

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
namespace Spheral {
  template class PiGaussianKernel< Dim< %(ndim)s >  >;
}
"""
