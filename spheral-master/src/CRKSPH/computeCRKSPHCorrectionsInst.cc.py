text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "computeCRKSPHCorrections.cc"

namespace Spheral {
template void computeCRKSPHCorrections(const FieldList<Dim<%(ndim)s> , Dim<%(ndim)s> ::Scalar>& m0,
                                       const FieldList<Dim<%(ndim)s> , Dim<%(ndim)s> ::Vector>& m1,
                                       const FieldList<Dim<%(ndim)s> , Dim<%(ndim)s> ::SymTensor>& m2,
                                       const FieldList<Dim<%(ndim)s> , Dim<%(ndim)s> ::ThirdRankTensor>& m3,
                                       const FieldList<Dim<%(ndim)s> , Dim<%(ndim)s> ::FourthRankTensor>& m4,
                                       const FieldList<Dim<%(ndim)s> , Dim<%(ndim)s> ::Vector>& gradm0,
                                       const FieldList<Dim<%(ndim)s> , Dim<%(ndim)s> ::Tensor>& gradm1,
                                       const FieldList<Dim<%(ndim)s> , Dim<%(ndim)s> ::ThirdRankTensor>& gradm2,
                                       const FieldList<Dim<%(ndim)s> , Dim<%(ndim)s> ::FourthRankTensor>& gradm3,
                                       const FieldList<Dim<%(ndim)s> , Dim<%(ndim)s> ::FifthRankTensor>& gradm4,
                                       const FieldList<Dim<%(ndim)s> , Dim<%(ndim)s> ::SymTensor>& H,
                                       const CRKOrder correctionOrder,
                                       FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>&,
                                       FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>&,
                                       FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor>&,
                                       FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>&,
                                       FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor>&,
                                       FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::ThirdRankTensor>&);

template void computeZerothCRKSPHCorrections(const FieldList<Dim<%(ndim)s> , Dim<%(ndim)s> ::Scalar>& m0,
                                             const FieldList<Dim<%(ndim)s> , Dim<%(ndim)s> ::Vector>& gradm0,
                                             FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>&,
                                             FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>&);

template void computeLinearCRKSPHCorrections(const FieldList<Dim<%(ndim)s> , Dim<%(ndim)s> ::Scalar>& m0,
                                             const FieldList<Dim<%(ndim)s> , Dim<%(ndim)s> ::Vector>& m1,
                                             const FieldList<Dim<%(ndim)s> , Dim<%(ndim)s> ::SymTensor>& m2,
                                             const FieldList<Dim<%(ndim)s> , Dim<%(ndim)s> ::Vector>& gradm0,
                                             const FieldList<Dim<%(ndim)s> , Dim<%(ndim)s> ::Tensor>& gradm1,
                                             const FieldList<Dim<%(ndim)s> , Dim<%(ndim)s> ::ThirdRankTensor>& gradm2,
                                             const FieldList<Dim<%(ndim)s> , Dim<%(ndim)s> ::SymTensor>& H,
                                             FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>&,
                                             FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>&,
                                             FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>&,
                                             FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor>&);

template void computeQuadraticCRKSPHCorrections(const FieldList<Dim<%(ndim)s> , Dim<%(ndim)s> ::Scalar>& m0,
                                                const FieldList<Dim<%(ndim)s> , Dim<%(ndim)s> ::Vector>& m1,
                                                const FieldList<Dim<%(ndim)s> , Dim<%(ndim)s> ::SymTensor>& m2,
                                                const FieldList<Dim<%(ndim)s> , Dim<%(ndim)s> ::ThirdRankTensor>& m3,
                                                const FieldList<Dim<%(ndim)s> , Dim<%(ndim)s> ::FourthRankTensor>& m4,
                                                const FieldList<Dim<%(ndim)s> , Dim<%(ndim)s> ::Vector>& gradm0,
                                                const FieldList<Dim<%(ndim)s> , Dim<%(ndim)s> ::Tensor>& gradm1,
                                                const FieldList<Dim<%(ndim)s> , Dim<%(ndim)s> ::ThirdRankTensor>& gradm2,
                                                const FieldList<Dim<%(ndim)s> , Dim<%(ndim)s> ::FourthRankTensor>& gradm3,
                                                const FieldList<Dim<%(ndim)s> , Dim<%(ndim)s> ::FifthRankTensor>& gradm4,
                                                const FieldList<Dim<%(ndim)s> , Dim<%(ndim)s> ::SymTensor>& H,
                                                FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>&,
                                                FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>&,
                                                FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor>&,
                                                FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>&,
                                                FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor>&,
                                                FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::ThirdRankTensor>&);
}

"""
