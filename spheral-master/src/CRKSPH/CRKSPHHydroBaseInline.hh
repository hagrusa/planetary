namespace Spheral {

//------------------------------------------------------------------------------
// Choose whether we want to sum for mass density, or integrate the continuity
// equation.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
MassDensityType
CRKSPHHydroBase<Dimension>::densityUpdate() const {
  return mDensityUpdate;
}

template<typename Dimension>
inline
void
CRKSPHHydroBase<Dimension>::
densityUpdate(const MassDensityType type) {
  mDensityUpdate = type;
}

//------------------------------------------------------------------------------
// Choose how we want to update the H tensor.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
HEvolutionType
CRKSPHHydroBase<Dimension>::HEvolution() const {
  return mHEvolution;
}

template<typename Dimension>
inline
void
CRKSPHHydroBase<Dimension>::
HEvolution(const HEvolutionType type) {
  mHEvolution = type;
}

//------------------------------------------------------------------------------
// Choose which order we use for the CRK Corrections
//------------------------------------------------------------------------------
template<typename Dimension>
inline
CRKOrder
CRKSPHHydroBase<Dimension>::correctionOrder() const {
  return mCorrectionOrder;
}

template<typename Dimension>
inline
void
CRKSPHHydroBase<Dimension>::
correctionOrder(const CRKOrder order) {
  mCorrectionOrder = order;
}

//------------------------------------------------------------------------------
// Choose which volume weighting to use.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
CRKVolumeType
CRKSPHHydroBase<Dimension>::volumeType() const {
  return mVolumeType;
}

template<typename Dimension>
inline
void
CRKSPHHydroBase<Dimension>::
volumeType(const CRKVolumeType x) {
  mVolumeType = x;
}

//------------------------------------------------------------------------------
// Access the flag determining if we're using the compatible energy evolution 
// algorithm.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
CRKSPHHydroBase<Dimension>::compatibleEnergyEvolution() const {
  return mCompatibleEnergyEvolution;
}

template<typename Dimension>
inline
void
CRKSPHHydroBase<Dimension>::compatibleEnergyEvolution(const bool val) {
  mCompatibleEnergyEvolution = val;
}

//------------------------------------------------------------------------------
// Access the flag determining if we're evolving total or specific energy
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
CRKSPHHydroBase<Dimension>::evolveTotalEnergy() const {
  return mEvolveTotalEnergy;
}

template<typename Dimension>
inline
void
CRKSPHHydroBase<Dimension>::evolveTotalEnergy(const bool val) {
  mEvolveTotalEnergy = val;
}

//------------------------------------------------------------------------------
// Access the flag determining if we're using the XSPH algorithm.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
CRKSPHHydroBase<Dimension>::XSPH() const {
  return mXSPH;
}

template<typename Dimension>
inline
void
CRKSPHHydroBase<Dimension>::XSPH(const bool val) {
  mXSPH = val;
}

//------------------------------------------------------------------------------
// The object defining how smoothing scales are evolved.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const SmoothingScaleBase<Dimension>&
CRKSPHHydroBase<Dimension>::
smoothingScaleMethod() const {
  return mSmoothingScaleMethod;
}

//------------------------------------------------------------------------------
// Fraction of the centroidal filtering to apply.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
CRKSPHHydroBase<Dimension>::filter() const {
  return mfilter;
}

template<typename Dimension>
inline
void
CRKSPHHydroBase<Dimension>::filter(const double val) {
  // VERIFY(val >= 0.0 and val <= 1.0);
  mfilter = val;
}

//------------------------------------------------------------------------------
// Parameter to determine the magnitude of the tensile small scale correction.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
CRKSPHHydroBase<Dimension>::
epsilonTensile() const {
  return mEpsTensile;
}

template<typename Dimension>
void
CRKSPHHydroBase<Dimension>::
epsilonTensile(const typename Dimension::Scalar val) {
  mEpsTensile = val;
}

//------------------------------------------------------------------------------
// Parameter to set the exponent used in the tensile small scale correction.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
CRKSPHHydroBase<Dimension>::
nTensile() const {
  return mnTensile;
}

template<typename Dimension>
inline
void
CRKSPHHydroBase<Dimension>::
nTensile(const typename Dimension::Scalar val) {
  mnTensile = val;
}
    
//------------------------------------------------------------------------------
// Per node correction limits.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
CRKSPHHydroBase<Dimension>::
correctionMin() const {
  return mCorrectionMin;
}

template<typename Dimension>
inline
void
CRKSPHHydroBase<Dimension>::
correctionMin(const double val) {
  mCorrectionMin = val;
}
    
template<typename Dimension>
inline
double
CRKSPHHydroBase<Dimension>::
correctionMax() const {
  return mCorrectionMax;
}

template<typename Dimension>
inline
void
CRKSPHHydroBase<Dimension>::
correctionMax(const double val) {
  mCorrectionMax = val;
}
    
//------------------------------------------------------------------------------
// The internal state field lists.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const FieldList<Dimension, int>&
CRKSPHHydroBase<Dimension>::
timeStepMask() const {
  return mTimeStepMask;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
CRKSPHHydroBase<Dimension>::
pressure() const {
  return mPressure;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
CRKSPHHydroBase<Dimension>::
soundSpeed() const {
  return mSoundSpeed;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
CRKSPHHydroBase<Dimension>::
specificThermalEnergy0() const {
  return mSpecificThermalEnergy0;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
CRKSPHHydroBase<Dimension>::
entropy() const {
  return mEntropy;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::SymTensor>&
CRKSPHHydroBase<Dimension>::
Hideal() const {
  return mHideal;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
CRKSPHHydroBase<Dimension>::
maxViscousPressure() const {
  return mMaxViscousPressure;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
CRKSPHHydroBase<Dimension>::
effectiveViscousPressure() const {
  return mEffViscousPressure;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
CRKSPHHydroBase<Dimension>::
viscousWork() const {
  return mViscousWork;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
CRKSPHHydroBase<Dimension>::
weightedNeighborSum() const {
  return mWeightedNeighborSum;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::SymTensor>&
CRKSPHHydroBase<Dimension>::
massSecondMoment() const {
  return mMassSecondMoment;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
CRKSPHHydroBase<Dimension>::
volume() const {
  return mVolume;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Vector>&
CRKSPHHydroBase<Dimension>::
massDensityGradient() const {
  return mMassDensityGradient;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Vector>&
CRKSPHHydroBase<Dimension>::
XSPHDeltaV() const {
  return mXSPHDeltaV;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Vector>&
CRKSPHHydroBase<Dimension>::
DxDt() const {
  return mDxDt;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Vector>&
CRKSPHHydroBase<Dimension>::
DvDt() const {
  return mDvDt;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
CRKSPHHydroBase<Dimension>::
DmassDensityDt() const {
  return mDmassDensityDt;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
CRKSPHHydroBase<Dimension>::
DspecificThermalEnergyDt() const {
  return mDspecificThermalEnergyDt;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::SymTensor>&
CRKSPHHydroBase<Dimension>::
DHDt() const {
  return mDHDt;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Tensor>&
CRKSPHHydroBase<Dimension>::
DvDx() const {
  return mDvDx;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Tensor>&
CRKSPHHydroBase<Dimension>::
internalDvDx() const {
  return mInternalDvDx;
}

template<typename Dimension>
inline
const FieldList<Dimension, std::vector<typename Dimension::Vector> >&
CRKSPHHydroBase<Dimension>::
pairAccelerations() const {
  return mPairAccelerations;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Vector>&
CRKSPHHydroBase<Dimension>::
deltaCentroid() const {
  return mDeltaCentroid;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
CRKSPHHydroBase<Dimension>::
A() const {
  return mA;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Vector>&
CRKSPHHydroBase<Dimension>::
B() const {
  return mB;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Tensor>&
CRKSPHHydroBase<Dimension>::
C() const {
  return mC;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Vector>&
CRKSPHHydroBase<Dimension>::
gradA() const {
  return mGradA;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Tensor>&
CRKSPHHydroBase<Dimension>::
gradB() const {
  return mGradB;
}
  
template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::ThirdRankTensor>&
CRKSPHHydroBase<Dimension>::
gradC() const {
  return mGradC;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
CRKSPHHydroBase<Dimension>::
m0() const {
  return mM0;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Vector>&
CRKSPHHydroBase<Dimension>::
m1() const {
  return mM1;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::SymTensor>&
CRKSPHHydroBase<Dimension>::
m2() const {
  return mM2;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::ThirdRankTensor>&
CRKSPHHydroBase<Dimension>::
m3() const {
  return mM3;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::FourthRankTensor>&
CRKSPHHydroBase<Dimension>::
m4() const {
  return mM4;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Vector>&
CRKSPHHydroBase<Dimension>::
gradm0() const {
  return mGradm0;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Tensor>&
CRKSPHHydroBase<Dimension>::
gradm1() const {
  return mGradm1;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::ThirdRankTensor>&
CRKSPHHydroBase<Dimension>::
gradm2() const {
  return mGradm2;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::FourthRankTensor>&
CRKSPHHydroBase<Dimension>::
gradm3() const {
  return mGradm3;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::FifthRankTensor>&
CRKSPHHydroBase<Dimension>::
gradm4() const {
  return mGradm4;
}

template<typename Dimension>
inline
const FieldList<Dimension, int>&
CRKSPHHydroBase<Dimension>::
surfacePoint() const {
  return mSurfacePoint;
}

template<typename Dimension>
inline
const FieldList<Dimension, int>&
CRKSPHHydroBase<Dimension>::
voidPoint() const {
  return mVoidPoint;
}

template<typename Dimension>
inline
const FieldList<Dimension, std::vector<typename Dimension::Vector>>&
CRKSPHHydroBase<Dimension>::
etaVoidPoints() const {
  return mEtaVoidPoints;
}

template<typename Dimension>
inline
const CRKSPHVoidBoundary<Dimension>&
CRKSPHHydroBase<Dimension>::
voidBoundary() const {
  return mVoidBoundary;
}

}
