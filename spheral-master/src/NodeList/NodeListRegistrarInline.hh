#include "NodeList/NodeList.hh"
#include "Field/FieldBase.hh"

#include <algorithm>
#include <sstream>

namespace Spheral {

//------------------------------------------------------------------------------
// Get the instance.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
NodeListRegistrar<Dimension>&
NodeListRegistrar<Dimension>::
instance() {
  if (mInstancePtr == 0) mInstancePtr = new NodeListRegistrar;
  CHECK(mInstancePtr != 0);
  return *mInstancePtr;
}

//------------------------------------------------------------------------------
// The number of NodeLists currently registered.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
NodeListRegistrar<Dimension>::
numNodeLists() const {
  return mNodeLists.size();
}

template<typename Dimension>
inline
int
NodeListRegistrar<Dimension>::
numFluidNodeLists() const {
  return mFluidNodeLists.size();
}

//------------------------------------------------------------------------------
// Iterators over the registered NodeLists.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename NodeListRegistrar<Dimension>::iterator
NodeListRegistrar<Dimension>::
begin() {
  return mNodeLists.begin();
}

template<typename Dimension>
inline
typename NodeListRegistrar<Dimension>::iterator
NodeListRegistrar<Dimension>::
end() {
  return mNodeLists.end();
}

template<typename Dimension>
inline
typename NodeListRegistrar<Dimension>::const_iterator
NodeListRegistrar<Dimension>::
begin() const {
  return mNodeLists.begin();
}

template<typename Dimension>
inline
typename NodeListRegistrar<Dimension>::const_iterator
NodeListRegistrar<Dimension>::
end() const {
  return mNodeLists.end();
}

// FluidNodeLists.
template<typename Dimension>
inline
typename NodeListRegistrar<Dimension>::fluid_iterator
NodeListRegistrar<Dimension>::
fluidBegin() {
  return mFluidNodeLists.begin();
}

template<typename Dimension>
inline
typename NodeListRegistrar<Dimension>::fluid_iterator
NodeListRegistrar<Dimension>::
fluidEnd() {
  return mFluidNodeLists.end();
}

template<typename Dimension>
inline
typename NodeListRegistrar<Dimension>::const_fluid_iterator
NodeListRegistrar<Dimension>::
fluidBegin() const {
  return mFluidNodeLists.begin();
}

template<typename Dimension>
inline
typename NodeListRegistrar<Dimension>::const_fluid_iterator
NodeListRegistrar<Dimension>::
fluidEnd() const {
  return mFluidNodeLists.end();
}

//------------------------------------------------------------------------------
// Generic method to find the proper place in a sequence to insert a Field
// or NodeList.
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename IteratorType, typename ThingyType>
inline
IteratorType
NodeListRegistrar<Dimension>::
findInsertionPoint(const ThingyType& thingy,
                   const IteratorType begin,
                   const IteratorType end) const {

  // If the input iterator sequence is empty, then the answer is easy!
  const int containerSize = std::distance(begin, end);
  if (containerSize == 0) return end;

  // Get the sequence of NodeLists represented by the input.
  std::vector<NodeList<Dimension>*> nodeListPtrs;
  nodeListPtrs.reserve(containerSize);
  for (IteratorType itr = begin; itr != end; ++itr) {
    NodeList<Dimension>* nodeListPtr = getNodeListPtr(*itr);
    CHECK(itr == begin or (nodeListPtr->name() > getNodeListPtr(*(itr - 1))->name()));
    nodeListPtrs.push_back(nodeListPtr);
  }
  CHECK(nodeListPtrs.size() == containerSize);

  // Now we can find where the specified thingy should go.
  NodeList<Dimension>* nodeListPtr = getNodeListPtr(thingy);
  typename std::vector<NodeList<Dimension>*>::iterator orderItr = std::upper_bound(nodeListPtrs.begin(),
                                                                                              nodeListPtrs.end(),
                                                                                              nodeListPtr,
                                                                                              NodeListComparator());
  const int displacement = std::distance(nodeListPtrs.begin(), orderItr);
  CHECK(displacement >= 0 && displacement <= containerSize);
  IteratorType result = begin + displacement;
  ENSURE(result >= begin && result <= end);
  return result;
}

}
