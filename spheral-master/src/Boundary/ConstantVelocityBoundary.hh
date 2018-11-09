//---------------------------------Spheral++----------------------------------//
// ConstantVelocityBoundary -- A boundary condition to enforce a constant 
// velocity on a given set of nodes.
//
// This boundary is very specialized -- it explicitly works on only one 
// NodeList.
//
// Created by JMO, Mon Sep  9 23:39:39 PDT 2002
//----------------------------------------------------------------------------//
#ifndef ConstantVelocityBoundary_HH
#define ConstantVelocityBoundary_HH

#include "Boundary.hh"
#include "DataOutput/registerWithRestart.hh"

#include <vector>

namespace Spheral {

template<typename Dimension> class NodeList;
template<typename Dimension, typename DataType> class Field;
template<typename Dimension, typename DataType> class FieldList;
template<typename Dimension> class DataBase;
class FileIO;

template<typename Dimension>
class ConstantVelocityBoundary: public Boundary<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::ThirdRankTensor ThirdRankTensor;

  // Constructors and destructors.
  ConstantVelocityBoundary(const NodeList<Dimension>& nodeList,
                           const std::vector<int>& nodeIndices);
  virtual ~ConstantVelocityBoundary();

  //**********************************************************************
  // All Boundary conditions must provide the following methods:
  // Use the given NodeList's neighbor object to select the ghost nodes.
  virtual void setGhostNodes(NodeList<Dimension>& nodeList);

  // For the computed set of ghost nodes, set the positions and H's.
  virtual void updateGhostNodes(NodeList<Dimension>& nodeList);

  // Apply the boundary condition to the given Field.
  virtual void applyGhostBoundary(Field<Dimension, int>& field) const;
  virtual void applyGhostBoundary(Field<Dimension, Scalar>& field) const;
  virtual void applyGhostBoundary(Field<Dimension, Vector>& field) const;
  virtual void applyGhostBoundary(Field<Dimension, Tensor>& field) const;
  virtual void applyGhostBoundary(Field<Dimension, SymTensor>& field) const;
  virtual void applyGhostBoundary(Field<Dimension, ThirdRankTensor>& field) const;

  // Find any internal nodes that are in violation of this Boundary.
  virtual void setViolationNodes(NodeList<Dimension>& nodeList);

  // For the computed set of nodes in violation of the boundary, bring them
  // back into compliance (for the positions and H's.)
  virtual void updateViolationNodes(NodeList<Dimension>& nodeList);

  // Apply the boundary condition to the violation node values in the given Field.
  virtual void enforceBoundary(Field<Dimension, int>& field) const;
  virtual void enforceBoundary(Field<Dimension, Scalar>& field) const;
  virtual void enforceBoundary(Field<Dimension, Vector>& field) const;
  virtual void enforceBoundary(Field<Dimension, Tensor>& field) const;
  virtual void enforceBoundary(Field<Dimension, SymTensor>& field) const;
  virtual void enforceBoundary(Field<Dimension, ThirdRankTensor>& field) const;
  //**********************************************************************

  // Allow read only access to the node indices and their forced velocities.
  const NodeList<Dimension>& nodeList() const;
  std::vector<int> nodeIndices() const;
  std::vector<Vector> velocityCondition() const;

  // Determine if the boundary is in a "valid", ready to use state.
  virtual bool valid() const;

  //******************************************************************************
  // Restart methods.
  virtual std::string label() const { return "ConstantVelocityBoundary"; }
  virtual void dumpState(FileIO& file, const std::string& pathName) const;
  virtual void restoreState(const FileIO& file, const std::string& pathName);
  //******************************************************************************

protected:
  //--------------------------- Protected Interface ---------------------------//

private:
  //--------------------------- Private Interface ---------------------------//
  const NodeList<Dimension>* mNodeListPtr;
  Field<Dimension, int> mNodes;
  Field<Dimension, Vector> mVelocity;

  // The restart registration.
  RestartRegistrationType mRestart;
};

}

#include "ConstantVelocityBoundaryInline.hh"

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class ConstantVelocityBoundary;
}

#endif
