//---------------------------------Spheral++----------------------------------//
// StateBase -- The base class for State and StateDerivatives method, providing
// the common methods for registering and storing Fields, as well as serving
// up FieldLists based on these Fields.
//
// Important usage notes!
// Copying StateBase objects (either by copy constructor or assignment) always
// copies by reference!  In other words, even if you have a StateBase object
// that is using internal storage to hold it's Fields, copying or assigning 
// another StateBase object with this one will only create new references to the
// original's Fields.  You have to explicitly call "copyFields()" if you want
// to do the expensive copy.
// 
// Created by JMO, Wed Aug 25 22:23:35 2004
//----------------------------------------------------------------------------//
#ifndef __Spheral_StateBase_hh__
#define __Spheral_StateBase_hh__

#include <string>
#include <utility>
#include <memory>
#include <vector>
#include <map>
#include <list>
#include <set>

namespace Spheral {

// Forward declaration.
template<typename Dimension> class NodeList;
template<typename Dimension> class FieldBase;
template<typename Dimension> class FieldListBase;
template<typename Dimension, typename DataType> class Field;
template<typename Dimension, typename DataType> class FieldList;
template<typename Dimension> class ConnectivityMap;
template<typename Dimension> class Mesh;

template<typename Dimension>
class StateBase {

public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Vector3d Vector3d;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Spheral::ConnectivityMap<Dimension> ConnectivityMapType;
  typedef typename Spheral::Mesh<Dimension> MeshType;

  typedef std::shared_ptr<ConnectivityMapType> ConnectivityMapPtr;
  typedef std::shared_ptr<MeshType> MeshPtr;

  typedef std::string KeyType;
  typedef typename FieldBase<Dimension>::FieldName FieldName;

  // Constructors, destructor.
  StateBase();
  StateBase(const StateBase& rhs);
  virtual ~StateBase();

  // Assignment.
  StateBase& operator=(const StateBase& rhs);

  // Test if two StateBases have equivalent fields.
  virtual bool operator==(const StateBase& rhs) const;

  // Test if the specified Field or key is currently registered.
  bool registered(const KeyType& key) const;
  bool registered(const FieldBase<Dimension>& field) const;
  bool registered(const FieldListBase<Dimension>& fieldList) const;
  bool fieldNameRegistered(const FieldName& fieldName) const;

  // Enroll a Field.
  virtual void enroll(FieldBase<Dimension>& field);
  virtual void enroll(std::shared_ptr<FieldBase<Dimension>>& fieldPtr);

  // Enroll a FieldList.
  virtual void enroll(FieldListBase<Dimension>& fieldList);

  // Enroll an externally held Mesh.
  void enrollMesh(MeshPtr meshPtr);

  // Enroll a std::vector<Scalar>.
  void enroll(const FieldName& key, std::vector<Scalar>& vec);

  // Return the field for the given key.
  template<typename Value>
  Field<Dimension, Value>& field(const KeyType& key,
                                             const Value& dummy) const;

  // Return FieldLists constructed from all registered Fields with the given name.
  template<typename Value>
  FieldList<Dimension, Value> fields(const std::string& name, 
                                                 const Value& dummy) const;

  // Return all the fields of the given Value.
  template<typename Value>
  std::vector<Field<Dimension, Value>*> allFields(const Value& dummy) const;

  // Return the vector<Scalar>.
  std::vector<Scalar>& array(const FieldName& key);

  // Return the complete set of keys registered.
  std::vector<KeyType> keys() const;

  // Return the set of known field names (unencoded from our internal mangling
  // convention with the NodeList name).
  std::vector<FieldName> fieldKeys() const;

  // A state object can carry around a reference to a ConnectivityMap.
  void enrollConnectivityMap(ConnectivityMapPtr connectivityMapPtr);
  const ConnectivityMapType& connectivityMap() const;

  // We also provide support for registering a mesh (though only one per StateBase).
  bool meshRegistered() const;
  const MeshType& mesh() const;
  MeshType& mesh();

  // Set the Fields equal to those in another State object.
  void assign(const StateBase<Dimension>& rhs);

  // Force the StateBase to create new internally owned copies of all state.
  virtual void copyState();

  // Construct the lookup key for the given field.
  static KeyType key(const FieldBase<Dimension>& field);

  // Construct the lookup key for the given FieldList.
  static KeyType key(const FieldListBase<Dimension>& fieldList);

  // Encode our underlying convention for combining Field and NodeList names into a 
  // single Key.
  static KeyType buildFieldKey(const KeyType& fieldName, const KeyType& nodeListName);
  static void splitFieldKey(const KeyType& key, KeyType& fieldName, KeyType& nodeListName);

protected:
  //--------------------------- Protected Interface ---------------------------//
  typedef std::map<KeyType, FieldBase<Dimension>*> StorageType;
  typedef std::list<std::shared_ptr<FieldBase<Dimension> > > CacheType;
  typedef std::map<KeyType, std::vector<Scalar>*> VectorStorageType;
  typedef std::list<std::shared_ptr<std::vector<Scalar> > > VectorCacheType;

  // Protected data.
  StorageType mStorage;
  CacheType mCache;
  VectorStorageType mVectorStorage;
  VectorCacheType mVectorCache;
  std::set<const NodeList<Dimension>*> mNodeListPtrs;
  ConnectivityMapPtr mConnectivityMapPtr;
  MeshPtr mMeshPtr;
};

}

#include "StateBaseInline.hh"

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class StateBase;
}

#endif

