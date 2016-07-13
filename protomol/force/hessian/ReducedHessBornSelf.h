/*  -*- c++ -*-  */
#ifndef REDUCEDHESSBORNSELF_H
#define REDUCEDHESSBORNSELF_H

#include <protomol/topology/ExclusionTable.h>
#include <protomol/type/Matrix3By3.h>

#include "protomol/type/Real.h"
#include "protomol/type/Vector3D.h"

namespace ProtoMol {
  class GenericTopology;

  class ReducedHessBornSelf {
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // New methods of class ReducedHessBornSelf
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    Matrix3By3 operator()(Real distSquared,
                          const Vector3D &diff,
                          const GenericTopology *topo,
                          int atom1, int atom2, 
                          int bornSwitch, Real dielecConst,
                          ExclusionClass excl) const;


  };
}
#endif
