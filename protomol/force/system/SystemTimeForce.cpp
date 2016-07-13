#include <protomol/force/system/SystemTimeForce.h>
#include <protomol/type/Vector3DBlock.h>

#include "protomol/force/Force.h"
#include "protomol/force/TimeForce.h"
#include "protomol/force/system/SystemForce.h"

namespace ProtoMol {
class GenericTopology;
class ScalarStructure;
}  // namespace ProtoMol

using namespace std;
using namespace ProtoMol;
//____ SystemTimeForce

SystemTimeForce::SystemTimeForce(Force *actualForce) :
  TimeForce(actualForce)
{}

void SystemTimeForce::evaluate(const GenericTopology *topo,
                               const Vector3DBlock *positions,
                               Vector3DBlock *forces,
                               ScalarStructure *energies) {
  preprocess(positions->size());
  (dynamic_cast<SystemForce *>(myActualForce))->evaluate(topo,
    positions,
    forces,
    energies);
  postprocess(topo, forces, energies);
}

void SystemTimeForce::parallelEvaluate(const GenericTopology *topo,
                                       const Vector3DBlock *positions,
                                       Vector3DBlock *forces,
                                       ScalarStructure *energies) {
  preprocess(positions->size());
  (dynamic_cast<SystemForce *>(myActualForce))->parallelEvaluate(topo,
    positions,
    forces,
    energies);
  postprocess(topo, forces, energies);
}

