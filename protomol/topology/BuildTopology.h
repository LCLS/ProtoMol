#ifndef BUILD_TOPOLOGY_H
#define BUILD_TOPOLOGY_H

#include <protomol/topology/CoulombSCPISMParameterTable.h>
#include <protomol/topology/ExclusionType.h>

namespace ProtoMol {
  class GenericTopology;
  class PAR;
  class PSF;
struct CoulombSCPISMParameterTable;

  void buildExclusionTable(GenericTopology *topo,
                           const ExclusionType &exclusionType);
  void buildTopology(GenericTopology *topo, const PSF &psf,
    const PAR &par, bool dihedralMultPSF, CoulombSCPISMParameterTable *SCPISMParameters);
  void buildMoleculeTable(GenericTopology *topo);
}

#endif // BUILD_TOPOLOGY_H
