#ifndef BUILD_TOPOLOGY_FROM_XML_H
#define BUILD_TOPOLOGY_FROM_XML_H

#include <protomol/topology/ExclusionType.h>
#include <protomol/type/Vector3DBlock.h>
#include <protomol/type/PDB.h>

#include <vector>
#include <string>
#include <set>

namespace ProtoMol {
  class GenericTopology;

  //build topo from XML file
  void buildTopologyFromXML(GenericTopology *topo,
                                Vector3DBlock &pos, Vector3DBlock &vel,
                                const string &fname,
                                std::vector<PDB::Atom>);

  //parse functions
  //bool parse_iparams(function &func, void *ft, void *ip, ostringstream &os);

  //atomic radius from lookup, choose from orig. and Bowman sets
  //double atom_radius( std::string atom_type, int set );

}

#endif // BUILD_TOPOLOGY_FROM_XML_H
