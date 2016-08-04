#include <protomol/topology/BuildTopologyFromXML.h>

#include <protomol/topology/BuildTopology.h>

#include <protomol/base/Exception.h>
#include <protomol/base/SystemUtilities.h>
#include <protomol/type/String.h>
#include <protomol/topology/TopologyUtilities.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/topology/LennardJonesParameterTable.h>
#include <vector>

#include <sstream>
#include <iomanip>

// GROMACS headers
#if defined(HAVE_GROMACS)
extern "C" {
#include <gromacs/tpxio.h>
}
#endif

using namespace std;
using namespace ProtoMol;
using namespace ProtoMol::Report;

// use GROMACS exclusions?
#define GROMACSEXCL

// Switch for which Van der Waal radius table to use for GB.
// 0 - Amber default
// 1 - Greg Bowman's modified
#define RADIUS_TABLE 0

// fudge for NETBEANS highlighting
// #define HAVE_GROMACS

// main build topology
void ProtoMol::buildTopologyFromXML(GenericTopology *topo, Vector3DBlock &pos,
                                    Vector3DBlock &vel, const string &fname,
                                    const string &posfname,
                                    std::vector<PDB::Atom> atoms) {
  // define versions of TPR file
  // Version 4.5.3 has tpx_version=73 and includes gb_radius in the tpr file
  //enum {GB_RADII_IN_TPR = 73};

#if defined(HAVE_GROMACS)
  report << "XML input: Atoms length " << atoms.size() << endr;
  
  // Check if XML exists
  if (!ProtoMol::SystemUtilities::exists(fname)) report << error << "XML Missing: " << fname << endr;
  
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // First, generate the array of atomtypes
  // Each time a new atom comes up, we need to check if it is
  // already in the vector....
  // NOTE:  this may take a while for large systems; however, it will cut
  // down on the size of the atomTypes vector, and therefore, the amount
  // access time in the back end.
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  topo->atoms.clear();
  topo->atomTypes.clear();
  topo->bonds.clear();
  topo->angles.clear();
  topo->dihedrals.clear();
  topo->impropers.clear();
  
  // Ryckert-Belleman
  topo->rb_dihedrals.clear();
  
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Get the atom types
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //create set for atom types
  std::set<std::string> atomtypes;
  
  //and residues
  std::set<std::string> residues;
  
  //loop over input atom data
  for(int i=0; i<atoms.size(); i++){
    report << "Element " << atoms[i].elementName << endr;
    atomtypes.insert(atoms[i].elementName);
    residues.insert(atoms[i].residueName);
  }
  
  report << "there are " << atomtypes.size() << " atom types" << endr;
  report << "there are " << residues.size() << " residue types" << endr;
  
  // resize atomtypes
  topo->atomTypes.resize(atomtypes.size());
  
  // save size for tests
  const unsigned atypesize = atomtypes.size();
  
  /*// get atomtypes
  t_atomtypes *atomtypes = &(mtop.atomtypes);
  
  // resize atomtypes
  topo->atomTypes.resize(atomtypes->nr);
  
  
  // loop and print data
  for (unsigned i = 0; i < atypesize; i++) {
    report << debug(810) << "atomtype[" << i << "]={radius="
    << atomtypes->radius[i]
    << ", volume=" << atomtypes->vol[i]
    << ", gb_radius=" << atomtypes->gb_radius[i]
    << ", surftens=" << atomtypes->surftens[i]
    << ", atomnumber=" << atomtypes->atomnumber[i]
    << ", S_hct=" << atomtypes->S_hct[i] << ")}" << endr;
  }*/
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Get the atoms
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  // used to get atom names to look up van der waal radii for GB
  //vector<string> atomTypeNames;
  
  // loop over all atoms in the Gromacs topology
  // molecule type topology data
  //for (int mt = 0; mt < mtop.nmoltype; mt++) {
  //  gmx_moltype_t *molt = &mtop.moltype[mt];
  //  t_atoms *atoms = &(molt->atoms);
  //  t_atom *atom = atoms->atom;
    
  // loop over atoms in molecule
  for (int i = 0; i < atoms.size(); i++) {
    // get gromacs data
    
    // atom type, test valid
    unsigned atype = atypesize;
    set<string>::iterator it;
    if ((it=std::find(atomtypes.begin(), atomtypes.end(), atoms[i].elementName)) != atomtypes.end())
    {
      // Element in set.
      //report << "Index " << distance(atomtypes.begin(), it) << endr;
      atype = distance(atomtypes.begin(), it);
    }
    
    if (atype >= atypesize) THROW("Undefined atom type.");
    
    // point to types
    AtomType *tempatomtype = &(topo->atomTypes[atype]);
    
    // create temporary atom
    Atom tempatom;
    
    // Update atom type if 'name' not initialized?
    if (!tempatomtype->name.length()) {
      // type not set
      // generate type name from first atom char and index
      stringstream ss;
      ss << atype;
      string str = string(atoms[i].elementName).substr((size_t)0, (size_t)1);
      
      tempatomtype->name = str + ss.str();
      tempatomtype->mass = 1.0;//atom[i].m; ####TODO get mass from XML
      tempatomtype->charge = 0.001;//atom[i].q; ####TODO get charge from XML
      // ####just take first char for now.
      tempatomtype->symbolName = str;
      
      // get radius, test file version
      // if (file_version >= GB_RADII_IN_TPR) {
      // version 4.5?
      //    tempatomtype->vdwR = atomtypes->gb_radius[atype];
      // } else {
      // }
      
    } else ;
    // ####TODO check new data is consistent? i.e. check new atom mass
    // is the same.
    
    // report types
    report << debug(810) << "Atom type " << tempatomtype->name << ", "
    << tempatomtype->mass << ", "
    << tempatomtype->charge << ", " << tempatomtype->symbolName
    << ", " << tempatomtype->vdwR << endr;
    
    
    // First, we need to find the index. (an integer corresponding
    // to the type of the atom
    tempatom.name = string(atoms[i].elementName);
    tempatom.type = atype;
    tempatom.residue_name = string(atoms[i].residueName);
    tempatom.residue_seq = atoms[i].residueNum;
    // Now, the scaled charge.  This is straightforward.
    tempatom.scaledCharge = tempatomtype->charge * Constant::SQRTCOULOMBCONSTANT;
    tempatom.scaledMass = tempatomtype->mass;
    
    //####Need this? atomTypeNames.push_back(string(*(atoms->atomtype[i])));
    
    // report atoms
    report << debug(810) << "Atom " << tempatom.name << ", " <<
    tempatom.type << ", " << tempatom.residue_name <<
    ", " << tempatom.residue_seq << ", " << tempatom.scaledCharge <<
    ", " << tempatom.scaledMass << endr;
    
    // Now we need the size of the group for heavy atom ordering
    // We need to parse the name for any H's then any numbers following
    // First, if the atom is an H then this is 0
    if (tempatom.name == "H") tempatom.hvyAtom = 0;
    else {
      // Otherwise, we need to parse..
      // Initialize to 1
      tempatom.hvyAtom = 1;
      for (unsigned pos = 0; pos < tempatom.name.size(); ++pos)
        if (tempatom.name[pos] == 'H') {
          string number = "";
          while (isdigit(tempatom.name[++pos]))
            number += tempatom.name[pos];
          
          if (number == "") number = "1"; // never entered loop, default is 1
          tempatom.hvyAtom += atoi(number.c_str());
        }
    }
    // C/C++ starts at 0, where PSF/PDB at 1
    tempatom.atomNum = atoms[i].elementNum - 1;//i; // ####atom->number - 1;*/
    // Also the molecule - using residue sequence for now
    topo->atoms.push_back(tempatom);
  }
  //}
  
  // calculate the # of degrees of freedom, if there are any bond constraints
  // they will be subtracted later by ModifierShake
  topo->degreesOfFreedom = 3 * topo->atoms.size() - 3;
  
  report << plain << "D.O.F. = " << topo->degreesOfFreedom << endr;


  
#endif
    //end here until viable force field produced
  report << error << "Force field incomplete for XML." << endr;
}

