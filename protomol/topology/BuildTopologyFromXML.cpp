#include <protomol/topology/BuildTopologyFromXML.h>

#include <protomol/topology/BuildTopology.h>

#include <protomol/base/Exception.h>
#include <protomol/base/SystemUtilities.h>
#include <protomol/type/String.h>
#include <protomol/topology/TopologyUtilities.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/topology/LennardJonesParameterTable.h>
#include <protomol/tinyxml2/tinyxml2.h>

#include <vector>
#include <sstream>
#include <iomanip>

#include <iostream>
#include <fstream>

// GROMACS headers
#if defined(HAVE_GROMACS)
extern "C" {
#include <gromacs/tpxio.h>
}
#endif

using namespace std;
using namespace ProtoMol;
using namespace ProtoMol::Report;
using namespace tinyxml2;

//~~~~structs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//masses with index
typedef struct mass_index {
  // members
  int index;
  float mass;
  
  mass_index(int indx, float ms){
    index = indx;
    mass = ms;
  }
  
  bool operator==(const int& l) const
  {
    return l == index;
  }
} mass_index;

//charge/epsilon/sigma with index
typedef struct electrostatic_index {
  // members
  int index;
  float charge, epsilon, sigma;
  
  electrostatic_index(int indx, float ch, float ep, float si){
    index = indx;
    charge = ch; epsilon = ep; sigma = si;
  }
  
  bool operator==(const int& l) const
  {
    return l == index;
  }
} electrostatic_index;
//~~~~End structs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
  // Read the XML file
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  /*string line;
  ifstream myfile (fname.c_str());
  if (myfile.is_open())
  {
    while ( getline (myfile,line) )
    {
      cout << line << '\n';
    }
    myfile.close();
  }*/
  //create XML object and read file into it
  XMLDocument doc;
  doc.LoadFile( fname.c_str() );
  
  if(doc.ErrorID()) report << error << "XML File error opening " << fname << endr;
  
  //~~~~get particle mass information~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //create vector for masses
  std::vector<mass_index> masses;

  //load element containing 'particle'
  tinyxml2::XMLElement *levelElement = doc.FirstChildElement("forcefield")->FirstChildElement("particles");
  
  //loop through it
  for (tinyxml2::XMLElement* child = levelElement->FirstChildElement(); child != NULL; child = child->NextSiblingElement())
  {
    // do something with each child element
    report << debug(800) << "Mass of particle " << child->Attribute( "mass" ) << ", " << child->Attribute( "index" ) << endr;
    const float xmlmass = atof(child->Attribute( "mass" ));
    const int xmlindex = atoi(child->Attribute( "index" ));
    masses.push_back(mass_index(xmlindex, xmlmass));
  }
  
  //test no error
  if(doc.ErrorID()) report << error << "XML File parsing masses error!" << endr;
  
  //~~~~get particle electrostatic information~~~~~~~~~~~~~~~~~~~~~~~~~
  //create vector for electrostatics
  std::vector<electrostatic_index> electrostatics;
  
  //load element containing 'particle'
  tinyxml2::XMLElement *levelElemente = doc.FirstChildElement("forcefield");
  
  //loop through it
  for (tinyxml2::XMLElement* child = levelElemente->FirstChildElement(); child != NULL; child = child->NextSiblingElement()){
    
    //find force of type NonbondedForce
    if(strcmp(child->Name(), "force") == 0 && strcmp(child->Attribute( "type" ), "NonbondedForce") == 0){
      report << "XML name " << child->Name() << ", " << child->Attribute( "type" ) << endr;
      
      //find children of particles, and partical count
      tinyxml2::XMLElement *forcenonbonded = child->FirstChildElement("particles");
      
      const int fnbcount = atoi(forcenonbonded->Attribute( "count" ));
      
      //report << "Count " << forcenonbonded->Attribute( "count" ) << endr;
      
      //loop through it
      for (tinyxml2::XMLElement* fnbchild = forcenonbonded->FirstChildElement(); fnbchild != NULL; fnbchild = fnbchild->NextSiblingElement()){
        
        //report << "Electrostatics " << fnbcount << endr;
        
        //save data
        electrostatics.push_back(electrostatic_index(atoi(fnbchild->Attribute( "index" )), atof(fnbchild->Attribute( "charge" )), atof(fnbchild->Attribute( "epsilon" )), atof(fnbchild->Attribute( "sigma" ))));
      }
                                                                                            
      //test right number
      if(electrostatics.size() != fnbcount){
        report << error << "Number of electrostatics wrong " << electrostatics.size() << ". " << fnbcount << endr;
      }
                                                                                            
      //go as there is only one set
      break;
    }
  }
  
  //test no error
  if(doc.ErrorID()) report << error << "XML File parsing electrostatics error!" << endr;
  
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Get the atom types
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //create set for atom types
  std::set<std::string> atomtypes;
  
  //and residues
  std::set<std::string> residues;
  
  //loop over input atom data
  for(int i=0; i<atoms.size(); i++){
    //report << "Element " << atoms[i].elementName << endr;
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
      //report << "setting type " << i << endr;
      // type not set
      // generate type name from first atom char and index
      stringstream ss;
      ss << atype;
      string str = string(atoms[i].elementName).substr((size_t)0, (size_t)1);
      
      tempatomtype->name = str + ss.str();
      
      //get mass, if available
      vector<mass_index>::iterator itm;
      if ((itm=std::find(masses.begin(), masses.end(), atoms[i].elementNum - 1)) != masses.end())
      {
        // Element in vector.
        //report << "Mass " << (*itm).mass << endr;
        tempatomtype->mass = (*itm).mass;
      }else{
        tempatomtype->mass = 0.0;//atom[i].m; ####TODO get mass from XML
        THROW("Mass of atom undefined.");
      }
      
      //get electrostatics, if available
      vector<electrostatic_index>::iterator ite;
      if ((ite=std::find(electrostatics.begin(), electrostatics.end(), atoms[i].elementNum - 1)) != electrostatics.end())
      {
        // Element in vector.
        //report << "Electrostatics " << (*ite).charge << endr;
        tempatomtype->charge = (*ite).charge;
        tempatomtype->sigma = (*ite).sigma * Constant::NM_ANGSTROM;
        tempatomtype->epsilon = (*ite).epsilon * Constant::KJ_KCAL;
      }else{
        tempatomtype->charge = 0.000;//atom[i].q; ####TODO get charge from XML
        THROW("Electrostatics of atom undefined.");
      }
      //~~~~~~~~
      
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

  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Get the forces
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Find the parameters from TPR
  int ignoredBonds = 0;   // preset ignored bonds
  int ignoredAngles = 0;  // and angles

  // ~~~~Bonds~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //load element containing 'forcefield'
  tinyxml2::XMLElement *levelElementb = doc.FirstChildElement("forcefield");
  
  //loop through it to find force, HarmonicBondForce
  for (tinyxml2::XMLElement* child = levelElementb->FirstChildElement(); child != NULL; child = child->NextSiblingElement()){
    
    //find force of type HarmonicBondForce
    if(strcmp(child->Name(), "force") == 0 && strcmp(child->Attribute( "type" ), "HarmonicBondForce") == 0){
      //report << "XML name " << child->Name() << ", " << child->Attribute( "type" ) << endr;
      
      //find children of force, HarmonicBondForce
      tinyxml2::XMLElement *forceharmonic = child->FirstChildElement("bonds");
      
      //number in file
      const int fhcount = atoi(forceharmonic->Attribute( "count" ));
      
      report << "Harmonic count " << fhcount << endr;
      
      //loop through it
      for (tinyxml2::XMLElement* fnbchild = forceharmonic->FirstChildElement(); fnbchild != NULL; fnbchild = fnbchild->NextSiblingElement()){
        
        report << debug(910) << "Harmonic k " << fnbchild->Attribute( "k" ) << endr;
        
        //create bond
        Bond tempbond;
        tempbond.restLength = atof(fnbchild->Attribute( "length" )) * Constant::NM_ANGSTROM;
        tempbond.springConstant = atof(fnbchild->Attribute( "k" )) * Constant::KJ_KCAL * Constant::ANGSTROM_NM * Constant::ANGSTROM_NM * 0.5;
        tempbond.atom1 = atof(fnbchild->Attribute( "particle1" )) - 1;
        tempbond.atom2 = atof(fnbchild->Attribute( "particle2" )) - 1;
        topo->bonds.push_back(tempbond);
        
        // populate the vector of bonds maintained at each atom
        topo->atoms[tempbond.atom1].mybonds.push_back((topo->bonds.size()) - 1);
        topo->atoms[tempbond.atom2].mybonds.push_back((topo->bonds.size()) - 1);
        
        if (!tempbond.springConstant) ignoredBonds++;

      }
      
      //test right number
      if(topo->bonds.size() != fhcount){
        report << error << "Number of bonds wrong " << topo->bonds.size() << ". " << fhcount << endr;
      }
      
      //go as there is only one set
      break;
    }
  }
  
  //test no error
  if(doc.ErrorID()) report << error << "XML File parsing Harmonic Bonds error!" << endr;

  // ~~~~Angles~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  
#endif
    //end here until viable force field produced
  report << error << "Force field incomplete for XML." << endr;
}

