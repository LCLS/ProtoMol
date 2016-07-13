#include "OutputDihedrals.h"

#include <protomol/ProtoMolApp.h>
#include <protomol/output/OutputCache.h>
#include <protomol/topology/GenericTopology.h>
#include <iomanip>
#include <set>

#include "protomol/base/Exception.h"
#include "protomol/base/Report.h"
#include "protomol/config/ConstraintValueType.h"
#include "protomol/config/Parameter.h"
#include "protomol/config/Value.h"
#include "protomol/output/Output.h"
#include "protomol/type/Real.h"

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;

using std::string;
using std::vector;
using std::set;
using std::setw;
using std::endl;
using std::flush;
using std::stringstream;
using std::setprecision;
using std::setiosflags;
using std::resetiosflags;
using std::ofstream;
using std::ifstream;

  //________________________________________________________ Output
  const string  OutputDihedrals::keyword("dihedralsFile");

  OutputDihedrals::OutputDihedrals(): Output(),
              myFilename(""),
				      myDihedralIndex(-1),
				      myDihedralsSetfile(""){}
  
  OutputDihedrals::OutputDihedrals(const string& filename, 
                                    int freq, 
                                    int index, 
                                    std::string dsetfile):
            Output(freq),
            myFilename(filename),
            myDihedralIndex(index),
            myDihedralsSetfile(dsetfile){
    }

  OutputDihedrals::~OutputDihedrals(){
  }

  void OutputDihedrals::doInitialize(){
    
    //test inline index
    if(myDihedralIndex < 0 || myDihedralIndex >= static_cast<int>(app->topology->dihedrals.size())){
      report << error << "[OutputDihedrals::doInitialize] Dihedral index "<<myDihedralIndex
             <<" out of range [0,"<<app->topology->dihedrals.size()-1<<"]."<<endr;
    }
            
    //load file if available
    if(myDihedralsSetfile.size()){
        
      ifstream dihedralsSetinput(string(myDihedralsSetfile).c_str(), std::ios::in);
      if(!dihedralsSetinput)
        report << error <<" Could not open \'"<<myDihedralsSetfile<<"\' for "<<getId()<<"."<<endr;
      int tempdihedral = 1;
      while (dihedralsSetinput >> tempdihedral)
        myDihedrals.push_back(tempdihedral);
    }else{
      //or just load index
      myDihedrals.push_back(myDihedralIndex);
    }

    file.open(myFilename.c_str(), ios::out);
    if (!file){ 
      THROWS("Failed to open dihedrals file '" << myFilename << "'");
    }else{
      file << "Dihedral index, Dihedral value." << endl;
    }

  }

  // The run function outputs the dihedral values
  void OutputDihedrals::doRun(long){

    //get the dihedral values from the chache
    vector<Real> dihedrals = app->outputCache.getDihedralPhis(myDihedrals);
    
    //output data to file
    for(unsigned i=0; i<dihedrals.size(); i++){
      file << myDihedrals[i] << " " << dihedrals[i] << endl;
    }

  }

  void OutputDihedrals::doFinalize(long){
    file.close();
  }
  
    
  Output* OutputDihedrals::doMake(const std::vector<Value>& values) const{
    return (new OutputDihedrals(values[0],values[1],values[2],values[3]));
  }

  void OutputDihedrals::getParameters(std::vector<Parameter> &parameter) const{
    parameter.push_back(Parameter(keyword, Value(myFilename, ConstraintValueType::NotEmpty())));
    Output::getParameters(parameter);
    parameter.push_back(Parameter("dihedralsIndex",Value(myDihedralIndex,ConstraintValueType::NotNegative()),1));
    parameter.push_back(Parameter("dihedralsSetfile",Value(myDihedralsSetfile),""));
  }
