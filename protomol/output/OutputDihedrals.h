/*  -*- c++ -*-  */
#ifndef OUTPUTDIHEDRALS_H
#define OUTPUTDIHEDRALS_H

#include <protomol/output/Output.h>
#include <protomol/type/Vector3DBlock.h>
#include <algorithm>
#include <set>
#include <string>
#include <vector>

namespace ProtoMol {
class Value;
struct Parameter;
}  // namespace ProtoMol

#ifdef HAVE_LIBFAH
#include <cbang/os/File.h>
#else
#include <fstream>
#endif

namespace ProtoMol {

  //________________________________________________________ OutputDihedrals
  /** 
      Writes the dihedral values to a file at given freqeuncy.
  */
  class OutputDihedrals : public Output {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    OutputDihedrals();
    OutputDihedrals(const std::string& filename, int freq, int index, std::string dsetfile);
    ~OutputDihedrals();
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class OutputDihedrals
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //  From class OutputFile
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    void doRun(long step);
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //  From class Output
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    void doInitialize();
    void doFinalize(long step);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Makeable
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    Output *doMake(const std::vector<Value>& values) const;
    std::string getIdNoAlias() const{ return keyword;}

    // Returns the identification string

    //virtual unsigned int getParameterSize() const {return 4;}
    //virtual bool adjustWithDefaultParameters(std::vector<Value>& values, const Configuration* config) const;
    void getParameters(std::vector<Parameter> &parameter) const;
  private:
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    static const std::string keyword;

  private:
    string myFilename;
    int myDihedralIndex;
    std::string myDihedralsSetfile;
    std::vector< int > myDihedrals;
    
  protected:
#ifdef HAVE_LIBFAH
    cb::File file;
#else
    std::ofstream file;
#endif

  };
}
#endif
