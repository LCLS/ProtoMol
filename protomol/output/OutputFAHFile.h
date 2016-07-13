/*  -*- c++ -*-  */
#ifndef PROTOMOL_OUTPUT_FAH_FILE_H
#define PROTOMOL_OUTPUT_FAH_FILE_H

#include <protomol/io/File.h>
#include <protomol/output/Output.h>
#include <string>
#include <vector>

namespace ProtoMol {
  class Configuration;
class Value;
struct Parameter;

  class OutputFAHFile : public Output, public File {
  public:
    static const std::string keyword;

  private:
    std::string filename;

  public:
    OutputFAHFile();
    OutputFAHFile(const std::string &filename, int freq);

    //   From class Output
  private:
    Output *doMake(const std::vector<Value> &values) const;
    void doInitialize();
    void doRun(long step);
    void doFinalize(long);
    bool isIdDefined(const Configuration *config) const;
    bool addDoKeyword() const {return false;}

    //  From class Makeabl
  public:
    std::string getIdNoAlias() const {return keyword;}
    void getParameters(std::vector<Parameter> &) const;
  };
}

#endif //  PROTOMOL_OUTPUT_FAH_FILE_H
