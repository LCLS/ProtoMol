/*  -*- c++ -*-  */
#ifndef PROTOMOL_OUTPUT_SCREEN_H
#define PROTOMOL_OUTPUT_SCREEN_H

#include <protomol/output/Output.h>
#include <string>
#include <vector>

#include "protomol/type/Real.h"

namespace ProtoMol {
  class Configuration;
class Value;
struct Parameter;

  class OutputScreen : public Output {
  public:
    static const std::string keyword;

  private:
    std::string unit;
    Real factor;

  public:
    OutputScreen();
    OutputScreen(int freq);

    //   From class Output
  private:
    Output *doMake(const std::vector<Value> &values) const;
    void doInitialize();
    void doRun(long step);
    void doFinalize(long) {}
    bool isIdDefined(const Configuration *config) const;
    bool addDoKeyword() const {return false;}

    //  From class Makeabl
  public:
    std::string getIdNoAlias() const {return keyword;}
    void getParameters(std::vector<Parameter> &) const;
  };
}
#endif //  PROTOMOL_OUTPUT_SCREEN_H
