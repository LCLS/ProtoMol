/*  -*- c++ -*-  */
#ifndef PROTOMOL_OUTPUT_FINAL_XYZ_VEL_H
#define PROTOMOL_OUTPUT_FINAL_XYZ_VEL_H

#include <protomol/output/Output.h>
#include <string>
#include <vector>

namespace ProtoMol {
  class Configuration;
class Value;
struct Parameter;

  class OutputFinalXYZVel : public Output {
  public:
    static const std::string keyword;

  private:
    std::string filename;

  public:
    OutputFinalXYZVel();
    OutputFinalXYZVel(const std::string &filename);

    //   From class Output
  private:
    Output *doMake(const std::vector<Value> &values) const;
    void doInitialize() {};
    void doRun(long) {};
    void doFinalize(long step);

    //  From class Makeabl
  public:
    std::string getIdNoAlias() const {return keyword;}
    void getParameters(std::vector<Parameter> &parameter) const;
  };
}
#endif //  PROTOMOL_OUTPUT_FINAL_XYZ_VEL_H
