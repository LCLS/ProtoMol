/*  -*- c++ -*-  */
#ifndef PROTOMOL_OUTPUT_DCD_TRAJECTORY_FORCE_H
#define PROTOMOL_OUTPUT_DCD_TRAJECTORY_FORCE_H

#include <protomol/output/Output.h>
#include <string>
#include <vector>

namespace ProtoMol {
  class DCDTrajectoryWriter;
class Configuration;
class Value;
struct Parameter;

  class OutputDCDTrajectoryForces : public Output {
  public:
    static const std::string keyword;

  private:
    DCDTrajectoryWriter *dCD;
    bool minimalImage;

  public:
    OutputDCDTrajectoryForces();
    OutputDCDTrajectoryForces(const std::string &filename, int freq, bool minimal);
    virtual ~OutputDCDTrajectoryForces();

    // From class Output
  private:
    Output *doMake(const std::vector<Value> &values) const;
    void doInitialize();
    void doRun(long step);
    void doFinalize(long step);

    // From class Makeable
  public:
    std::string getIdNoAlias() const {return keyword;}
    void getParameters(std::vector<Parameter> &parameter) const;
    bool adjustWithDefaultParameters(std::vector<Value> &values,
                                     const Configuration *config) const;
  };
}
#endif // PROTOMOL_OUTPUT_DCD_TRAJECTORY_FORCE_H
