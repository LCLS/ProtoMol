/*  -*- c++ -*-  */
#ifndef PROTOMOL_OUTPUT_DCD_TRAJECTORY_H
#define PROTOMOL_OUTPUT_DCD_TRAJECTORY_H

#include <protomol/output/Output.h>
#include <protomol/type/Vector3DBlock.h>
#include <string>
#include <vector>

namespace ProtoMol {
  class DCDTrajectoryWriter;
class Configuration;
class Value;
struct Parameter;

  class OutputDCDTrajectory : public Output {
  public:
    static const std::string keyword;

  private:
    DCDTrajectoryWriter *dCD;
    bool minimalImage;
    int frameOffset;
    string filename;
    
    //checkpoint correction at first output
    bool firstWrite;

    //additions to cache data
    std::vector<Vector3DBlock> cachedCoords;
    int cachesize;
    unsigned int cacheoffset;

  public:
    OutputDCDTrajectory();
    OutputDCDTrajectory(const std::string &filename, int freq,
                            bool minimal, int frameoffs, int cachesz);
    virtual ~OutputDCDTrajectory();

    //   From class Output
  private:
    Output *doMake(const std::vector<Value> &values) const;
    void doInitialize();
    void doRun(long step);
    void doFinalize(long step);

    //  From class Makeabl
  public:
    std::string getIdNoAlias() const {return keyword;}
    void getParameters(std::vector<Parameter> &parameter) const;
    bool adjustWithDefaultParameters(std::vector<Value> &values,
                                     const Configuration *config) const;
  };
}
#endif //  PROTOMOL_OUTPUT_DCD_TRAJECTORY_H
