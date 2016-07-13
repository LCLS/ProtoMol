/*  -*- c++ -*-  */
#ifndef PROTOMOL_OUTPUT_CHECKPIINT_H
#define PROTOMOL_OUTPUT_CHECKPIINT_H

#include <protomol/base/StringUtilities.h>
#include <protomol/base/Timer.h>
#include <protomol/output/Output.h>
#include <string>
#include <vector>

namespace ProtoMol {
  class Configuration;
class Value;
struct Parameter;

  class OutputCheckpoint : public Output {
  public:
    static const std::string keyword;

  protected:
    int current;
    std::string name;
    std::string posBase, velBase;

  public:
    OutputCheckpoint() : current(0) {}
    OutputCheckpoint(const std::string &name, int freq, int start,
                      const std::string &posbase, const std::string &velbase);

  private:
    void WritePositions(long step);
    void WriteVelocities(long step);
    void WriteConfig(long step);

  public:
    void doIt(long step);

  private:
    //  From class Output
    Output *doMake(const std::vector<Value> &values) const;
    void doInitialize();
    void doRun(long step);
    void doFinalize(long) {}
    bool isIdDefined(const Configuration *config) const;
    bool addDoKeyword() const {return false;}

  public:
    //  From class Makeable
    std::string getIdNoAlias() const {return keyword;}
    void getParameters(std::vector<Parameter> &) const;
    bool adjustWithDefaultParameters(std::vector<Value> &values,
                                     const Configuration *config) const;
  };
}

#endif //  PROTOMOL_OUTPUT_CHECKPIINT_H
