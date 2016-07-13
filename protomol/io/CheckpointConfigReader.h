/*  -*- c++ -*-  */
#ifndef CHECKPOINTCONFIGREADER_H
#define CHECKPOINTREADERREADER_H

#include <protomol/base/Random.h>
#include <protomol/config/Configuration.h>
#include <protomol/integrator/Integrator.h>
#include <protomol/io/Reader.h>
#include <iosfwd>
#include <string>

namespace ProtoMol {
  //_________________________________________________________CheckpointConfigReader
  /*
   */
class Configuration;
class Integrator;
class Random;

  class CheckpointConfigReader : public Reader {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors (both default here), assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    CheckpointConfigReader();
    CheckpointConfigReader(const std::string &filename);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Reader
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual bool tryFormat();
    virtual bool read() { return !file.fail(); }

    bool readBase( Configuration& conf, Random &rand );
    bool readIntegrator( Integrator* integ );
  };
}
#endif /* CHECKPOINTCONFIGREADER_H */
