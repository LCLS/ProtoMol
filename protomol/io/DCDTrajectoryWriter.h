/*  -*- c++ -*-  */
#ifndef DCDTRAJECTORYWRITER_H
#define DCDTRAJECTORYWRITER_H

#include <protomol/base/SystemUtilities.h>
#include <protomol/io/Writer.h>
#include <protomol/type/TypeSelection.h>
#include <protomol/type/XYZ.h>
#include <iosfwd>
#include <string>
#include <vector>

#include "protomol/type/Real.h"

namespace ProtoMol {
  //____DCDTrajectoryWriter

  /**
   * Writes DCD trajectories and updates the number of coordinate sets
   * after each write, no need to know the final number of sets.
   */
class Vector3DBlock;
struct XYZ;

  class DCDTrajectoryWriter : public Writer {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Typedef
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    typedef TypeSelection::Int<4>::type int32;
    typedef TypeSelection::Float<4>::type float4;
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors (both default here), assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    explicit DCDTrajectoryWriter(Real timestep = 1.0, unsigned int firststep =
                                   1,
                                 bool isLittleEndian = ISLITTLEENDIAN);
    explicit DCDTrajectoryWriter(const std::string &filename, Real timestep =
                                   1.0, unsigned int firststep = 1,
                                 bool isLittleEndian = ISLITTLEENDIAN);
    explicit DCDTrajectoryWriter(const char *filename, Real timestep = 1.0,
                                 unsigned int firststep = 1,
                                 bool isLittleEndian = ISLITTLEENDIAN);
    explicit DCDTrajectoryWriter(std::ios::openmode mode,
                                 const std::string &filename, Real timestep =
                                   1.0, unsigned int firststep = 1,
                                 bool isLittleEndian = ISLITTLEENDIAN);
    explicit DCDTrajectoryWriter(std::ios::openmode mode, const char *filename,
                                 Real timestep = 1.0, unsigned int firststep =
                                   1,
                                 bool isLittleEndian = ISLITTLEENDIAN);
    virtual ~DCDTrajectoryWriter();
    
    //special create for append to DCD
    explicit DCDTrajectoryWriter(std::ios::openmode mode, int frameoffs,
                                 const std::string &filename, Real timestep =
                                   1.0, unsigned int firststep = 1,
                                 bool isLittleEndian = ISLITTLEENDIAN);
    
    // Need this implementation, otherwise const char* will bee converted to 
    // bool or int ...

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class DCD
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    bool openWith(Real timestep, unsigned int firststep = 1,
                  bool isLittleEndian = ISLITTLEENDIAN);
    bool openWith(const std::string &filename, Real timestep,
                  unsigned int firststep = 1,
                  bool isLittleEndian = ISLITTLEENDIAN);
    bool openWith(const char *filename, Real timestep,
                  unsigned int firststep = 1,
                  bool isLittleEndian = ISLITTLEENDIAN);

    bool write(const Vector3DBlock &coords);
    bool write(const std::vector<Vector3DBlock> &cachedCoords);

    void setLittleEndian(bool littleEndian);
    void setTimestep(Real timestep);
    void setFirststep(unsigned int firststep);

    bool reopen(unsigned int numAtoms, unsigned int numberOfSets);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Friends
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    friend DCDTrajectoryWriter &operator<<(DCDTrajectoryWriter &dcdWriter,
                                           const Vector3DBlock &coords);

    friend DCDTrajectoryWriter &operator<<(DCDTrajectoryWriter &dcdWriter,
                                           const XYZ &xyz);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    int myFrameOffset;
    bool myIsLittleEndian;
    unsigned int myFirstStep;
    Real myTimeStep;
    bool firstWrite;
    std::vector<float4> myX;
    std::vector<float4> myY;
    std::vector<float4> myZ;
    
  };

  //____INLINES
}
#endif /* DCDTRAJECTORYWRITER_H */
