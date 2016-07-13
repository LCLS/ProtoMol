#include <protomol/io/Writer.h>
#include <iostream>

#include "protomol/io/File.h"

using namespace std;
using namespace ProtoMol;

//____ Writer
Writer::Writer() : File(ios::out | ios::trunc) {}
Writer::Writer(const string &filename) :
  File(ios::out | ios::trunc, filename) {}
Writer::Writer(ios::openmode mode) : File(ios::out | mode) {}
Writer::Writer(ios::openmode mode, const string &filename) :
  File(ios::out | mode, filename) {}
