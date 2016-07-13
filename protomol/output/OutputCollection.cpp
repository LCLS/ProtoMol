#include <protomol/ProtoMolApp.h>
#include <protomol/base/Exception.h>
#include <protomol/base/MathUtilities.h>
#include <protomol/output/Output.h>
#include <protomol/output/OutputCollection.h>
#include <algorithm>
#include <string>

#include "protomol/base/PMConstants.h"
#include "protomol/output/OutputCache.h"

using namespace ProtoMol;


OutputCollection::~OutputCollection() {
  for (iterator i = begin(); i != end(); i++) delete (*i);
}


void OutputCollection::initialize(const ProtoMolApp *app) {
  this->app = app;
  for (iterator i = begin(); i != end(); i++) (*i)->initialize(app);
}


bool OutputCollection::run(long step) {
  bool outputRan = false;

  app->outputCache.uncache();
  for (iterator i = begin(); i != end(); ++i)
    outputRan |= (*i)->run(step);

  return outputRan;
}


void OutputCollection::finalize(long step) {
  app->outputCache.uncache();
  for (iterator i = begin(); i != end(); i++) (*i)->finalize(step);
}


void OutputCollection::adoptOutput(Output *output) {
  if (!output) THROW("null pointer");
  outputList.push_back(output);
}


long OutputCollection::getNext() const {
  long next = Constant::MAX_LONG;
  for (const_iterator i = begin(); i != end(); i++)
    next = min((*i)->getNext(), next);

  return next;
}
