#ifndef MODULE_H
#define MODULE_H

#include <protomol/base/Exception.h>
#include <set>
#include <string>

namespace ProtoMol {
  class ForceFactory;
  class GenericTopology;
  class ProtoMolApp;

  class Module {
  public:
    typedef std::set<std::string> module_deps_t;

    Module() {}
    virtual ~Module() {}

    virtual const std::string getName() const = 0;
    virtual int getPriority() const;
    virtual const std::string getHelp() const;
    virtual void getDependencies(module_deps_t &deps) const;

    virtual void configure(ProtoMolApp *app);
    virtual void read(ProtoMolApp *app);
    virtual void registerForces(ProtoMolApp *app);
    virtual void postBuild(ProtoMolApp *app);
    virtual void addModifiers(ProtoMolApp *app);
  protected:
    virtual void init(ProtoMolApp *app);

    friend class ModuleManager;
  };

  struct moduleLess {
    bool operator()(const Module *m1, const Module *m2) const;
  };
}

#endif // MODULE_H
