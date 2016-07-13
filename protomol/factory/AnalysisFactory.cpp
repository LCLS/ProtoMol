#include <protomol/analysis/AnalysisCollection.h>
#include <protomol/base/StringUtilities.h>
#include <protomol/config/Configuration.h>
#include <protomol/factory/AnalysisFactory.h>
#include <protomol/factory/HelpTextFactory.h>
#include <stddef.h>
#include <map>
#include <ostream>
#include <set>
#include <utility>

#include "protomol/analysis/Analysis.h"
#include "protomol/base/Exception.h"
#include "protomol/base/Factory.h"
#include "protomol/config/Parameter.h"
#include "protomol/config/Value.h"
#include "protomol/type/SimpleTypes.h"

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;

void AnalysisFactory::registerAllExemplarsConfiguration(Configuration *config) const {
	for( pointers_t::const_iterator itr = pointers.begin(); itr != pointers.end(); ++itr ) {
		const Analysis *prototype = ( *itr );
		vector<Parameter> parameter;
		prototype->getParameters(parameter);

		for( unsigned int i = 0; i < parameter.size(); i++ ) {
			config->registerKeyword(parameter[i].keyword, parameter[i].defaultValue);
			if( !parameter[i].text.empty()) {
				config->setText(parameter[i].keyword, parameter[i].text);
			}
		}

		if( prototype->addDoKeyword()) {
			config->registerKeyword("do" + prototype->getId(), Value(true, Value::undefined));
		}
	}
}

Analysis *AnalysisFactory::make(const string &id, const vector<Value> &values) const {
	// Find prototype
	const Analysis *prototype = getPrototype(id);

	if( prototype == NULL ) {
		THROWS(" Could not find any match for '" << id << "' in " << Analysis::scope << "Factory. Possible outputs are:\n" << *this);
	}

	// Make
	Analysis *newObj = prototype->make(values);
	if( !newObj ) { THROW("Could not make output from prototype"); }

	// Adjust external alias
	newObj->setAlias(id);
	return newObj;
}

AnalysisCollection* AnalysisFactory::makeCollection(const Configuration *config) const {
	AnalysisCollection *res = new AnalysisCollection();

	for( Configuration::const_iterator i = config->begin(); i != config->end(); ++i ) {
		if(( *i ).second.valid()) {
			const Analysis *prototype = getPrototype(( *i ).first);

			if( prototype != NULL ) {
				if( prototype->isIdDefined(config)) {
					vector<Parameter> parameter;
					prototype->getParameters(parameter);
					vector<Value> values(parameter.size());

					for( unsigned int k = 0; k < parameter.size(); k++ ) {
						if( config->valid(parameter[k].keyword) && !( parameter[k].keyword.empty())) {
							values[k].set(( *config )[parameter[k].keyword]);
						} else if( parameter[k].defaultValue.valid()) {
							values[k].set(parameter[k].defaultValue);
						} else {
							values[k] = parameter[k].value;
							values[k].clear();
						}
					}

					if( !prototype->checkParameters(values)) {
						prototype->adjustWithDefaultParameters(values, config);
					}

					prototype->assertParameters(values);

					res->adoptAnalysis(make(( *i ).first, values));
				}
			}
		}
	}

	return res;
}

void AnalysisFactory::registerHelpText() const {
	for( exemplars_t::const_iterator i = exemplars.begin(); i != exemplars.end(); ++i ) {
		HelpText helpText;
		i->second->getParameters(helpText.parameters);
		helpText.id = i->second->getIdNoAlias();
		helpText.text = i->second->getText();
		helpText.scope = i->second->getScope();
		if( i->second->addDoKeyword()) { helpText.parameters.push_back(Parameter("do" + i->second->getId(), Value(true, Value::undefined), Text("flag to switch on/off the output"))); }
		HelpTextFactory::registerExemplar(i->second->getId(), helpText);

		string txt = string("a parameter of " + i->second->getId());
		for( unsigned int j = 0; j < helpText.parameters.size(); j++ ) {
			helpText.text = txt;
			if( !equalNocase(i->second->getId(), helpText.parameters[j].keyword)) {
				if( !helpText.parameters[j].text.empty()) {
					helpText.text = txt + ", " + helpText.parameters[j].text;
				}
				HelpTextFactory::registerExemplar(helpText.parameters[j].keyword, helpText);
			}
		}

		HelpText alias;
		alias.text = "alias for \'" + i->second->getId() + "\'";
		alias.scope = i->second->getScope();
		for( exemplars_t::const_iterator j = aliasExemplars.begin(); j != aliasExemplars.end(); ++j ) {
			if( j->second->getIdNoAlias() == i->second->getIdNoAlias()) {
				alias.id = j->first;
				HelpTextFactory::registerExemplar(alias.id, alias);
			}
		}
	}
}
