/* -*- c++ -*- */
#ifndef ANALYSIS_FACTORY_H
#define ANALYSIS_FACTORY_H

#include <protomol/analysis/Analysis.h>
#include <protomol/base/Factory.h>
#include <protomol/config/Value.h>
#include <string>
#include <vector>

namespace ProtoMol {
	class AnalysisCollection;
class Configuration;
class Value;

	class AnalysisFactory : public Factory<Analysis> {
		public:
			virtual void registerHelpText () const;

		public:
			void registerAllExemplarsConfiguration (Configuration *config) const;
			Analysis *make (const std::string &id, const std::vector<Value> &values) const;
			AnalysisCollection *makeCollection (const Configuration *config) const;
	};
}

#endif	// ANALYSIS_FACTORY_H
