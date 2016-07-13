#ifndef PROTOMOL_ANALYSIS_ANALYSIS_H
#define PROTOMOL_ANALYSIS_ANALYSIS_H

#include <protomol/base/Makeable.h>
#include <string>
#include <vector>

namespace ProtoMol {
	class ProtoMolApp;
	class Configuration;
	class Value;
	struct Parameter;

	class Analysis : public Makeable<Analysis> {
		public:
			static const std::string scope;
			const ProtoMolApp *app;

		protected:
			bool onOutput;
			bool mShouldStop;

		public:
			Analysis(bool output = 1);

			// / To initialize the object, before the simulation starts.
			virtual void initialize (const ProtoMolApp *app);

			// / Called at each step (e.g., printing total energy on the screen),
			// / takes care of the output frequency.  Returns true if it ran.
			virtual bool run (long step);

			// Called after running to determine if the simulation should be stopped
			bool shouldStop();

			// / At the end of the simulation (e.g., writing final positions), and
			// / calls first run() to ensure that run is called for the last
			// / step, if needed.
			virtual void finalize (long step);

			// / Factory method to create a complete output object from its prototy
			virtual Analysis *make (const std::vector<Value> &values) const;

			// / Should return true if the concrete object is defined/specified in
			// / Configuration by the user. Normally if gedId() has a valid valu
			// / in Configuration.
			virtual bool isIdDefined (const Configuration *config) const;

			// / Defines if the output object supports do<getId()> to enable or disabl
			// / the output.
			virtual bool addDoKeyword() const { return true; }

			bool isOnOutput() const { return onOutput; }

		private:
			virtual void doInitialize () = 0;
			virtual void doRun (long step) = 0;
			virtual void doFinalize(long step) {}

		public:
			virtual std::string getScope() const { return scope; }
			void getParameters (std::vector<Parameter> &parameter) const;
			bool adjustWithDefaultParameters (std::vector<Value> &values, const Configuration *config) const;
	};
}
#endif	//  PROTOMOL_ANALYSIS_ANALYSIS_H
