#ifndef PROTOMOL_ANALYSIS_DIHEDRAL_H
#define PROTOMOL_ANALYSIS_DIHEDRAL_H

#include <protomol/analysis/Analysis.h>
#include <protomol/base/StringUtilities.h>
#include <map>
#include <string>
#include <vector>

#include "protomol/type/Real.h"

namespace ProtoMol {
	class Configuration;
class Value;
struct Parameter;

	class AnalysisDihedral : public Analysis {
		public:
			static const std::string keyword;

		private:
			const std::vector<int> mIndex;
			const Real mPsiMin, mPsiMax;
			const Real mPhiMin, mPhiMax;

			struct Dihedral {
				int a, b, c, d;
				Dihedral() : a(0), b(0), c(0), d(0) {}
			};
			std::vector<Dihedral> mPhiDihedral, mPsiDihedral;

			std::vector<int> residues_alpha_c, residues_phi_n, residues_psi_c;

		public:
			AnalysisDihedral();
			AnalysisDihedral(const std::vector<int> index, const Real psiMin, const Real psiMax, const Real phiMin, const Real phiMax);

		public:
			void doIt (long step);

		private:
			Analysis *doMake (const std::vector<Value> &values) const;
			void doInitialize ();
			void doRun (long step);
			void doFinalize(long) {}
			bool isIdDefined (const Configuration *config) const;
			bool addDoKeyword() const { return false; }
			const Real computeDihedral (const int a1, const int a2, const int a3, const int a4) const;

		public:
			std::string getIdNoAlias() const { return keyword; }
			void getParameters (std::vector<Parameter> &) const;
			bool adjustWithDefaultParameters (std::vector<Value> &values, const Configuration *config) const;
	};
}

#endif	// PROTOMOL_ANALYSIS_DIHEDRAL_H
