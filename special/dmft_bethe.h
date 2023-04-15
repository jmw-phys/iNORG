
#pragma once

/*
coded by Jia-Ming Wang (jmw@ruc.edu.cn, RUC, China) date 2022
*/

//! not finished yet 

#include "specs.h"
#include "prmtr.h"
#include "model.h"
#include "green.h"
#include "hyberr.h"
#include "bath.h"
#include "impurity.h"
#include "norg.h"

// the bethe lattice model to study the multi-orbital DMFT calculation.

class DMFT {
public:
	const MyMpi& mm;
	const Prmtr& p;
	const Model& mdl;
	Int iter_cnt;
	VecReal gloc_err;
	
	Real n_1;				// number of electrons on d orbital
	Real n_2;				// number of electrons on p orbital
	// Real ded;				// delta epsd to adjust n_1
	// Real dep;				// delta epsp to adjust n_2

	// The iteration variable.
	ImGreen g_loc;			// if in mode = 0, the DMFT iteration on G (default).
	ImGreen se;				// if in mode = 1, the DMFT iteration on self-energy.
	bool imp_backup;		// to tell if you has the impurity model's back up.
private:
	bool converged() const {
		const Real dev = DEV(gloc_err);
		if(gloc_err[gloc_err.size()-2]<1.E-5 && gloc_err[gloc_err.size()-2]<gloc_err[gloc_err.size()-1]) return true;
		if(gloc_err[gloc_err.size()-1]<1.E-6) return true;
		if (dev > 1.E-4) { return false; }
		else if (dev > 1.E-10) {
			for_Int(i, 1, gloc_err.size()) {
				if (gloc_err[0] > gloc_err[i]) { return false; }
			}
			return true;
		}
		else{ return true; }
	}
	void append_gloc_err(Real err) {
		for_Int(i, 1, gloc_err.size()) {
			gloc_err[i - 1] = gloc_err[i];
		}
		gloc_err[gloc_err.size() - 1] = err;
	}
	void print_log(const Str& lbl, std::ostream& os = std::cout) const {
		using namespace std;
		os << iofmt("sci");
		os << setw(4) << iter_cnt
			<< "  " << setw(w_Real) << gloc_err[gloc_err.size() - 1]
			<< "  " << setw(w_Real) << n_1
			<< "  " << setw(w_Real) << n_2
			<< "  " << setw(w_Real) << ded
			<< "  " << setw(w_Real) << dep
			// << "  " << setw(w_Real) << mdl.epsd
			// << "  " << setw(w_Real) << mdl.epsp
			<< "  " << present()
			<< "  " << lbl << endl;
	}
	void log(const Str& lbl) {
		if (mm) {
			print_log(lbl);
			OFS ofs(p.ofx + ".log.txt", std::ios::app);
			print_log(lbl, ofs);
			ofs.close();
		}
	}
	void find_gfloc(const Green& seloc, Green& gfloc) const;
	ImGreen find_weiss(ImGreen& g_loc) const;
	ImGreen find_hb(const ImGreen& g0eff_inv) const;
	
	ReGreen find_lattice_se(const ReGreen& g_loc) const;
	
	ImGreen find_imp_se(const ImGreen& g_imp, const ImGreen& hb_imp) const;
	
	bool check_gloc(const Str& file);
	bool check_seloc(const Str& file);
	void save_the_backup(Bath& bth, NORG& solver, Int iter_cnt = 999);
	void do_you_also_have_these_backup(Bath& bth, NORG& solver, bool if_backuped=false);

	// cala by self-energy.
	ImGreen find_gloc_by_se(const ImGreen& se_i) const;

	void find_gloc_by_se(Green& g_loc, const Green& se_i, Int kind) const;

	// (Deprecated) returen all the eta = 0.01 + a[i] * w, for green function.
	void print_gloc_by_sefe_with_different_eta(const Green& se_i, VecReal a) const;
	// (Deprecated)
	void print_gloc_by_gimp_with_different_eta(const Green& gimp, VecReal a) const;
	
	ImGreen find_hb_by_se(const ImGreen& se_i) const;
	

public:
	DMFT(const MyMpi& mm_i, const Prmtr& prmtr_i, const Model& mdl_i);
	// For mode = 1, DMFT iteration by self-energy, for mode = 0, DMFT iteration by g_imp.
	DMFT(const MyMpi& mm_i, const Prmtr& prmtr_i, const Model& mdl_i, const Int mode);
};
