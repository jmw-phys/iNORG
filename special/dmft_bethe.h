
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
#include "occler.h"

// the bethe lattice model to study the multi-orbital DMFT calculation.

class DMFT {
// ------------------------------- since the 2023.04.15 -------------------------------
	const MyMpi& mm;
	Prmtr& p;

// bethe lattice part:
	Real bethe_u;				// intra-orbital coupling U
	Real bethe_u12;				// inter-orbital coupling U12
	Real bethe_mu;				// Chemical potential /mu
	VecReal bethe_t;			// hopping strength t

	
	VEC<ImGreen> se_input;		// Pulay mixing
    VEC<ImGreen> res_past;		// Pulay mixing: residual 

	ImGreen g0_loc() const;		// the local lattice Green's function with no interaction.
	void set_parameter();		// set the parameters for the DMFT calculation.

public:
	Int iter_cnt;
	VecReal gloc_err;
	
	VecReal n_eles;				// number of electrons on d orbital

	// The iteration variable.
	ImGreen g_loc;			// if in mode = 0, the DMFT iteration on G (default).
	ImGreen se;				// if in mode = 1, the DMFT iteration on self-energy.
	bool imp_backup;		// to tell if you has the impurity model's back up.

private:
	void auto_nooc(Str mode, const Impurity &imp);

	bool converged() const;

	void append_gloc_err(Real err) {
		for_Int(i, 1, gloc_err.size()) {
			gloc_err[i - 1] = gloc_err[i];
		}
		gloc_err[gloc_err.size() - 1] = err;
	}

	void print_log(const Str& lbl, std::ostream& os = std::cout) const {
		using namespace std;
		Str temp; for_Int(i, 0, p.npartical.size()) {if (i == 0) temp += STR(p.npartical[i]); else /*if (i % 2 == 0)*/ temp += "-" + STR(p.npartical[i]);}
		
		os << iofmt();
		os << setw(4) << iter_cnt;
		// os << setw(8) << p.U;
		os << iofmt("sci");
		os << "  " << setw(w_Real) << gloc_err[gloc_err.size() - 1];
		os << setw(4 + 2 * p.nband) << temp;
		for_Int(i, 0, n_eles.size()) {
			os << "  " << setw(w_Real) << n_eles[i];
		}
		os << "  " << present();
		os << "  " << lbl << endl;
	}

	void log(const Str& lbl) {
		if (mm) {
			print_log(lbl);
			OFS ofs("log.txt", std::ios::app);
			print_log(lbl, ofs);
			ofs.close();
		}
	}

	ImGreen find_hb(const ImGreen& g0eff_inv) const;
	ImGreen find_imp_se(const ImGreen& g_imp, const ImGreen& hb_imp) const;

	bool check_gloc(const Str& file);
	bool check_seloc(const Str& file);
	// void save_the_backup(Bath& bth, NORG& solver, Int iter_cnt = 999);
	// void do_you_also_have_these_backup(Bath& bth, NORG& solver, bool if_backuped=false);

	ImGreen find_gloc_by_se(const ImGreen& se_i) const;
	ReGreen find_gloc_by_se(const ReGreen& se_i) const;

	ImGreen find_hb_by_se(const ImGreen& se_i) const;


	// Calculate Green's function using parity degeneracy
	void get_gimp_evenodd(NORG & norg,Green & gfimp) const;

	
    void pulay_mixing(const ImGreen& seimp);
	Real dot(const ImGreen& r1, const ImGreen& r2) {
		Real t(0.);
		for_Int(i, 0, r1.g.size()) {
			for_Int(row, 0, r1.g[i].nrows()) {
				for_Int(col, 0, r1.g[i].ncols()) {
					t += r1[i][row][col].real() * r2[i][row][col].real() + r1[i][row][col].imag() * r2[i][row][col].imag();
				}
			}
		}
		return t;
	}
public:
	// For mode = 1, DMFT iteration by self-energy, for mode = 0, DMFT iteration by g_imp.
	DMFT(const MyMpi& mm_i, Prmtr& prmtr_i, const Int mode);
};
