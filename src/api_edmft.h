#pragma once

/*
code developed and maintained by (jmw@ruc.edu.cn, RUC, China) date 2022 - 2024
*/

#include "specs.h"
#include "prmtr.h"
#include "bath.h"
// #include "model.h"
#include "green.h"
#include "impurity.h"
#include "norg.h"
// #include "asnci.h"
#include "occler.h"



class APIedmft{
	const MyMpi& mm;					// parameters
	Prmtr& p;							// parameters
	Int num_orbital;					// The number of orbiatals.
	Int num_nondegenerate;				// The number of nondegenerate orbiatal.
	Real num_omg;						// Number of positive imaginary frequencies used for bath fitting.
	MatCmplx imfrq_hybrid_function;		// Imaginary frequencies hybrid function, nrows():number of nondegenerate orbiatal.
	VecInt or_deg_idx;					// orbitals  degenerate idx.
	VecReal solver_eimp_data;			// Impurity energy level for orbitals, correspondingly mean Impurity energy.
	// ctqmcdata solver_ctqmc_data;
	Real Uc, Jz, mu;
	Int nband;							
	Int norbs;
	Int ful_pcl_sch;
	Int iter_count;
	// VecInt artificial_symm;

	
	// NORG coding console
	Str mode;
	Real fit_range;
	VecInt restrain, fit_nbaths;
	VecReal weight_nooc, weight_freze;

	// NORG test part
	//Int nimp;
	VecReal bathose, bathhop;

	// Iteration info
	Real sig_err;
	VecReal n_eles, fit_err;

public:
	VEC<VecReal> t_ose;						// hopping integral for all sites
	VEC<VecReal> t_hyb;						// H_0 onset energy
	VecReal muvec;
private:
	void update(const Str& file);
	void edmft_back_up(const Str& status);

	bool if_lock(const Str file) const;

	void read_eDMFT(const Str& file);

	// NORG choose_cauculation_style(Str mode, Impurity &imp);
	// void fitting();

/*
	void test_for_fitting(const Bath& bth, const ImGreen& hby_i, Int num = 666);
	ImGreen fix_se(const ImGreen& se) const;
*/

void read_norg_setting(
	const std::string& filename,
	std::vector<double>& Ed,
	std::vector<int>& Deg,
	double& J,
	std::string& CoulombF,
	double& beta,
	double& U,
	std::vector<double>& noc1,
	std::vector<double>& noc2,
	std::vector<int>& restrain,
	std::vector<int>& fit_nbaths);

	void process_occupation_sets(
		VEC<VecReal>& occnums,
		VEC<VecReal>& occweights,
		const VecInt& row_sizes,
		const VecReal& occnum_lst);

	void auto_nooc(Str mode, const Impurity &imp);

	void print_log(const Str& lbl, std::ostream& os = std::cout) const;

	void log(const Str& lbl) {
		if (mm) {
			print_log(lbl);
			OFS ofs("log.norg", std::ios::app);
			print_log(lbl, ofs);
			ofs.close();
		}
	}

public:
	APIedmft(const MyMpi& mm_i, Prmtr& p, const Str& file = empty_str);

};