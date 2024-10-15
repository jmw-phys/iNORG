#pragma once
/*
coded by Jia-Ming Wang (jmw@ruc.edu.cn, RUC, China)
date 2022 - 2023
*/

#include "specs.h"
#include "prmtr.h"
#include "bath.h"
#include "model.h"
#include "green.h"
#include "impurity.h"
#include "norg.h"
// #include "asnci.h"
#include "occler.h"

#include <chrono>
#include <thread>
#include <cstdio>



// A class for transform the data of ZEN.
// - Preprocessor. (APIzen)
//   - read hybrid function
//   - do the fitting
//   - Temporarily store the build matrix information required by NORG

class APIzen{
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

	
	// NORG coding console
	Str mode;
	VecInt restrain, distribute;
	Real weight_nooc, weight_freze;

	// NORG test part
	//Int nimp;
	VecReal bathose, bathhop;

public:
	VEC<VecReal> t_ose;						// hopping integral for all sites
	VEC<VecReal> t_hyb;						// H_0 onset energy
	VecReal muvec;
	Int dmft_cnt;
private:
	void update(const Str& file);

	bool if_lock(const Str file) const;

	void read_ZEN(const Str& file);

	NORG choose_cauculation_style(Str mode, Impurity &imp);
	// void fitting();

	void test_for_fitting(const Bath& bth, const ImGreen& hby_i, Int num = 666);
/*
	ImGreen fix_se(const ImGreen& se) const;
*/
	void seimp_fixer(ImGreen& seimp_in);

	void read_norg_setting( const std::string& filename, std::string& CoulombF, double& U, double& J, std::vector<int>& restrain, std::vector<int>& distribute );

	void auto_nooc(Str mode, const Impurity &imp);
public:
	APIzen(const MyMpi& mm_i, Prmtr& p, const Str& file = empty_str);

};