#pragma once

/*
coded by Jia-Ming Wang (jmw@ruc.edu.cn, RUC, China) date 2022 - 2023
*/

#include "specs.h"

class Prmtr {
private:
	const Int np;					// number of mpi processes
public:
	// model related, model parameters
	Str project;					// project name
	Int nband;						// number of band in lattice model.
	Int norbs;						// number of spin-orbitals in lattice model.
	Real bandw;						// The bandwidth of this model.
	
	// square lattice, model parameters
	mutable Real U;					// Hubbard interaction U
	mutable Real jz;				// hund interaction Jz
	mutable Real Uprm;				// U^' as U - 2 * Jz
	mutable Real mu;				// chemical potential
	mutable Real degel;				// Degenerate energy levels
	VecReal	eimp;					// the impurity energy (spinless orbital).
	VecReal	bsr;					// the vector of bath sum rule.
	VecReal bethe_t;				// hopping strength t for bethe lattice
	//ReGreen
	Real freq_upp;					// upper bound of real frequency, 1.5 * (band upper bound) suggested
	Real freq_low;					// lower bound of real frequency, 1.5 * (band lower bound) suggested
	Real dlt_freq;					// delta frequency, freq = n * dlt_freq, freq[i] = freq_low + dlt_freq * i;
	Real eta_freq;					// eta of omega + i * eta
	Int nfreq;						// number of real frequencies = 1 + Int_ROUND((freq_upp - freq_low) / dlt_freq)	// dmft related
    VecCmplx Re_z;		

	//ImGreen		
	mutable Real beta;						// unit_omg = PI/beta;
	mutable Real unit_omg;					// unit imaginary frequency, omg_n = (2 n + 1) unit_omg, 0.01 or 0.02 suggested for zero temperature
	mutable Real max_omg;					// imaginary frequency cutoff, 4 * (half bandwidth) suggested
	mutable Real nmesh;						// imaginary frequency number form Zen.
	mutable Int num_omg;					// number of positive imaginary frequencies
	mutable VecCmplx Im_z;

	//fitting related		
	mutable Real fit_max_omg;				// imaginary frequency upper bound used for bath fitting
	mutable Int fit_num_omg;				// number of positive imaginary frequencies used for bath fitting
	mutable Real fit_pow;					// 1 / (omg_n + fit_rsd) ^ fit_pow is the weight used in bath fitting
	mutable Real fit_rsd;					// 1 / (omg_n + fit_rsd) ^ fit_pow is the weight used in bath fitting

	//iter related
	bool imp_backup;				// imp backup
	Int iter_max;					// DMFT max number of iterations
	Int gauss_n_max;				// max gauss_n for Green function integration
	Int gauss_n_min;				// min gauss_n for Green function integration

	// New prmtr for NORG solver.
	Int norg_sets;					// number of norg sets (spinless orbital).
	Int iter_max_norg;				// the max NORG iteration times.
	bool if_norg_imp;				// if norg rotate NO with impurity orbitals, if true, the imp_orb will rotate, treated same as bath orbital.
	Int if_norg_degenerate;			// set 0: the norg running  in the non-degenerate state; set 1: the norg running for finding the degenerate state.
	mutable VecInt nI2B;			// the number of bath sites for each impurity.
	mutable VecInt nO2sets;			// The number of orbital to each sets.
	mutable VecInt npartical;		// number of particals.
	mutable Int nbath;				// number of bath sites, must be an integer multiple of 4 
	mutable Int norbit;				// number of NORG orbital(imp + bath).
	mutable Int n_rot_orb;			// rotation orbital number.

	mutable MatInt control_divs;	// to set the number of division and the shortcut restratin.
	mutable VecInt templet_restrain;// to set the restrain templet one for all distribute;
	mutable VecInt templet_control;	// to set the control templet one for all distribute;
	mutable VecInt stage2_restrain;	// to set the stage2 restrain one for all distribute;
	mutable VecInt stage2_control;	// to set the stage2 control one for all distribute;
	mutable Int ndiv;				// the the divsion's nubmer.
	mutable VEC<MatReal> rotationU;	// save the rotation matrix for convience.

	
	mutable Str nooc_mode;

	//--------------------------------------------------------special for the hhd function(arXiv:2209.14178v1)-----------------------------------------------------------------
	Real alpha;						// the parameter for the hhd function(arXiv:2209.14178v1) to control the strength of the interaction in two wide band.
	Real delta;						// the parameter for the hhd function(arXiv:2209.14178v1) to control the strength of the interaction in two wide band.


	//---------------------------------------------------------------------------------------------------------------------------------------

	
	// derived
	Str ofx;				// output filename prefix

private:
	void set_inert_values();
	void set_values();
	void print(std::ostream &os, const Str &var, const Str &val, const Str &comment) const {
		using namespace std;
		Str true_var = var == "\"\"" ? "" : var;
		os << rght_justify(true_var, 16) << (true_var == "" ? "   " : " = ") << left_justify(val, w_Real)
			<< "    # " + comment << endl;
	}
	void derive();

	VEC<MatReal> uormat_initialize() const {
		VEC<MatReal> uormat_i;
		for_Int(i, 0, nO2sets.size()) {
			MatReal temp(dmat(nO2sets[i], 1.));
			uormat_i.push_back(std::move(temp));
		}
		return uormat_i;
	}


public:
	Prmtr(const MyMpi& mm);
	void after_modify_prmtr() const;
	void according_nppso(const VecInt& nppsos) const;
	void according_controler(const Vec<VecInt>& controler, const VecInt& or_deg) const;
	void derive_ImGreen() const;
	void recalc_partical_number() const;
	void print(std::ostream &os = std::cout) const;
	void change_the_norg_restrain_and_div(VecInt new_restrain, VecInt new_control) const;

	// Imag frequency
	Real Imomg(Int n) const { return imag(Im_z[n]); }
	Cmplx Imz(Int n) const { return Im_z[n]; }
	// Real frequency
	Real Reomg(Int n) const { return real(Re_z[n]); }
	Cmplx Rez(Int n) const { return Re_z[n]; }
};
