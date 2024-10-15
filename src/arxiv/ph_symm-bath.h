#pragma once

/*
coded by Rong-Qiang He (rqhe@ruc.edu.cn, RUC, China) date 2021-02-19
code developed and maintained by (jmw@ruc.edu.cn, RUC, China) date 2023
*/

#include "specs.h"
#include "prmtr.h"
#include "green.h"
#include "ph-symm_hyberr.h"

class Bath {
public:
	const MyMpi& mm;
	const Prmtr& p;			// parameters
	Int nb;					// number of bath sites
	UURandomSFMT uur;
	VecReal ose;			// on-site energies for bath sites
	VecReal hop;			// hopping coefficients between impurity site and bath sites

	VEC<VecReal> vec_ose,vec_hop;
	MatReal info;			// save the print out the NAV5(nmin, err, err_crv, err_reg, /*err_bsr,*/ a_norm)

private:
	void regularize_ose_hop() {
		slctsort(ose, hop);
		// after a unitary transformation, hop can always be 
		// non-negative when there is only one impurity site
		hop = ABS(hop);
	}
	
	void init_ose_hop() {
		uur(ose);
		ose -= 0.5;
		ose *= 4.;
		uur(hop);
		hop -= 0.5;
		hop *= 4.;
		regularize_ose_hop();
	}

	VecReal next_initial_fitting_parameters(const VecReal& a0, const Int& ntry_fine, Int& itry);
	
	void perturb(VecReal& a, Real amplitude, Real exponent) {
		for_Int(i, 0, a.size()) {
			Real ra = uur() - 0.5;
			Real re = uur() - 0.5;
			a[i] = a[i] * (1. + amplitude * ra) + SIGN(std::pow(10., (8 * ABS(re) - exponent)), a[i]);
		}
	}

public:
	Bath(const MyMpi& mm_i, const Prmtr& prmtr_i);
	void number_bath_fit(const ImGreen& hb_i, Int iter, Int mode = 1);
	
	void init_vec_ose_hop();

	void write_ose_hop(Int iter_cnt = -1, const Str& bath_name = empty_str) const;
	void read_ose_hop();

	std::tuple<Real, VecReal, Int> bath_fit_number_contest(const VecReal& a0, Int nb, const ImGreen&hb_i, Int orb_i,Int mode);
	MatReal find_hop() const;


	void bath_fit(const ImGreen& hb_i, Int iter);
};