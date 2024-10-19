#pragma once

/*
code developed and maintained by (jmw@ruc.edu.cn, RUC, China) date 2022 - 2024
*/

#include "specs.h"
#include "prmtr.h"
#include "green.h"
#include "hyberr.h"

class Bath {
public:
	const MyMpi& mm;
	const Prmtr& p;				// parameters
	const Int npart;			// number of fitting parts
	VecInt ni;					// number of impurity orbital to be fitted in each part
	VecInt nb;					// number of bath orbital for fitting in each part
	UURandomSFMT uur;
	Vec<Vec<MatReal>> fvb; 		// fvb[nparts][the:0||eps:1][imp_idx][bath_idx]
	MatReal info;				// save the print out the NAV5(nmin, err, err_crv, err_reg, a_norm)
private:
	//Block diagonalising hybridisation functions and separating them
	VEC<ImGreen> generate_hb(const ImGreen& hb) const;
	VEC<ImGreen> generate_hb(const ImGreen& hb, const VecInt or_deg) const;
	
	Vec<MatReal> generate_bs();

	VecReal number_bath_fit_part(const ImGreen& vhb_i, const MatReal& vbs_i, Int idx);

	std::tuple<Real, VecReal, Int> bath_fit_number_contest(const VecReal& a0, const ImGreen& vhb_i, const MatReal& vbs_i, Int idx);

	VecReal next_initial_fitting_parameters(const VecReal& a0, const Int& ntry_fine, Int& itry);

	//using Box-Muller method to generate normal distribution random number
	VecReal perturb(const VecReal& a_i, Real scale) {
		VecReal a(a_i);
		for_Int(i, 0, a.size()) {
			Real u_1 = uur();
			Real u_2 = uur();
			Real z_0 = SQRT(-2 * LOG(u_1)) * cos(2 * pi_Real * u_2);
			//Real z_1 = SQRT(-2 * LOG(u_1)) * sin(2 * pi_Real * u_2);
			a[i] += scale * z_0;
		}
		return a;
	}
public:
	Bath(const MyMpi& mm_i, const Prmtr& prmtr_i);

	void number_bath_fit(const ImGreen& hb_i, const VecInt or_deg);

	void fvb_init();

	VecReal fvb2arr(const Int idx);

	Vec<Vec<MatReal>> arr2fvb(const Vec<VecReal>& arr);
	
	// Roughly making the impurity model symmetric
	void fix_fvb_symmetry();

	void bth_write_fvb(Int iter_cnt = -1, const Str bath_name = empty_str) const;

	void bth_read_fvb(const Str& file);
	
	MatReal find_hop() const;
};

