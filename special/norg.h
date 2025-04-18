#pragma once

/*
coded by Jia-Ming Wang (jmw@ruc.edu.cn, RUC, China) date 2021 - 2023
*/

#include "specs.h"
#include "prmtr.h"
#include "densitymatrix.h"
#include "nocspace.h"
#include "state.h"
#include "crrltfun.h"
#include "crrvec.h"
#include "krylov.h"
#include "impurity.h"


// impurity model
class NORG {
	// typedef VEC<VecInt> Tab;
	//*****************iteration***************
	Int iter_norg_cnt;						// count the NORG optimize iteration times.
	Int norg_stab_cnt;						// count the NORG stability iteration times.
	VecReal occnum_pre;						// occupation number record before one optimize iteration.
	VecReal occnum_lst;						// occupation number record after one optimize iteration.
	Real	groune_pre;						// ground state energy record before one optimize iteration.
	Real energy_err;						// ground state energy Error: error between the two NORG iterations
	Real occupation_err;					// occupation Error			: error between the two NORG iterations
	Real correctionerr;						// correction Error			: error between two correction vector modify
	mutable VecCmplx checkg_pre;			// correction vec check point before one optimize iteration.
	mutable VecCmplx checkg_lst;			// correction vec check point after one optimize iteration.

public:
	const MyMpi &mm;						// parameters
	const Prmtr &p;							// parameters
	VEC<MatReal> uormat;					// unitary_orbital_rotation_matrix

	Impdata impH;
	Real	groune_lst;						// ground state energy record after one optimize iteration.
	MatReal&  final_ground_state;
	VecReal  occnum;
	ImGreen impgreen;

	mutable NocSpace scsp;					// The main space.
	mutable DensityMat oneedm;				// The main space's density matrix.
	
	VecReal nparticle;						// The impurity's particle number.
private:
	/*// ! abandon
	// It assume that we already have the hopint from the Impurity class, but still have not rotated it yet.
	VecReal set_row_primeter_byimpH(const VEC<MatReal>& uormat_i, const MatReal& impH_i);

	// C^+_i C^+_j C_k C_l h_inter from [i][l][j][k] to [alpha][eta][beta][gamma]
	// It assume that we already have the hopint and h_inter from the Impurity class, with four Fermi interaction.
	VecReal set_row_primeter_byimpH(const VEC<MatReal>& uormat_i, const Impdata& impH_i);
	*/

	// Already have the hopint and h_inter from the Impurity class, with four Fermi interaction.
	// Speed up by Operator calss. C^+_i C^+_j C_k C_l h_inter from [i][l][j][k] to [alpha][eta][beta][gamma]
	// void set_row_primeter_byimpH(const VEC<MatReal>& uormat_i, const Impdata& impH_i, std::map<Int, Real>& oper_i);
	void set_row_primeter_byimpH(const VEC<MatReal>& uormat_i, const Impdata& impH_i, VecReal& oper_i);

	// only change for first norg_set of the first div.
	VecInt nppso(const VecInt &a, Int positon)
	{
		VecInt nppso_i(a);
		Int i = abs(positon);
		if (positon > 0) nppso_i[(i - 1) * 2] += 1;
		if (positon < 0) nppso_i[(i - 1) * 2] -= 1;
		return nppso_i;
	}

// ! NOT pass the test.
	VecInt nppso_all(const VecInt &a, Int positon)
	{
		VecInt nppso_i(a);
		Int i = abs(positon);
		if (positon > 0) nppso_i[(i - 1)] += 1;
		if (positon < 0) nppso_i[(i - 1)] -= 1;
		return nppso_i;
	}

	bool converged();

	bool green_modify_converged() const;
	bool green_modify_converged(Real correctionerr_i) const;

	void readmatrix(MatReal& m, const Str& file);
/*
	// using the correct vector to modify the shortcut space.
	void upgrade_space(NocSpace& scsp_i, NocSpace& scsp_ipl, VecReal& state_pl, NocSpace& scsp_imi, VecReal& state_mi);
*/

	VEC<MatReal> uormat_initialize();

	void show_the_nozero_number_of_tabel();
/*

	void update_by_correct_vec() const;

	void update_by_correct_vec(const Int orbital, const VecReal& ground_state_temp, const Real& ge,  const VEC<MatReal>& uormat_i) const;

	void rotation_the_space(const VEC<MatReal>& uormat_i) const;
*/
public:
	NORG(const MyMpi& mm_i, const Prmtr& prmtr_i);
	NORG(const MyMpi& mm_i, const Prmtr& prmtr_i, const Tab& table);
	// NORG(const MyMpi& mm_i, const Prmtr& prmtr_i, VecInt nparticals);
	NORG(const MyMpi& mm_i, const Prmtr& prmtr_i, Str tab_name);

	// NORG(const MyMpi& mm_i, const Impurity& imp_i, const Prmtr& prmtr_i);
	
	// for multi-number space.
	void up_date_h0_to_solve(const Impdata& impH, const VecReal sub_energy = Vec<Real>()); 

	void up_date_h0_to_solve(const Impdata& impH, const Int mode);
	

/*
	// In this method we calculate the correction vector using the Krylov-space approach to modify the U.
	void modify_by_krylov_correct_vec();

	// (Deprecated) In this method we calculate the correction vector by the conjugategradient method.
	void modify_by_correct_vector() ;

	// (Deprecated)
	void writematrix(const MatReal& m, const Str U_name, Int iter_norg_cnt) const;
*/

//---------------------------------------calculate the physical operator---------------------------------

	Real sz_imp_sz_bath(const Int imp_postition, const VecReal& vgs_i);

	// To check the NO-interactions check.
	Real return_nointeractions_ground_state_energy(const MatReal& h0_i) const;
//--------------------------------------- for the Green function---------------------------------
	
	// void get_gimp_by_krylov_CV_modify(Green& imp_i) const;
	// // Only use for test the validity.
	// void get_gimp_by_krylov_CV_modify(ReGreen& imp_i) const;

	// Int get_gimp_with_possible_degeneracy(ImGreen& imp_i, Int iter_cont = 999);

	// void get_gimp_by_krylov(const ImGreen& imp_i);

	// void find_g(Green &g) ;
	
	// Int find_g_with_possible_degeneracy(VEC<ReGreen> &g, Int kind = 0) ;

	// give the impurity green from krylov correction vec.
	void get_g_by_KCV(Green& imp_i);
	// give the impurity green matrix at onece from krylov correction vec.
	void get_g_by_KCV_spup(Green& imp_i);

	// git the impurity green from continue fraction.
	void get_g_by_CF(Green& imp_i);

	void get_gimp(Green& imp_i);

	void get_gimp(Green& imp_i, VecInt or_deg);

	void get_gimp_eigpairs(Green& imp_i);

	void get_gimp_eigpairs(Green& imp_i, VecInt or_deg);

	void get_gimp_all(Green& imp_i);

	void find_g0(Green& imp_i);

	// To calc the excitation spectrum for quasiparticle of holon-doublon bound states.
	void get_gimp_hdQPs(Green& imp_i);
//--------------------------------------- for the io---------------------------------
	void write_norg_info(Int iter_cnt) const;
	
	//void write_occupation_info() const;// ! Abandon

	void write_state_info(Int iter_cnt) const;

	MatReal see_MatReal(const VEC<MatReal>& uormat_i) const {
		MatReal transform_uormat(dmat(p.norbit, 1.));
		Int counter(0);
		for (const auto& uormat_ii : uormat_i) {
			for_Int(i, 0, uormat_ii.nrows()) {
				for_Int(j, 0, uormat_ii.ncols()) {
					transform_uormat[i + counter][j + counter] = uormat_ii[i][j];
				}
			}
			counter += uormat_ii.nrows();
		}
		return transform_uormat;
	}

	VecReal write_impurtiy_occupation(Int iter_cnt = -1, const Str& phy_name = empty_str) const;

	void PIO_occweight(VecReal occnum_lst) const {
		MatReal occnum, occweight;
		occnum = occnum_lst.mat(p.norg_sets, p.n_rot_orb/p.norg_sets);
		occweight = occnum;
		for_Int(i, 0, p.norg_sets) for_Int(j, 0, occnum.ncols()) occweight[i][j] = MIN(occweight[i][j],1 - occnum[i][j]) < 1e-14 ? 0 : MIN(occweight[i][j],1 - occnum[i][j]);
		Str nppso = scsp.nppso_str();
		if (mm) PIO(NAV4(p.if_norg_imp, nppso, occnum, occweight));
	}

	// print the Double Occupancy matrix
	MatReal print_DO(const DensityMat& dm) const {
		MatReal docc = dm.find_imp_double_occupancy();

		// if (mm) PIO(NAV1(docc));
		return docc;
	}

	// For speed up reason, we need to save(read) the nature transform matrix(NTR).
	bool check_NTR() {
		IFS checker("ru" + scsp.nppso_str() + ".bi");
		return checker.good();
	}
	MatReal save_NTR();
	MatReal read_NTR();
};

