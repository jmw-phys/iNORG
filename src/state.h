#pragma once

/*
code developed and maintained by (jmw@ruc.edu.cn, RUC, China) date 2022 - 2024
*/

#include "specs.h"
#include "prmtr.h"
#include "nocspace.h"	// cause now This class was only for the NORG imp solver, at the Shortcut space.

/// <summary>
/// In this class will find all the possible state(StateStatistics) for one gived Combined division (ComDiv).
/// impurity model
/// At the Shortcut space(NocSpace) we set first of two divisons as impurity and the active orbital(bath site).
/// For now it's support for find two fermion act on one state.
/// </summary>
class StateStatistics 
{
	typedef VEC<std::tuple<Int, Int, MatInt>> hopdata;
	typedef VEC<std::pair<std::array<int, 4>, MatInt>> furfrm;
private:
	static const BnmlFast bf;
public:
	const NocSpace& space;								// parameters
	const Int div_idx;										// The divison idx at all Combined division space.
	const Int idx;											// The idx in one divison.
	const MatInt occ_n;										// occupy number in each division.
	// const MatInt sit_mat;								// spin - orbits number in each division.
	// const VecInt sit_n;									// spin - orbits number in each division.
	const ComDivs cfg;										// The configuration of the state.
	// const hopdata hopup;		// For the upper spin-orbits, a set of occupy number in each division
	// const hopdata hopdw;		// For the down spin-orbits, a set of occupy number in each division
	//ComDivs cfgdw;										// The configuration of the state.

	// const VecOnb cfgup;										// Here only for better reading reason.
	// const VecOnb cfgdw;										// Here only for better reading reason.

	// for the electronic statistics of the state.
	Vec< VEC<Int> > filled_spinless;

	// VEC<VecInt> off_diagonal_term;	// a set of valuable term with the annihilation, creation orbit's positon, and Colum idx(H_i).
private:
	ComDivs cfgdetail();
	// To statistic the occupation situation and filling the filled_spinless matrix.
	void statistics();



	// fermion anticommutativity, for two orbitals changed.
	//Int anticommu(const Int& pi,const Int& pj);

	// find the div's idx according to the "newdiv".
	Int find_newdiv_idx(MatInt newdiv);

	template<typename T>
	constexpr const T& myclamp(const T& v, const T& lo, const T& hi) {
		return (v < lo) ? lo : ((hi < v) ? hi : v);
	}

public:
	//state(Int idxD, VecInt bases_i) :idx(idxD), bases(bases_i) {};
	StateStatistics(const Int& h_i, const Int& comdiv_idx, const NocSpace& s_i);

	// The out put<0>[i]: the annihilation divsion position; <1>[i]: the creation dision positon; <2>[i]: The new divs.
	hopdata divocchop_ingroup(const Int& ComDiv, Idx sets_n);
	
	// [0]~[3] i-j-k-l orbit's position; (*this).second: The new divs.
	furfrm divs_change_fourFermi(const Int& ComDiv, Idx sets_i, Idx sets_j);

	// The out put<0>[i]: the annihilation divsion position; <1>[i]: the creation dision positon; <2>[i]: The new divs.
	VEC<MatInt> interation_soc_hop(const Int& ComDiv);

	// Output all off-diagonal H matrix term in one row.(For all sites hoping?)
	// [i][0]:annihilation orbit's position; [i][1]:creation orbit's positon; [i][2]:Colum idx(i); [i][3]:sign(fermion anticommutativity)
	// PS.:The matrix Max size: (N*N-(N-N_{imp})-N_{imp}*N_{imp},2), and it's spinless.
	VEC<VecInt> find_off_diagonal_term(const hopdata &hopup, const hopdata &hopdw);

	// Output all off-diagonal H0 matrix term with one spinless orbit.
	// [i][0]:annihilation orbit's position; [i][1]:creation orbit's positon; [i][2]:Colum idx(i); [i][3]:sign(fermion anticommutativity)
	VEC<VecInt> find_each_spiless_group_off_diagonal_term(const hopdata &hopspi, const Int sets_n);

	// Output all off-diagonal H0 matrix and all H_{ijkl} term with one spinless orbit.
	//  [0]~[3] i-j-k-l orbit's position; [4]:Colum idx(i);[5]:sign(fermion anticommutativity)
	VEC<std::array<int, 6>> find_off_diagonal_term_fourFermi(const furfrm &furfrm, const Int sets_i, const Int sets_j);
	
	// Output all off-diagonal SOC spinless orb term.
	VEC<VecInt> off_diagonal_soc_term(const VEC<MatInt> &hop_soc);

	Mat<Char> show_cfg();

	VecBool cfg2vecbool();
	
	VecBool cfgex2vecbool(Int ex_pos) ;
};
