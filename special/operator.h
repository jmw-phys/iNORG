#pragma once

/*
coded by Jia-Ming Wang (jmw@ruc.edu.cn, RUC, China) date 2022 - 2023
*/

#include "specs.h"
#include "prmtr.h"
#include "nocspace.h"
#include "state.h"

#include <iostream>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <mpi.h>

// // impurity model
// typedef Vec<VEC<Int>> Tab;  REPLACE by VEC<Int> @2024-05-17
// typedef VEC<Int> Tab;  //! REPLACE by VEC<Int> @2024-05-17

// ! not finished.（may not need）
// class TAB { }; 



class Operator {
public:
	const MyMpi& mm;				// parameters
	const Prmtr& p;					// parameters
	const NocSpace& scsp;			// NocSpace
	const Idx dim;
    // const VecReal& coefficient;      // coefficient for all the H's terms

	VEC<Int> table;						// The movement by femi operator contain.
	Real groundstate_energy;		// The ground state energy on this shortcut restratin.

	
	MatReal ground_state;			// if degeneracy happened, we only get one of them.
	// std::map<Int, Real> oper_value;	// The Hamiltonian operator in the shortcut space.
	VecReal oper_value;				// The Hamiltonian operator in the shortcut space.
private:

	template<typename T>
	constexpr const T& myclamp(const T& v, const T& lo, const T& hi) {
		return (v < lo) ? lo : ((hi < v) ? hi : v);
	}

public:
	// Expand the Shortcut space under Number of particles(NumberSpa).
	Operator(const MyMpi& mm_i, const Prmtr& prmtr_i, const NocSpace& s_i);
	// Operator(const MyMpi& mm_i, const Prmtr& prmtr_i, const Tab &tab, const NocSpace& s_i);
	Operator(const MyMpi& mm_i, const Prmtr& prmtr_i, const NocSpace& s_i, Str tab_name);

	// [i][0]row Idx; [i][1]colum Idx;[i][2] for the position for the hopint.
	Tab find_fullH_idx();

	void write_the_Tab(const std::string& filepath)  const;

	VEC<int> read_the_Tab(const std::string& filepath) const;

	void clear_TAB() { VEC<int>().swap(table); }

	// By using the find_h_idx() to find the idx first then speed up the procedure.
	SparseMatReal find_hmlt(const Tab h_idx) const;


	// Using Lanczos algorithm to find several lowest eigenpairs of a Hermitian matrix
	MatReal lowest_eigpairs(const Idx n, bool if_need_fast = true, Int wish_nev = 1);

	//Compute the Single particle excited state(sn_prtcl_ex_state) 
	//here only support for the div[0] 
	//$$\left|f_{0}\right\rangle=A\left|\psi_{0}\right\rangle$$
	VecReal sn_prtcl_ex_state(const Int ex_orbital_position, const VecReal ground_state, const Int crtann)const;

	// void save_the_Tab(Tab& tab, Str name) const;

	// Tab read_the_Tab(Str name) const;

	// 1 Ln: N; 2 Ln: Sz; 3 Ln: P;
	MatReal local_multiplets_state(const MatReal& g_state) const;
};