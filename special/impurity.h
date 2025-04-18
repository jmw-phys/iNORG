#pragma once

/*
coded by Rong-Qiang He (rqhe@ruc.edu.cn, RUC, China)
date 2021-02-25
*/

#include "specs.h"
#include "prmtr.h"
#include "green.h"
#include "bath.h"
#include "model.h"

// a basis state (i.e., a configuration) is represented by an integer.
// the lower 16 bits represent a spin-down configuration
// the higher 16 bits represent a spin-up configuration

// tensor V_{ijkl} = i * n^3 + j * n^2 + k * n + l
typedef std::pair<MatReal, VecReal> Impdata;
class Impurity {
private:
	const MyMpi& mm;			// parameters
	const Prmtr& p;				// parameters
	const Bath& bth;

    const Int ni;				//number of impurity
	const Int ns;				// number of sites,ns=ni+nb
	const Int nb;				// number of bath sites
	MatReal h0;					// hopping factors
	VecReal pos_imp;			// position of imp site
	VecReal imp_lvl;			// impurity energy level
public:
	Impdata impH;				// The impurity H construction data

private:
	//hopping  factors;when bath parameters is unusual,we just need to modify bath.hop() 
	void set_factor();
	//interaction  factors
	VecReal set_interaction();
	//---------------------------------------special for the hhd function(arXiv:2209.14178v1)------------------------------------------------
	VecReal set_3band_interaction_withalpha();
	void modify_Impdata_for_half_fill_hhd(Impdata& impH_i);
	//---------------------------------------------------------------------------------------------------------------------------------------

public:
	Impurity(const MyMpi& mm_i, const Prmtr& prmtr_i, const Bath& bth_i, const Str& file = empty_str);
	Impurity(const MyMpi& mm_i, const Prmtr& prmtr_i, const Bath& bth_i, const VecInt or_deg);
	void find_g0(Green& g0) const;
	void find_all_g0(Green& g0) const;
	
	void find_hb(Green& hb) const;
	
	MatReal find_hop_for_test() const;

	void update(Str mode = empty_str);
	
	// From U(n_up n_dn) to the U(n_up - 0.5)(n_dw - 0.5)	@23.04.19
	void modify_Impdata_for_half_fill(Impdata& impH_i);

	// void write_H0info(const Bath &b) const;
	void write_H0info(const Bath &b, Int ndeg = -1, Int iter_cnt = -1) const;

	// void save() const;

	// void read(const Str& file) const;
};
