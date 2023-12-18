/*
coded by Jia-Ming Wang (jmw@ruc.edu.cn, RUC, China) date 2022 - 2023
*/

#include "nocspace.h"
#define nip p.norg_sets
#define i2b p.nI2B
#define nob p.norbit

/*
// the NocSpace, is start frome the shortcutspace.
NocSpace::NocSpace(const MyMpi& mm_i, const Prmtr& prmtr_i, const Int& NumberSpa) :
	mm(mm_i), p(prmtr_i), ndivs(prmtr_i.ndiv), nspa(NumberSpa), sit_mat(p.norg_sets, prmtr_i.ndiv, 0),
	hopint(nob, nob, 0.),
	control_divs(p.norg_sets + 1, prmtr_i.ndiv, 0), orbital_divcnt(prmtr_i.ndiv, 0), dim(0)
{
	set_control();
	find_combined_number_subspaces();
	//WRN("Finish finding the subspace and the Dim is:" + NAV(dim));
}
*/

NocSpace::NocSpace(const MyMpi& mm_i, const Prmtr& prmtr_i, const Int& NumberSpa) :
 	mm(mm_i), p(prmtr_i), ndivs(prmtr_i.ndiv), nspa(NumberSpa), sit_mat(p.norg_sets, prmtr_i.ndiv, 0),
	hopint(nob, nob, 0.), coefficient(int(std::pow(nob,4) + hopint.size() + 1), 0.),
 	control_divs(p.norg_sets + 1, prmtr_i.ndiv, 0), orbital_divcnt(prmtr_i.ndiv, 0), dim(0)
{
 	set_control();
 	find_combined_number_subspaces();
}

// nppso mean: number of partical per spin orbital.
NocSpace::NocSpace(const MyMpi& mm_i, const Prmtr& prmtr_i, const VecInt& nppso_i) :
	mm(mm_i), p(prmtr_i), ndivs(prmtr_i.ndiv), nspa(SUM(nppso_i)), sit_mat(p.norg_sets, prmtr_i.ndiv, 0),
	hopint(nob, nob, 0.), coefficient(int(std::pow(nob,4) + hopint.size() + 1), 0.), nppso(nppso_i),
	control_divs(p.norg_sets + 1, prmtr_i.ndiv, 0), orbital_divcnt(prmtr_i.ndiv, 0), dim(0)
{
	set_control();
	// find_combined_number_subspaces(1);
	// find_all_noc_subspaces();
	// find_all_noc_subspaces_by_row();

	// find_all_noc_subspaces();// ! tested for the YBCO-CDMFT version.
	find_thought_noc_subspaces();
	

	if(mm) PIO("Begin find_combined_number_subspaces()"+ NAV(dim)+"   "+present());
}

// nppso mean: number of partical per spin orbital.
NocSpace::NocSpace(const MyMpi& mm_i, const Prmtr& prmtr_i, const VecInt& nppso_i, Str tab_name) :
	mm(mm_i), p(prmtr_i), ndivs(prmtr_i.ndiv), nspa(SUM(nppso_i)), sit_mat(p.norg_sets, prmtr_i.ndiv, 0),
	hopint(nob, nob, 0.), coefficient(int(std::pow(nob,4) + hopint.size() + 1), 0.), nppso(nppso_i),
	control_divs(p.norg_sets + 1, prmtr_i.ndiv, 0), orbital_divcnt(prmtr_i.ndiv, 0), dim(read_the_Tab(tab_name))
{
	set_control();
	if(mm) WRN("Begin find_combined_number_subspaces()"+ NAV(dim));
}

// using the table to construce the NOC space.
NocSpace::NocSpace(const MyMpi& mm_i, const Prmtr& prmtr_i, const VecInt& nppso_i, const Tab& tab) :
	mm(mm_i), p(prmtr_i), ndivs(prmtr_i.ndiv), nspa(SUM(nppso_i)), sit_mat(p.norg_sets, prmtr_i.ndiv, 0),
	hopint(nob, nob, 0.), coefficient(int(std::pow(nob,4) + hopint.size() + 1), 0.), nppso(nppso_i),
	control_divs(p.norg_sets + 1, prmtr_i.ndiv, 0), orbital_divcnt(prmtr_i.ndiv, 0), dim(tab[0].size())
{
	set_control();
	if(mm) WRN("Begin find_combined_number_subspaces()"+ NAV(dim));
}

//----------------------------- judge function -----------------------------

// here when return true, mean satisfy the PHSs rule.
bool NocSpace::check_correlated_column(const Int& col_pos, const VecInt& div_colsum) const {
	Int change_nhole(0), change_nelec(0);
	if (p.if_norg_imp) {
		change_nhole = orbital_divcnt[col_pos] - div_colsum[col_pos];
		change_nelec = div_colsum[div_colsum.size() - col_pos - 1];
		if (change_nhole + change_nelec <= control_divs[0][div_colsum.size() - col_pos - 1]) return true;
	}
	else {
		change_nhole = orbital_divcnt[col_pos] - div_colsum[col_pos];
		change_nelec = div_colsum[div_colsum.size() - col_pos];
		if (change_nhole + change_nelec <= control_divs[0][div_colsum.size() - col_pos]) return true;
	}
	return false;
}


// here when return true, mean satisfy the PHSs rule.
bool NocSpace::check_if_PHSs(const VecInt& div_colsum) const {
	/*
	Int change_nhole(0), change_nelec(0);
	if (p.if_norg_imp) {
		for_Int(col_idx, col_pos, div_colsum.size() / 2) {
			change_nhole += orbital_divcnt[col_idx] - div_colsum[col_idx];
			change_nelec += div_colsum[div_colsum.size() - 1 - col_idx];
		}
		if (change_nhole + change_nelec <= control_divs[0][div_colsum.size() - col_pos - 1]) return true;
	}
	else {
		for_Int(col_idx, col_pos, div_colsum.size() / 2) {
			change_nhole += orbital_divcnt[col_idx] - div_colsum[col_idx];
			change_nelec += div_colsum[div_colsum.size() - col_idx];
		}
		if (change_nhole + change_nelec <= control_divs[0][div_colsum.size() - col_pos]) return true;
	}
	// if (mm) WRN(NAV2(col_pos, div_colsum));
	// if (mm) WRN(NAV2(change_nhole, change_nelec));
	*/
	bool thought_NOOC(true);
	if (p.if_norg_imp) {
		ERR("Since this part was not in used, FOR here I have not test!");
		return true;
	}
	else {
		// In this part of code we consider the PHSs rule with out the "*" terms.
		// VecInt sum_hole_chages(div_colsum.size() / 2 - 1, 0), sum_elec_chages(div_colsum.size() / 2 - 1, 0);
		VecInt temp_hole_chages(div_colsum.size() / 2 - 1, 0), temp_elec_chages(div_colsum.size() / 2 - 1, 0);

		Int length = div_colsum.size() / 2 - 1;
		for_Int(col_idx, 0, length) {
			temp_hole_chages[col_idx] = orbital_divcnt[col_idx + 1] - div_colsum[col_idx + 1];
			temp_elec_chages[col_idx] = div_colsum[div_colsum.size() - col_idx - 1];
		}
		// for_Int(col_idx, 0, length) {
		// 	sum_hole_chages[col_idx] = SUM_0toX(temp_hole_chages.reverse(), col_idx);
		// 	sum_elec_chages[col_idx] = SUM_0toX(temp_elec_chages.reverse(), col_idx);
		// }
		// if((SUM(temp_hole_chages) + SUM(temp_elec_chages)) == 1) thought_NOOC = true;
		if((SUM(temp_hole_chages) + SUM(temp_elec_chages)) <= control_divs[0][div_colsum.size() - 2]) return true; // Only suit for the 8 divisions.

		for_Int(col_pos, 0, length) if((temp_hole_chages[col_pos] + temp_elec_chages[col_pos]) > control_divs[0][div_colsum.size() - col_pos - 1]) return false;

		for_Int(col_pos, 0, length)
			if (SUM(temp_hole_chages.truncate(col_pos, length)) + SUM(temp_elec_chages.truncate(col_pos, length)) > 0 \
				&& (SUM_0toX(temp_hole_chages, col_pos) != 0 || SUM_0toX(temp_elec_chages, col_pos) != 0)) return false;
		
/* // way two
		MatInt temp_div_orbital_elec_and_hole = orbital_divcnt.mat(2,orbital_divcnt.size() / 2).tr().truncate_row(1,orbital_divcnt.size() / 2).tr();
		VecInt sum_div_orbital_elec_and_hole(temp_div_orbital_elec_and_hole.ncols(), 0);
		for_Int(i, 0, length) sum_div_orbital_elec_and_hole[i] = temp_div_orbital_elec_and_hole[0][i] + temp_div_orbital_elec_and_hole[1][i];
		for_Int(col_pos, 0, length)
			if (SUM_0toX(sum_div_orbital_elec_and_hole.reverse(), col_pos) > SUM_0toX(sum_div_orbital_elec_and_hole, col_pos+1) \
				&& sum_hole_chages[col_pos] + sum_elec_chages[col_pos] < SUM_0toX(temp_hole_chages, col_pos+1) + SUM_0toX(temp_elec_chages, col_pos+1)) return false;
		// for_Int(col_pos, 0, length)
		// 	if (sum_hole_chages[col_pos] + sum_elec_chages[col_pos] <= control_divs[0][div_colsum.size() - col_pos - 1] \
		// 		&& SUM_0toX(temp_hole_chages, col_pos) == 0 && SUM_0toX(temp_elec_chages, col_pos) == 0) return true;
*/
	}
	return thought_NOOC;
}


Vec<MatInt> NocSpace::multi_judger_with_return(VEC<VEC<int> >& s, const VEC<VEC<int> >& a) const {
	VEC<MatInt> split_spilss_divs;
	
	if (p.if_norg_imp) {
		VecPartition div_dim(mm.np(), mm.id(), s.size()); VecInt splited_judge(div_dim.len(), 0), judge(s.size());
		for_Int(i, div_dim.bgn(), div_dim.end()) {
			MatInt spilss_div(read_from_col_lable(s[i], a).mat(control_divs.ncols(), control_divs.nrows() - 1).tr());
			if (ifin_NocSpace(spilss_div, nppso)) splited_judge[i - div_dim.bgn()] = 1;
		}
		judge = mm.Allgatherv(splited_judge, div_dim);

		for_Int(i, 0, s.size()) if(judge[i]) {
			MatInt spilss_div(read_from_col_lable(s[i], a).mat(control_divs.ncols(), control_divs.nrows() - 1).tr());
			if (ifin_NocSpace(spilss_div, nppso)) split_spilss_divs.push_back(spilss_div);
		}
	} else for_Int(i, 0, s.size()) {
		MatInt spilss_div(read_from_col_lable(s[i], a).mat(control_divs.ncols() - 2, control_divs.nrows() - 1).tr());
		if (suit_NOOC(spilss_div, nppso)) { //!! here need to chage the restraint.
			// if (mm) WRN(NAV(spilss_div));
			VecInt left_div_size(p.norg_sets), left_n(p.norg_sets);
			VEC<VEC<VecInt>> left_ii(p.norg_sets);

			// for the two "*" terms.
			for_Int(j, 0, p.norg_sets) {
				Int n(nppso[j] - SUM(spilss_div[j])), sit_i(sit_mat[j][0]), sit_b(sit_mat[j][ndivs / 2]); left_n[j] = n;
				Int possible = n < MIN(sit_i, sit_b) ? n+1 : n+1 - MAX(n-sit_b,0) - MAX(n-sit_i,0);	left_div_size[j] =  possible;
				VecInt a_beg(2, 0);
				a_beg = possible == n+1 ? Vec{0, n}:  Vec{n - sit_b, sit_b};
				for_Int(k, 0, possible) {left_ii[j].push_back(a_beg); a_beg[0]++; a_beg[1]--;}
				// if(mm) WRN(NAV(left_ii[j].size()));
			}
			
			if (!p.if_norg_imp) for (const auto& x : cart_product(left_ii)) {
				MatInt unnooc(p.norg_sets, 2);
				for_Int(j, 0, x.size()) unnooc[j] = x[j];
				MatInt div(concat(concat(unnooc.tr()[0],spilss_div.tr().truncate_row(0,spilss_div.ncols()/2).vec()), \
						   concat(unnooc.tr()[1],spilss_div.tr().truncate_row(spilss_div.ncols()/2,spilss_div.ncols()).vec())).\
						   mat(sit_mat.ncols(),sit_mat.nrows()).tr());				
				split_spilss_divs.push_back(div);
			}
		}
	}

	Vec<MatInt> out = Vec(split_spilss_divs);
	return out;
}

VecInt NocSpace::multi_judger(const VEC<VEC<int> >& s, const VEC<VEC<int> >& a) const{
	VecPartition div_dim(mm.np(), mm.id(), s.size());
	// MatInt div_dim()
	VecInt splited_out(div_dim.len(), 0);
	// VecInt temp_splited_out()
	for_Int(i, div_dim.bgn(), div_dim.end()){
		// VEC<Int> x(s[i]);
		// VecInt spilss_div_v(read_from_col_lable(x,a));
		// MatInt spilss_div_tr(spilss_div_v.mat(control_divs.ncols(), control_divs.nrows() - 1));
		// MatInt spilss_div(spilss_div_tr.tr());
		MatInt spilss_div(read_from_col_lable(s[i],a).mat(control_divs.ncols(), control_divs.nrows() - 1).tr());
		if (ifin_NocSpace(spilss_div, nppso)) splited_out[i-div_dim.bgn()] = 1;
		// if (ifin_NocSpace_judge_by_nppso(spilss_div, nppso)) splited_out[i-div_dim.bgn()] = 1;
	}
	VecInt out(s.size());
	out = mm.Allgatherv(splited_out, div_dim);
	return out;
}

VecInt NocSpace::multi_judger_by_row(const VEC<VEC<int> >& s, const VEC<VEC<int> >& a) const{
	VecPartition div_dim(mm.np(), mm.id(), s.size());
	VecInt splited_out(div_dim.len(), 0);
	for_Int(i, div_dim.bgn(), div_dim.end()){
		MatInt spilss_div(read_from_col_lable(s[i],a).mat( control_divs.nrows() - 1, control_divs.ncols()));
		if (ifin_NocSpace(spilss_div, nppso)) splited_out[i-div_dim.bgn()] = 1;
	}
	VecInt out(s.size());
	out = mm.Allgatherv(splited_out, div_dim);
	return out;
}

bool NocSpace::ifin_NocSpace(MatInt& spilss_div) const
{
	if (SUM(spilss_div) == nspa) {
		for_Int(i, 1, control_divs.nrows()) for_Int(j, 0, ndivs) if (spilss_div[i - 1][j] < 0 || spilss_div[i - 1][j] > control_divs[i][j]) return false;
		VecInt spilsdiv_countvec(ndivs, 0);
		for_Int(j, 0, spilss_div.ncols())for_Int(i, 0, spilss_div.nrows()) spilsdiv_countvec[j] += spilss_div[i][j];
		if (spilsdiv_countvec[1] >= orbital_divcnt[1] + control_divs[0][1] && \
			spilsdiv_countvec[1] <= orbital_divcnt[1] && \
			spilsdiv_countvec[2] >= orbital_divcnt[2] + control_divs[0][2] && \
			spilsdiv_countvec[2] <= orbital_divcnt[2] && \
			spilsdiv_countvec[4] >= 0 && \
			spilsdiv_countvec[4] <= control_divs[0][4] && \
			spilsdiv_countvec[5] >= 0 && \
			spilsdiv_countvec[5] <= control_divs[0][5]) return true;
	}
	return false;
}

bool NocSpace::ifin_NocSpace(MatInt& spilss_div, const VecInt& nppso) const
{
	if(ndivs%2) ERR("ndivs is not a even number!");
	for_Int(i, 0, p.norg_sets)	if(SUM(spilss_div[i]) != nppso[i]) return false;
	if (SUM(spilss_div) == nspa) {
		for_Int(i, 1, control_divs.nrows()) for_Int(j, 0, ndivs) if (spilss_div[i - 1][j] < 0 || spilss_div[i - 1][j] > control_divs[i][j]) return false;
		VecInt spilsdiv_countvec(ndivs, 0);
		for_Int(j, 0, spilss_div.ncols()) for_Int(i, 0, spilss_div.nrows()) spilsdiv_countvec[j] += spilss_div[i][j];

		Int start = p.if_norg_imp ? 0 : 1;
		for_Int(i, start, ndivs) if(!check_each_column(i, spilsdiv_countvec)) return false;
		if(p.nooc_mode == STR("nooc")) return true;
		else if(p.nooc_mode == STR("cpnooc")) {for_Int(i, start, ndivs/2 - 1) if(!check_correlated_column(i, spilsdiv_countvec)) return false;}		// Here we set middle of one were "nooc", and others were "cnooc".
		else if(p.nooc_mode == STR("cnooc")) {for_Int(i, start, ndivs/2) if(!check_correlated_column(i, spilsdiv_countvec)) return false;}
		else if (p.nooc_mode == STR("phess")) {return check_if_PHSs(spilsdiv_countvec);}
		else ERR("nooc_mode in put was wrong!");
		return true;
	}
	else return false;
}

// In here we only consider the impurity not rotate.
bool NocSpace::suit_NOOC(MatInt& spilss_div, const VecInt& nppso) const
{
	if(ndivs%2) ERR("ndivs is not a even number!");
	for_Int(i, 0, p.norg_sets)	if( sit_mat[i][0] + sit_mat[i][ndivs/2] < nppso[i] - SUM(spilss_div[i]) || nppso[i] - SUM(spilss_div[i]) < 0) return false;

	// if(mm) WRN(NAV(spilss_div));
	MatInt tmp_cntrl_divs = concat(control_divs.tr().truncate_row(1, ndivs/2).vec(), \
	control_divs.tr().truncate_row(ndivs/2+1, ndivs).vec()).mat(spilss_div.ncols(),spilss_div.nrows()+1).tr();
	// if(mm) WRN(NAV(tmp_cntrl_divs));
	for_Int(i, 1, tmp_cntrl_divs.nrows()) for_Int(j, 0, ndivs - 2) if (spilss_div[i - 1][j] < 0 || spilss_div[i - 1][j] > tmp_cntrl_divs[i][j]) return false;
	VecInt spils_c(spilss_div.ncols(), 0), zero(1,0);
	for_Int(j, 0, spilss_div.ncols()) for_Int(i, 0, spilss_div.nrows()) spils_c[j] += spilss_div[i][j];
	VecInt countvec = concat(concat(zero,spils_c.truncate(0, spils_c.size()/2)),concat(zero,spils_c.truncate(spils_c.size()/2, spils_c.size())));
	// if(mm) WRN(NAV(countvec));
	Int countvec_length = countvec.size();
	for_Int(i, 1, countvec_length) if (!check_each_column(i, countvec)) return false;
	if (p.nooc_mode == STR("nooc")) return true;
	else if (p.nooc_mode == STR("cpnooc")) { for_Int(i, 1, countvec_length / 2 - 1) if (!check_correlated_column(i, countvec)) return false; }
	else if (p.nooc_mode == STR("cnooc")) { for_Int(i, 1, countvec_length / 2) if (!check_correlated_column(i, countvec)) return false; }
	else if (p.nooc_mode == STR("phess")) {return check_if_PHSs(countvec);}
	// else if (p.nooc_mode == STR("phess")) {if (mm && check_if_PHSs(countvec) == true) WRN(NAV2(spilss_div, countvec)); return check_if_PHSs(countvec);}
	// else ERR("nooc_mode in put was wrong!");
	
	// if (mm) WRN(NAV2(spilss_div, countvec));
	return true;
	// else ERR("nooc_mode in put was wrong!");
}

/*
bool NocSpace::ifin_NocSpace(const VecBool new_cfig, const VecInt& nppso) const
{
	MatInt spilss_div(sit_mat.nrows(), sit_mat.ncols(), 0);
	{//chnage new_cfig to spilss_div
		VecInt sit_vec(sit_mat.vec());
		for_Int(i, 0, sit_vec.size()) for_Int(j, SUM_0toX(sit_vec, i), SUM_0toX(sit_vec, i+1)) spilss_div.vec()[i] += new_cfig[j] ? 1 : 0;
	}


	if(ndivs%2) ERR("ndivs is not a even number!");
	for_Int(i, 0, p.norg_sets)	if( sit_mat[i][0] + sit_mat[i][ndivs/2] < nppso[i] - SUM(spilss_div[i]) || nppso[i] - SUM(spilss_div[i]) < 0) return false;

	// if(mm) WRN(NAV(spilss_div));
	MatInt tmp_cntrl_divs = concat(control_divs.tr().truncate_row(1, ndivs/2).vec(), \
	control_divs.tr().truncate_row(ndivs/2+1, ndivs).vec()).mat(spilss_div.ncols(),spilss_div.nrows()+1).tr();
	// if(mm) WRN(NAV(tmp_cntrl_divs));
	for_Int(i, 1, tmp_cntrl_divs.nrows()) for_Int(j, 0, ndivs - 2) if (spilss_div[i - 1][j] < 0 || spilss_div[i - 1][j] > tmp_cntrl_divs[i][j]) return false;
	VecInt spils_c(spilss_div.ncols(), 0), zero(1,0);
	for_Int(j, 0, spilss_div.ncols()) for_Int(i, 0, spilss_div.nrows()) spils_c[j] += spilss_div[i][j];
	VecInt countvec = concat(concat(zero,spils_c.truncate(0, spils_c.size()/2)),concat(zero,spils_c.truncate(spils_c.size()/2, spils_c.size())));
	// if(mm) WRN(NAV(countvec));
	for_Int(i, 1, ndivs) if(!check_each_column(i, countvec)) return false;
	if(p.nooc_mode == STR("nooc")) return true;
	else if(p.nooc_mode == STR("cpnooc")) {for_Int(i, 1, ndivs/2 - 1) if(!check_correlated_column(i, countvec)) return false;}
	else if(p.nooc_mode == STR("cnooc")) {for_Int(i, 1, ndivs/2) if(!check_correlated_column(i, countvec)) return false;}
	else if (p.nooc_mode == STR("phess")) {return check_if_PHSs(countvec);}
	else ERR("nooc_mode in put was wrong!");
	// else ERR("nooc_mode in put was wrong!");
}
*/

bool NocSpace::ifin_NocSpace_more_strict(MatInt& spilss_div, const VecInt& nppso) const
{
	for_Int(i, 0, p.norg_sets)	if(SUM(spilss_div[i]) != nppso[i]) return false;
	if (SUM(spilss_div) == nspa) {
		for_Int(i, 1, control_divs.nrows()) for_Int(j, 0, ndivs) if (spilss_div[i - 1][j] < 0 || spilss_div[i - 1][j] > control_divs[i][j]) return false;
		VecInt spilsdiv_countvec(ndivs, 0);
		for_Int(j, 0, spilss_div.ncols())for_Int(i, 0, spilss_div.nrows()) spilsdiv_countvec[j] += spilss_div[i][j];
		/*
		if (spilsdiv_countvec[1] >= orbital_divcnt[1] + control_divs[0][1] && \
			spilsdiv_countvec[1] <= orbital_divcnt[1] && \
			spilsdiv_countvec[2] >= orbital_divcnt[2] + control_divs[0][2] && \
			spilsdiv_countvec[2] <= orbital_divcnt[2] && \
			spilsdiv_countvec[4] >= 0 && \
			spilsdiv_countvec[4] <= control_divs[0][4] && \
			spilsdiv_countvec[5] >= 0 && \
			spilsdiv_countvec[5] <= control_divs[0][5] && \
			orbital_divcnt[1] - spilsdiv_countvec[1] + spilsdiv_countvec[5] <= control_divs[0][5]) return true;
		*/
		for_Int(i, 1, ndivs) if(check_each_column(i, spilsdiv_countvec)) return false;
	//  && orbital_divcnt[1] - spilsdiv_countvec[1] + spilsdiv_countvec[5] <= control_divs[0][5]
		return true;
	}
	else return false;
}

bool NocSpace::ifin_NocSpace_judge_by_nppso(const MatInt& spilss_div, const VecInt& nppso) const
{
	for_Int(i, 0, p.norg_sets)	if(SUM(spilss_div[i]) != nppso[i]) return false;

	return true;
}

bool NocSpace::if_div_in_restraint(const VecInt& restraint, const Int position, const Int max, const Int now) const
{
	if(p.if_norg_imp){
		if (position < ndivs / 2 && now >= (max + restraint[position])) return true;
		if (position >= ndivs / 2 && now <= restraint[position]) return true;
	} else {
		if (position == 0 || position == ndivs / 2) return true;
		if (position < ndivs / 2 && now >= (max + restraint[position])) return true;
		if (position > ndivs / 2 && now <= restraint[position]) return true;
	}
	return false;
}

bool NocSpace::if_col_divs_in_restraint(const Int& restraint, const VEC<Int>& divcol_i, Int col_idx) const
{
	VecInt divcol(divcol_i);	
	if(p.if_norg_imp){
		if(col_idx < ndivs / 2 && SUM(divcol) >= (orbital_divcnt[col_idx] + restraint)) return true;
		if(col_idx >= ndivs / 2 && SUM(divcol) <= restraint ) return true;
	} else {
		if (col_idx == 0 || col_idx == ndivs / 2) return true;
		if(col_idx < ndivs / 2 && SUM(divcol) >= (orbital_divcnt[col_idx] + restraint)) return true;
		if(col_idx > ndivs / 2 && SUM(divcol) <= restraint ) return true;
	}
	return false;
}

bool NocSpace::if_row_divs_in_restraint(const Int& restraint, const VEC<Int>& divrow_i, VecInt count_sit_mat) const
{
	if(SUM(Vec(divrow_i)) != restraint) return false;
	// VecInt countvec(sit_mat.ncols(), 0);
	for_Int(col_idx, 1, divrow_i.size()){
		// if(!(col_idx < ndivs / 2 && col_idx >= (divrow_i[col_idx] + control_divs[0][col_idx]))) return false;
		// if(!(col_idx > ndivs / 2 && col_idx <= control_divs[0][col_idx] )) return false;
		// if (col_pos == ndivs / 2) return true;
		if (col_idx < ndivs / 2 && divrow_i[col_idx] < count_sit_mat[col_idx] + control_divs[0][col_idx] ) return false;
		if (col_idx > ndivs / 2 && divrow_i[col_idx] > control_divs[0][col_idx]) return false;
	}
	return true;
}


//----------------------------- space function before judge -----------------------------


void NocSpace::find_all_possible_state_by_col(VEC<VEC<Int> >& a, VEC<VEC<Int> >& s) const
{
	VEC<VEC<Int> > a_lable;
	Int counter(0);
	for_Int(col, 0, control_divs.ncols())
	{
		VEC<VEC<Int>> temp_a, a_rol_temp;
		for_Int(row, 1, control_divs.nrows())
		{
			VEC<Int> one_div;
			for_Int(spl, 0, control_divs[row][col] + 1)
			{
				if (if_div_in_restraint(control_divs[0], col, control_divs[row][col], spl))
					one_div.push_back(spl);
			}
			a_rol_temp.push_back(one_div);
		}
		temp_a = cart_product(a_rol_temp);
		for (const auto &one_divs : temp_a)	{
			if (if_col_divs_in_restraint(control_divs[0][col], one_divs, col))
				a.push_back(one_divs);
		}
		VEC<Int> a_lable_i;
		for_Int(i, counter, a.size()){
			a_lable_i.push_back(i);
			counter++;
		}
		a_lable.push_back(a_lable_i);
		// if(mm) WRN(NAV2( col,a.size()));
	}

	s = cart_product_monitor_col(a_lable, a);
	// s = cart_product(a_lable);
}

void NocSpace::find_all_possible_state_by_row(VEC<VEC<Int> >& a, VEC<VEC<Int> >& s) const
{
	VEC<VEC<Int> > a_lable;
	Int counter(0);
	for_Int(row, 1, control_divs.nrows())
	{
		VEC<VEC<Int>> temp_a, a_row_temp;
		for_Int(col, 0, control_divs.ncols())
		{
			VEC<Int> one_div;
			for_Int(spl, 0, control_divs[row][col] + 1)
			{
				if (if_div_in_restraint(control_divs[0], col, control_divs[row][col], spl))
					one_div.push_back(spl);
			}
			a_row_temp.push_back(one_div);
		}
		temp_a = cart_product(a_row_temp);
		for (const auto &one_divs : temp_a)	{
			if (if_row_divs_in_restraint(nppso[row-1], one_divs, sit_mat[row-1])) a.push_back(one_divs);
		}
		VEC<Int> a_lable_i;
		for_Int(i, counter, a.size()){
			a_lable_i.push_back(i);
			counter++;
		}
		a_lable.push_back(a_lable_i);
		// if(mm) WRN(NAV2( row,a.size()));
	}

	s = cart_product_monitor_row(a_lable, a);
	// s = cart_product(a_lable);
}

void NocSpace::find_all_possible_state_by_nooc(VEC<VEC<Int> >& a, VEC<VEC<Int> >& s) const
{
	VEC<VEC<Int> > a_lable;
	Int counter(0);
	for_Int(col, 0, control_divs.ncols()) if (p.if_norg_imp || (col != 0 && col != ndivs / 2))
	{
		VEC<VEC<Int>> temp_a, a_rol_temp;
		for_Int(row, 1, control_divs.nrows())
		{
			VEC<Int> one_div;
			for_Int(spl, 0, control_divs[row][col] + 1)
			{
				if (if_div_in_restraint(control_divs[0], col, control_divs[row][col], spl))
					one_div.push_back(spl);
			}
			a_rol_temp.push_back(one_div);
		}
		temp_a = cart_product(a_rol_temp);
		for (const auto &one_divs : temp_a)	{
			if (if_col_divs_in_restraint(control_divs[0][col], one_divs, col))
				{a.push_back(one_divs);
				// if(mm)WRN(NAV(Vec(one_divs)));
				}
		}
		VEC<Int> a_lable_i;
		for_Int(i, counter, a.size()){
			a_lable_i.push_back(i);
			counter++;
		}
		a_lable.push_back(a_lable_i);
		// if(mm) WRN(NAV2( col,a.size()));
	}

	// s = cart_product_monitor_col(a_lable, a);
	s = cart_product(a_lable);
}

void NocSpace::find_all_possible_state(VEC<VEC<Int> >& a, VEC<VEC<Int> >& s) const {
	VEC<VEC<Int> > a_lable;
	Int counter(0);
	for_Int(col, 0, control_divs.ncols()) {
		VEC<VEC<Int>> temp_a;
		VEC<VEC<Int> > a_rol_temp;
		for_Int(row, 1, control_divs.nrows()) {
			VEC<Int> one_div;
			for_Int(spl, 0, control_divs[row][col] + 1) {
				if (if_div_in_restraint(control_divs[0], col, control_divs[row][col], spl))
					one_div.push_back(spl);
			}
			a_rol_temp.push_back(one_div);
		}
		temp_a = cart_product(a_rol_temp);
		for (const auto& one_divs : temp_a) {
			if (if_col_divs_in_restraint(control_divs[0][col], one_divs, col))
				a.push_back(one_divs);
		}
		VEC<Int> a_lable_i;
		for_Int(i, counter, a.size()) {
			a_lable_i.push_back(i);
			counter++;
		}
		a_lable.push_back(a_lable_i);
	}
	s = cart_product(a_lable);
}



//----------------------------- space function after judge -----------------------------

void NocSpace::find_thought_noc_subspaces()
{
	Idx length_ECNS;						// length for each combined number subspaces(ECNS).
	idx_div.push_back(dim);
	VEC<VEC<Int> > a, s;

	if(p.if_norg_imp) find_all_possible_state_by_col(a, s);
	else find_all_possible_state_by_nooc(a, s);
	// if(mm) WRN(NAV(s.size()));
	Vec<MatInt> spilss_divs = multi_judger_with_return(s, a);
	// if(mm) WRN(NAV(spilss_divs.size()));
	for_Int(k, 0, spilss_divs.size()) {
		const MatInt & spilss_div(spilss_divs[k]);
		div.push_back(spilss_div);
		divs_to_idx.insert(std::pair<std::string, Int>(spilss_div.vec().string(), dim));
		length_ECNS = 1;
		for_Int(i, 1, control_divs.nrows()) for_Int(j, 0, ndivs) length_ECNS = length_ECNS * bf(control_divs[i][j], spilss_div[i - 1][j]);
		dim += length_ECNS;
		idx_div.push_back(dim);
	}
}

void NocSpace::find_all_noc_subspaces()
{
	Idx length_ECNS;						// length for each combined number subspaces(ECNS).
	idx_div.push_back(dim);
	VEC<VEC<Int> > a, s;


	find_all_possible_state_by_col(a, s);
	// if(mm) WRN(NAV(s.size()));
	for (const auto& x : s) {
		VecInt spilss_div_v(read_from_col_lable(x,a));
		MatInt spilss_div_tr(spilss_div_v.mat(control_divs.ncols(), control_divs.nrows() - 1));
		MatInt spilss_div(spilss_div_tr.tr());
		if (ifin_NocSpace(spilss_div, nppso)) {
			// if (mm) WRN(NAV(spilss_div));
			div.push_back(spilss_div);
			divs_to_idx.insert(std::pair<std::string, Int>(spilss_div.vec().string(), dim));
			length_ECNS = 1;
			for_Int(i, 1, control_divs.nrows()) for_Int(j, 0, ndivs) length_ECNS = length_ECNS * bf(control_divs[i][j], spilss_div[i - 1][j]);
			dim += length_ECNS;
			idx_div.push_back(dim);
		}
	}
			// WRN(NAV(dim));
}

// find all the possible state speed up with MPI.
void NocSpace::find_all_noc_subspaces_multi() {
	Idx length_ECNS;						// length for each combined number subspaces(ECNS).
	idx_div.push_back(dim);
	VEC<VEC<Int> > a, s;

	find_all_possible_state_by_col(a, s);
	if(mm) PIO(NAV(s.size())+"   "+present());
	IFS ifs(STR("judger"+nppso_str()+ ".bdat"));
	VecInt judger_out(s.size());
	if(!ifs) judger_out = multi_judger(s, a);
	else biread(ifs, CharP(judger_out.p()), judger_out.szof());
	// for (const auto &x : out)
	for_Int(ii, 0, s.size()) if(judger_out[ii]) {
		VEC<Int>& x = s[ii];
		VecInt spilss_div_v(read_from_col_lable(x,a));
		MatInt spilss_div_tr(spilss_div_v.mat(control_divs.ncols(), control_divs.nrows() - 1));
		MatInt spilss_div(spilss_div_tr.tr());
		// if (ifin_NocSpace(spilss_div, nppso)) 
		div.push_back(spilss_div);
		divs_to_idx.insert(std::pair<std::string, Int>(spilss_div.vec().string(), dim));
		length_ECNS = 1;
		for_Int(i, 1, control_divs.nrows()) for_Int(j, 0, ndivs) length_ECNS = length_ECNS * bf(control_divs[i][j], spilss_div[i - 1][j]);
		dim += length_ECNS;
		idx_div.push_back(dim);
	}
}

void NocSpace::find_all_noc_subspaces_by_row()
{
	Idx length_ECNS;						// length for each combined number subspaces(ECNS).
	idx_div.push_back(dim);
	VEC<VEC<Int> > a, s;

	find_all_possible_state_by_row(a, s);
	if(mm) PIO(NAV(s.size())+"   "+present());
	IFS ifs(STR("judger"+nppso_str()+ ".bdat"));
	VecInt judger_out(s.size());
	if(!ifs) judger_out = multi_judger_by_row(s, a);
	else biread(ifs, CharP(judger_out.p()), judger_out.szof());
	// for (const auto &x : out)
	for_Int(ii, 0, s.size()) if(judger_out[ii]) {
		MatInt spilss_div(read_from_col_lable(s[ii],a).mat( control_divs.nrows() - 1, control_divs.ncols()));
		div.push_back(spilss_div);
		divs_to_idx.insert(std::pair<std::string, Int>(spilss_div.vec().string(), dim));
		length_ECNS = 1;
		for_Int(i, 1, control_divs.nrows()) for_Int(j, 0, ndivs) length_ECNS = length_ECNS * bf(control_divs[i][j], spilss_div[i - 1][j]);
		dim += length_ECNS;
		idx_div.push_back(dim);
	}
	// if(mm) {
	// 	OFS ofs;	ofs.open("judger"+nppso_str()+ ".bdat");
	// 	biwrite(ofs, CharP(judger_out.p()), judger_out.szof());
	// 	ofs.close();
	// }
}

//----------------------------- (old/stable)space function -----------------------------

// find all the combined number subspaces OR find all the combined number subspaces with special partical control.
void NocSpace::find_combined_number_subspaces(const Int mode)
{
	Idx length_ECNS;			// length for each combined number subspaces(ECNS).
	idx_div.push_back(dim);
	VEC<VEC<Int> > a;
	for_Int(row, 1, control_divs.nrows()) {
		for_Int(col, 0, control_divs.ncols()) {
			VEC<Int> one_div;
			for_Int(spl, 0, control_divs[row][col] + 1) {
				if(if_div_in_restraint(control_divs[0],col, control_divs[row][col],spl)) one_div.push_back(spl);
			}
			a.push_back(one_div);
		}
	}
	Idx total_posibile(1);
	for (const auto& u : a) total_posibile *= u.size();

	if(mode == 1) {
		for_Idx(x, 0, total_posibile) {
			VecInt spilss_div_v(free_div_base_decode(x, a));
			MatInt spilss_div(spilss_div_v.mat(control_divs.nrows() - 1, control_divs.ncols()));
			if (ifin_NocSpace(spilss_div, nppso)) {
				div.push_back(spilss_div);	length_ECNS = 1;
				for_Int(i, 1, control_divs.nrows()) for_Int(j, 0, ndivs) length_ECNS = length_ECNS * bf(control_divs[i][j], spilss_div[i - 1][j]);
				dim += length_ECNS;			idx_div.push_back(dim);
			}
		}
	}

	if(mode == 0) {
		for_Idx(x, 0, total_posibile) {
			VecInt spilss_div_v(free_div_base_decode(x, a));
			MatInt spilss_div(spilss_div_v.mat(control_divs.nrows() - 1, control_divs.ncols()));
			if (ifin_NocSpace(spilss_div)) {
				div.push_back(spilss_div);	length_ECNS = 1;
				for_Int(i, 1, control_divs.nrows()) for_Int(j, 0, ndivs) length_ECNS = length_ECNS * bf(control_divs[i][j], spilss_div[i - 1][j]);
				dim += length_ECNS;			idx_div.push_back(dim);
			}
		}
	}
}

//----------------------------- basic function -----------------------------

void NocSpace::set_control()
{
	control_divs = p.control_divs;
	sit_mat = control_divs.truncate_row(1, control_divs.nrows());
	for_Int(j, 0, control_divs.ncols())for_Int(i, 1, control_divs.nrows()) orbital_divcnt[j] += control_divs[i][j];
	// sit_mat.reset(control_divs.truncate_row(1, control_divs.nrows()));
}

VEC<VEC<Int> > NocSpace::cart_product_monitor_col (const VEC<VEC<Int> >& v, const VEC<VEC<Int> >& a)const
{
    VEC<VEC<Int> > s = {{}};
    for (const auto& u : v) {
		VEC<VEC<Int> > r;
		for (const auto& x : s) {
			// Int sum(0);
			// for (const auto x_i : x) sum += SUM(VecInt(a[x_i]));
			// if (sum <= nspa) {
			VecInt sum(ndivs,0); bool judge(true);
			for (const auto x_i : x) sum += VecInt(a[x_i]);
			for_Int(i, 0, ndivs)	if(sum[i] > nppso[i]) judge = false;
			// if (s[0].size() == sit_mat.size()) for_Int(i, 0, ndivs)	if(sum[i] != nppso[i]) judge = false;
			if (judge) {
				for (const auto y : u) {
					r.push_back(x);
					r.back().push_back(y);
				}
			}
		}
		s = move(r);
    }
    return s;
}

VEC<VEC<Int> > NocSpace::cart_product_monitor_row (const VEC<VEC<Int> >& v, const VEC<VEC<Int> >& a)const
{
    VEC<VEC<Int> > s = {{}};
    for (const auto& u : v) {
        VEC<VEC<Int> > r;
        for (const auto& x : s) {
			bool judge(true); 
			VecInt m_row(ndivs, 0), count_sit_mat(ndivs, 0);
			for_Int(i, 0, x.size()) {count_sit_mat += sit_mat[i]; m_row += Vec(a[x[i]]);}
			if (if_row_divs_in_restraint(SUM(nppso.truncate(0, x.size())), m_row.stdvec(), count_sit_mat))
				for (const auto y : u) {
					r.push_back(x);
					r.back().push_back(y);
				}
		}
        s = move(r);
    }
    return s;
}

VecInt NocSpace::free_div_base_decode(Int idx, VEC<VEC<Int> > v) const {
	VecInt base(v.size() + 1, 1);
	for_Int(i, 1, base.size()) base[i] = base[i - 1] * v[i - 1].size();

	VecInt rep(base.size() - 1, 0);
	for_Int(i, 0, rep.size()) rep[i] = v[i][(idx % base[i + 1]) / base[i]];
	if (idx >= base[base.size() - 1]) ERR(STR("free div base decode error."));
	return rep;
}

VecInt NocSpace::read_from_col_lable(const VEC<Int> x, const VEC<VEC<Int> > a) const {
	VecInt ren;
	for (const auto& lable : x) {
		VecInt ren_i(ren);
		VecInt col(a[lable]);
		ren.reset(concat(ren_i, col));
	}
	return ren;
}

Int NocSpace::wherein_NocSpace(const Int& h_i)const
{
	Int comdiv_idx(0);
	for_Int(j, 0, idx_div.size()) {
		if (h_i < idx_div[j]) {
			comdiv_idx = j - 1;
			break;
		}
		else if (h_i >= idx_div[idx_div.size() - 1]) comdiv_idx = idx_div.size() - 1;
	}
	return comdiv_idx;
}

void NocSpace::print(std::ostream& os) const {
#define nocspace_print(var, comment) print(os, NAME(var), STR(var), comment)

	using namespace std;
	Str NOOC = p.nooc_mode;
	Str nppsos = nppso_str();
	Int norbits = SUM(sit_mat);
	Str rotation_imp = p.if_norg_imp ? "Yes" : "NO";


	os << "// NORG setting" << endl;

	// nocspace_print(ndivs, "The amount of divisons's number. ");
	nocspace_print(dim, "the dimension of the Shortcut space.");
	nocspace_print(nppsos, "The number of partical per spin orbital. ");
	nocspace_print(nspa, "The amount of partical's number.");
	nocspace_print(norbits, "The amount of orbital(imp+bath) number.");
	nocspace_print(NOOC, "Correlation nature orbital occupation constraint.");
	nocspace_print(rotation_imp, "if rotation impurity orbital is used. ");
	// nocspace_print(ndivs, "The amount of divisons's number. ");
	nocspace_print(control_divs, "to set the number of division and the shortcut restratin.");
	// nocspace_print(h0, "transformed hopping integral");
	// nocspace_print(mu, "-mu");
	nocspace_print(p.U, "The Hubbard term U.");
	nocspace_print(p.Uprm, "The U^' term");
	nocspace_print(p.jz, "The hund coupling");
	nocspace_print(p.mu, "The hund coupling");
    nocspace_print(p.fit_max_omg, "The fitting max omg");
    nocspace_print(p.fit_num_omg, "The max number of fitting points.");
	// u_hbd, p.U12, p.U12, p.U12

	os << "// prmtr print end  " << present() << endl;

#undef nocspace_print
}
