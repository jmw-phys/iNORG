/*
coded by Jia-Ming Wang (jmw@ruc.edu.cn, RUC, China) date 2022 - 2023
*/

#include "operator.h"

using namespace std;


Operator::Operator(const MyMpi& mm_i, const Prmtr& prmtr_i, const NocSpace& s_i):
	mm(mm_i), p(prmtr_i), scsp(s_i), table(find_fullH_idx()), dim(s_i.dim), oper_value(pow(p.norbit, 4) + pow(p.norbit, 2) + 1)
{
	// for (const auto &i : table[2]) oper_value.insert(std::make_pair(ABS(i), 0.));
	// if (mm) WRN(NAV2(table[2].size(), oper_value.size()));
}

Operator::Operator(const MyMpi& mm_i, const Prmtr& prmtr_i,const Tab &tab, const VecReal& coefficient_i):
	mm(mm_i), p(prmtr_i), scsp(NocSpace(mm_i, prmtr_i)), table(tab), dim(tab[0].size()), oper_value(pow(p.norbit, 4) + pow(p.norbit, 2) + 1)
{
}

Operator::Operator(const MyMpi& mm_i, const Prmtr& prmtr_i, const NocSpace& s_i,const Tab &per_table):
	mm(mm_i), p(prmtr_i), scsp(s_i), table(per_table), dim(s_i.dim), oper_value(pow(p.norbit, 4) + pow(p.norbit, 2) + 1)
{
}

Operator::Operator(const MyMpi& mm_i, const Prmtr& prmtr_i, const NocSpace& s_i, Str tab_name):
	mm(mm_i), p(prmtr_i), scsp(s_i), table(read_the_Tab(tab_name)), dim(s_i.dim), oper_value(pow(p.norbit, 4) + pow(p.norbit, 2) + 1)
{
}

// ! (Abandon), this function can fully represent by the "find_fullH_idx()" function.
Tab Operator::find_h_idx(){
	clock_t t_find_hmlt_table;
	t_find_hmlt_table = clock(); // TIME_BGN("find_hmlt_table" + NAV(mm.id()), t_find_hmlt_table);
	VecPartition row_H(mm.np(), mm.id(), dim);
	Tab h_idxs(3);
	MatInt mat_hop_pos(scsp.hopint.nrows(),scsp.hopint.ncols());
	for_Int(i, 0, mat_hop_pos.nrows()) for_Int(j, 0, mat_hop_pos.ncols()) mat_hop_pos[i][j] = i * mat_hop_pos.ncols() + j;
	Int h_hbd_idx(mat_hop_pos.size()	+ 1);
	Int h_orb_ud_idx(mat_hop_pos.size() + 2);
	Int h_orb_uu_idx(mat_hop_pos.size() + 3);
	Int h_orb_dd_idx(mat_hop_pos.size() + 4);

	// Int h_orb_j_idx(mat_hop_pos.size() + 5);
	for_Int(h_i, row_H.bgn(), row_H.end()) {
		// To save as sparse matrix, [0]: row number;[1]: colum number;[2]: idx.
		VecInt h_idx(3, 0);
		Int sparse_idx(h_i - row_H.bgn());
		//WRN("wherein_NocSpace" + NAV(h_i - scsp.idx_div[scsp.wherein_NocSpace(h_i)]));
		StateStatistics a(h_i, scsp.wherein_NocSpace(h_i), scsp);
		for (const auto &x : a.filled_spinless) {
			for (const auto &i : x) {
				h_idx = { sparse_idx, h_i, mat_hop_pos[i][i] + 1 };
				for_Int(pos, 0, 3) h_idxs[pos].push_back(h_idx[pos]);
			}
		}


		#define ndiv scsp.ndivs
		{ // normal form of interation.
		// add the U.
		for_Int(i, 0, p.nband) if( a.cfg.cf[(i * 2) * ndiv].isocc(0) && a.cfg.cf[(i * 2 + 1) * ndiv].isocc(0) ){
			h_idx = { sparse_idx, h_i, h_hbd_idx};
			for_Int(pos, 0, 3) h_idxs[pos].push_back(h_idx[pos]);
		}

		// add the up-down term.
		for_Int(i, 0, p.nband) {
			for_Int(j, 0, p.nband) if(i != j && a.cfg.cf[(i * 2) * ndiv].isocc(0) && a.cfg.cf[(j * 2 + 1) * ndiv].isocc(0) ){
				h_idx = { sparse_idx, h_i, h_orb_ud_idx};
				for_Int(pos, 0, 3) h_idxs[pos].push_back(h_idx[pos]);
			}
		}

		// add the up-up term.
		for_Int(i, 0, p.nband) {
			for_Int(j, 0, p.nband) if(i != j && a.cfg.cf[(i * 2) * ndiv].isocc(0) && a.cfg.cf[(j * 2) * ndiv].isocc(0) ){
				h_idx = { sparse_idx, h_i, h_orb_uu_idx};
				for_Int(pos, 0, 3) h_idxs[pos].push_back(h_idx[pos]);
			}
		}

		// add the down-down term.
		for_Int(i, 0, p.nband) {
			for_Int(j, 0, p.nband) if(i != j && a.cfg.cf[(i * 2 + 1) * ndiv].isocc(0) && a.cfg.cf[(j * 2 + 1) * ndiv].isocc(0) ){
				h_idx = { sparse_idx, h_i, h_orb_dd_idx};
				for_Int(pos, 0, 3) h_idxs[pos].push_back(h_idx[pos]);
			}
		}
		}
		// // add the U.
		// for_Int(i, 0, p.nband) {
		// 	if((Real(a.cfg.cf[(i * 2) * ndiv][0]) - 0.5) * (Real(a.cfg.cf[(i * 2 + 1) * ndiv][0]) - 0.5) > 0 ) h_idx = { sparse_idx, h_i, h_hbd_idx };
		// 	else h_idx = { sparse_idx, h_i, -h_hbd_idx };
		// 	for_Int(pos, 0, 3) h_idxs[pos].push_back(h_idx[pos]);
		// }

		// // add the up-down term.
		// for_Int(i, 0, p.nband) {
		// 	for_Int(j, 0, p.nband) {
		// 		if (i != j && (Real(a.cfg.cf[(i * 2) * ndiv][0]) - 0.5) * (Real(a.cfg.cf[(j * 2 + 1) * ndiv][0]) - 0.5) > 0) h_idx = { sparse_idx, h_i, h_hbd_idx };
		// 		// if (i != j && a.cfg.cf[(i * 2) * ndiv].isocc(0) && a.cfg.cf[(j * 2 + 1) * ndiv].isocc(0)) h_idx = { sparse_idx, h_i, h_hbd_idx };
		// 		else h_idx = {sparse_idx, h_i, -h_orb_ud_idx};
		// 		for_Int(pos, 0, 3) h_idxs[pos].push_back(h_idx[pos]);
		// 	}
		// }

		// // add the up-up term.
		// for_Int(i, 0, p.nband) {
		// 	for_Int(j, 0, p.nband) {
		// 	if(i != j && (Real(a.cfg.cf[(i * 2) * ndiv][0]) - 0.5) * (Real(a.cfg.cf[(j * 2) * ndiv][0]) - 0.5) > 0) h_idx = { sparse_idx, h_i, h_orb_uu_idx};
		// 	else h_idx = { sparse_idx, h_i, -h_orb_uu_idx};
		// 		for_Int(pos, 0, 3) h_idxs[pos].push_back(h_idx[pos]);
		// 	}
		// }

		// // add the down-down term.
		// for_Int(i, 0, p.nband) {
		// 	for_Int(j, 0, p.nband) {
		// 	if(i != j && (Real(a.cfg.cf[(i * 2 + 1) * ndiv][0]) - 0.5) * (Real(a.cfg.cf[(j * 2 + 1) * ndiv][0]) - 0.5) > 0) h_idx = { sparse_idx, h_i, h_orb_dd_idx};
		// 	else h_idx = { sparse_idx, h_i, -h_orb_dd_idx};
		// 		for_Int(pos, 0, 3) h_idxs[pos].push_back(h_idx[pos]);
		// 	}
		// }
 
		// for_Int(i, 0, p.nband) {
		// 	if((Real(a.cfg.cf[0][i]) - 0.5) * (Real(a.cfg.cf[0][p.nband + i]) - 0.5) > 0 ) h_idx = { sparse_idx, h_i, h_hbd_idx };
		// 	else h_idx = { sparse_idx, h_i, -h_hbd_idx };
		// 	for_Int(pos, 0, 3) h_idxs[pos].push_back(h_idx[pos]);
		// }

		/* 
		#define ndiv scsp.ndivs

		for_Int(i, 0, p.norg_sets)
		{
			for_Int(j, i, p.norg_sets)
			{
				if (i != j)
				{
					if (i / 2 == j / 2 && i + 1 == j)
					{
						if (a.cfg.cf[i * ndiv][0] && a.cfg.cf[j * ndiv][0])
							h_idx = {sparse_idx, h_i, h_hbd_idx};
						// if (a.cfg.cf[i * ndiv][0] != a.cfg.cf[j * ndiv][0])
						// 	h_idx = {sparse_idx, h_i, -h_hbd_idx};
						for_Int(pos, 0, 3) h_idxs[pos].push_back(h_idx[pos]);
					}
					if (i / 2 != j / 2 && i % 2 == 1 && j % 2 == 0)
					{
						if (a.cfg.cf[i * ndiv][0] && a.cfg.cf[j * ndiv][0])
							h_idx = {sparse_idx, h_i, h_orb_ud_idx};
						// if (a.cfg.cf[i * ndiv][0] != a.cfg.cf[j * ndiv][0])
						// 	h_idx = {sparse_idx, h_i, -h_orb_ud_idx};
						for_Int(pos, 0, 3) h_idxs[pos].push_back(h_idx[pos]);
					}
					if (i / 2 != j / 2 && i % 2 == 0 && j % 2 == 1)
					{
						if (a.cfg.cf[i * ndiv][0] && a.cfg.cf[j * ndiv][0])
							h_idx = {sparse_idx, h_i, h_orb_ud_idx};
						// if (a.cfg.cf[i * ndiv][0] != a.cfg.cf[j * ndiv][0])
						// 	h_idx = {sparse_idx, h_i, -h_orb_ud_idx};
						for_Int(pos, 0, 3) h_idxs[pos].push_back(h_idx[pos]);
					}
					if (i / 2 != j / 2 && i % 2 == 1 && j % 2 == 1)
					{
						if (a.cfg.cf[i * ndiv][0] && a.cfg.cf[j * ndiv][0])
							h_idx = {sparse_idx, h_i, h_orb_uu_idx};
						// if (a.cfg.cf[i * ndiv][0] != a.cfg.cf[j * ndiv][0])
						// 	h_idx = {sparse_idx, h_i, -h_orb_uu_idx};
						for_Int(pos, 0, 3) h_idxs[pos].push_back(h_idx[pos]);
					}
					if (i / 2 != j / 2 && i % 2 == 0 && j % 2 == 0)
					{
						if (a.cfg.cf[i * ndiv][0] && a.cfg.cf[j * ndiv][0])
							h_idx = {sparse_idx, h_i, h_orb_dd_idx};
						// if (a.cfg.cf[i * ndiv][0] != a.cfg.cf[j * ndiv][0])
						// 	h_idx = {sparse_idx, h_i, -h_orb_dd_idx};
						for_Int(pos, 0, 3) h_idxs[pos].push_back(h_idx[pos]);
					}
				}
			}
		}
 		*/

		for_Int(idx_sets, 0, p.norg_sets)
		{
			// off_diagonal_term
			// i[0]:annihilation orbit's position; i[1]:creation orbit's positon; i[2]:Colum idx(i); i[3]:sign(anticommutativity)
			VEC<VecInt> off_dt_next(a.find_each_spiless_group_off_diagonal_term(a.divocchop_ingroup(a.div_idx, idx_sets), idx_sets));
			for (const auto &i : off_dt_next){
				h_idx = {sparse_idx, i[2], i[3] * (mat_hop_pos[i[1]][i[0]] + 1)};
				for_Int(pos, 0, 3) h_idxs[pos].push_back(h_idx[pos]);
			}
		}

		// ! The SOC term.(not finished.)
		// {
		// 	// off_diagonal_term for SOC
		// 	// i[0]:annihilation orbit's position; i[1]:creation orbit's positon; i[2]:Colum idx(i); i[3]:sign(anticommutativity)
		// 	VEC<VecInt> off_d_interation(a.off_diagonal_soc_term(a.interation_soc_hop(a.div_idx)));
		// 	for (const auto &i : off_d_interation){
		// 		h_idx = {sparse_idx, i[2], i[3] * (h_orb_j_idx)};
		// 		for_Int(pos, 0, 3) h_idxs[pos].push_back(h_idx[pos]);
		// 	}
		// }
	}
	// TIME_END("t_find_hmlt_table" + NAV(mm.id()), t_find_hmlt_table);
	for_Int(jj, 0, 3) h_idxs[jj].shrink_to_fit();
	return std::move(h_idxs);
}


Tab Operator::find_fullH_idx()
{
	clock_t t_find_hmlt_table;
	t_find_hmlt_table = clock(); // TIME_BGN("find_hmlt_table" + NAV(mm.id()), t_find_hmlt_table);
	VecPartition row_H(mm.np(), mm.id(), dim);
	Tab h_idxs(3);
	MatInt mat_hop_pos(scsp.hopint.nrows(),scsp.hopint.ncols());
	for_Int(i, 0, mat_hop_pos.nrows()) for_Int(j, 0, mat_hop_pos.ncols()) mat_hop_pos[i][j] = i * mat_hop_pos.ncols() + j;
	Mat<MatInt> tensor_u(scsp.hopint.nrows(), scsp.hopint.ncols(), mat_hop_pos);
	for_Int(i, 0, tensor_u.nrows()) for_Int(j, 0, tensor_u.ncols()) for_Int(k, 0, mat_hop_pos.nrows()) for_Int(l, 0, mat_hop_pos.ncols())
		tensor_u[i][j][k][l] = i * std::pow(p.norbit, 3) + j * std::pow(p.norbit, 2) + k * std::pow(p.norbit, 1) + l;
		// tensor_u[i][j][k][l] = (i * tensor_u.ncols() + j) * mat_hop_pos.size() + k * mat_hop_pos.ncols() + l;
	for_Int(h_i, row_H.bgn(), row_H.end()) {
		// To save as sparse matrix, [0]: row number;[1]: colum number;[2]: idx.
		VecInt h_idx(3, 0);
		Int sparse_idx(h_i - row_H.bgn());
		//WRN("wherein_NocSpace" + NAV(h_i - scsp.idx_div[scsp.wherein_NocSpace(h_i)]));

		StateStatistics a(h_i, scsp.wherein_NocSpace(h_i), scsp);
		
		//! Diagonal term. n_i
		for (const auto &x : a.filled_spinless) {
			for (const auto &i : x) {
				h_idx = { sparse_idx, h_i, mat_hop_pos[i][i] + 1 };
				for_Int(pos, 0, 3) h_idxs[pos].push_back(h_idx[pos]);
			}
		}
/*
		//! Diagonal term. n_i n_j
		for_Int(b1, 0, p.nband) {
			if( a.cfg.cf[(b1 * 2) * ndiv].isocc(0) && a.cfg.cf[(b1 * 2 + 1) * ndiv].isocc(0)){
				Int interact_pos = mat_hop_pos.size() + 1 + tensor_u[SUM_0toX(p.nO2sets, (b1 * 2))][SUM_0toX(p.nO2sets, (b1 * 2) + 1)][SUM_0toX(p.nO2sets, (b1 * 2) + 1)][SUM_0toX(p.nO2sets, (b1 * 2))];
				h_idx = { sparse_idx, h_i, interact_pos };
				for_Int(pos, 0, 3) h_idxs[pos].push_back(h_idx[pos]);
			}
			for_Int(b2, 0, p.nband) {
				if (b1 != b2) {
					if ((a.cfg.cf[(b1 * 2) * ndiv].isocc(0) && a.cfg.cf[(b2 * 2 + 1) * ndiv].isocc(0))) {
						Int interact_pos = mat_hop_pos.size() + 1 + tensor_u[SUM_0toX(p.nO2sets, (b1 * 2))][SUM_0toX(p.nO2sets, (b2 * 2) + 1)][SUM_0toX(p.nO2sets, (b2 * 2) + 1)][SUM_0toX(p.nO2sets, (b1 * 2))];
						h_idx = { sparse_idx, h_i, interact_pos };
						for_Int(pos, 0, 3) h_idxs[pos].push_back(h_idx[pos]);
					}
					if ((a.cfg.cf[(b1 * 2) * ndiv].isocc(0) && a.cfg.cf[(b2 * 2) * ndiv].isocc(0))) {
						Int interact_pos = mat_hop_pos.size() + 1 + tensor_u[SUM_0toX(p.nO2sets, (b1 * 2))][SUM_0toX(p.nO2sets, (b2 * 2))][SUM_0toX(p.nO2sets, (b2 * 2))][SUM_0toX(p.nO2sets, (b1 * 2))];
						h_idx = { sparse_idx, h_i, interact_pos };
						for_Int(pos, 0, 3) h_idxs[pos].push_back(h_idx[pos]);
					}
					if ((a.cfg.cf[(b1 * 2 + 1) * ndiv].isocc(0) && a.cfg.cf[(b2 * 2 + 1) * ndiv].isocc(0))) {
						Int interact_pos = mat_hop_pos.size() + 1 + tensor_u[SUM_0toX(p.nO2sets, (b1 * 2 + 1))][SUM_0toX(p.nO2sets, (b2 * 2) + 1)][SUM_0toX(p.nO2sets, (b2 * 2) + 1)][SUM_0toX(p.nO2sets, (b1 * 2 + 1))];
						h_idx = { sparse_idx, h_i, interact_pos };
						for_Int(pos, 0, 3) h_idxs[pos].push_back(h_idx[pos]);
					}
				}
			}
		}
*/
		//! Diagonal term. n_i n_j
		for_Int(set_i, 0, p.norg_sets) {
			for_Int(set_j, 0, p.norg_sets) {
				VEC<Int> N_veci(a.filled_spinless[set_i]), N_vecj(a.filled_spinless[set_j]);
				for (const auto& N_i : N_veci) for (const auto& N_j : N_vecj) {
					// if (N_i != N_j) {
						h_idx = { sparse_idx, h_i, int(mat_hop_pos.size() + 1 + tensor_u[N_i][N_j][N_j][N_i]) }; for_Int(pos, 0, 3) h_idxs[pos].push_back(h_idx[pos]);
						// ? if this term(↓) is necessary? only can mute this term when the interaction happen in different sets.
						// h_idx = { sparse_idx, h_i,-int(mat_hop_pos.size() + 1 + tensor_u[N_i][N_j][N_i][N_j]) }; for_Int(pos, 0, 3) h_idxs[pos].push_back(h_idx[pos]);
					// }
				}
			}
		}
		//! off-Diagonal term. C\power-index{i|†}C\index{j}|| C\power-index{i|†}C\power-index{j|†}C\power-index{k}C\index{l}
		for_Int(sets_i, 0, p.norg_sets) {
			{
				// off_diagonal_term Two-Fermi
				// i[0]:annihilation orbit's position; i[1]:creation orbit's positon; i[2]:Colum idx(i); i[3]:sign(anticommutativity)
				VEC<VecInt> off_dt_next(a.find_each_spiless_group_off_diagonal_term(a.divocchop_ingroup(a.div_idx, sets_i), sets_i));
				for (const auto &i : off_dt_next){
					h_idx = {sparse_idx, i[2], i[3] * (mat_hop_pos[i[1]][i[0]] + 1)};
					for_Int(pos, 0, 3) h_idxs[pos].push_back(h_idx[pos]);

					//! Diagonal+off-Diagonal term. n_sC_e^+C_f
					if (p.if_norg_imp) for (const auto& x : a.filled_spinless) {
						for (const auto& N_s : x) if (N_s != i[0] && N_s != i[1]) {
							//																		[i   ][j  ][k  ][l   ]
							h_idx = { sparse_idx, i[2],-i[3] * int(mat_hop_pos.size() + 1 + tensor_u[N_s][i[1]][N_s][i[0]]) }; for_Int(pos, 0, 3) h_idxs[pos].push_back(h_idx[pos]);
							h_idx = { sparse_idx, i[2], i[3] * int(mat_hop_pos.size() + 1 + tensor_u[N_s][i[1]][i[0]][N_s]) }; for_Int(pos, 0, 3) h_idxs[pos].push_back(h_idx[pos]);
							h_idx = { sparse_idx, i[2], i[3] * int(mat_hop_pos.size() + 1 + tensor_u[i[1]][N_s][N_s][i[0]]) }; for_Int(pos, 0, 3) h_idxs[pos].push_back(h_idx[pos]);
							h_idx = { sparse_idx, i[2],-i[3] * int(mat_hop_pos.size() + 1 + tensor_u[i[1]][N_s][i[0]][N_s]) }; for_Int(pos, 0, 3) h_idxs[pos].push_back(h_idx[pos]);
						}
					}
				}
			}

			if (p.if_norg_imp) {				
				// for_Int(sets_j, 0, p.norg_sets) {
				for_Int(sets_j, 0, p.norg_sets) if (sets_j != sets_i) {
					// off_diagonal_term Four-Fermi
					// [0]~[3] i-j-k-l orbit's position(C^+_i C^+_j C_k C_l); [4]:Colum idx(i);[5]:sign(fermion anticommutativity)
					VEC<array<int, 6>> off_dt_next(a.find_off_diagonal_term_fourFermi(a.divs_change_fourFermi(a.div_idx, sets_i, sets_j), sets_i, sets_j));
					for (const auto& i : off_dt_next) {
						h_idx = { sparse_idx, i[4], i[5] * int(mat_hop_pos.size() + 1 + tensor_u[i[0]][i[1]][i[2]][i[3]]) };
						for_Int(pos, 0, 3) h_idxs[pos].push_back(h_idx[pos]);
					}
				}
			}
		}
	}
	// TIME_END("t_find_hmlt_table" + NAV(mm.id()), t_find_hmlt_table);
	for_Int(jj, 0, 3) h_idxs[jj].shrink_to_fit();
	return std::move(h_idxs);
}

SparseMatReal Operator::find_hmlt(const Tab h_idx) const
{
	// clock_t t_find_hmlt;
	// TIME_BGN("find_hmlt" + NAV(mm.id()), t_find_hmlt);
	VecPartition row_H(mm.np(), mm.id(), dim);
	SparseMatReal hmlt_splited(row_H.len(), dim, mm);
	Real diagonal(0.);	Bool flag(false);
	for_Idx(pos, 0, h_idx[2].size())
	{
		Int coeff_comm = h_idx[2][pos] >= 0 ? 1 : -1;
		if (h_idx[1][pos] == h_idx[0][pos] + row_H.bgn()) {diagonal += coeff_comm * oper_value[(abs(h_idx[2][pos]))];flag=true;}
		else {
			if (flag) {hmlt_splited.addelement(diagonal, h_idx[1][pos-1], h_idx[0][pos-1]); diagonal = 0.;flag=false;}
			hmlt_splited.addelement(coeff_comm * oper_value[(abs(h_idx[2][pos]))], h_idx[1][pos], h_idx[0][pos]);
		}
	}
	// TIME_END("find_hmlt" + NAV(mm.id()), t_find_hmlt);
	hmlt_splited.shrink_to_fit();
	hmlt_splited.fix();
	return hmlt_splited;
}


// Try to find the Ground state.
MatReal Operator::lowest_eigpairs(const Idx n, bool if_need_fast, Int wish_nev)
//	n      ( input) is the dimension of the matrix
//	wish_nev    ( input) is the expected maximum eigenstates to be found
//	actual number of eigenpairs is returned
//	nev_eff is the number of eigenvectors corresponding to complete degenerate levels reached
{
	UURandomSFMT uur;
	VecReal eval(wish_nev, 0.);				// is the eigenvalues
	// VEC < VecReal> eigenvec_i;			// a VEC to contain all the eigenvector.
	MatReal eigenvec_i(wish_nev, n);		// a VEC to contain all the eigenvector.
	VecInt ev_dgcy(wish_nev, 0.);			// cotain all the degeneracy eigenvector's number for each eigenvalue.
	// if (mm)PIO("find_hmlt BEGIN :::");
	SparseMatReal sep_hmltoperator = find_hmlt(table);
	// if(mm) WRN("finished fiding the hmlt")
//#ifdef _ASSERTION_
//	if (mm.np() == 1) {
//		WRN("test begin:::");
//		if (sep_hmltoperator.if_hermitian() == false)ERR("This sparse matrix is not Hermitian matrix!!!!");
//		WRN("test end  !!!");
//	}
//#endif
	VecReal inital_state(n, 0.); uur(inital_state); inital_state -= 0.5;
	VecInt krylov_space_size = lanczos(eval, eigenvec_i, ev_dgcy, n, wish_nev, sep_hmltoperator, inital_state, mm, if_need_fast, 1000);
	if(mm) {cout <<"PIO: krylov_space_size = ";cout <<iofmt("sci"); for_Int(i, 0, krylov_space_size.size()) cout << eval[i] << ","<< krylov_space_size[i] << "; "; cout << std::endl;}
	// if(mm) std::cout << "The eigenvalue" << iofmt("sci") << eval << std::endl;
	groundstate_energy = eval[0];
	if(ev_dgcy[0] != p.degel) {p.degel = ev_dgcy[0]; if(mm) WRN(NAV2(p.degel,ev_dgcy.mat(1,ev_dgcy.size())))}
	// MatReal eigenvec(eigenvec_i.size(),n);
	// for_Int(i, 0, eigenvec_i.size()) eigenvec[i] = eigenvec_i[i];
	// ground_state = eigenvec_i[0];
	ground_state.reset(eigenvec_i.truncate_row(0, p.degel));
	// if(mm) WRN("finished lanczos.")
#ifdef _ASSERTION_
	//WRN("TEST_Lanczos is right:::" + NAV(mm.np()));
	//VecReal test_a(eigenvec[0]);
	//VecReal test_b(sep_hmltoperator * test_a);
	//WRN("TEST_Lanczos[0] state" + NAV3(test_a.isnormalized(), eval[0], test_b.avg_abs_elem_diff(eval[0] * test_a)));
#endif
	// if(mm) WRN(NAV(groundstate_energy));
	return eigenvec_i;
}


VecReal Operator::sn_prtcl_ex_state(const Int imp_div, const VecReal ground_state, const Int crtann) const
{
	VecPartition row_H(mm.np(), mm.id(), dim);
	VecReal ex_state_part(row_H.len(), 0.);
	for_Int(h_i, row_H.bgn(), row_H.end())
	{
		Int subscsp(scsp.wherein_NocSpace(h_i));
		ComDivs cfg(h_i - scsp.idx_div[subscsp], (scsp.div[subscsp]), (scsp.sit_mat), true);
		if (crtann == +1)if (cfg.cf[imp_div * scsp.ndivs].isuno(0))ex_state_part[h_i - row_H.bgn()] = 1.;
		if (crtann == -1)if (cfg.cf[imp_div * scsp.ndivs].isocc(0))ex_state_part[h_i - row_H.bgn()] = 1.;
	}
	ex_state_part *= ground_state.truncate(row_H.bgn(), row_H.end());
	VecReal ex_state(dim, 0.);
	ex_state = mm.Allgatherv(ex_state_part, row_H);
	return std::move(ex_state);
}

//------------------------------------------------------------------ io ------------------------------------------------------------------

void Operator::save_the_Tab(Tab& tab, Str name) const{
	Int size_temp(tab[0].size());
	Int size = mm.Allreduce(size_temp);
	VecInt v_size_i(1); v_size_i[0] = tab[0].size();
	VecPartition split_v_size(mm.np(), mm.id(), mm.np());
	VecInt v_size = mm.Allgatherv(v_size_i, split_v_size);
	if(mm) {// write the Tab's size's info
		OFS ofs;	ofs.open(name + ".inf");
		ofs << setw(9) << "dim" << setw(p_Real) << dim << endl;
		ofs << setw(9) << "size" << setw(p_Real) << size << endl;
		for_Int(i, 0, mm.np())	{
			ofs << setw(9) << "size_np"+STR(i) << setw(p_Real) << v_size[i] << endl;
		}
		ofs.close();
	}

	OFS ofs;	ofs.open(name + ".bdat");
	for_Int(i, 0, tab.size()) {
		VecPartition split_table(mm.np(), mm.id(), size, v_size);
		VecInt temp_tabi = mm.Gatherv(Vec(tab[i]), split_table);
		if(mm) {
			biwrite(ofs, CharP(temp_tabi.p()), temp_tabi.szof());
			
			// if(mm) WRN(NAV2(temp_tabi.size(),temp_tabi.truncate(size-100,size)));
		}		
	}
	ofs.close();
}

void Operator::write_the_multiTab(Str name) const{
	Int size(table[0].size());

	OFS ofs;	ofs.open(name + to_string(mm.id()) + ".txt");
	ofs << "\t" << setw(w_Real) << "row_idx";
	ofs << "\t" << setw(w_Real) << "col_idx";
	ofs << "\t" << setw(w_Real) << "hmt_idx";
	ofs << endl; 
	for_Int(i, 0, size) 
	{
		ofs << iofmt("sci");
		ofs << "\t" << setw(w_Real) << table[0][i];
		ofs << "\t" << setw(w_Real) << table[1][i];
		ofs << "\t" << setw(w_Real) << table[2][i];
		ofs << endl; 
	}
	ofs.close();
}

Tab Operator::read_the_Tab(Str name) const{
	
	// WRN("Here is fine"+NAV(name));
	VecInt v_size(mm.np(), 0);
	Int size(-1);	Tab tab(3);
	{
		IFS ifs(STR(name + ".inf"));	Str strr;
		while(1) {// read the Tab's size's info
			ifs >> strr;
			if(strr == "size")	ifs >> size;
			for_Int(i, 0, mm.np()) if(strr == "size_np"+STR(i))	ifs >> v_size[i];
			if (!ifs) break;
		}
	}
	// if(mm) WRN("Here is fine"+NAV(size));
	{
		IFS ifs(STR(name + ".bdat"));
		// VecPartition split_tab(mm.np(), mm.id(), size);		
		VecPartition split_tab(mm.np(), mm.id(), size, v_size);
		for_Int(i, 0, tab.size()) {
			VecInt temp(size); biread(ifs, CharP(temp.p()), temp.szof());
			// if(mm) WRN(NAV(temp.truncate(size-100,size)));
			tab[i] = temp.truncate(split_tab.bgn(),split_tab.end()).stdvec();
			// SWAP(tab[i], temp.truncate(split_tab.bgn(),split_tab.end()).stdvec());
		}
	}
	return tab;
}

// ! only suit for the two orbital cases.
MatReal Operator::local_multiplets_state(const MatReal& g_state) const {
	MatReal local_multiplets_state(3, 6, 0.);
	VecReal gs_p = g_state[0]; gs_p.normalize();

	//set the first line for N;
	// {0, 1, 2,..., 2*n}
	local_multiplets_state[0] = VecReal{ 0,  1,  2,  2,  3,  4 };
	//set the second line for Sz;
	// {0, 1/2, 2/2,..., n/2}
	local_multiplets_state[1] = VecReal{ 0,0.5,  0,  1,0.5,  0 };

	for_Int(h_i, 0, dim) {
		StateStatistics a(h_i, scsp.wherein_NocSpace(h_i), scsp);
		VecInt imp(4, 0);// the last one left for check the right.
		for_Int(i, 0, p.norbs) imp[i] = a.cfg.cf[(i)*ndiv].isocc(0) ? 1 : 0;
		switch (SUM(imp)) {
		case 0:
			local_multiplets_state[2][0] += gs_p[h_i] * gs_p[h_i];
			break;
		case 1:
			local_multiplets_state[2][1] += gs_p[h_i] * gs_p[h_i];
			break;
		case 2:
			if (SUM(imp.mat(2, 2).tr()[0]) == 1) local_multiplets_state[2][2] += gs_p[h_i] * gs_p[h_i];
			else local_multiplets_state[2][3] += gs_p[h_i] * gs_p[h_i];
			break;
		case 3:
			local_multiplets_state[2][4] += gs_p[h_i] * gs_p[h_i];
			break;
		case 4:
			local_multiplets_state[2][5] += gs_p[h_i] * gs_p[h_i];
			break;
		default:
			cout << "Input is not within the specified range" << endl;
		}
	}
	return local_multiplets_state;
}
