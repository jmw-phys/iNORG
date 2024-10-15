/*
code developed and maintained by (jmw@ruc.edu.cn, RUC, China) date 2023.03.03
*/

// adaptive natural sampling configuration interaction

#include "asnci.h"

using namespace std;
// ! NOW this class need to modify since need rotation the whole impurity.

Asnci::Asnci(const NORG& norg, Idx trncat_size):
    dim(trncat_size), mm(norg.mm), p(norg.p), impH(norg.impH), mayhop(find_mayhop()),
    nosp(norg.scsp), groundE(norg.groune_lst), 
    // core_dim(Int(trncat_size / Int(mayhop.size() )))
    core_dim(trncat_size)
{
    // inital = git_nci(norg);
    trncat = truncation(git_nci(norg.final_ground_state[0]));
}

Asnci::Asnci(const NORG& norg, Idx trncat_size, Int ex_pos):
    dim(trncat_size), mm(norg.mm), p(norg.p), impH(norg.impH), mayhop(find_mayhop()),
    nosp(norg.scsp), groundE(norg.groune_lst), 
    // core_dim(Int(trncat_size / Int(mayhop.size() / 4.0)))
    core_dim(trncat_size),
    trncat(git_nci_no_rank(norg.final_ground_state[0], ex_pos))
{
    // !Nci nci = git_nci(norg.final_ground_state[0], ex_pos);
    // !trncat = truncation(nci);
    // trncat = git_nci_no_rank(norg.final_ground_state[0], ex_pos);

    dim = trncat.first.size();

    // add the map



    if(mm) WRN(NAV2(dim,cfig_idx.size()))
    // // for_Int(i, 0, dim) cfig_idx.insert(pair<array<UInt,6>, Int>(trncat.first[i], i));
    for_Int(i, 0, dim) {
        if(!cfig_idx.insert(make_pair(trncat.first[i], i)).second){
            
            if(mm) WRN(NAV3(show_string(trncat.first[i]), cfig_idx.at(trncat.first[i]), i));
            
        }
    }
    // // for_Int(i, 0, dim) cfig_idx[trncat.first[i]] = i;


    if(mm) WRN(NAV(cfig_idx.size()))
    find_fullH_idx();

    // trncat = truncation(git_nci(norg.final_ground_state[0], ex_pos));
}

NORG Asnci::get_norg(Tab table, Int mode) {
    NORG norg(mm, p, table);
    norg.up_date_h0_to_solve(impH, mode);
    return norg;    
}


void Asnci::asnci_gimp(Green& imp_i, Int ex_pos)
{
    Int crtorann = ex_pos > 0 ? 1 : -1, pos(ABS(ex_pos) - 1);
    Operator oper(mm, p, find_fullH_idx(), coefficient);
    CrrltFun temp_green(mm, p, oper, find_ex_state(), crtorann);
    if (imp_i.type_info() == STR("ImGreen"))
    {
        ImGreen green_function(1, p);
        if (ex_pos > 0) temp_green.find_gf_greater(groundE, green_function);
        if (ex_pos < 0) temp_green.find_gf_lesser(groundE, green_function);
        // for_Int(n, 0, green_function.nomgs) imp_i[n][pos][pos] += green_function[n][0][0];
    }
    if(imp_i.type_info() == STR("ReGreen")) {
        ReGreen green_function(1, p);
        if(ex_pos > 0) temp_green.find_gf_greater(groundE, green_function);
        if(ex_pos < 0) temp_green.find_gf_lesser(groundE, green_function);
        // for_Int(n, 0, green_function.nomgs) imp_i[n][pos][pos] += green_function[n][0][0];
    }
}


//------------------------------------------------------------------ private ------------------------------------------------------------------

VEC<Int> Asnci::find_mayhop() {
    VEC<Int>    mayhop_i;
    if(impH.first.size() == 0) ERR("hop_h.size() = 0");
    // for_Int(i, 0, hop_h.size())  if(ABS(hop_h.vec()[i]) > 1e-6) mayhop_i.push_back(i);
    for_Int(i, 0, impH.first.size())  if(ABS(impH.first.vec()[i]) > 0.0) mayhop_i.push_back(i);
    mayhop_i.shrink_to_fit();
    if(mm) WRN(NAV2(impH.first, mayhop_i.size()));
    return move(mayhop_i);
}

VecBool Asnci::ints2vecbool(std::array<UInt,6> nums) {
    VecBool cfig(p.norbit, false);
    for_Int(i, 0, p.norbs){
        UInt num = nums[i];
        VecBool result(p.nI2B[i]+1, false);
        for_Int(j, 0, result.size()) {
            result[j] = (num >> j) & 1;
        }
        for_Int(j, 0, result.size()) cfig[result.size() * i + j] = result[j];
    }
    return cfig;
}

VecBool Asnci::back2nospace(const VecBool& new_cfig, const Int change_pos){
    Mat<bool> cfig(new_cfig.mat(nosp.sit_mat.nrows(), p.nI2B[0]+1));
    Idx pos(ABS(change_pos) - 1);
    // cfig[pos][0] ^= 1;
    if(change_pos > 0 && !cfig[pos][0]) cfig[pos][0] = true;
    else if(change_pos < 0 && cfig[pos][0]) cfig[pos][0] = false;
    else cfig.reset();
    return cfig.vec();
}

std::array<UInt,6> Asnci::vecbool2ints(const VecBool &vec) {
    std::array<UInt,6> nums;
    for_Int(i, 0, nums.size()){
        UInt result = 0;
        Int_for(j, (SUM_0toX(p.nI2B,i)+i), (SUM_0toX(p.nI2B,i+1)+i+1)){
            result <<= 1;
            if (vec[j]) result |= 1;
        }
        nums[i] = result;
    }
    return std::move(nums);
}

VecBool Asnci::change_cfg(const VecBool& cfg, Int pos) {
    VecBool cfg_new(cfg);
    Int crt(pos / impH.first.ncols()), ann(pos % impH.first.ncols());
    if((!(!cfg[crt]) && cfg[ann])) {
        Str cf; cf.reserve(cfg.size());
        for (const auto &b : cfg) cf += b ? '1':'0';
        ERR("Some thing wrong with the cfg :" + NAV3(cf, crt, ann));
    }
    cfg_new[crt] = true; cfg_new[ann] = false;
    return cfg_new;
}

Mat<Char> Asnci::show_string(const std::array<UInt,6>& cifg) {
    VecBool cfg(ints2vecbool(cifg));
    Vec<Char> out(cfg.size());
    for_Int(i,0,cfg.size()) out[i] = cfg[i] ? '1':'0';
    // if(mm) WRN(NAV(out.mat(p.norbs,p.nI2B[0]+1)));
    return out.mat(p.norbs,p.nI2B[0]+1);
}

// get the norg ground state to Nci's space
Nci Asnci::git_nci(const VecReal& ground_state) {
    Nci                         natural_cfg;
    VEC<array<UInt,6>>&         cfigs(natural_cfg.first);
    VEC<Real>&                  ranks(natural_cfg.second);

    VecReal grndste_norm = SQR(ground_state);
    VecIdx groundstate_idx(sort_indexes(grndste_norm));
    for_Int(i, 0, core_dim) {
        // if(grndste_norm[i] > 1e-5){
        {
            const Idx& ci_idx = groundstate_idx[i];
            StateStatistics cig(ci_idx, nosp.wherein_NocSpace(ci_idx), nosp);

            // {// open to test the ground state.            
            //     if(mm) cig.string();
            //     if(mm) WRN(NAV4(i,ci_idx,grndste_norm[ci_idx],ex_pos));
            //     if(mm) show_string(cig.cfg2nums());
            // }
            cfigs.push_back(vecbool2ints(cig.cfg2vecbool()));
            ranks.push_back(ground_state[ci_idx]);
        }
    }
    expand(natural_cfg); cfigs.shrink_to_fit(); ranks.shrink_to_fit();
    return natural_cfg;
}

// get the norg ground state's ex-state to Nci's space
Nci Asnci::git_nci(const VecReal& ground_state, const Int ex_pos) {
    Nci                         natural_cfg;
    VEC<array<UInt,6>>&         cfigs(natural_cfg.first);
    VEC<Real>&                  ranks(natural_cfg.second);

    VecReal grndste_norm = SQR(ground_state);
    VecIdx groundstate_idx(sort_indexes(grndste_norm));

    for_Int(i, 0, core_dim) {
        const Idx& ci_idx = groundstate_idx[i];
        StateStatistics cig(ci_idx, nosp.wherein_NocSpace(ci_idx), nosp);

        if(check_ifinex(cig, ex_pos)){
            cfigs.push_back(vecbool2ints(cig.cfgex2vecbool(ex_pos)));
            
            // {// test for the cfig we in put.
                // if(mm) cig.string();
                // if(mm) WRN(NAV(ex_pos));
                // if(mm) show_string(cig.cfg2ex2nums(ex_pos));
            // }
        }
    }
    for_Int(i, 0, core_dim) ranks.push_back(cfi2rank(cfigs[i], cfigs));

    {
        VecReal probability(core_dim); for_Int(i, 0, core_dim) probability[i] = grndste_norm[groundstate_idx[i]];
        if(mm) WRN(NAV(probability.truncate(0,100)));

        VecReal norm = SQR(Vec(ranks));
        VecIdx ranks_idx(sort_indexes(norm));
        VecReal ranks_sort(core_dim); for_Int(i, 0, core_dim) ranks_sort[i] = ranks[ranks_idx[i]];
        if(mm) WRN(NAV(ranks_sort.truncate(0,100)));
    }
    if(mm) WRN("Done here, finish the git_nci() clc, need to expand"+NAV2(cfigs.size(),ranks.size()))
    // expand(natural_cfg); cfigs.shrink_to_fit(); ranks.shrink_to_fit();
    return natural_cfg;
}

// get the  FULL norg ground state's ex-state
Nci Asnci::git_nci_no_rank(const VecReal& ground_state, const Int ex_pos) {
    Nci                         natural_cfg;
    VEC<array<UInt,6>>&         cfigs(natural_cfg.first);
    VEC<Real>&                  coefs(natural_cfg.second);

    VecReal grndste_norm = SQR(ground_state);
    VecIdx groundstate_idx(sort_indexes(grndste_norm));

    for_Int(i, 0, ground_state.size()) {
        const Idx ci_idx = groundstate_idx[i];
        StateStatistics cig(ci_idx, nosp.wherein_NocSpace(ci_idx), nosp);

        if(check_ifinex(cig, ex_pos)){
            cfigs.push_back(vecbool2ints(cig.cfgex2vecbool(ex_pos)));
            // {// test for the cfig we in put.
                // if(mm) cig.string();
                // if(mm) WRN(NAV(ex_pos));
                // if(mm) WRN(NAV3(ci_idx, cig.show_cfg(), show_string(vecbool2ints(cig.cfgex2vecbool(ex_pos)))))
            // }
            coefs.push_back(ground_state[ci_idx]);
        }
    }

    if(mm) WRN("Done here, finish the git_nci() clc, need to expand"+NAV2(cfigs.size(),coefs.size()))
    if(mm) WRN(NAV(coefs[core_dim - 1]));
    // expand_no_rank(natural_cfg, groundstate_idx, ex_pos); cfigs.shrink_to_fit(); coefs.shrink_to_fit();
    return natural_cfg;
}


void Asnci::expand(Nci& natural_cfgs) {
    VEC<array<UInt,6>>&        cfigs(natural_cfgs.first);
    VEC<Real>&                 ranks(natural_cfgs.second);

    Vec<VecBool> cfigs_core(cfigs.size());
    for_Idx(i, 0, cfigs.size()){
        cfigs_core[i] = ints2vecbool(cfigs[i]);
    }
    for_Idx(i, 0, cfigs_core.size()){
        for(const auto j : mayhop) if(judge(cfigs_core[i], j)) {
            VecBool new_cfig(change_cfg(cfigs_core[i], j));
            Real rank = cfi2rank(new_cfig, cfigs_core);
            if(rank > 0.){
                ranks.push_back(cfi2rank(new_cfig, cfigs_core));
                cfigs.push_back(vecbool2ints(new_cfig));
            }


            {// test
                // VecBool a(80, false);
                // VecInt b(10, 2);
                // array<UInt,10> c;
                // if(mm) WRN(NAV3((a.szof()), sizeof(b.szof()), sizeof(c)));
                // if(mm) show_string(cfigs[i]);
                // if(mm) {Int orbit_idx(j/(p.nI2B[0]+1)%p.norbs),ctr(j/(p.norbs*(p.nI2B[0]+1))%(p.nI2B[0]+1)), ann(j%(p.norbs*(p.nI2B[0]+1))%(p.nI2B[0]+1));
                //         WRN(NAV5(i, SQR(ranks[i]), orbit_idx, ctr, ann));}
                // if(mm) show_string(vecbool2ints(new_cfig));
                // if(mm) show_string(vecbool2ints(ints2vecbool(cfigs[i])));
                // if(mm) WRN(NAV3(j,j/(p.nI2B[0]+1), j%(p.nI2B[0]+1)));
                // if(mm) show_string(vecbool2ints(new_cfig));
                // if(mm) show_string(vecbool2ints(ints2vecbool(vecbool2ints(new_cfig))));
            }
        }
    }
    if(mm) WRN(NAV(cfigs.size()));
}

void Asnci::expand_no_rank(Nci& natural_cfgs, const VecIdx& prob_idx, const Int& ex_pos) {
    VEC<array<UInt,6>>&        cfigs(natural_cfgs.first);
    VEC<Real>&                 coefs(natural_cfgs.second);

    Vec<VecBool> cfigs_core(core_dim);
    for_Idx(i, 0, core_dim){
        cfigs_core[i] = ints2vecbool(cfigs[i]);
    }
    for_Idx(i, 0, core_dim){
        for(const auto j : mayhop) if(judge(cfigs_core[i], j)) {
            // if(mm) WRN(NAV(show_string(vecbool2ints(cfigs_core[i]))));
            VecBool new_cfig(change_cfg(cfigs_core[i], j));

            // if(mm) WRN(NAV(show_string(vecbool2ints(new_cfig))));

            // if(back2nospace(new_cfig, -ex_pos).size() > 1){
                // if(mm) WRN(NAV(show_string(vecbool2ints(back2nospace(new_cfig, -ex_pos)))));}
            // else if(mm) WRN("we skip this cfig.");

            VecBool new_cfig_innorg (back2nospace(new_cfig, -ex_pos));
            if (new_cfig_innorg.size() > 1 && !nosp.ifin_NocSpace(new_cfig_innorg, nosp.nppso))
            {
                // if(mm) WRN(NAV(show_string(vecbool2ints(new_cfig))));
                cfigs.push_back(vecbool2ints(new_cfig));
                coefs.push_back(0.);
            }
        }
    }
    if(mm) WRN(NAV2(cfigs.size(), coefs.size()));
}

Nci Asnci::truncation(Nci inital) {
    VEC<array<UInt,6>>  cfigs;
    VEC<Real>           ranks;

    VecReal norm = SQR(Vec(inital.second));
    VecIdx groundstate_idx(sort_indexes(norm));
    dim = dim < inital.first.size() ? dim : inital.first.size();
    if(mm) WRN(NAV(dim));
    for_Int(cunt, 0, dim){// add the map
        const Idx& ci_idx = groundstate_idx[cunt];
        cfig_idx.insert(pair<array<UInt,6>, Int>(inital.first[ci_idx], ci_idx));
        cfigs.push_back((inital.first[ci_idx]));
        ranks.push_back((inital.second[ci_idx]));
    }
    cfigs.shrink_to_fit(); ranks.shrink_to_fit();

    // if(mm) WRN(NAV(Vec(ranks).truncate(0,1024)));
    return pair(cfigs, ranks);
}

Real Asnci::cfi2rank(const VecBool& alpha, const Vec<VecBool>& beta) {
    Real rank(0.);
    for_Idx(i, 0, beta.size()){
        Real upper = hamilton_value(alpha, beta[i]);
        rank += upper / (groundE - hamilton_value(alpha));
    }
    return rank;
}

Real Asnci::cfi2rank(const array<UInt,6>& alpha_number, const VEC<array<UInt,6>>& betas_number) {
    Real rank(0.);

    VecBool alpha(ints2vecbool(alpha_number));
    Vec<VecBool> cfigs_core(betas_number.size());
    for_Idx(i, 0, betas_number.size()) cfigs_core[i] = ints2vecbool(betas_number[i]);

    rank = cfi2rank(alpha, cfigs_core);
    return rank;
}

// ?! The H_a,b need to consider the Anticommutativity？
Real Asnci::hamilton_value(const VecBool& alpha, const VecBool& beta_i) {
    Real value(0.), uz(p.U), jz = (p.jz);
    Idx nimp(p.norbs), nband(p.nband), norb(p.norbit);
    if(beta_i.size() == 0 || (alpha == beta_i)) {
        // Vec<Char> a_cfig(alpha.size()); for_Int(i, 0, a_cfig.size()) a_cfig[i] = alpha[i];
        VecBool a_cfig(alpha); 

        // add the on site energy.
        for_Int(i, 0, norb) {
                if(a_cfig[i]) value += impH.first[i][i];
            }

        VecBool impcfig = a_cfig.mat(nband, p.nI2B[0] + 1).tr()[0];

        // add the U.
        for_Int(i, 0, nband) if (impcfig[2 * i] && impcfig[2 * i + 1]) value += uz;

        // add the up-down term.
        for_Int(i, 0, nband)
            for_Int(j, 0, nband) if (i != j && impcfig[2 * i] && impcfig[2 * j + 1])
                value += uz - 2 * jz;

        // add the up-up term.
        for_Int(i, 0, nband)
            for_Int(j, 0, nband) if (i != j && impcfig[2 * i] && impcfig[2 * j])
                value += uz - 3 * jz;

        // add the down-down term.
        for_Int(i, 0, nband)
            for_Int(j, 0, nband) if (i != j && impcfig[2 * i + 1] && impcfig[2 * j + 1])
                value += uz - 3 * jz;

        return value;
    }
    else {
        Int crt(-1), ann(-1);
        for_Int(i, 0, norb) if(alpha[i] ^ beta_i[i])  beta_i[i] ? ann = i : crt = i;
        if(crt >= 0 && ann >= 0) value = crt > ann ? impH.first[crt][ann] : - impH.first[crt][ann];
        else {
            // if(mm) {
            //     WRN("some thing wrong in here!");
            //     show_string(vecbool2ints(alpha)); show_string(vecbool2ints(beta_i));}
            // ERR("some thing wrong in here!"+NAV2(alpha.string(),beta_i.string()));
            }
        return value;
    }
}


Tab Asnci::find_fullH_idx()
{
	clock_t t_find_hmlt_table;
	t_find_hmlt_table = clock(); 
    // TIME_BGN("find_hmlt_table" + NAV(mm.id()), t_find_hmlt_table);
	VecPartition row_H(mm.np(), mm.id(), dim);
	Tab h_idxs(3);
	MatInt mat_hop_pos(impH.first.nrows(), impH.first.ncols());
	for_Int(i, 0, mat_hop_pos.nrows()) for_Int(j, 0, mat_hop_pos.ncols()) mat_hop_pos[i][j] = i * mat_hop_pos.ncols() + j;
	Int h_hbd_idx(mat_hop_pos.size()	+ 1);
	Int h_orb_ud_idx(mat_hop_pos.size() + 2);
	Int h_orb_uu_idx(mat_hop_pos.size() + 3);
	Int h_orb_dd_idx(mat_hop_pos.size() + 4);

	// Int h_orb_j_idx(mat_hop_pos.size() + 5);

    const VEC<array<UInt,6>>&   cfigs(trncat.first);
    const VEC<Real>&                  ranks(trncat.second);

    Real uz(p.U), jz = (p.jz);
    Idx nimp(p.norbs), nband(p.nband), norb(p.norbit);
    
    Idx cut_H_need(0);
	for_Int(h_i, row_H.bgn(), row_H.end()) {
        VecBool a_cfig(Vec(ints2vecbool(cfigs[h_i])));
		
		// To save as sparse matrix, [0]: row number;[1]: colum number;[2]: idx.
		VecInt h_idx(3, 0);
		Int sparse_idx(h_i - row_H.bgn());		
    // ! For the diagonal term:
        // add the on site energy.
        for_Int(i, 0, norb) {
                if(a_cfig[i]) {
                    h_idx = { sparse_idx, h_i, mat_hop_pos[i][i] + 1 };
                    for_Int(pos, 0, 3) h_idxs[pos].push_back(h_idx[pos]);
                }
            }

        // add the interation.
        VecBool impcfig = a_cfig.mat(nband, p.nI2B[0] + 1).tr()[0];

        // add the U.
        for_Int(i, 0, nband) if (impcfig[2 * i] && impcfig[2 * i + 1]) {
			h_idx = { sparse_idx, h_i, h_hbd_idx};
			for_Int(pos, 0, 3) h_idxs[pos].push_back(h_idx[pos]);
		}

        // add the up-down term.
        for_Int(i, 0, nband)
            for_Int(j, 0, nband) if (i != j && impcfig[2 * i] && impcfig[2 * j + 1]){
				h_idx = { sparse_idx, h_i, h_orb_ud_idx};
				for_Int(pos, 0, 3) h_idxs[pos].push_back(h_idx[pos]);
			}

        // add the up-up term.
        for_Int(i, 0, nband)
            for_Int(j, 0, nband) if (i != j && impcfig[2 * i] && impcfig[2 * j]){
				h_idx = { sparse_idx, h_i, h_orb_uu_idx};
				for_Int(pos, 0, 3) h_idxs[pos].push_back(h_idx[pos]);
			}

        // add the down-down term.
        for_Int(i, 0, nband)
            for_Int(j, 0, nband) if (i != j && impcfig[2 * i + 1] && impcfig[2 * j + 1]){
				h_idx = { sparse_idx, h_i, h_orb_dd_idx};
				for_Int(pos, 0, 3) h_idxs[pos].push_back(h_idx[pos]);
			}
    // ! off diagonal term:
        for(const auto j : mayhop) if(judge(a_cfig, j)) {
            VecBool b_cfig(change_cfg(a_cfig, j));
            Int crt(-1), ann(-1);
            for_Int(i, 0, norb) if(a_cfig[i] != b_cfig[i]) b_cfig[i] ? ann = i : crt = i;
            Int Anticommutativity = crt > ann ? 1 : - 1;
            array<UInt,6>   nums(vecbool2ints(b_cfig));

            if(crt >= 0 && ann >= 0) {
                auto it = cfig_idx.find(nums);
                if (it != cfig_idx.end()) {
                    	h_idx = {sparse_idx, it->second, Anticommutativity * (mat_hop_pos[crt][ann] + 1)};
                    	for_Int(pos, 0, 3) h_idxs[pos].push_back(h_idx[pos]);
                } else {
                    cut_H_need++;
                }
            }
        }
	}
    if(mm) WRN(NAV(cut_H_need));
	// TIME_END("t_find_hmlt_table" + NAV(mm.id()), t_find_hmlt_table);
	for_Int(jj, 0, 3) h_idxs[jj].shrink_to_fit();
	return std::move(h_idxs);
}