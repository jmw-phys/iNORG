#pragma once
/*
code developed and maintained by (jmw@ruc.edu.cn, RUC, China) date 2023.03.03
*/

#include <bitset>
#include <algorithm>
#include <vector>
#include <string>
#include"norg.h"

// adaptive sampling natural configuration interaction (Asnci)
// natural configuration interaction (Nci)

// typedef pair<VEC<Str>,VEC<Real>> Nci;
// typedef pair<VEC<__uint128_t>,VEC<Real>> Nci;
// typedef std::pair<VEC<UInt>,VEC<Real>> Nci;
// typedef std::pair<MatInt,VecReal> Nci;


typedef std::pair<VEC<std::array<UInt,6>>,VEC<Real>> Nci;
typedef Vec<VEC<Int>> Tab;

// ! NOW this class need to modify since need rotation the whole impurity.
class Asnci 
{
	const MyMpi& mm;				// parameters
	const Prmtr& p;					// parameters
	const NocSpace& nosp;			// parameters
    const Real& groundE;            // The final ground state energy
    const Impdata impH;             // The impurity structure
    const VecReal coefficient;      // coefficient for all the H's terms

    VEC<Int> mayhop;                // the position of orbital hop may happen
    
public:
    const Idx core_dim;             // The core space size
    Idx dim;                        // The truncated space size
    Nci trncat;                     // The truncation NCI
    // The cfig's idx
    std::map<std::array<UInt,6>, Int> cfig_idx;

private:
    inline VecIdx sort_indexes(const VecReal& vec) {
        VecIdx idx(vec.size());
        std::iota(idx.begin(), idx.end(), 0);

        std::sort(idx.begin(), idx.end(),
                [&vec](int i, int j) { return vec[i] > vec[j]; });

        return idx;
    }

    VEC<Int> find_mayhop();
    
    Nci git_nci(const VecReal& ground_state);

    Nci git_nci(const VecReal& ground_state, const Int ex_pos);

    Nci git_nci_no_rank(const VecReal& ground_state, const Int ex_pos);
    
    void expand(Nci& natural_cfgs);

    VecBool back2nospace(const VecBool& new_cfig, const Int change_pos);

    void expand_no_rank(Nci& natural_cfgs, const VecIdx& prob_idx, const Int& ex_pos);

    // Str to_binary_string(unsigned long num) {
    //     using namespace std;
    //     bitset<sizeof(unsigned long) * 8> binary(num);
    //     return binary.to_string();
    // }

    VecBool ints2vecbool(std::array<UInt,6> nums);

	std::array<UInt,6> vecbool2ints(const VecBool &vec);

/*
    bool judge(Str& cfg_str) {
        for_Idx(i, 0, cfg_str.size()) if((cfg_str[i] != '1' && cfg_str[i] != '0')) return false;
        return true;
    }

    bool judge(const Str& cfg_str,const Int& pos) {
        Int crt(pos/hop_h.ncols()), ann(pos%hop_h.ncols());
        if(cfg_str[ann]!='1') return false;
        if(cfg_str[crt]!='0') return false;
        return true;
    }
*/

    bool judge(const VecBool& cfg, const Int& pos) {
        MatReal hop_h = impH.first;
        Int crt(pos / hop_h.ncols()), ann(pos % hop_h.ncols());
        if((cfg[ann]) && (!cfg[crt]))return true;
        return false;
    }

    bool check_ifinex(const StateStatistics& cig, const Int& ex_pos) {
        // bool flag(true);
        Idx pos(ABS(ex_pos) - 1);
        if(ex_pos > 0 && cig.cfg.cf[pos].isuno(0)) return true;
        if(ex_pos < 0 && cig.cfg.cf[pos].isocc(0)) return true;
        return false;
    }


    VecBool change_cfg(const VecBool& cfg, Int pos);

    // The hamilton form:
    Real hamilton_value(const VecBool& alpha, const VecBool& beta_i = VecBool());

    Real cfi2rank(const VecBool& alpha, const Vec<VecBool>& beta);

    Real cfi2rank(const std::array<UInt,6>& alpha_number, const VEC<std::array<UInt,6>>& betas_number);

    Nci truncation(Nci inital);

    VecReal find_ex_state() {return VecReal(trncat.second);}

    Mat<Char> show_string(const std::array<UInt,6>& cifg);
public:

// Asnci: mode = 0: assume NO converged;  
//        mode = 1: NO converg with adaptive sampling configuration interaction.
// Asnci(const NORG& norg, Idx trncat_size, const Int mode = 0);
Asnci(const NORG& norg, Idx trncat_size);

// to use the ASNCI here we suppose the we alreadly find the ground state from NORG.
Asnci(const NORG& norg, Idx trncat_size, Int ex_pos);

// out put the NORG class with the table
NORG get_norg(Tab table, Int mode = 0);

// out put the Table
Tab find_fullH_idx();

void asnci_gimp(Green& imp_i, Int pos);

};