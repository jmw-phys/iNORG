#pragma once
/*
coded by Jia-Ming Wang (jmw@ruc.edu.cn, RUC, China) date 2023.03.03
*/


#include <bitset>
#include <string>
#include"norg.h"

// adaptive sampling natural configuration interaction (Asnci)
// natural configuration interaction (Nci)

// typedef pair<VEC<Str>,VEC<Real>> Nci;
// typedef pair<VEC<__uint128_t>,VEC<Real>> Nci;
// typedef std::pair<VEC<UInt>,VEC<Real>> Nci;
// typedef std::pair<MatInt,VecReal> Nci;


typedef std::pair<VEC<std::array<UInt,10>>,VEC<Real>> Nci;
typedef Vec<VEC<Int>> Tab;

class Asnci 
{
	const MyMpi& mm;				// parameters
	const Prmtr& p;					// parameters
	const NocSpace& nosp;			// parameters
    const Real& groundE;            // The final ground state energy

    VEC<Int> mayhop;                // the position of orbital hop may happen
    
public:
    const MatReal hop_h;            // The hopping H(H_0)
    const Idx core_dim;             // The core space size
    const Idx dim;                  // The truncated space size
    Nci trncat;                     // The truncation NCI
    // The cfig's idx
    std::map<std::array<UInt,10>, Int> cfig_idx;

private:

    VEC<Int> find_mayhop();
    
    Nci git_nci(const VecReal& ground_state);

    Nci git_nci(const VecReal& ground_state, const Int ex_pos);
    
    void expand(Nci& natural_cfgs);

    Str to_binary_string(unsigned long num) {
        using namespace std;
        bitset<sizeof(unsigned long) * 8> binary(num);
        return binary.to_string();
    }

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

    bool check_ifinex(const StateStatistics& cig, const Int& ex_pos) {
        // bool flag(true);
        Int pos(ABS(ex_pos));
        if(ex_pos > 0 && cig.occ_n.tr()[0][pos] == 0) true;
        if(ex_pos < 0 && cig.occ_n.tr()[0][pos] == 1) true;
        return false;
    }

    Str change_cfg_str(const Str& cfg_str, Int pos) {
        Int crt(pos/hop_h.ncols()), ann(pos%hop_h.ncols());
        Str temp_c = cfg_str;
        temp_c[ann] = (temp_c[ann]=='1') ? '0' : 'x';
        temp_c[crt] = (temp_c[crt]=='0') ? '1' : 'x';
        return temp_c;
    }

    // The hamilton form:
    Real hamilton_value(const Str alpha, const Str beta_i = Str());

    Real cfi2rank(Str alpha, Vec<Str> beta) {
        Real rank(0.);
        for_Idx(i, 0, beta.size()){
            Real upper = hamilton_value(alpha, beta[i]);
            rank += upper / (groundE - hamilton_value(alpha));
        }
        return rank;
    }

    Nci truncation(Nci inital);

    Tab find_table(Str inter_type);

    VecReal find_ex_state() {return VecReal(trncat.second); }
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
Tab find_h_idx();

void asnci_gimp(Green& imp_i, Int pos);

};