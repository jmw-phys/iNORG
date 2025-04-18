#include "impurity.h"

Impurity::Impurity(const MyMpi &mm_i, const Prmtr &prmtr_i, const Bath &bth_i, const Str& file)
    : mm(mm_i), p(prmtr_i), bth(bth_i), nb(p.nbath), ni(p.norbs), ns(p.norbit), pos_imp(p.norbs), h0(p.norbit, p.norbit, 0.), imp_lvl(p.norbs)
{
    // if (!file.empty()) read(file);
    for_Int(i, 0, ni) imp_lvl[i] = p.eimp[i] - p.mu;
}

Impurity::Impurity(const MyMpi &mm_i, const Prmtr &prmtr_i, const Bath &bth_i, const VecInt or_deg)
    : mm(mm_i), p(prmtr_i), bth(bth_i), nb(p.nbath), ni(p.norbs), ns(p.norbit), pos_imp(p.norbs), h0(p.norbit, p.norbit, 0.), imp_lvl(p.norbs, 0.)
{
    // if (!file.empty()) read(file);
    // VecInt ordeg(concat(or_deg,or_deg).mat(2, p.nband).tr().vec());
    VecInt ordeg(or_deg);
    
    VecReal deg_lvl(MAX(ordeg), 0.);
    for_Int(i, 0, ni) deg_lvl[ordeg[i] - 1] += p.eimp[i] - p.mu;
    for_Int(i, 0, MAX(ordeg)) {
        Int cnt(0);
        for_Int(j, 0, ni) if(i == ordeg[j] - 1) cnt++;
        deg_lvl[i] = deg_lvl[i] / Real(cnt);
    }
    for_Int(i, 0, ni) imp_lvl[i] = deg_lvl[ordeg[i] - 1];
}

using namespace std;

// rely on imp model's frame
void Impurity::find_g0(Green &g0) const {

    MatCmplx Z(ns, ns);
    for_Int(n, 0, g0.nomgs) {
        Z = g0.z(n);
        MatCmplx g0_z = matinvlu(Z - cmplx(h0));
        for_Int(i, 0, p.nband) {
            g0[n][i][i] = g0_z[pos_imp[2 * i]][pos_imp[2 * i]];
        }
    }
}

// rely on imp model's frame
void Impurity::find_all_g0(Green &g0) const {

    MatCmplx Z(ns, ns);
    for_Int(n, 0, g0.nomgs) {
        Z = g0.z(n);
        MatCmplx g0_z = matinvlu(Z - cmplx(h0));
        for_Int(i, 0, p.norbit) {
            g0[n][i][i] = g0_z[i][i];
        }
    }
}

//rely on imp model's frame
void Impurity::find_hb(Green &hb) const {
    for_Int(i,0,p.nband){
        const Int nb_i = p.nI2B[2*i];
        for_Int(n,0,hb.nomgs){
            const VecCmplx Z(nb_i, hb.z(n));
            VecCmplx S = INV(Z - cmplx(bth.vec_ose[i]));
            VecCmplx V = cmplx(bth.vec_hop[i]);
            hb[n][i][i] = -SUM(V * S * V.co());
            //hb[n][i][i] = DOT(cmplx(bth.vec_hop[i]), S * cmplx(bth.vec_hop[i]));
        }
    }
}

MatReal Impurity::find_hop_for_test() const
{
    VecReal hop({ -0.312546, 0.159346, -0.0619612, -0.312546, 0.159346 });
    VecReal ose({ -0.474855,0.0667285, 0,           0.474855, -0.0667285 });
    MatReal h0(0, 0, 0.);
    for_Int(i, 0, p.nband) {
        const Int nb_i = p.nI2B[i * 2];
        MatReal h0_i(1 + nb_i, 1 + nb_i, 0.);
        for_Int(j, 0, nb_i) {
            h0_i[0][j + 1] = hop[j];
            h0_i[j + 1][0] = hop[j];
            h0_i[j + 1][j + 1] = ose[j];
        }
        h0_i.reset(direct_sum(h0_i, h0_i));
        h0.reset(direct_sum(h0, h0_i));
    }
    return h0;
}

void Impurity::update(Str mode) {
    if (mode.empty()) {
        set_factor();
        impH = std::make_pair(h0, set_interaction());
        // modify_Impdata_for_half_fill(impH);
    }
    else if (mode == "behte") {
        for_Int(i, 0, ni) imp_lvl[i] = p.eimp[i] - p.mu;
        set_factor();
        impH = std::make_pair(h0, set_interaction());
        modify_Impdata_for_half_fill(impH);
    }
    else if (mode == "behte_alpha") {
        for_Int(i, 0, ni) imp_lvl[i] = p.eimp[i] - p.mu;
        set_factor();
        impH = std::make_pair(h0, set_3band_interaction_withalpha());
        // modify_Impdata_for_half_fill(impH);
        modify_Impdata_for_half_fill_hhd(impH);
    }
    // if(mm) WRN(NAV(h0));
}

void Impurity::write_H0info(const Bath &b, Int ndeg, Int iter_cnt) const {
    OFS ofs;
    if(iter_cnt < 0) ofs.open("h0.txt");
    if(iter_cnt > 0) ofs.open(iox + "zic" + prefill0(iter_cnt, 3) +".h0.txt");
    using namespace std;
    if(ndeg > 0) for_Int(i, 0, ndeg)	{
        ofs << "degeneracy band "<< i+1  << "nmin: " << b.info[i][0] << " err: " << b.info[i][1] << " err_crv: " << b.info[i][2] << " err_regE: " << b.info[i][3] << " err_regV: " << b.info[i][4] << " err_bsr: " << b.info[i][5] <<" norm: " << b.info[i][6]<< "  " << endl;
    }
    else for_Int(i, 0, p.nband) {
        ofs << "band "<< i+1  << "nmin: " << b.info[i][0] << " err: " << b.info[i][1] << " err_crv: " << b.info[i][2] << " err_regE: " << b.info[i][3] << " err_regV: " << b.info[i][4] << " err_bsr: " << b.info[i][5] <<" norm: " << b.info[i][6]<< "  " << endl;
    }
    for_Int(i, 0, p.nband)	{
        ofs << "band "<< i+1 << endl;
        Int begin(i*2 * (p.nI2B[i*2] + 1)), end((i*2 + 1) * (p.nI2B[i*2] + 1));
        ofs << iofmt() << h0.truncate(begin, begin, end, end) << endl;
    }
}

//---------------------------------------------Private function---------------------------------------------


void Impurity::set_factor() {
    // set hyb part and bath part
    h0 = bth.find_hop();
    
    // // h0 = find_hop_for_test();
    // set imp part
    MatReal h0loc(ni,ni,0.);
    for_Int(i, 0, ni) h0loc[i][i] = imp_lvl[i];
    
    // find i_th imp in which site
    Int site = 0;
    for_Int(i, 0, ni) {
        pos_imp[i] = site;
        site += p.nI2B[i] + 1;
    }

    // set h0loc
    for_Int(i, 0, ni) {
        for_Int(j, 0, ni) {
            h0[pos_imp[i]][pos_imp[j]] = h0loc[i][j];
        }
    }
}

// four-fermion operator terms for C^+_i C^+_j C_k C_l; 
// C^+_i C^+_j C_k C_l h_inter from [i][l][j][k] to [alpha][eta][beta][gamma]
VecReal Impurity::set_interaction() {
    Int n = p.norbit;
    VecReal interaction(std::pow(n, 4), 0.);
    // Mat<MatReal> imp_interact(p.norbs, p.norbs, MatReal(p.norbs, p.norbs, 0.));
    MatReal imp_dd_interact(p.norbs, p.norbs,  0.);


    for_Int(b1, 0, p.nband) {// NO double counting term for impurity
        imp_dd_interact[2 * b1][2 * b1 + 1] = p.U;
        for_Int(b2, b1, p.nband) if (b1 != b2) { // same spin orientation
            imp_dd_interact[2 * b1][2 * b2] = p.Uprm - p.jz;
            imp_dd_interact[2 * b1 + 1][2 * b2 + 1] = p.Uprm - p.jz;
        }
        for_Int(b2, 0, p.nband) if (b1 != b2) {
            imp_dd_interact[2 * b1][2 * b2 + 1] = p.Uprm;
        }
    }

/* 
    for_Int(b1, 0, p.nband) {// with double counting term for impurity, so we need to divide by 2
        imp_dd_interact[2 * b1][2 * b1 + 1] = p.U/2.0;
        imp_dd_interact[2 * b1 + 1][2 * b1] = p.U/2.0;
        for_Int(b2, 0, p.nband) if (b1 != b2) { // same spin orientation
            imp_dd_interact[2 * b1][2 * b2] = (p.Uprm - p.jz)/2.0;
            imp_dd_interact[2 * b1 + 1][2 * b2 + 1] = (p.Uprm - p.jz)/2.0;
        }
        for_Int(b2, 0, p.nband) if (b1 != b2) {
            imp_dd_interact[2 * b1][2 * b2 + 1] = p.Uprm/2.0;
            imp_dd_interact[2 * b1 + 1][2 * b2] = p.Uprm/2.0;
        }
    }
*/
    
    if(mm) WRN(NAV(imp_dd_interact));
    for_Int(N_i, 0, p.norbs)for_Int(N_j, 0, p.norbs) {
        // if (imp_dd_interact[N_i][N_j] != 0)
            interaction[SUM_0toX(p.nO2sets, N_i) * std::pow(n, 3) + SUM_0toX(p.nO2sets, N_j) * std::pow(n, 2) + SUM_0toX(p.nO2sets, N_j) * std::pow(n, 1) + SUM_0toX(p.nO2sets, N_i)] = imp_dd_interact[N_i][N_j];

    }
    // if (mm) WRN(NAV5(p.U, p.Uprm ,p.nO2sets, interaction.size(), SUM(interaction[0][0])));
    return interaction;
}

void Impurity::modify_Impdata_for_half_fill(Impdata& impH_i){
    //! modify the non-interacting part of impurity Hamiltonian
	MatReal& h0 = impH_i.first;
	Int norb2set = p.nO2sets[0];
	for_Int(i, 0, p.norg_sets) {
        if      (p.nband == 1)h0[i * norb2set][i * norb2set] -= p.U * 0.5;
        else if (p.nband == 2)h0[i * norb2set][i * norb2set] -= p.U * 0.5 + p.Uprm * 0.5 + (p.Uprm - p.jz) * 0.5;
        else if (p.nband == 3)h0[i * norb2set][i * norb2set] -= p.U * 0.5 + 2 * (p.Uprm * 0.5) + 2 * (p.Uprm - p.jz) * 0.5;
        else if (p.nband  > 3)h0[i * norb2set][i * norb2set] -= p.U * 0.5 + (p.nband-1) * (p.Uprm * 0.5) + (p.nband-1) * (p.Uprm - p.jz) * 0.5;
        else ERR("impurity input is wrong");
	}
}

//--------------------------------------------------------special for the hhd function(arXiv:2209.14178v1)-----------------------------------------------------------------


// four-fermion operator terms for C^+_i C^+_j C_k C_l; 
// C^+_i C^+_j C_k C_l h_inter from [i][l][j][k] to [alpha][eta][beta][gamma]
VecReal Impurity::set_3band_interaction_withalpha() {
    Int n = p.norbit;
    VecReal interaction(std::pow(n, 4), 0.);
    if (p.nband != 3) ERR("this function is only for 3-band model");
    // Mat<MatReal> imp_interact(p.norbs, p.norbs, MatReal(p.norbs, p.norbs, 0.));
    MatReal imp_dd_interact(p.norbs, p.norbs,  0.);

    for_Int(b1, 0, p.nband) {// NO double counting term for impurity
        imp_dd_interact[2 * b1][2 * b1 + 1] = p.U;
        for_Int(b2, b1, p.nband) if (b1 != b2) { // same spin orientation
            imp_dd_interact[2 * b1][2 * b2] = p.Uprm - p.jz;
            imp_dd_interact[2 * b1 + 1][2 * b2 + 1] = p.Uprm - p.jz;
        }
        for_Int(b2, 0, p.nband) if (b1 != b2) {
            imp_dd_interact[2 * b1][2 * b2 + 1] = p.Uprm;
        }
    // setting the interation between two wide bands.
        for_Int(b2, b1, p.nband) if (b1 == 0 && b2 == 1) { // same spin orientation
            imp_dd_interact[2 * b1][2 * b2] *= p.alpha;
            imp_dd_interact[2 * b1 + 1][2 * b2 + 1] *= p.alpha;
        }
        for_Int(b2, 0, p.nband) if (b1 == 0 && b2 == 1 || b2 == 0 && b1 == 1 ) {
            imp_dd_interact[2 * b1][2 * b2 + 1] *= p.alpha;
        }
    }



    if(mm) WRN(NAV(imp_dd_interact));
    for_Int(N_i, 0, p.norbs)for_Int(N_j, 0, p.norbs) {
        // if (imp_dd_interact[N_i][N_j] != 0)
            interaction[SUM_0toX(p.nO2sets, N_i) * std::pow(n, 3) + SUM_0toX(p.nO2sets, N_j) * std::pow(n, 2) + SUM_0toX(p.nO2sets, N_j) * std::pow(n, 1) + SUM_0toX(p.nO2sets, N_i)] = imp_dd_interact[N_i][N_j];

    }
    // if (mm) WRN(NAV5(p.U, p.Uprm ,p.nO2sets, interaction.size(), SUM(interaction[0][0])));
    return interaction;
}

void Impurity::modify_Impdata_for_half_fill_hhd(Impdata& impH_i){
    //! modify the non-interacting part of impurity Hamiltonian
	MatReal& h0 = impH_i.first;
	Int norb2set = p.nO2sets[0];
	for_Int(i, 0, p.norg_sets) {
        h0[i * norb2set][i * norb2set] -= p.U * 0.5;
        // interorbital part
        // special chage for the 3 band hhd case. (arXiv:2209.14178v1)
        Real iner_part_per_band = ((p.nband - 1) * (p.Uprm * 0.5) + (p.nband - 1) * (p.Uprm - p.jz) * 0.5) / (p.nband - 1);
        if (i/2 < 2) h0[i * norb2set][i * norb2set] -= (iner_part_per_band * p.alpha) + iner_part_per_band * (p.nband - 2);
        if (i/2 >= 2) h0[i * norb2set][i * norb2set] -= iner_part_per_band * (p.nband - 1);
	}
}



//---------------------------------------------------------------------------------------------------------------------------------------