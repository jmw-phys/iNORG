/*
coded by Jia-Ming Wang (jmw@ruc.edu.cn, RUC, China) date 2022 - 2023
*/

#include "occler.h"
using namespace std;

Occler::Occler(const MyMpi& mm_i, Prmtr& prmtr_i):
    mm(mm_i),// nparticals(prmtr_i.npartical), // p(prmtr_i), 
    sub_energy({1., 0., 0.}), tp(prmtr_i)
{
    // initial_uormat();
}


// old version
NORG Occler::find_ground_state_partical(const Impdata &impH_i){
    Int counter(0);
    while(1){
        if(mm) std::cout << "The " << ++counter << "-th NORG begin" << std::endl;	// norg counter
        NORG a(mm, tp);                          IFS ifs_a("ru"+nppso_str(a.scsp.nppso)+".bi"); 
        if(ifs_a) {for_Int(i, 0, a.uormat.size()) biread(ifs_a, CharP(a.uormat[i].p()), a.uormat[i].szof());}
        a.up_date_h0_to_solve(impH_i);            sub_energy[0] = a.groune_lst;
        if(mm){
            OFS ofs_a; ofs_a.open("ru"+nppso_str(a.scsp.nppso)+".bi"); 
            for_Int(i, 0, a.uormat.size()) biwrite(ofs_a, CharP(a.uormat[i].p()), a.uormat[i].szof());}

        tp.templet_control[1]--;                 tp.templet_control[tp.ndiv-1]++;  tp.after_modify_prmtr(); tp.recalc_partical_number();
            if(mm) std::cout << "The " << ++counter << "-th NORG begin" << std::endl;	// norg counter
            NORG a_m(mm, tp);                    IFS ifs_m("ru"+nppso_str(a_m.scsp.nppso)+".bi"); 
            if(ifs_m) {for_Int(i, 0, a_m.uormat.size()) biread(ifs_m, CharP(a_m.uormat[i].p()), a_m.uormat[i].szof());}
            a_m.up_date_h0_to_solve(impH_i);      sub_energy[1] = a_m.groune_lst;
        if(mm){
            OFS ofs_m; ofs_m.open("ru"+nppso_str(a_m.scsp.nppso)+".bi"); 
            for_Int(i, 0, a_m.uormat.size()) biwrite(ofs_m, CharP(a_m.uormat[i].p()), a_m.uormat[i].szof());}
            
            if(tp.templet_control[1] == 0 && sub_energy[1] < sub_energy[0])  return a_m;
            tp.templet_control[1]++;             tp.templet_control[tp.ndiv-1]--;  tp.after_modify_prmtr(); tp.recalc_partical_number();

        tp.templet_control[1]++;     tp.templet_control[tp.ndiv-1]--;  tp.after_modify_prmtr(); tp.recalc_partical_number();
            if(mm) std::cout << "The " << ++counter << "-th NORG begin" << std::endl;	// norg counter
            NORG a_p(mm, tp);                    IFS ifs_p("ru"+nppso_str(a_p.scsp.nppso)+".bi"); 
            if(ifs_p) {for_Int(i, 0, a_p.uormat.size()) biread(ifs_p, CharP(a_p.uormat[i].p()), a_p.uormat[i].szof());}
            a_p.up_date_h0_to_solve(impH_i);          sub_energy[2] = a_p.groune_lst;
        if(mm){
            OFS ofs_p; ofs_p.open("ru"+nppso_str(a_p.scsp.nppso)+".bi"); 
            for_Int(i, 0, a_p.uormat.size()) biwrite(ofs_p, CharP(a_p.uormat[i].p()), a_p.uormat[i].szof());}

            if(tp.templet_control[tp.ndiv-1] == 0 && sub_energy[2] < sub_energy[0])  return a_p;
            tp.templet_control[1]--;     tp.templet_control[tp.ndiv-1]++;  tp.after_modify_prmtr(); tp.recalc_partical_number();

        Int check = if_ground_state();
        if (check == 0) return a;
        // if (check == 1 && counter == 3) return a_p;
        // if (check == 2 && counter == 3) return a_m;
        if (check == 3) {
            if(mm) WRN("the occler is not converged.");
            return a;
        }
    }
    ERR("There some thing wrong in Occler.cpp!")
}



// ! tested.
NORG Occler::find_ground_state_partical(const Impdata &impH_i, const VecInt& or_deg)
{
    Int counter_norg(0);
    std::map<std::string, Real> sub_energy_data;
    tp.nooc_mode = STR("cnooc");
    while(1){
            Int counter(0);
            VEC<VecInt> nppsos = list_all_posible_nppsos(tp.npartical, or_deg);
            sub_energy.reset(nppsos.size(),0.); sub_energy = 0.;
        for(const auto& nppso: nppsos) {
            if(MIN(nppso)==0 || MAX(nppso)==tp.nO2sets[0]) break;
            if (mm) std::cout << "The " << ++counter_norg << "-th NORG begin" << std::endl; // norg counter

            tp.according_nppso(tp.npartical = nppso);
            auto iter = sub_energy_data.find(nppso.string());
            if (iter != sub_energy_data.end()) {
                sub_energy[counter] = sub_energy_data[nppso.string()];
                if(mm) PIO("The subspace: "+nppso_str(nppso)+", and energy: "+STR(sub_energy[counter]));
            }
            else {
                std::streambuf* old = std::cout.rdbuf();
                std::cout.rdbuf(NULL); // redirect std::cout to nowhere
                
                NORG a(mm, tp);
                // IFS ifs_a("ru" + nppso_str(a.scsp.nppso) + ".bi");
                // if (ifs_a) for_Int(i, 0, a.uormat.size()) biread(ifs_a, CharP(a.uormat[i].p()), a.uormat[i].szof());
                std::cout.rdbuf(old);
                a.read_NTR();
                if(IFS("ru" + nppso_str(a.scsp.nppso) + ".bi")) if(mm) WRN("find it in the private space");
                a.up_date_h0_to_solve(impH_i, sub_energy.truncate(0, counter)); sub_energy[counter] = a.groune_lst;
                sub_energy_data[nppso.string()] = a.groune_lst;
                if(mm) PIO("The subspace: "+nppso_str(nppso)+", and energy: "+STR(sub_energy[counter]));
            }

            if (mm) {
                // OFS ofs_a;
                // ofs_a.open("ru" + nppso_str(a.scsp.nppso) + ".bi"); 
                // for_Int(i, 0, a.uormat.size()) biwrite(ofs_a, CharP(a.uormat[i].p()), a.uormat[i].szof());
            }
            // if(mm) WRN(NAV4(counter_norg, counter, a.groune_lst, a.scsp.nppso.mat(1,nppsos[0].size())));
            counter++;
        }

        if(sub_energy.idx_min() == 0) {

            tp.according_nppso(tp.npartical = nppsos[0]);
            NORG a(mm, tp);

            // IFS ifs_a("ru" + nppso_str(a.scsp.nppso) + ".bi");
            // if (ifs_a) for_Int(i, 0, a.uormat.size()) biread(ifs_a, CharP(a.uormat[i].p()), a.uormat[i].szof());
            a.up_date_h0_to_solve(impH_i);
            return a;
        } else {
            // nppsos.clear();
            // nppsos = list_all_posible_nppsos(nppsos[sub_energy.idx_min()], or_deg);
            // if(mm) WRN(NAV(sub_energy.stdvec().at(MIN(sub_energy))));
            // tp.npartical = nppsos[sub_energy.stdvec().at(MIN(sub_energy))];
            // if(mm) WRN(NAV2(sub_energy.idx_min(),nppsos[sub_energy.idx_min()].mat(1,10)));
            tp.npartical = nppsos[sub_energy.idx_min()];
            // break;
        }
    }
}


VecInt Occler::find_gs_nppso(const Impdata &impH_i, const VecInt& or_deg)
{
    Int counter_norg(0);
    VEC<MatReal> u_temp;
    while(1){
            Int counter(0);
            VEC<VecInt> nppsos = list_all_posible_nppsos(tp.npartical, or_deg);
            sub_energy.reset(nppsos.size(),0.); sub_energy = 0.;
        for(const auto& nppso: nppsos) {
            if(MIN(nppso)==0 || MAX(nppso)==tp.nO2sets[0]) break;
            if (mm) std::cout << "The " << ++counter_norg << "-th NORG begin" << std::endl; // norg counter

            tp.according_nppso(tp.npartical = nppso);
            NORG a(mm, tp);
            // IFS ifs_a("ru" + nppso_str(a.scsp.nppso) + ".bi");
            // if (ifs_a) for_Int(i, 0, a.uormat.size()) biread(ifs_a, CharP(a.uormat[i].p()), a.uormat[i].szof());
            a.up_date_h0_to_solve(impH_i, sub_energy.truncate(0, counter)); sub_energy[counter] = a.groune_lst;
            if (mm) {
                // OFS ofs_a;
                // ofs_a.open("ru" + nppso_str(a.scsp.nppso) + ".bi"); 
                // for_Int(i, 0, a.uormat.size()) biwrite(ofs_a, CharP(a.uormat[i].p()), a.uormat[i].szof());
            }
            // if(mm) WRN(NAV4(counter_norg, counter, a.groune_lst, a.scsp.nppso.mat(1,nppsos[0].size())));
            counter++;
        }

        if (sub_energy.idx_min() == 0) {
            return nppsos[0];
        }
        else {
            tp.npartical = nppsos[sub_energy.idx_min()];
        }
    }
}

//-----------------------------------------private-----------------------------------------


// void Occler::initial_uormat(){
// 	for_Int(i, 0, tp.nI2B.size()) {
// 		MatReal temp(dmat(tp.nI2B[i], 1.));
// 		uormat.push_back(std::move(temp));
// 	}
// }

Int Occler::if_ground_state() {
    if(sub_energy[0] <= sub_energy[1] && sub_energy[0] <= sub_energy[2]) return 0;
    else
    {
        if (sub_energy[1] >= sub_energy[0] && sub_energy[0] > sub_energy[2]) {
            tp.templet_control[1]++;     tp.templet_control[tp.ndiv-1]--;  tp.after_modify_prmtr(); tp.recalc_partical_number();
            return 1;
        }
        if ( sub_energy[1] < sub_energy[0] && sub_energy[0] <= sub_energy[2]) {
            tp.templet_control[1]--;     tp.templet_control[tp.ndiv-1]++;  tp.after_modify_prmtr(); tp.recalc_partical_number();
            return 2;
        };
        if ( sub_energy[1] > sub_energy[0] && sub_energy[0] > sub_energy[2])
            ERR("The form of energy is no longer dominated by the quadratic form, and the Occler class needs to be changed.");
    }
    return 3;
}

// bool Occler::if_ground_state(VecReal sub_energy) {
// }

VEC<VecInt> Occler::list_all_posible_nppsos(const VecInt& nppso_i, const VecInt& or_deg) const {
    if(nppso_i.size() != or_deg.size()) ERR("The or_deg's number should same as the nppso.")
	VEC<Int> idx(MAX(or_deg),0); Int cter(0);
	// for_Int(i, 0, or_deg.size()*2) if(cter < or_deg[i/2]) idx[cter++] = i; // for spin inversion symmetry
	for_Int(i, 0, or_deg.size()) if(cter < or_deg[i]) idx[cter++] = i;           // list all the possible
    VEC<VEC<Int>> nppsos_ii(MAX(or_deg)), nppsos_idx;
    for_Int(i, 0, idx.size()){
        Int idx_m(nppso_i[idx[i]]-1), idx_p(nppso_i[idx[i]]+1);
        nppsos_ii[i].push_back(nppso_i[idx[i]]);
        if(idx_m > 0) nppsos_ii[i].push_back(idx_m);
        if(idx_p > 0) nppsos_ii[i].push_back(idx_p);
        if(MIN(VecInt(nppsos_ii[i])) < 0) ERR("nppsos_ii's element should not has the minus.")
    }
    // if(mm) WRN(NAV2(nppsos_ii[0].size(),nppsos_ii[1].size()));
    nppsos_idx = cart_product(nppsos_ii);
    VEC<VecInt> nppsos(nppsos_idx.size(), tp.npartical);
    // for_Int(i, 0, nppsos_idx.size()) for_Int(j, 0, or_deg.size()*2) nppsos[i][j] = nppsos_idx[i][or_deg[j/2] - 1]; // for spin inversion symmetry
    for_Int(i, 0, nppsos_idx.size()) for_Int(j, 0, or_deg.size()) nppsos[i][j] = nppsos_idx[i][or_deg[j] - 1];        // list all the possible
    return nppsos;
}