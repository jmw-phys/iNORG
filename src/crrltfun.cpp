/*
code developed and maintained by (jmw@ruc.edu.cn, RUC, China) date 2022 - 2024
*/

#include "crrltfun.h"
using namespace std;

CrrltFun::CrrltFun(const MyMpi& mm_i, const Prmtr& prmtr_i, const NocSpace& old_nosp_i, const NocSpace& main_nosp, const VecReal& vgs_i, const Int position)
    :Operator(mm_i, prmtr_i, main_nosp),
    old_nosp(old_nosp_i), new_nosp(main_nosp), crtorann(main_nosp.nspa - old_nosp.nspa),
    ex_state(project_uplwer_parical_space(vgs_i, crtorann, position))
{}

CrrltFun::CrrltFun(const MyMpi& mm_i, const Prmtr& prmtr_i, const NocSpace& old_nosp_i, const NocSpace& main_nosp, const Operator& oper, const VecReal& vgs_i, const Int position)
    :Operator(oper),
    old_nosp(old_nosp_i), new_nosp(main_nosp), crtorann(main_nosp.nspa - old_nosp.nspa),
    ex_state(project_uplwer_parical_space(vgs_i, crtorann, position))
{}

CrrltFun::CrrltFun(const MyMpi& mm_i, const Prmtr& prmtr_i, const NocSpace& old_nosp_i, const NocSpace& main_nosp, const Operator& oper, const VecReal& vgs_i, const Int pos_in_set, const Int pos_in_div)
    :Operator(oper),
    old_nosp(old_nosp_i), new_nosp(main_nosp), crtorann(main_nosp.nspa - old_nosp.nspa),
    ex_state(project_uplwer_parical_space(vgs_i, crtorann, pos_in_set, pos_in_div))
{}

CrrltFun::CrrltFun(const MyMpi& mm_i, const Prmtr& prmtr_i, const Operator& oper, const VecReal ex_state_i, const Int crtorann_i)
    :Operator(oper),
    old_nosp(NocSpace(mm_i, prmtr_i)), new_nosp(NocSpace(mm_i, prmtr_i)), crtorann(crtorann_i),
    ex_state(ex_state_i)
{}


// ?CrrltFun::CrrltFun(const NocSpace& old_nosp_i, const Operator& main_op, const VecReal& vgs_i, const Int position = 0)
// ?    :Operator(main_op),
// ?    old_nosp(old_nosp_i), new_nosp(main_op.scsp), crtorann(main_op.scsp.nspa - old_nosp.nspa),
// ?    ex_state(project_uplwer_parical_space(vgs_i, crtorann, position))
// ?{}

// (Deactivate) Only for test
ImGreen CrrltFun::find_density_density_correlation_function(const Real &ge0) {

    Real upper_fraction(mm.Allreduce(DOT(ex_state, ex_state)));
    Vectrdgnl trdgnl(find_trdgnl(ex_state, -1));
    // if (mm) WRN("Finish the find_trdgnl");
    VecReal zerod(get<0>(trdgnl).size(), 0.);
    VecReal zeros(get<1>(trdgnl).size(), 0.);
    VecCmplx a(cmplx(get<0>(trdgnl), zerod));
    VecCmplx b(cmplx(get<1>(trdgnl), zeros));
    Int n(a.size());
    ImGreen green(1, p);
    for_Int(w, 0, green.nomgs) {
        //Cmplx z(cmplx(ge0, green.omg(w)));
        Cmplx z = green.z(w) + ge0;
        Cmplx c(z + a[n - 2] - b[n - 2] * b[n - 2] / (z + a[n - 1]));
        for (Int i = n; i >= 3; i--) c = z + a[i - 3] - b[i - 3] * b[i - 3] / c;
        Cmplx under_fraction(c);
        Cmplx gaz(upper_fraction / under_fraction);
        green[w][0][0] = gaz;
    }
    return green;
}


void CrrltFun::find_gf(const Real& ge0, Green &g0, bool is_greater) 
{
    // Control sign: 1 for greater, -1 for lesser
    const Int sign = is_greater ? 1 : -1;
    
    Real upper_fraction(mm.Allreduce(DOT(ex_state, ex_state)));
    // if (mm) WRN(NAV(upper_fraction));
    //VECtrdgnl trdignl(find_trdgnl_first(ex_state));
    VEC<Real> ltd;	        // diagonal elements 
    VEC<Real> lt_sd;	    // sub-diagonal elements
    const SparseMatReal sep_h = find_hmlt_V2(table);
    
    // Initialize vectors
    VecReal v0(ex_state);
    v0 *= INV(SQRT(mm.Allreduce(DOT(v0, v0))));
    VEC<VecReal> v0_saved;
    v0_saved.push_back(v0);
    VecReal v1(sep_h * v0);
    Real a_i(0.), b_i(0.);
    a_i = mm.Allreduce(DOT(v1, v0));
    ltd.push_back(a_i);

    // Initial iterations based on Green function type
    const Int init_iters = g0.type_info() == STR("ImGreen") ? 300 : 3000;
    for_Int(i, 0, init_iters) {
        find_trdgnl_one_step(v0_saved, v0, v1, a_i, b_i, sep_h);
        ltd.push_back(a_i); 
        lt_sd.push_back(b_i); 
        v0_saved.push_back(v0);
        // if(v0_saved.size() > 10) v0_saved.erase(v0_saved.begin());
    }

    // Calculate Green function
    Cmplx zero(0.,0.);
    VecCmplx green_pre(g0.nomgs, zero);
    VecCmplx green_error(g0.nomgs, zero);
    
    for_Int(iter_time, 0, 500) {
        find_trdgnl_one_step(v0_saved, v0, v1, a_i, b_i, sep_h);
        ltd.push_back(a_i); 
        lt_sd.push_back(b_i); 
        v0_saved.push_back(v0);
        // if(v0_saved.size() > 10) v0_saved.erase(v0_saved.begin());
        VecReal a(ltd), b(lt_sd);
        Int n(a.size());

        for_Int(w, 0, g0.nomgs) {
            // Use sign parameter to control the signs
            Cmplx z = g0.z(w) + sign * ge0;
            Cmplx c(z - sign * a[n - 2] - b[n - 2] * b[n - 2] / (z - sign * a[n - 1]));
            for (Int i = n; i >= 3; i--) c = z - sign * a[i - 3] - b[i - 3] * b[i - 3] / c;
            Cmplx gaz = upper_fraction / c;
            green_error[w] = gaz - green_pre[w];
            green_pre[w] = gaz;
        }

        // Convergence check
        if(iter_time % 50 == 0 && iter_time != 0){
            if(compare_error(ltd[ltd.size() - 1], ltd[ltd.size() - 2]) < 1E-6) break;
            VecCmplx check_err = g0.type_info() == STR("ImGreen") ? 
                green_error.truncate(0,int(p.fit_num_omg/2)) : green_error;
            if(mm) WRN(NAV(ABS(SUM(check_err))/check_err.size()));
            if(ABS(SUM(check_err))/check_err.size() < 1E-14) break;
        }
    }

    // Update results
    for_Int(w, 0, g0.nomgs) g0[w][0][0] += green_pre[w];
    if (mm) PIO("The size of a and b in " + STR(is_greater ? "greater" : "lesser") + 
                ":" + NAV2(ltd.size(), lt_sd.size()));
}

void CrrltFun::find_gf_greater(const Real& ge0, Green &g0) {
    find_gf(ge0, g0, true);
}

void CrrltFun::find_gf_lesser(const Real& ge0, Green &g0) {
    find_gf(ge0, g0, false);
}

/*
// ! only suit for the two orbital cases.
VecReal CrrltFun::find_hd_exstate(Int position, Int ex_idx) {
    VEC<VecBool> ex_VecBools;
    {
        ex_VecBools.push_back(VecBool({false, false}));// 0  , true, true
        ex_VecBools.push_back(VecBool({false,  true}));// 1  , true, true
        ex_VecBools.push_back(VecBool({ true, false}));// 2  , true, true
        ex_VecBools.push_back(VecBool({ true,  true}));// 3  , true, true
        //for special case
        ex_VecBools.push_back(VecBool());// 4
    }
#define ndiv new_nosp.ndivs
    if (p.norbs > 4) ERR(NAV(p.norbs));
    VecPartition vp(mm.np(), mm.id(), new_nosp.dim);
    
    VecReal hdex_state(ex_state);
    VecBool imp(4, false), imp_t(2, false);// the last one left for check the right.
    for_Int(h_i, 0, vp.len()) {
        StateStatistics a(vp.bgn() + h_i, new_nosp.wherein_NocSpace(vp.bgn() + h_i), new_nosp);
        for_Int(i, 0, p.norbs) imp[i] = a.cfg.cf[(i)*ndiv].isocc(0) ? true : false;
        bool check(true);
        if (position == 0) {
            imp_t = imp.truncate(2,4);
            if (ex_idx < 4) check = (imp_t == ex_VecBools[ex_idx]) && imp.truncate(0,2) == VecBool({ true,  true}) ? true : false;
            else if (ex_idx == 4) check = ((imp_t == ex_VecBools[1]) || (imp_t == ex_VecBools[2])) && imp.truncate(0,2) == VecBool({ true,  true}) ? true : false;                            // #4 = |1> + |2>
            else if (ex_idx == 5) check = ((imp_t == ex_VecBools[0]) || (imp_t == ex_VecBools[1]) || (imp_t == ex_VecBools[2])) && imp.truncate(0,2) == VecBool({ true,  true}) ? true : false; // #5 = |0> + |1> + |2> 
            else if (ex_idx == 6) check = imp[0] && (!imp[1]) && (!imp[2]) && imp[3];// kondo peak
        }
        if (position == 2) {
            imp_t = imp.truncate(0,2);
            if (ex_idx < 4) check = (imp_t == ex_VecBools[ex_idx]) && imp.truncate(2,4) == VecBool({ true,  true}) ? true : false;
            else if (ex_idx == 4) check = ((imp_t == ex_VecBools[1]) || (imp_t == ex_VecBools[2])) && imp.truncate(2,4) == VecBool({ true,  true}) ? true : false;                            // #4 = |1> + |2>
            else if (ex_idx == 5) check = ((imp_t == ex_VecBools[0]) || (imp_t == ex_VecBools[1]) || (imp_t == ex_VecBools[2])) && imp.truncate(2,4) == VecBool({ true,  true}) ? true : false; // #5 = |0> + |1> + |2> 
        }
        hdex_state[h_i] = check ? hdex_state[h_i] : 0.;
    }
#undef ndiv

    return hdex_state;
}
*/


// ! only suit for the three orbital cases.
VecReal CrrltFun::find_hd_exstate(Int position, Int ex_idx) {
    VEC<VecBool> ex_VecBools;
    {
        //       |00> |01> |10> |11>
        //  |00>  #0   #1   #2   #3
        //  |01>  #4   #5   #6   #7
        //  |10>  #8   #9  #10  #11
        //  |11> #12  #13  #14  #15
        ex_VecBools.push_back(VecBool({false, false, false, false, true, true}));// 0
        ex_VecBools.push_back(VecBool({false, false, false,  true, true, true}));// 1
        ex_VecBools.push_back(VecBool({false, false,  true, false, true, true}));// 2
        ex_VecBools.push_back(VecBool({false, false,  true,  true, true, true}));// 3
        ex_VecBools.push_back(VecBool({false,  true, false, false, true, true}));// 4
        ex_VecBools.push_back(VecBool({false,  true, false,  true, true, true}));// 5
        ex_VecBools.push_back(VecBool({false,  true,  true, false, true, true}));// 6
        ex_VecBools.push_back(VecBool({false,  true,  true,  true, true, true}));// 7
        ex_VecBools.push_back(VecBool({ true, false, false, false, true, true}));// 8
        ex_VecBools.push_back(VecBool({ true, false, false,  true, true, true}));// 9
        ex_VecBools.push_back(VecBool({ true, false,  true, false, true, true}));// 10
        ex_VecBools.push_back(VecBool({ true, false,  true,  true, true, true}));// 11
        ex_VecBools.push_back(VecBool({ true,  true, false, false, true, true}));// 12
        ex_VecBools.push_back(VecBool({ true,  true, false,  true, true, true}));// 13
        ex_VecBools.push_back(VecBool({ true,  true,  true, false, true, true}));// 14
        ex_VecBools.push_back(VecBool({ true,  true,  true,  true, true, true}));// 15
        //for special case
        ex_VecBools.push_back(VecBool()); // 16
    }
#define ndiv new_nosp.ndivs
    if (p.norbs != 6) ERR(NAV(p.norbs));
    VecPartition vp(mm.np(), mm.id(), new_nosp.dim);
    
    VecReal hdex_state(ex_state);
    VecBool imp(6, false);// the last one left for check the right.
    for_Int(h_i, 0, vp.len()) {
        StateStatistics a(vp.bgn() + h_i, new_nosp.wherein_NocSpace(vp.bgn() + h_i), new_nosp);
        for_Int(i, 0, p.norbs) imp[i] = a.cfg.cf[(i)*ndiv].isocc(0) ? true : false;
        bool check(true);
        // if (position == 0) check = !(imp[1]) && imp[0] && imp[3] && !(imp[2]);
        // if      (position == 0) check = !(imp[2]) && !(imp[3]) && imp[0] && imp[1];
        if (position == 0)      check = (imp[0]) && (!imp[1]) && (!imp[4]) && (imp[5]); //kondo peak
        else if (position == 2) check = (!imp[0]) && (!imp[1]) && (imp[2]) && (imp[3]);//hd peak
        else if (ex_idx < 16 && position == 4) check = imp == ex_VecBools[ex_idx] ? true : false;
        else if (ex_idx == 16 && position == 4) check = ((!imp[2]) && (!imp[3]) && (imp[4]) && (imp[5])) || ((!imp[0]) && (!imp[1]) && (imp[4]) && (imp[5])); // #16 h*d+*hd
        else if (ex_idx == 17 && position == 4) check = ((imp == ex_VecBools[1]) || (imp == ex_VecBools[2]) || (imp == ex_VecBools[4]) || (imp == ex_VecBools[8])) ? true : false; // #17
        else if (ex_idx == 18 && position == 4) check = ((imp == ex_VecBools[5]) || (imp == ex_VecBools[6]) || (imp == ex_VecBools[9]) || (imp == ex_VecBools[10])) ? true : false; // #18
        else if (ex_idx == 19 && position == 4) check = ((imp == ex_VecBools[7]) || (imp == ex_VecBools[11]) || (imp == ex_VecBools[13]) || (imp == ex_VecBools[14])) ? true : false; // #19
        else if (ex_idx == 20 && position == 4) check = ((imp == ex_VecBools[3]) || (imp == ex_VecBools[12])) ? true : false; // hdd+dhd
        else if (ex_idx == 21 && position == 4) check = ((imp == ex_VecBools[1]) || (imp == ex_VecBools[2]) || (imp == ex_VecBools[4]) || (imp == ex_VecBools[8]) || (imp == ex_VecBools[5]) || (imp == ex_VecBools[6]) || (imp == ex_VecBools[9]) || (imp == ex_VecBools[10])) ? true : false; // #17 || #18
        // else if (position == 4) check = !(imp[0]) && !(imp[1]) && !(imp[2]) && !(imp[3]) && imp[4] && imp[5];
        // else if (position == 4) check = !(imp[0]) && !(imp[1]) &&  (imp[2]) &&  (imp[3]) && imp[4] && imp[5];
        // else if (position == 4) check = (!(imp[0]) && !(imp[1]) && imp[4] && imp[5]) || (!(imp[2]) && !(imp[3]) && imp[4] && imp[5]);
        // else if (position == 4) check = (!(imp[0]) && !(imp[1]) && imp[4] && imp[5] && (imp[2]) && imp[3]) || (!(imp[0]) && !(imp[1]) && imp[4] && imp[5] && !(imp[3]) && !(imp[2]));
        else {ERR("The position is wrong!" + NAV(position));}
        // if(mm && check) WRN(NAV(imp));
        hdex_state[h_i] = check ? hdex_state[h_i] : 0.;
    }
    // if (position == 2) 	for_Int(h_i, vp.bgn(), vp.end()) {
    //     StateStatistics a(h_i, new_nosp.wherein_NocSpace(h_i), new_nosp);
    //     for_Int(i, 0, p.norbs) imp[i] = a.cfg.cf[(i)*ndiv].isocc(0) ? true : false;
    //     bool check = !(imp[0]) && !(imp[1]) && imp[2] && imp[3];
    //     if(mm && check) WRN(NAV(imp));
    //     hdex_state[h_i] = check ? hdex_state[h_i] : 0.;
    // }
#undef ndiv

    return hdex_state;
}

void CrrltFun::find_hd_greater(const Real& ge0, Green &g0, Int position, Int ex_idx) 
{
    VecReal hd_exstate(find_hd_exstate(position, ex_idx));
    Real upper_fraction(mm.Allreduce(DOT(hd_exstate, hd_exstate)));
    // if (mm) WRN(NAV(upper_fraction));
    VEC<Real> ltd;	        // diagonal elements 
    VEC<Real> lt_sd;	    // sub-diagonal elements
    const SparseMatReal sep_h = find_hmlt_V2(table);
    VecReal v0(hd_exstate);
    v0 *= INV(SQRT(mm.Allreduce(DOT(v0, v0))));
    VEC<VecReal> v0_saved;
    v0_saved.push_back(v0);
    VecReal v1(sep_h * v0);
    Real a_i(0.), b_i(0.);
    a_i = mm.Allreduce(DOT(v1, v0));
    ltd.push_back(a_i);
    if(g0.type_info() == STR("ImGreen")){
        for_Int(i, 0, 200) {
        find_trdgnl_one_step(v0_saved, v0, v1, a_i, b_i, sep_h);
        ltd.push_back(a_i); lt_sd.push_back(b_i); v0_saved.push_back(v0);
        // if(v0_saved.size() > 10) v0_saved.erase(v0_saved.begin());
        }
    }
    if(g0.type_info() == STR("ReGreen")){
        for_Int(i, 0, 3000) {
        find_trdgnl_one_step(v0_saved, v0, v1, a_i, b_i, sep_h);
        ltd.push_back(a_i); lt_sd.push_back(b_i); v0_saved.push_back(v0);
        // if(v0_saved.size() > 10) v0_saved.erase(v0_saved.begin());
        }
    }
    // ImGreen last_green(1, p);
    Cmplx zero(0.,0.);
    VecCmplx green_pre(g0.nomgs, zero);
    VecCmplx green_error(g0.nomgs, zero);
    Int iter_time(0);
    while (true) {
        find_trdgnl_one_step(v0_saved, v0, v1, a_i, b_i, sep_h);
        ltd.push_back(a_i); lt_sd.push_back(b_i);
        VecReal a(ltd), b(lt_sd);
        Int n(a.size());
        //ImGreen green_i(1, p);
        for_Int(w, 0, g0.nomgs) {
            Cmplx z = g0.z(w) + ge0;
            Cmplx c(z - a[n - 2] - b[n - 2] * b[n - 2] / (z - a[n - 1]));
            for (Int i = n; i >= 3; i--) c = z - a[i - 3] - b[i - 3] * b[i - 3] / c;
            //Cmplx under_fraction(c);
            Cmplx gaz = upper_fraction / c;
            green_error[w] = gaz - green_pre[w];
            green_pre[w] = gaz;
        }
        if(iter_time % 50 == 0 && iter_time != 0){
            if(compare_error(ltd[ltd.size() - 1], ltd[ltd.size() - 2]) < 1E-6) break;
            VecCmplx check_err = g0.type_info() == STR("ImGreen") ? green_error.truncate(0,int(p.fit_num_omg/2)): green_error;
            if(mm) WRN(NAV(ABS(SUM(check_err))/check_err.size()));
            if(ABS(SUM(check_err))/check_err.size() < 1E-14) break;
        }
        iter_time++;
        // if (mm && ltd.size() > 200 && g0.type_info() == STR("ImGreen")) PIO("The size of a and b in greaer:" + NAV3(ltd.size(), lt_sd.size(), SUM(green_error)));
    }
    for_Int(w, 0, g0.nomgs) g0[w][0][0] += green_pre[w];
    if (mm) PIO("The size of a and b in greaer:" + NAV2(ltd.size(), lt_sd.size()));
    // return last_green;
}

// void CrrltFun::find_gf_greater(const Real& ge0, Green &g0, Int kind) 
// {
//     Real upper_fraction(DOT(ex_state, ex_state));
//     //VECtrdgnl trdignl(find_trdgnl_first(ex_state));
//     VEC<Real> ltd;	        // diagonal elements 
//     VEC<Real> lt_sd;	    // sub-diagonal elements
//     const SparseMatReal sep_h = find_hmlt_V2(table);
//     Cmplx zero(0.,0.);
//     VecCmplx green_pre(g0.nomgs, zero);
//     VecCmplx green_error(g0.nomgs, zero);
//     if(kind == 0){
//         VecReal v0(ex_state);
//         v0.normalize();
//         VecReal v0_saved(v0);
//         VecReal v1(sep_h * v0);
//         Real a_i(0.), b_i(0.);
//         a_i = DOT(v1, v0);
//         ltd.push_back(a_i);
//         for_Int(i, 0, 3000) {
//             find_trdgnl_one_step(v0_saved, v0, v1, a_i, b_i, sep_h);
//             ltd.push_back(a_i); lt_sd.push_back(b_i);
//         }
//         // ImGreen last_green(1, p);
//         while (true) {
//             find_trdgnl_one_step(v0_saved, v0, v1, a_i, b_i, sep_h);
//             ltd.push_back(a_i); lt_sd.push_back(b_i);
//             VecReal a(ltd), b(lt_sd);
//             Int n(a.size());
//             //ImGreen green_i(1, p);
//             for_Int(w, 0, g0.nomgs) {
//                 Cmplx z = g0.z(w, kind) + ge0;
//                 Cmplx c(z - a[n - 2] - b[n - 2] * b[n - 2] / (z - a[n - 1]));
//                 for (Int i = n; i >= 3; i--) c = z - a[i - 3] - b[i - 3] * b[i - 3] / c;
//                 //Cmplx under_fraction(c);
//                 Cmplx gaz = upper_fraction / c;
//                 green_error[w] = gaz - green_pre[w];
//                 green_pre[w] = gaz;
//             }
//             if (ABS(SUM(green_error)) < 1.E-10 * g0.nomgs) break;
//             a_saved.reset(a); b_saved.reset(b);
//         }
//     }
//     if(kind > 0){
//         VecReal a(a_saved), b(b_saved);
//         Int n(a.size());
//         for_Int(w, 0, g0.nomgs) {
//             Cmplx z = g0.z(w, kind) + ge0;
//             Cmplx c(z - a[n - 2] - b[n - 2] * b[n - 2] / (z - a[n - 1]));
//             for (Int i = n; i >= 3; i--) c = z - a[i - 3] - b[i - 3] * b[i - 3] / c;
//             //Cmplx under_fraction(c);
//             Cmplx gaz = upper_fraction / c;
//             green_pre[w] = gaz;
//         }
//     }
//     for_Int(w, 0, g0.nomgs) g0[w][0][0] += green_pre[w];
//     if (mm && kind == 0) PIO("The size of a and b in greaer:" + NAV2(a_saved.size(), b_saved.size()));
//     // return last_green;
// }

// void CrrltFun::find_gf_lesser(const Real& ge0, Green &g0, Int kind) 
// {
//     Real upper_fraction(DOT(ex_state, ex_state));
//     //VECtrdgnl trdignl(find_trdgnl_first(ex_state));
//     VEC<Real> ltd;	        // diagonal elements 
//     VEC<Real> lt_sd;	    // sub-diagonal elements
//     const SparseMatReal sep_h = find_hmlt_V2(table);
//     Cmplx zero(0.,0.);
//     VecCmplx green_pre(g0.nomgs, zero);
//     VecCmplx green_error(g0.nomgs, zero);
//     if(kind == 0){
//         VecReal v0(ex_state);
//         v0.normalize();
//         VecReal v0_saved(v0);
//         VecReal v1(sep_h * v0);
//         Real a_i(0.), b_i(0.);
//         a_i = DOT(v1, v0);
//         ltd.push_back(a_i);
//         for_Int(i, 0, 3000) {
//             find_trdgnl_one_step(v0_saved, v0, v1, a_i, b_i, sep_h);
//             ltd.push_back(a_i); lt_sd.push_back(b_i);
//         }
//         // ImGreen last_green(1, p);
//         while (true) {
//             find_trdgnl_one_step(v0_saved, v0, v1, a_i, b_i, sep_h);
//             ltd.push_back(a_i); lt_sd.push_back(b_i);
//             VecReal a(ltd), b(lt_sd);
//             Int n(a.size());
//             //ImGreen green_i(1, p);
//             for_Int(w, 0, g0.nomgs) {
//                 Cmplx z = g0.z(w, kind) - ge0;
//                 Cmplx c(z + a[n - 2] - b[n - 2] * b[n - 2] / (z + a[n - 1]));
//                 for (Int i = n; i >= 3; i--) c = z + a[i - 3] - b[i - 3] * b[i - 3] / c;
//                 //Cmplx under_fraction(c);
//                 Cmplx gaz = upper_fraction / c;
//                 green_error[w] = gaz - green_pre[w];
//                 green_pre[w] = gaz;
//             }
//             if (ABS(SUM(green_error)) < 1.E-10 * g0.nomgs) break;
//             a_saved.reset(a); b_saved.reset(b);
//         }
//     }
//     if(kind > 0){
//         VecReal a(a_saved), b(b_saved);
//         Int n(a.size());
//         for_Int(w, 0, g0.nomgs) {
//             Cmplx z = g0.z(w, kind) - ge0;
//             Cmplx c(z + a[n - 2] - b[n - 2] * b[n - 2] / (z + a[n - 1]));
//             for (Int i = n; i >= 3; i--) c = z + a[i - 3] - b[i - 3] * b[i - 3] / c;
//             //Cmplx under_fraction(c);
//             Cmplx gaz = upper_fraction / c;
//             green_pre[w] = gaz;
//         }
//     }
//     for_Int(w, 0, g0.nomgs) g0[w][0][0] += green_pre[w];
//     if (mm && kind == 0) PIO("The size of a and b in lesser:" + NAV2(a_saved.size(), b_saved.size()));
//     // return last_green;
// }

//---------------------------------------------Private function---------------------------------------------

// (Deactivate) 
CrrltFun::Vectrdgnl CrrltFun::find_trdgnl(const VecReal& initial_vector, const Int crtann) {
    // for tridiagonal matrix
    VEC<Real> ltd;	    // diagonal elements 
    VEC<Real> lt_sd;	// sub-diagonal elements
    SparseMatReal sep_h = find_hmlt_V2(table);
    VecReal v0(initial_vector);
    v0.normalize();
    VecReal v0_saved(v0);
    //DBG(NAV(v0_saved.truncate(0, 1000)));
    VecReal v1(sep_h * v0);
    Idx k = 0;
    Real a(0.), b(0.);
    while (true)
    {
        a = mm.Allreduce(DOT(v1, v0));
        ltd.push_back(a);
        if (k >= 30) if (k > v0_saved.size() - 20) break;
        k++;
        v1 -= a * v0;
        v1 -= mm.Allreduce(DOT(v1, v0_saved)) * v0_saved;
        b = v1.norm();
        lt_sd.push_back(b);
        v1 *= (1 / b);
        SWAP(v0, v1);
        v1 *= -b;
        v1 += sep_h * v0;
    }
    VecReal va_i(ltd), vb_i(lt_sd);
#ifdef _ASSERTION_
    if (va_i.size() - 1 != vb_i.size())ERR("Wrong with count orbit idx: " + NAV2(va_i.size(), vb_i.size()));
#endif
    return std::make_tuple(va_i, vb_i);
}

void CrrltFun::find_trdgnl_one_step(const VEC<VecReal>& initial_vector, VecReal& v0, VecReal& v1, Real&a, Real& b, const SparseMatReal& sep_h){
    v1 -= a * v0;
    b = SQRT(mm.Allreduce(DOT(v1, v1)));
    // v1 -= mm.Allreduce(DOT(v1, initial_vector)) * initial_vector;
    for_Int(i, 0, initial_vector.size() - 1) {
        v1 -= mm.Allreduce(DOT(v1, initial_vector[i])) * initial_vector[i];
        // WRN("For the iter" + NAV(i));
    }
    v1 *= INV(SQRT(mm.Allreduce(DOT(v1, v1))));
    SWAP(v0, v1);
    v1 *= -b;
    v1 += sep_h * v0;
    a = mm.Allreduce(DOT(v1, v0));
}

//--------------------------------------------- ex_state function part ---------------------------------------------
void CrrltFun::add_ex_state_part_in_rotation(const VecReal &initial_vector, VecReal& ex_state_part, const Int& set_n, const Int orb_before_rot, const Int& h_i) const {
    const Int subnosp(old_nosp.wherein_NocSpace(h_i));

    StateStatistics a(h_i, subnosp, old_nosp);

    VecReal rotation_coefficients = p.rotationU[set_n].tr()[orb_before_rot];
    for_Int(pos, 0, rotation_coefficients.size()) {
        MatInt new_nospdiv = old_nosp.div[subnosp];
        Int div_idx_in_one_set(0);
        for_Int(i, 1, p.ndiv + 1) if (SUM_0toX(old_nosp.sit_mat[set_n], i) > pos) {div_idx_in_one_set = i - 1; break; }
        if (crtorann == -1) --new_nospdiv[set_n][div_idx_in_one_set];
        if (crtorann == +1) ++new_nospdiv[set_n][div_idx_in_one_set];
        if (new_nosp.ifin_NocSpace(new_nospdiv, new_nosp.nppso)) {
            VecOnb exd_cf = a.cfg.cf;
            Int orbit_pos_in_div = pos - SUM_0toX(old_nosp.sit_mat[set_n], div_idx_in_one_set);
            if (if_in_this_orbital(exd_cf, crtorann, set_n * new_nosp.ndivs + div_idx_in_one_set, orbit_pos_in_div)) {
                if (crtorann == -1) exd_cf[set_n * new_nosp.ndivs + div_idx_in_one_set] = exd_cf[set_n * new_nosp.ndivs + div_idx_in_one_set].ann(orbit_pos_in_div);
                if (crtorann == +1) exd_cf[set_n * new_nosp.ndivs + div_idx_in_one_set] = exd_cf[set_n * new_nosp.ndivs + div_idx_in_one_set].crt(orbit_pos_in_div);
                const ComDivs b(exd_cf, new_nospdiv, new_nosp.sit_mat);
                Int begin_idx(-1);
                begin_idx = new_nosp.divs_to_idx.at(new_nospdiv.vec().string());
                Int sign  = a.cfg.sgn(-1, pos + SUM_0toX(p.nO2sets, set_n));
                if (begin_idx == -1) ERR("wrong with ex_state" + NAV2(a.cfg.ne, new_nospdiv));
                ex_state_part[begin_idx + b.idx] += sign * rotation_coefficients[pos] * initial_vector[h_i];
            }
        }
    }
}

VecReal CrrltFun::project_uplwer_parical_space(const VecReal& initial_vector, const Int crtann, const Int norg_set, const Int orbit_pos_in_div) const
// By using the C^+ or C on a spinless orbital(for the spinless reason ONE normal orbital has two spinless orbits).
// the norg_set only suppose for the even number.
{
    // clock_t t_find_Newstate; TIME_BGN("find_Newstate" + NAV(mm.id()), t_find_Newstate);
    VecReal ex_state_part(new_nosp.dim, 0.);
    VecPartition row_H(mm.np(), mm.id(), old_nosp.dim);
    if (p.if_norg_imp) for_Int(h_i, row_H.bgn(), row_H.end()) {
        add_ex_state_part_in_rotation(initial_vector, ex_state_part, norg_set, orbit_pos_in_div, h_i);
    }
    else for_Int(h_i, row_H.bgn(), row_H.end()) {
        const Int subnosp(old_nosp.wherein_NocSpace(h_i));
        MatInt new_nospdiv = old_nosp.div[subnosp];
        if (crtann == -1) --new_nospdiv[norg_set][0];
        if (crtann == +1) ++new_nospdiv[norg_set][0];
        // if (new_nosp.ifin_NocSpace(new_nospdiv)) {
        if (new_nosp.ifin_NocSpace(new_nospdiv, new_nosp.nppso)) {
            const ComDivs group(h_i - old_nosp.idx_div[subnosp], (old_nosp.div[subnosp]), (old_nosp.sit_mat), true);
            VecOnb exd_cf = group.cf;
            if (if_in_this_orbital(exd_cf, crtann, norg_set * new_nosp.ndivs, orbit_pos_in_div)) {
                if (crtann == -1) exd_cf[norg_set * new_nosp.ndivs] = exd_cf[norg_set * new_nosp.ndivs].ann(orbit_pos_in_div);
                if (crtann == +1) exd_cf[norg_set * new_nosp.ndivs] = exd_cf[norg_set * new_nosp.ndivs].crt(orbit_pos_in_div);
                const ComDivs b(exd_cf, new_nospdiv, old_nosp.sit_mat);
                Int begin_idx(-1);
                begin_idx = new_nosp.divs_to_idx.at(new_nospdiv.vec().string());
                if (begin_idx == -1) ERR("wrong with ex_state" + NAV2(group.ne, new_nospdiv));
                ex_state_part[begin_idx + b.idx] = exd_cf[norg_set * new_nosp.ndivs].sgn(orbit_pos_in_div) * initial_vector[h_i]; //! Here the "sgn(orbit_pos_in_div)" is useful only on the ↑ spin orbital.
            }
        }
    }
    VecReal ex_state(mm.Allreduce(ex_state_part));
    // TIME_END("find_Newstate" + NAV(mm.id()), t_find_Newstate);
    // if(mm) WRN(NAV4(initial_vector.truncate(0,1000), ex_state.truncate(0,1000), SUM(initial_vector), SUM(ex_state)));
    VecPartition vp(mm.np(), mm.id(), new_nosp.dim);
    VecReal ex_state_p(vp.len());
    for_Int(h_i, 0, vp.len()) ex_state_p[h_i] = ex_state[vp.bgn() + h_i];
    return ex_state_p;
}

bool CrrltFun::if_in_this_orbital(const VecOnb &exd_cf, const Int crtann, const Int div_in_one_row_pos, const Int orbit_pos_in_div) const {
    if (crtann == -1) return exd_cf[div_in_one_row_pos].isocc(orbit_pos_in_div);
    if (crtann == +1) return exd_cf[div_in_one_row_pos].isuno(orbit_pos_in_div);
    ERR("Some thing wrong with this function, if_in_this_orbital ()!")
}
