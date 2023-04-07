/*
coded by Jia-Ming Wang (jmw@ruc.edu.cn, RUC, China) date 2021 - 2023
*/

#include "norg.h"
#define condition iter_norg_cnt < p.iter_max_norg  && !converged()

NORG::NORG(const MyMpi& mm_i, const Prmtr& prmtr_i) :
	mm(mm_i), p(prmtr_i), impgreen(prmtr_i.norbs, prmtr_i), uormat(uormat_initialize()),
	occupation_err(1.), energy_err(1.), correctionerr(1.), occnum_pre(prmtr_i.n_rot_orb, 0.),
	iter_norg_cnt(0), groune_pre(1.e99), groune_lst(0.), occnum_lst(prmtr_i.n_rot_orb, 0.),
	scsp(mm_i, prmtr_i, prmtr_i.npartical), oneedm(mm, prmtr_i, scsp),norg_stable_count(0)
{
	show_the_nozero_number_of_tabel();
}

// NORG::NORG(const MyMpi& mm_i, const Prmtr& prmtr_i, VecInt nparticals) :
// 	mm(mm_i), p(prmtr_i), impgreen(prmtr_i.norbs, prmtr_i), uormat(uormat_initialize()),
// 	occupation_err(1.), energy_err(1.), correctionerr(1.), occnum_pre(prmtr_i.n_rot_orb, 0.),
// 	iter_norg_cnt(0), groune_pre(1.e99), groune_lst(0.), occnum_lst(prmtr_i.n_rot_orb, 0.),
// 	scsp(mm_i, prmtr_i, nparticals), oneedm(mm, prmtr_i, scsp),norg_stable_count(0)
// {
// 	show_the_nozero_number_of_tabel();
// }

NORG::NORG(const MyMpi& mm_i, const Prmtr& prmtr_i, const Tab& table) :	
	mm(mm_i), p(prmtr_i), impgreen(prmtr_i.norbs, prmtr_i), uormat(uormat_initialize()),
	occupation_err(1.), energy_err(1.), correctionerr(1.), occnum_pre(prmtr_i.n_rot_orb, 0.),
	iter_norg_cnt(0), groune_pre(1.e99), groune_lst(0.), occnum_lst(prmtr_i.n_rot_orb, 0.),
	scsp(mm_i, prmtr_i, prmtr_i.npartical, table), oneedm(mm, prmtr_i, scsp, table),norg_stable_count(0)
{

}

NORG::NORG(const MyMpi& mm_i, const Prmtr& prmtr_i, Str tab_name) :
	mm(mm_i), p(prmtr_i), impgreen(prmtr_i.norbs, prmtr_i), uormat(uormat_initialize()),
	occupation_err(1.), energy_err(1.), correctionerr(1.), occnum_pre(prmtr_i.n_rot_orb, 0.),
	iter_norg_cnt(0), groune_pre(1.e99), groune_lst(0.), occnum_lst(prmtr_i.n_rot_orb, 0.),
	scsp(mm_i, prmtr_i, prmtr_i.npartical, tab_name), oneedm(mm, prmtr_i, scsp, tab_name),norg_stable_count(0)
{
	show_the_nozero_number_of_tabel();
}

void NORG::up_date_h0_to_solve(const Impdata& impH_i, const VecReal sub_energy) {
	if(mm) std::cout << std::endl;						// blank line
	impH = impH_i;
	// if (mm) PIO(NAV2(h0,scsp.dim));
	// scsp.coefficient = set_row_primeter_byimpH(uormat, impH_i);
	set_row_primeter_byimpH(uormat, impH_i, oneedm.oper_value);
	oneedm.update(); final_ground_state = oneedm.ground_state;
	// if(mm) PIO(NAV(oneedm.sum_off_diagonal()));
	groune_lst = oneedm.groundstate_energy;
	if (mm)	{
		scsp.print();
		// write_norg_info(iter_norg_cnt);
		// write_state_info(iter_norg_cnt);
	}

	if(sub_energy.size() != 0 && MIN(sub_energy) < groune_lst * (1 + 2e-3) )  {
		if(mm) std::cout <<iofmt("def")<< "The energy level is so hight, "<< NAV(groune_lst) << ".\nWhich reach to " \
		<<100*MIN(sub_energy)/groune_lst<<"%, so, the Ground state mustn't in this sub-space, and calc. will be stoped." <<"\n"<< std::endl;
		return ;}
	
	while (iter_norg_cnt < p.iter_max_norg && !converged()) {
		iter_norg_cnt++;
		if(mm) PIO("The iteration counting: " + NAV(iter_norg_cnt));
		VEC<MatReal> uormat_new(oneedm.find_unitary_orbital_rotation_matrix());
		// if(mm) WRN(NAV(see_MatReal(uormat_new)));
		for_Int(i, 0, uormat.size()) uormat[i] = uormat_new[i] * uormat[i];
		// if(mm) WRN(NAV(uormat[0]));
		// scsp.coefficient = set_row_primeter_byimpH(uormat, impH_i);
		set_row_primeter_byimpH(uormat, impH_i, oneedm.oper_value);	//if (mm) scsp.print();
		// if (mm)PIO("ground_state size" + NAV(oneedm.ground_state.size()));
		groune_pre = groune_lst;	occnum_pre = occnum_lst;
		oneedm.update(); final_ground_state = oneedm.ground_state;
		// if(mm) PIO(NAV(see_MatReal(oneedm.dm)));
		// if(mm) PIO(NAV(oneedm.sum_off_diagonal()));
		occnum_lst = VECVectoVec(oneedm.occupationnumber);
		groune_lst = oneedm.groundstate_energy;
		// if(mm) WRN(NAV(oneedm.dm[0]));
		if (mm) {
			// WRN(NAV(iter_norg_cnt));
			// for_Int(i, 0, uormat.size()) WRN(NAV(uormat_new[i]));
			// write_norg_info(iter_norg_cnt);
			// write_state_info(iter_norg_cnt);
		}
		occupation_err = SQRT(SUM(SQR(occnum_pre - occnum_lst)) / p.norbit);
		energy_err = 2 * (groune_pre - groune_lst) / (ABS(groune_pre) + ABS(groune_lst) + 1.);
		if (mm) {
			PIO(NAV2(energy_err, occupation_err));
			std::cout << "groundE_pre " << iofmt("sci") << groune_pre << std::setw(4) << "   groundE_lst " << groune_lst  << "  " << present() << std::endl;
		}
		if(mm) std::cout << std::endl;						// blank line
	}
	iter_norg_cnt = 0;		occupation_err = 1.;		energy_err = 1.;
	final_ground_state = oneedm.ground_state;	norg_stable_count = 0.; occnum = occnum_lst;
	/*--------------------------------print put--------------------------------*/
	MatReal occnum, occweight;
	occnum = occnum_lst.mat(p.norg_sets, p.n_rot_orb/p.norg_sets);
	occweight = occnum;
	for_Int(i, 0, p.norg_sets) for_Int(j, 0, occnum.ncols()) occweight[i][j] = MIN(occweight[i][j],1 - occnum[i][j]) < 1e-14 ? 0 : MIN(occweight[i][j],1 - occnum[i][j]);
	Str nppso = scsp.nppso_str();
	if (mm) WRN(NAV4(p.if_norg_imp, nppso, occnum, occweight));

	p.rotationU = uormat;
	// Finish the NORGiteration.
	// oneedm.save_the_Tab(oneedm.table, "mainspace");
}

// mode 0: for not rotation; mode 1: for rotation the orbitals; 
void NORG::up_date_h0_to_solve(const Impdata& impH_i, const Int mode) {
	if(mm) std::cout << std::endl;						// blank line
	// if (mm) WRN(NAV2(impH_i.first,scsp.dim));
	impH = impH_i;
	
	// {// test
	// 	if(mm) WRN(NAV2(scsp.coefficient.size(), present()));
	// 	VecReal a(set_row_primeter_byimpH(uormat, impH_i));
	// 	if(mm) WRN(NAV2( a.size(), present()));
	// }

	// scsp.coefficient = set_row_primeter_byimpH(uormat, impH_i);
	set_row_primeter_byimpH(uormat, impH_i, oneedm.oper_value);

	oneedm.update(); final_ground_state = oneedm.ground_state;
	// if(mm) PIO(NAV(oneedm.sum_off_diagonal()));
	groune_lst = oneedm.groundstate_energy;
	if (mm)	{
		scsp.print();
		// write_norg_info(iter_norg_cnt);
		// write_state_info(iter_norg_cnt);
	}
	
	while (iter_norg_cnt < p.iter_max_norg && !converged()) {
		iter_norg_cnt++;
		if(mm) PIO("The iteration counting: " + NAV(iter_norg_cnt));
		if(mode) {
			VEC<MatReal> uormat_new(oneedm.find_unitary_orbital_rotation_matrix());
			// if(mm) WRN(NAV(see_MatReal(uormat_new)));
			for_Int(i, 0, uormat.size()) uormat[i] = uormat_new[i] * uormat[i];}
		// if(mm) WRN(NAV(uormat[0]));
		// scsp.coefficient = set_row_primeter_byimpH(uormat, impH_i);	//if (mm) scsp.print();
		set_row_primeter_byimpH(uormat, impH_i, oneedm.oper_value);
		// if (mm)PIO("ground_state size" + NAV(oneedm.ground_state.size()));
		groune_pre = groune_lst;	occnum_pre = occnum_lst;
		oneedm.update(); final_ground_state = oneedm.ground_state;
		// if (mm) PIO(NAV(see_MatReal(oneedm.dm)));
		// if(mm) PIO(NAV(oneedm.sum_off_diagonal()));
		occnum_lst = VECVectoVec(oneedm.occupationnumber);				PIO_occweight(occnum_lst);
		groune_lst = oneedm.groundstate_energy;
		// if(mm) WRN(NAV(oneedm.dm[0]));
		if (mm) {
			// WRN(NAV(iter_norg_cnt));
			// for_Int(i, 0, uormat.size()) WRN(NAV(uormat_new[i]));
			// write_norg_info(iter_norg_cnt);
			// write_state_info(iter_norg_cnt);
		}
		occupation_err = SQRT(SUM(SQR(occnum_pre - occnum_lst)) / p.norbit);
		energy_err = 2 * (groune_pre - groune_lst) / (ABS(groune_pre) + ABS(groune_lst) + 1.);
		if (mm) {
			PIO(NAV2(energy_err, occupation_err));
			std::cout << "groundE_pre " << iofmt("sci") << groune_pre << std::setw(4) << "   groundE_lst " << groune_lst  << "  " << present() << std::endl;
		}
		if(mm) std::cout << std::endl;						// blank line
	}
	iter_norg_cnt = 0;		occupation_err = 1.;		energy_err = 1.;
	final_ground_state = oneedm.ground_state;	norg_stable_count = 0.; occnum = occnum_lst;
	//----------------------print put--------------------------------
	if(mm) std::cout << "---------- THE NORG CONVERGED ----------" << std::endl;;
	PIO_occweight(occnum);
	p.rotationU = uormat;
	write_impurtiy_occupation();
	// Finish the NORGiteration.
	// oneedm.save_the_Tab(oneedm.table, "mainspace");
}

void NORG::get_g_by_KCV(Green& imp_i)
{
	NocSpace scsp_1pone(mm, p, nppso(p.npartical, 1));
	NocSpace scsp_1mone(mm, p, nppso(p.npartical, -1));
	Operator n1pone(mm, p, scsp_1pone);
	Operator n1mone(mm, p, scsp_1mone);
	for_Int(i, 0, imp_i.g[0].nrows()) {
	// for_Int(i, 0, 3){
		// {Int i(0);
			// scsp_1pone.coefficient = set_row_primeter_byimpH(uormat, impH);
			// scsp_1mone.coefficient = set_row_primeter_byimpH(uormat, impH);
			set_row_primeter_byimpH(uormat, impH, n1pone.oper_value);
			set_row_primeter_byimpH(uormat, impH, n1mone.oper_value);

			Crrvec greaer(scsp, n1pone, 0, i, final_ground_state, groune_lst);
			Crrvec lesser(scsp, n1mone, 0, i, final_ground_state, groune_lst);
			
			{// test
				if(mm) WRN(NAV2(greaer.ex_state.norm(), lesser.ex_state.norm()));
			}

		for_Int(j, 0, imp_i.g[0].ncols()) {
		// for_Int(j, 0, 3){
			// {Int j(i);
			VecCmplx green_function, vec_z;
			if(imp_i.type_info() == STR("ImGreen")) vec_z = p.Im_z;
			else if(imp_i.type_info() == STR("ReGreen")) vec_z = p.Re_z;
			else ERR("The imp_i.type_info() was wrong");
			green_function  = greaer.krylov_space_for_green(0, j, vec_z);
			green_function += lesser.krylov_space_for_green(0, j, vec_z);
			for_Int(n, 0, imp_i.g.size()) imp_i[n][i][j] = green_function[n];
		}
	if(mm) std::cout << present() << std::endl;
	}
}

void NORG::get_g_by_KCV_spup(Green& imp_i)
{
	StdVecInt difference = {1, -1};
	for(const auto ii: difference)
	{
		NocSpace scsp_1(mm, p, nppso(p.npartical, ii));
		Operator opr_sub(mm, p, scsp_1);
		// scsp_1.coefficient = set_row_primeter_byimpH(uormat, impH);
		set_row_primeter_byimpH(uormat, impH, opr_sub.oper_value);
		Crrvec greaer(scsp, opr_sub, final_ground_state, groune_lst, imp_i);
	}
	if(mm) std::cout << present() << std::endl;
}

void NORG::get_gimp(Green& imp_i)
{
	for_Int(i, 0, p.nband) {
		StdVecInt difference = {(i+1), -(i+1)};
		for(const auto ii: difference)
		{
			NocSpace scsp_sub(mm, p, nppso(p.npartical, ii));
			Operator opr_sub(mm, p, scsp_sub);
			// scsp_sub.coefficient = set_row_primeter_byimpH(uormat, impH);
			set_row_primeter_byimpH(uormat, impH, opr_sub.oper_value);
			CrrltFun temp_green(mm, p, scsp, scsp_sub, opr_sub.table, final_ground_state, i * 2);
			if(imp_i.type_info() == STR("ImGreen")) {
				ImGreen green_function(1, p);
				if(ii > 0) temp_green.find_gf_greater(groune_lst, green_function);
				if(ii < 0) temp_green.find_gf_lesser(groune_lst, green_function);
				for_Int(n, 0, green_function.nomgs) imp_i[n][i][i] += green_function[n][0][0];
			}
			if(imp_i.type_info() == STR("ReGreen")) {
				ReGreen green_function(1, p);
				if(ii > 0) temp_green.find_gf_greater(groune_lst, green_function);
				if(ii < 0) temp_green.find_gf_lesser(groune_lst, green_function);
				for_Int(n, 0, green_function.nomgs) imp_i[n][i][i] += green_function[n][0][0];
			}
		}
		if (mm) PIO("finished the " + STR(i) + " find_g_norg   " + present());
	}
}

void NORG::get_gimp(Green& imp_i, VecInt or_deg)
{
	VecInt idx(MAX(or_deg),0); Int cter(0);
	for_Int(i, 0, or_deg.size()) if(cter < or_deg[i]) idx[cter++] = i; 
	for (int &i : idx) {
	// for_Int(i, 0, p.nband) {
		StdVecInt difference = {(i+1), -(i+1)};
		for(const auto ii: difference)
		{
			NocSpace scsp_sub(mm, p, nppso(p.npartical, ii));
			Operator opr_sub(mm, p, scsp_sub);
			// scsp_sub.coefficient = set_row_primeter_byimpH(uormat, impH);
			set_row_primeter_byimpH(uormat, impH, opr_sub.oper_value);
			CrrltFun temp_green(mm, p, scsp, scsp_sub, opr_sub.table, final_ground_state, i * 2);
			if(imp_i.type_info() == STR("ImGreen")) {
				ImGreen green_function(1, p);
				if(ii > 0) temp_green.find_gf_greater(groune_lst, green_function);
				if(ii < 0) temp_green.find_gf_lesser(groune_lst, green_function);
				for_Int(n, 0, green_function.nomgs) imp_i[n][i][i] += green_function[n][0][0];
			}
			if(imp_i.type_info() == STR("ReGreen")) {
				ReGreen green_function(1, p);
				if(ii > 0) temp_green.find_gf_greater(groune_lst, green_function);
				if(ii < 0) temp_green.find_gf_lesser(groune_lst, green_function);
				for_Int(n, 0, green_function.nomgs) imp_i[n][i][i] += green_function[n][0][0];
			}
		}
		if (mm) PIO("finished the " + STR(i) + " find_g_norg   " + present());
	}
	for_Int(i, 0, or_deg.size()) for_Int(n, 0, imp_i.nomgs) imp_i[n][i][i] = imp_i[n][idx[or_deg[i] - 1]][idx[or_deg[i] - 1]];
}

void NORG::readmatrix(MatReal& m, const Str& file)
{
	IFS ifs(file);
	if (!ifs) {
		if (mm) WRN(STR("file opening failed with ") + NAV(file))
	}
	else {
	for_Int(i, 0, m.nrows()) {
		for_Int(j, 0, m.ncols()) ifs >> m[i][j];
	}
	}
	ifs.close();
}

void NORG::get_g_by_CF(Green& imp_i)
{
	NocSpace scsp_1pone(mm, p, nppso(p.npartical, 1));
	NocSpace scsp_1mone(mm, p, nppso(p.npartical, -1));
	Operator n1pone(mm, p, scsp_1pone);
	Operator n1mone(mm, p, scsp_1mone);

	{//n=0 using the lanczos continue fraction.
		// scsp_1pone.coefficient = set_row_primeter_byimpH(uormat, impH);
		// scsp_1mone.coefficient = set_row_primeter_byimpH(uormat, impH);
		set_row_primeter_byimpH(uormat, impH, n1pone.oper_value);
		set_row_primeter_byimpH(uormat, impH, n1mone.oper_value);

		for_Int(i, 0, imp_i.g[0].nrows()){
			CrrltFun greaer(mm, p, scsp, scsp_1pone, n1pone.table, final_ground_state, 0 * 2, i);
			CrrltFun lesser(mm, p, scsp, scsp_1mone, n1mone.table, final_ground_state, 0 * 2, i);

			{// test
				if(mm) WRN(NAV2(greaer.ex_state.norm(), lesser.ex_state.norm()));
			}
			ImGreen green_function(1, p);
			greaer.find_gf_greater(groune_lst, green_function);
			lesser.find_gf_lesser(groune_lst, green_function);
			// green_function.write("continued_fraction-imp_green_function", 2 * i + 1);
			for_Int(n, 0, green_function.nomgs) imp_i[n][i][i] = green_function[n][0][0];
			if (mm) PIO("finished the find_g_norg");

		}
	}
}

VEC<MatReal> NORG::uormat_initialize()
{
	VEC<MatReal> uormat_i;
	if(p.if_norg_imp) {
			for_Int(i, 0, p.nO2sets.size()) {
			MatReal temp(dmat(p.nO2sets[i], 1.));
			uormat_i.push_back(std::move(temp));
		}
	} else {	
		for_Int(i, 0, p.nI2B.size()) {
			MatReal temp(dmat(p.nI2B[i], 1.));
			uormat_i.push_back(std::move(temp));
		}
	}
	return std::move(uormat_i);
}

bool NORG::converged() 
{
		const Real gseerr = ABS(energy_err);
		if (gseerr > 1.E-3) { return false; }
		else if (gseerr < 1.E-10) { return true; }
		else if (energy_err < 1.E-7) { 
			norg_stable_count++;
			if (norg_stable_count > 10) return true;
			}
		else if (iter_norg_cnt > 40) {
			const Real occerr = ABS(occupation_err);
			if (occerr < 1.E-8) {
				WRN("The iteration stoped since the occerr is less than 1.E-8.");
				return true;
			}
		}
		return false;
}

bool NORG::green_modify_converged() const
{
	const Real cverr = ABS(correctionerr);
	if (cverr < 1.E-5) { return true; }
	return false;
}

bool NORG::green_modify_converged(Real correctionerr_i) const
{
	const Real cverr = ABS(correctionerr_i);
	if (cverr < 1.E-5) { return true; }
	return false;
}

//------------------------------------------------------------------ spin-spin ------------------------------------------------------------------

Real NORG::sz_imp_sz_bath(const Int imp_postition, const VecReal& vgs_i)
{
	VecPartition row_H(mm.np(), mm.id(), scsp.dim);
	Real sz_imp_sz_bath_i(0.);
	// VecReal gstate_part = vgs_i.truncate(row_H.bgn(), row_H.end());
	for_Int(h_i, row_H.bgn(), row_H.end())
	{
		StateStatistics a(h_i, scsp.wherein_NocSpace(h_i), scsp);
		Int conuter = 0;
		// for imp up spin:
#define ndiv scsp.ndivs
		if(a.occ_n[imp_postition][0] == 1) {
			conuter += SUM(a.occ_n[imp_postition].truncate(1,5));
			conuter -= SUM(a.occ_n[imp_postition + 1].truncate(1,5));

		}
		// for imp down spin:
		if(a.occ_n[imp_postition + 1][0] == 1){
			conuter -= SUM(a.occ_n[imp_postition].truncate(1,5));
			conuter += SUM(a.occ_n[imp_postition + 1].truncate(1,5));
		}
		sz_imp_sz_bath_i += conuter * vgs_i[h_i] * vgs_i[h_i];
	}
	Real sz_imp_sz_bath = mm.Allreduce(sz_imp_sz_bath_i) / 4.;
	return sz_imp_sz_bath;
}

//------------------------------------------------------------------- private -------------------------------------------------------------------

void NORG::show_the_nozero_number_of_tabel() 
{
	// Real size_one(oneedm.table[2].size()), size_two(n1mone.table[2].size()), size_tree(n1pone.table[2].size());
	// LLInt size_of_main_t(mm.Allreduce(size_one)), size_of_main_n1mt(mm.Allreduce(size_two)), size_of_main_n1pt(mm.Allreduce(size_tree));
	// if(mm) PIO(NAV3(size_of_main_t, size_of_main_n1mt, size_of_main_n1pt));
	Real size_one(oneedm.table[2].size());
	LLInt size_of_main_t(mm.Allreduce(size_one));
	if(mm) PIO(NAV2(size_of_main_t, scsp.dim)+"   "+present());
}

MatReal NORG::save_transform_uormat(){
	MatReal transform_uormat(dmat(p.norbit, 1.));
	Int counter(0);
	if (!uormat.empty())
	{
		for (const auto& uormat_ii : uormat)
		{
			counter++;
			for_Int(i, 0, uormat_ii.nrows()) {
				for_Int(j, 0, uormat_ii.ncols()) {
					transform_uormat[i + counter][j + counter] = uormat_ii[i][j];
				}
			}
			counter += uormat_ii.nrows();
		}
	}
	return transform_uormat;
}

/*
VecReal NORG::set_row_primeter_byimpH(const VEC<MatReal>& uormat_i, const Impdata& impH_i)
{
	MatReal hopint(impH_i.first), transform_uormat(dmat(p.norbit, 1.));
	Int counter(0);

	if (!uormat_i.empty()) {
		for (const auto& uormat_ii : uormat_i) {
			if(p.if_norg_imp) counter = p.nO2sets[0] - p.nI2B[0];
			for_Int(i, 0, uormat_ii.nrows()) {
				for_Int(j, 0, uormat_ii.ncols()) {
					transform_uormat[i + counter][j + counter] = uormat_ii[i][j];
				}
			}
			counter += uormat_ii.nrows();
		}
	}

	// if (mm) WRN(NAV2(hopint, transform_uormat));
	hopint = transform_uormat * hopint * transform_uormat.ct();
	// Mat<MatReal> h_correcation = impH_i.second;
	// if(mm) WRN(present());
	VEC<Real> orbital_correction;
	if(p.if_norg_imp){
		for_Int(alpha, 0, p.norbit) {
			for_Int(eta, 0, p.norbit) {
				MatReal h_inter_bta_gmm(p.norbit, p.norbit, 0.);
				//? MatReal row_i(transform_uormat[alpha].mat(1, p.norbit)), col_l(transform_uormat[eta].mat(p.norbit, 1));
				//? MatReal temp_U_i_l = col_l * row_i;
				for_Int(i, 0, p.norbit) {
					for_Int(l, 0, p.norbit) {
						//// ADD(h_inter_bta_gmm, h_inter_bta_gmm,transform_uormat[alpha][i] * transform_uormat.ct()[l][eta] * (transform_uormat * impH_i.second[i][l] * transform_uormat.ct()));
						h_inter_bta_gmm += transform_uormat[alpha][i] * transform_uormat.ct()[l][eta] * (transform_uormat * impH_i.second[i][l] * transform_uormat.ct());
						//? h_inter_bta_gmm += temp_U_i_l[i][l] * (transform_uormat * impH_i.second[i][l] * transform_uormat.ct());
					}
				} 
				// h_correcation[alpha][eta] = h_inter_bta_gmm;
				for_Int(i_l, 0, p.norbit * p.norbit) orbital_correction.push_back(std::move(h_inter_bta_gmm.vec()[i_l]));
			}
		}
	}
	// if(mm) WRN(present()+NAV3(orbital_correction[6311], orbital_correction[6312], orbital_correction[6313]));

	VecReal coefficient_i(hopint.vec());
	VecReal check_point(1, 0.); coefficient_i.reset(concat(check_point, coefficient_i));

	return concat(coefficient_i, Vec(orbital_correction));
}

void NORG::set_row_primeter_byimpH(const VEC<MatReal>& uormat_i, const Impdata& impH_i, std::map<Int, Real>& oper_i)
{
	using namespace std;
	MatReal hopint(impH_i.first), transform_uormat(dmat(p.norbit, 1.));
	Int counter(0);
	if (uormat_i.empty() || oper_i.empty()) ERR("the uormat_i or oper_value is empty!!");
	
	if (!uormat_i.empty()) {
		for (const auto& uormat_ii : uormat_i) {
			if(!p.if_norg_imp) counter = p.nO2sets[0] - p.nI2B[0];
			for_Int(i, 0, uormat_ii.nrows()) {
				for_Int(j, 0, uormat_ii.ncols()) {
					transform_uormat[i + counter][j + counter] = uormat_ii[i][j];
				}
			}
			counter += uormat_ii.nrows();
		}
	}
	// if(mm) WRN(NAV(transform_uormat));
	hopint = transform_uormat * hopint * transform_uormat.ct();

	// if (mm) WRN(NAV(oper_i.size()));
	for (const auto& [key, value] : oper_i) {
		// if (mm) WRN(NAV2(key, value));
		if (key <= hopint.size()) oper_i[key] = hopint.vec()[key - 1];
		else if (key > hopint.size() && key <= std::pow(p.norbit, 4) + hopint.size()) {
			Int alpha, eta, beta, gamma, n(p.n_rot_orb), left(key - hopint.size() - 1);

			alpha 	= (left / (n * n)) / n;
			eta 	= (left / (n * n)) % n;
			beta 	= (left % (n * n)) / n;
			gamma 	= (left % (n * n)) % n;

			if(p.if_norg_imp) {
				Real v_aebg(0.);
				for_Int(i, 0, n) {
				for_Int(l, 0, n) {
				for_Int(j, 0, n) {
				for_Int(k, 0, n) {
					v_aebg += transform_uormat[alpha][i] * transform_uormat[eta][l] * \
					(transform_uormat[beta][j] * impH_i.second[i][l][k][j] * transform_uormat[gamma][k]);
				}
				}
				}
				}
				oper_i[key] = v_aebg;
			} else 
			{
				oper_i[key] = impH_i.second[alpha][eta][beta][gamma];
			}
		}
		else ERR("the key is out of range!!");
		// if (mm) WRN(NAV2(key, value));
	}
}
*/

void NORG::set_row_primeter_byimpH(const VEC<MatReal>& uormat_i, const Impdata& impH_i, VecReal& oper_i)
{
	using namespace std;
	MatReal hopint(impH_i.first), transform_uormat(dmat(p.norbit, 1.));
	Int counter(0), n(p.n_rot_orb);
	if (uormat_i.empty() || oper_i.size() == 0) ERR("the uormat_i or oper_value is empty!!");
	
	if (!uormat_i.empty()) {
		for (const auto& uormat_ii : uormat_i) {
			if(p.if_norg_imp == false) counter += p.nO2sets[0] - p.nI2B[0];
			for_Int(i, 0, uormat_ii.nrows()) {
				for_Int(j, 0, uormat_ii.ncols()) {
					transform_uormat[i + counter][j + counter] = uormat_ii[i][j];
				}
			}
			counter += uormat_ii.nrows();
		}
	}
	{// test the transform_uormat and density matrix.
		MatReal up0_dm(see_MatReal(oneedm.dm).truncate(0, 0, p.nO2sets[0], p.nO2sets[0])), up0_U(transform_uormat.truncate(0, 0, p.nO2sets[0], p.nO2sets[0]));
		if (mm) WRN(NAV2(up0_dm, up0_U));
	}
	hopint = transform_uormat * hopint * transform_uormat.ct();

	// alpha:a, beta=b, gamma=g, eta=e
	if (p.if_norg_imp) {
	VecReal ijke = (impH_i.second.mat(n * n * n, n) * transform_uormat.ct()).vec();
	VecReal ajke = (transform_uormat * ijke.mat(n, n * n * n)).vec();
	VecReal	keaj = ajke.mat(n * n, n * n).tr().vec();
	// VecReal geaj = (transform_uormat.ct().tr() * keaj.mat(n, n * n * n)).vec();
	// VecReal geab = (geaj.mat(n * n * n, n) * transform_uormat.tr()).vec();
	VecReal geaj = (transform_uormat * keaj.mat(n, n * n * n)).vec();
	VecReal geab = (geaj.mat(n * n * n, n) * transform_uormat.ct()).vec();
	VecReal abge = (geab.mat(n * n, n * n).tr()).vec();
	oper_i = concat(concat(VecReal{ 0. }, hopint.vec()), abge);
	} else { 
		oper_i = concat(concat(VecReal{ 0. }, hopint.vec()), impH_i.second);
	}
	// if(mm) WRN(NAV(transform_uormat));
}


//------------------------------------------------------------------ io ------------------------------------------------------------------

void NORG::write_norg_info(Int iter_cnt) const {
	using namespace std;
	Real sc = 0;
	for_Int(i, 0, final_ground_state.size()) sc -= SQR(final_ground_state[i]) * log(SQR(final_ground_state[i]));
	OFS ofs_app_scorrelation(tox + "slater_correlation.txt", std::ios::app);
	ofs_app_scorrelation << iofmt("sci");
	ofs_app_scorrelation << setw(4) << iter_cnt;
	ofs_app_scorrelation << "\t" << setw(w_Real) << sc;
	ofs_app_scorrelation << endl; ofs_app_scorrelation.close();

	OFS ofs_app_energy(tox + "ground_energy.txt", std::ios::app);
	ofs_app_energy << iofmt("sci");
	ofs_app_energy << setw(4) << iter_cnt;
	ofs_app_energy << "\t" << setw(w_Real) << groune_lst;
	ofs_app_energy << endl; ofs_app_energy.close();

	// Only show for the upper spin.
	if(iter_cnt != 0)for_Int(i, 0, p.nband) {
		OFS ofs_app_occupation(tox +"imp"+ STR(i*2) + "_bath_occupation_number.txt", std::ios::app);
		ofs_app_occupation << iofmt("sci");
		ofs_app_occupation << setw(4) << iter_cnt;
		for_Int(j, 0, p.nI2B[i*2]) { ofs_app_occupation << "  " << setw(w_Real) << occnum_lst[SUM_0toX(p.nI2B, i*2) + j]; }
		ofs_app_occupation << endl;
		ofs_app_occupation.close();
	}
	
	// Only show for the upper spin.
	if(iter_cnt != 0)for_Int(i, 0, p.nband) {
		VEC<VecReal> dmeigen = oneedm.check_dm_get_occupation_number();
		OFS ofs_app_occupation_all(tox +"imp"+ STR(i*2) + "_occupation_number.txt", std::ios::app);
		ofs_app_occupation_all << iofmt("sci");
		ofs_app_occupation_all << setw(4) << iter_cnt;
		for_Int(j, 0, dmeigen[i].size()) { ofs_app_occupation_all << "  " << setw(w_Real) << dmeigen[i][j]; }
		ofs_app_occupation_all << endl;
		ofs_app_occupation_all.close();
	}

}

void NORG::write_occupation_info() const {
	using namespace std;
	OFS ofs; ofs.open("nmat.txt");
	VEC<MatReal> dmtemp(oneedm.dm);
	VecReal counter(3);
	ofs << "#   < n_i >   data:"<< endl;
	for_Int(orb_i, 0, p.norbs)	{
		ofs << iofmt();
		std::string temp = (orb_i % 2) == 0 ? STR(Int(orb_i / 2) + 1) + "up" : STR(Int(orb_i / 2) + 1) + "dn";
		ofs << setw(6) << temp << setw(p_Real) << dmtemp[orb_i][0][0] << endl;
		(orb_i % 2) == 0 ? counter[0] += dmtemp[orb_i][0][0] : counter[1] += dmtemp[orb_i][0][0];
	}
	counter[2] = counter[0] + counter[1];
	ofs << setw(6) << "sup" << setw(p_Real) << counter[0] << endl;
	ofs << setw(6) << "sdn" << setw(p_Real) << counter[1] << endl;
	ofs << setw(6) << "sum" << setw(p_Real) << counter[2];
	ofs.close();
}

void NORG::write_state_info(Int iter_cnt) const {
	using namespace std;
	// OFS ofs_app_state(STR(iter_cnt) + "write_state_info.txt", ios::app);
	OFS ofs_app_state; ofs_app_state.open(STR(iter_cnt) + "write_state_info.txt");
	ofs_app_state << setw(4) << "iter_cnt";
	ofs_app_state << "\t" << setw(w_Real) << "state";
	ofs_app_state << "\t" << setw(w_Real) << "state_norm";
	ofs_app_state << "\t" << setw(w_Real) << "norm__sort";
	ofs_app_state << endl; 
	VecReal temp_norm_state = SQR(final_ground_state);
	VecReal temp_norm__sort = temp_norm_state; sort(temp_norm__sort);
	if(mm) WRN(NAV3(final_ground_state.norm(),temp_norm_state.norm(),temp_norm__sort.norm()));
	for_Int(i, 0, final_ground_state.size()) 
	{
		ofs_app_state << iofmt("sci");
		ofs_app_state << setw(4) << iter_cnt;
		ofs_app_state << "\t" << setw(w_Real) << final_ground_state[i];
		ofs_app_state << "\t" << setw(w_Real) << temp_norm_state[i];
		ofs_app_state << "\t" << setw(w_Real) << temp_norm__sort[final_ground_state.size()-1-i];
		ofs_app_state << endl; 
	}
	ofs_app_state.close();
}

void NORG::write_impurtiy_occupation() const {
	using namespace std;
	if (mm) {
		OFS ofs;
		ofs.open("nmat.txt");
		VecReal counter(3);
		MatReal dm_origional = see_MatReal(uormat).ct() * see_MatReal(oneedm.dm) * see_MatReal(uormat);
		// WRN(NAV(dm_origional));
		VecReal particals(dm_origional.diagonal().mat(p.norg_sets, p.nO2sets[0]).tr()[0]);

		ofs << "#   < n_i >   data:" << endl;
		for_Int(orb_i, 0, p.norbs) {
			ofs << iofmt();
			std::string temp = (orb_i % 2) == 0 ? STR(Int(orb_i / 2) + 1) + "up" : STR(Int(orb_i / 2) + 1) + "dn";
			ofs << setw(6) << temp << setw(p_Real) << particals[orb_i] << endl;
		}
		counter[0] = SUM(particals.mat(p.nband, 2).tr()[0]);
		counter[1] = SUM(particals.mat(p.nband, 2).tr()[1]);

		counter[2] = counter[0] + counter[1];
		ofs << setw(6) << "sup" << setw(p_Real) << counter[0] << endl;
		ofs << setw(6) << "sdn" << setw(p_Real) << counter[1] << endl;
		ofs << setw(6) << "sum" << setw(p_Real) << counter[2];
		ofs.close();

		std::cout << "\nnorg ground state energy: " << groune_lst << "  " << present() << std::endl;
		std::cout << std::endl;							// blank line
	}
}



Real NORG::return_nointeractions_ground_state_energy(const MatReal& h0_i) const {
	MatReal htemp = h0_i;
	VecReal vtemp(htemp.ncols(),0.);
	Real energy = 0.;
	heevr(htemp, vtemp);
	for_Int(i, 0, vtemp.size()) {
		if(vtemp[i]<0) energy += vtemp[i];
	}
	return energy;
}
