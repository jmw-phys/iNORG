/*
coded by Jia-Ming Wang (jmw@ruc.edu.cn, RUC, China) date 2022
*/
//! not finished yet 


#include "dmft_bethe.h"


DMFT::DMFT(const MyMpi& mm_i, const Prmtr& prmtr_i, const Model& mdl_i, const Int mode) :
	mm(mm_i), p(prmtr_i), mdl(mdl_i), iter_cnt(0), g_loc(p.norbs, p),
	gloc_err({ 1.e99,1.e98,1.e97 ,1.e96 ,1.e95}), se(p.norbs, p), imp_backup(p.imp_backup),
	//se_err({ 1.e99 }),
	n_1(0.), n_2(0.)
{
	// make random seed output together
	{ mm.barrier(); SLEEP(1); }
	log("initial");
	Model bethe(mm, p);
	Bath bth(mm, p);
	Impurity imp(mm, p, bth, mdl);
	NORG norg(mm, p);
	g_loc = mdl.g0_loc(); 
	if(mode == 0 ){
		do_you_also_have_these_backup(bth, norg, imp_backup);
		if (!check_gloc(bix + "000.mb.g_imp.txt"))  {
			if (mm && imp_backup) WRN("The initial g_imp wasn't exist, but don't worry, I will use impurity's info.");
			if (mm &&!imp_backup) WRN("Now we will start with the new case.");
		}
	}
	if (mode == 1){
		do_you_also_have_these_backup(bth, norg, imp_backup);
		if (!check_seloc(bix + "000.mb.se_imp.txt")) {
			if (mm && imp_backup) WRN("The initial se_imp wasn't exist, but don't worry, I will use impurity's info.");
			if (mm &&!imp_backup) WRN("Now we will start with the new case.");
		}
	}
	if (mm) bth.write_ose_hop(iter_cnt);						if (mm) g_loc.write("g_0loc", iter_cnt);
	while (iter_cnt < p.iter_max && !converged()) 
	{
		++iter_cnt;ImGreen hb(p.norbs, p);
		if(mode == 1) {hb = find_hb_by_se(se);					if (mm) hb.write("hb", iter_cnt);}
		if(mode == 0) {hb = find_hb(g_loc);						if (mm) hb.write("hb", iter_cnt);}
		if(!imp_backup) bth.number_bath_fit(hb,iter_cnt, 0);
		imp.update();											if (mm) bth.write_ose_hop(iter_cnt);
		ImGreen hb_imp(p.norbs, p);   	imp.find_hb(hb_imp); 	if (mm) hb_imp.write("hb_imp", iter_cnt);
		ImGreen g0_imp(p.norbs, p);   	imp.find_g0(g0_imp); 	if (mm) g0_imp.write("g0_imp", iter_cnt);
		norg.up_date_h0_to_solve(imp.h0); 
		norg.get_gimp(); ImGreen g_imp = norg.impgreen; 		if (mm) g_imp.write("g_imp", iter_cnt);
		// ImGreen gmb_imp_modifyed(p.norbs, p);	norg.get_gimp_by_krylov_CV_modify(gmb_imp_modifyed);	if (mm) gmb_imp_modifyed.write("gmb_imp_modifyed", p.iter_max_cv_norg);
		// ReGreen gre_imp_modifyed(p.norbs, p);	norg.get_gimp_by_krylov_CV_modify(gre_imp_modifyed);	if (mm) gre_imp_modifyed.write("gre_imp_modifyed", p.iter_max_cv_norg);
		// norg.modify_by_krylov_correct_vec(); 

		for_Int(n, 0, g_loc.nomgs) g_imp[n] -= cmplx(real(g_imp[n]));
		ImGreen se_imp = g0_imp.inverse() - g_imp.inverse();	if (mm) se_imp.write("se_imp", iter_cnt);

		// ImGreen g_imp_krylov(p.norbs, p); norg.get_gimp_by_krylov(g_imp_krylov); if (mm) g_imp_krylov.write("krylovG", iter_cnt);
		ImGreen g_imp_krylov_2(p.norbs, p); norg.get_g_by_krylov_CV(g_imp_krylov_2); if (mm) g_imp_krylov_2.write("krylovG_2side", iter_cnt);
		// norg.get_gimp_with_possible_degeneracy(g_imp, iter_cnt);if (mm) g_imp.write("g_imp", iter_cnt);

        // g_loc = 0.5 * g_loc + 0.5 * g_imp;
		// g_loc = 0.3 * g_loc + 0.7 * g_imp;

		if (mode == 0) {
			append_gloc_err(g_imp.error(g_loc));				log("gerr_update");
			g_loc = g_imp;
		}
		// obtain se from impurity model.
		if (mode == 1) {
			append_gloc_err(se.error(se_imp));					log("gerr_update");
			se = se_imp;	g_loc = find_gloc_by_se(se);
		}
		imp_backup = false;										if (mm) save_the_backup(bth, norg, iter_cnt);
	}
	ImGreen g_imp_mb(p.norbs, p);
	Int degeneracy(0), countk(0);
	degeneracy=norg.get_gimp_with_possible_degeneracy(g_imp_mb);	if (mm)	g_imp_mb.write("g_imp_mb");
	ImGreen g0_imp_mb(p.norbs, p);   	imp.find_g0(g0_imp_mb); 	if (mm) g0_imp_mb.write("g0_imp_mb");
	ImGreen se_mb = g0_imp_mb.inverse() - g_imp_mb.inverse();		if (mm) se_mb.write("se_loc_mb");
	if(mm) PIO(NAV(degeneracy));									if (mm) save_the_backup(bth, norg);
	VEC<ReGreen> g_imp_re_vec;
	norg.find_g_with_possible_degeneracy(g_imp_re_vec);
	for (const auto &g_imp_re : g_imp_re_vec){
		if (mm)	g_imp_re.write("g_imp_re_with_eta"+STR(p.a[countk]));	
		ReGreen g0_imp_re(p.norbs, p);imp.find_g0(g0_imp_re,countk);if (mm) g0_imp_re.write("g0_imp_re_with_eta"+STR(p.a[countk]));
		ReGreen se = g0_imp_re.inverse() - g_imp_re.inverse();		if (mm) se.write("se_loc_with_eta"+STR(p.a[countk]));
		ReGreen g_loc(p.norbs, p);find_gloc_by_se(g_loc,se,countk);	if (mm) g_loc.write("g_loc_with_eta"+STR(p.a[countk]));
		countk++;
	}
	// VecReal ai({0.0, 0.05, 0.1, 0.2});							if (mm) print_gloc_by_sefe_with_different_eta(se, ai);
	// VecReal aj({0.0, 0.05, 0.1, 0.2});							if (mm) print_gloc_by_gimp_with_different_eta(g_re, aj);
}


// -------------------------------------------------- bethe lattice --------------------------------------------------

ImGreen DMFT::find_weiss(ImGreen& g_loc) const
{
	ImGreen weiss(g_loc.norbs, p);
	for_Int(n, 0, g_loc.nomgs) {
		for_Int(i,0,p.norbs){
			weiss[n][i][i] = INV((g_loc.z(n) + p.bethe_mu - p.t[i] * p.t[i] * g_loc[n][i][i]));
		}
	}
}


void DMFT::print_gloc_by_sefe_with_different_eta(const Green& se_i, VecReal a) const{
	// ReGreen g_loc(p.norbs, p);
	// if (a.size() == 0) {
	// 	for_Int(n, 0, g_loc.nomgs) {
	// 		for_Int(i, 0, p.norbs) {
	// 			Cmplx z = g_loc.z(n) + p.bethe_mu - se_i[n][i][i];
	// 			Cmplx temp = SQRT(SQR(z) - 4. * SQR(p.t[i]));
	// 			if (imag(temp) < 0)	temp *= -1.;
	// 			g_loc[n][i][i] = (z - temp) / (2. * SQR(p.t[i]));
	// 		}
	// 	}
	// 	if(mm) g_loc.write("g_loc");
	// }
	// else {
	// 	for_Int(a_i, 0, a.size()) {
	// 		for_Int(n, 0, g_loc.nomgs) {
	// 			for_Int(i, 0, p.norbs) {
	// 				Cmplx z = g_loc.z(n) + I * a[a_i] * ABS(real(g_loc.z(n)))  + p.bethe_mu - se_i[n][i][i];
	// 				Cmplx temp = SQRT(SQR(z) - 4. * SQR(p.t[i]));
	// 				if (imag(temp) < 0) temp *= -1.;
	// 				g_loc[n][i][i] = (z - temp) / (2. * SQR(p.t[i]));
	// 			}
	// 		}
	// 	if(mm) g_loc.write("g_loc_with"+STR(a[a_i]));
	// 	}
	// }
}

void DMFT::print_gloc_by_gimp_with_different_eta(const Green& gimp, VecReal a) const{

}

void DMFT::find_gloc_by_se(Green& g_loc, const Green& se_i, Int kind) const{
	for_Int(n, 0, g_loc.nomgs) {
		for_Int(i,0,p.norbs){
			Cmplx z = g_loc.z(n, kind) + p.bethe_mu - se_i[n][i][i];
			// if(mm) WRN(NAV(g_loc.z(n, kind)));
			Cmplx temp = SQRT(SQR(z) - 4 * SQR(p.t[i]));
			if(imag(temp) < 0) temp *= -1;
			g_loc[n][i][i] = (z - temp)/ (2 * SQR(p.t[i]));
		}
	}
}

ImGreen DMFT::find_gloc_by_se(const ImGreen& se_i) const
{
	ImGreen gloc_temp(g_loc.norbs, p);
	for_Int(n, 0, gloc_temp.nomgs) {
		for_Int(i,0,p.norbs){
			Cmplx z = gloc_temp.z(n) + p.bethe_mu - se_i[n][i][i];
			Cmplx temp = SQRT(SQR(z) - 4 * SQR(p.t[i]));
			if(imag(temp) < 0) temp *= -1;
			// temp = real(temp) + I * ABS(imag(temp));
			gloc_temp[n][i][i] = (z - temp)/ (2 * SQR(p.t[i]));
		}
	}
	return gloc_temp;
}

ImGreen DMFT::find_hb(const ImGreen& g_loc) const
{
	ImGreen hb(p.norbs, p);
	for_Int(i,0,p.norbs){
		for_Int(n, 0, g_loc.nomgs){
			hb[n][i][i]=-p.t[i] * p.t[i] * g_loc[n][i][i];
			hb[n][i][i] -= real(hb[n][i][i]);
		}
	}
	return hb;
}

ImGreen DMFT::find_hb_by_se(const ImGreen& se_i) const
{
	ImGreen hb(p.norbs, p);
	for_Int(i,0,p.norbs){
		for_Int(n, 0, g_loc.nomgs){
			Cmplx z = g_loc.z(n) + p.bethe_mu - se_i[n][i][i];
			hb[n][i][i]= -z + g_loc.inverse()[n][i][i];
		}
	}
	return hb;
}

bool DMFT::check_gloc(const Str& file){
	IFS ifs(file);
	if(ifs && !file.empty()) {
		ImGreen g_loc_temp(p,file);
		g_loc = g_loc_temp;
		return true;
	}
	return false;
}

bool DMFT::check_seloc(const Str& file){
	IFS ifs(file);
	if(ifs && !file.empty()) {
		ImGreen g_loc_temp(p,file);
		ImGreen se_temp(p,file);
		se = se_temp;
		g_loc = find_gloc_by_se(se_temp);
		return true;
	}
	return false;
}

// Wrong?!
ReGreen DMFT::find_lattice_se(const ReGreen& g_loc_Re) const{
	ReGreen se(p.norbs, p);
	ReGreen g_loc_inv = g_loc_Re.inverse();
	for_Int(n, 0, g_loc_Re.nomgs) {
		for_Int(i,0,p.norbs){
			se[n][i][i] = g_loc_Re.z(n, 0.) + p.bethe_mu - p.t[i] * p.t[i] * g_loc_Re[n][i][i] - g_loc_inv[n][i][i];
		}
	}
	return se;
}

ImGreen DMFT::find_imp_se(const ImGreen& g_imp, const ImGreen& hb_imp) const{
	ImGreen se(p.norbs, p);
	// ImGreen g_loc_inv = g_loc.inverse();
	for_Int(n, 0, se.nomgs) {
		for_Int(i, 0, p.norbs){
			se[n][i][i] = g_imp.z(n) + p.bethe_mu + g_imp.inverse()[n][i][i] - hb_imp[n][i][i];
		}
	}
	return se;
}

void DMFT::save_the_backup(Bath& bth, NORG& solver, Int iter_cnt){
	{
		OFS ofs;
    	if (iter_cnt == 999) { ofs.open(iox + "output." + "ose"); }
		else ofs.open(bix + "zic" + prefill0(iter_cnt, 3) + "ose");
		for_Int(i, 0, p.norbs) biwrite(ofs, CharP(bth.vec_ose[i].p()), bth.vec_ose[i].szof());
	}{
		OFS ofs;
    	if (iter_cnt == 999) { ofs.open(iox + "output." + "hop"); }
		else ofs.open(bix + "zic" + prefill0(iter_cnt, 3) + "hop");
		for_Int(i, 0, p.norbs) biwrite(ofs, CharP(bth.vec_hop[i].p()), bth.vec_hop[i].szof());
	}{
		OFS ofs;
    	if (iter_cnt == 999) { ofs.open(iox + "output." + "rotuation_u"); }
		ofs.open(bix + "zic" + prefill0(iter_cnt, 3) + "rotuation_u");
		for_Int(i, 0, solver.uormat.size()) biwrite(ofs, CharP(solver.uormat[i].p()), solver.uormat[i].szof());
	}
}

void DMFT::do_you_also_have_these_backup(Bath& bth, NORG& solver, bool if_backuped){
	if(if_backuped) {
		{
			Str na_ose(bix + prefill0(0, 3) + "ose");
			IFS ifs(na_ose);
			for_Int(i, 0, p.norbs) biread(ifs, CharP(bth.vec_ose[i].p()), bth.vec_ose[i].szof());
		}{
			Str na_hop(bix + prefill0(0, 3) + "hop");
			IFS ifs(na_hop);
			for_Int(i, 0, p.norbs) biread(ifs, CharP(bth.vec_hop[i].p()), bth.vec_hop[i].szof());
		}{
			Str na_rotuation_u(bix + prefill0(0, 3) + "rotuation_u");
			IFS ifs(na_rotuation_u);
			for_Int(i, 0, solver.uormat.size()) biread(ifs, CharP(solver.uormat[i].p()), solver.uormat[i].szof());
		}
		for_Int(i, 0, p.norbs){
			slctsort(bth.vec_ose[i], bth.vec_hop[i]);
			bth.vec_hop[i] = ABS(bth.vec_hop[i]);
		}
	}
}


//------------------------------------------------------ (Deprecated) -----------------------------------------------------------
