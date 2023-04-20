/*
coded by Jia-Ming Wang (jmw@ruc.edu.cn, RUC, China) date 2022
*/
//! not finished yet 


#include "dmft_bethe.h"

DMFT::DMFT(const MyMpi& mm_i, const Prmtr& prmtr_i, const Int mode) :
	mm(mm_i), p(prmtr_i), iter_cnt(0), g_loc(p.nband, p), n_1(0.), n_2(0.),
	gloc_err({ 1.e99,1.e98,1.e97 ,1.e96 ,1.e95}), se(p.nband, p), imp_backup(p.imp_backup)
{
	// make random seed output together
	{ mm.barrier(); SLEEP(1); }
	log("initial");	set_parameter();
	if(mm) WRN("TO here is fine");
	Bath bth(mm, p);
	Impurity imp(mm, p, bth);
	NORG norg(mm, p);
	g_loc = g0_loc(); 				if (mm) g_loc.write("g_0loc", iter_cnt);
	// if (mode == 0) {
	// 	do_you_also_have_these_backup(bth, norg, imp_backup);
	// 	if (!check_gloc("000.mb.g_imp.txt"))  {
	// 		if (mm && imp_backup) WRN("The initial g_imp wasn't exist, but don't worry, I will use impurity's info.");
	// 		if (mm &&!imp_backup) WRN("Now we will start with the new case.");
	// 	}
	// }
	// if (mode == 1) {
	// 	do_you_also_have_these_backup(bth, norg, imp_backup);
	// 	if (!check_seloc("000.mb.se_imp.txt")) {
	// 		if (mm && imp_backup) WRN("The initial se_imp wasn't exist, but don't worry, I will use impurity's info.");
	// 		if (mm &&!imp_backup) WRN("Now we will start with the new case.");
	// 	}
	// }
	// if (mm) bth.write_ose_hop(iter_cnt);						if (mm) g_loc.write("g_0loc", iter_cnt);
	while (iter_cnt < p.iter_max && !converged()) 
	{
		++iter_cnt;ImGreen hb(p.norbs, p);
		if(mode == 1) {hb = find_hb_by_se(se);					if (mm) hb.write("hb", iter_cnt);}
		if(mode == 0) {hb = find_hb(g_loc);						if (mm) hb.write("hb", iter_cnt);}
		if (!imp_backup) bth.bath_fit(hb, iter_cnt);
		imp.update("bethe");											if (mm) bth.write_ose_hop(iter_cnt);
		ImGreen hb_imp(p.norbs, p);   	imp.find_hb(hb_imp); 	if (mm) hb_imp.write("hb_imp", iter_cnt);
		ImGreen g0_imp(p.norbs, p);   	imp.find_g0(g0_imp); 	if (mm) g0_imp.write("g0_imp", iter_cnt);
		norg.up_date_h0_to_solve(imp.impH, 1);

		ImGreen g0imp(p.nband, p);	imp.find_g0(g0imp);							if (mm)	g0imp.write("g0imp");
		ImGreen gfimp(p.nband, p);	norg.get_gimp(gfimp);						if (mm) gfimp.write("gfimp");
		ImGreen seimp(p.nband, p);	seimp = g0imp.inverse() - gfimp.inverse();	if (mm) seimp.write("seimp");

		if (mode == 0) {
			append_gloc_err(gfimp.error(g_loc));				log("gerr_update");
			g_loc = gfimp;
		}
		// obtain se from impurity model.
		if (mode == 1) {
			append_gloc_err(se.error(seimp));					log("gerr_update");
			se = seimp;	g_loc = find_gloc_by_se(se);
		}
		imp_backup = false;										if (mm) save_the_backup(bth, norg, iter_cnt);
	}
		ReGreen g0_imp_re(p.norbs, p);imp.find_g0(g0_imp_re);		if (mm) g0_imp_re.write("g0_imp_re");
		ReGreen gfimp_re(p.nband, p);	norg.get_gimp(gfimp_re);	if (mm) gfimp_re.write("gfimp");
		ReGreen se = g0_imp_re.inverse() - gfimp_re.inverse();		if (mm) se.write("se_loc_re");
		// ReGreen g_loc(p.norbs, p);find_gloc_by_se(g_loc,se);		if (mm) g_loc.write("g_loc_re");
}


// -------------------------------------------------- bethe lattice --------------------------------------------------

ImGreen DMFT::find_gloc_by_se(const ImGreen& se_i) const
{
	ImGreen gloc_temp(g_loc.norbs, p);
	for_Int(n, 0, gloc_temp.nomgs) {
		for_Int(i,0,p.norbs){
			Cmplx z = gloc_temp.z(n) + bethe_mu - se_i[n][i][i];
			Cmplx temp = SQRT(SQR(z) - 4 * SQR(bethe_t[i]));
			if(imag(temp) < 0) temp *= -1;
			// temp = real(temp) + I * ABS(imag(temp));
			gloc_temp[n][i][i] = (z - temp)/ (2 * SQR(bethe_t[i]));
		}
	}
	return gloc_temp;
}

ImGreen DMFT::find_hb(const ImGreen& g_loc) const
{
	ImGreen hb(p.norbs, p);
	for_Int(i,0,p.norbs){
		for_Int(n, 0, g_loc.nomgs){
			hb[n][i][i]=-bethe_t[i] * bethe_t[i] * g_loc[n][i][i];
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
			Cmplx z = g_loc.z(n) + bethe_mu - se_i[n][i][i];
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

ImGreen DMFT::find_imp_se(const ImGreen& g_imp, const ImGreen& hb_imp) const{
	ImGreen se(p.norbs, p);
	// ImGreen g_loc_inv = g_loc.inverse();
	for_Int(n, 0, se.nomgs) {
		for_Int(i, 0, p.norbs){
			se[n][i][i] = g_imp.z(n) + bethe_mu + g_imp.inverse()[n][i][i] - hb_imp[n][i][i];
		}
	}
	return se;
}

void DMFT::save_the_backup(Bath& bth, NORG& solver, Int iter_cnt){
	{
		OFS ofs;
    	if (iter_cnt == 999) { ofs.open( "output.ose"); }
		else ofs.open("zic" + prefill0(iter_cnt, 3) + "ose");
		for_Int(i, 0, p.norbs) biwrite(ofs, CharP(bth.vec_ose[i].p()), bth.vec_ose[i].szof());
	}{
		OFS ofs;
    	if (iter_cnt == 999) { ofs.open( "output.hop"); }
		else ofs.open("zic" + prefill0(iter_cnt, 3) + "hop");
		for_Int(i, 0, p.norbs) biwrite(ofs, CharP(bth.vec_hop[i].p()), bth.vec_hop[i].szof());
	}{
		OFS ofs;
    	if (iter_cnt == 999) { ofs.open( "output.rotuation_u"); }
		ofs.open("zic" + prefill0(iter_cnt, 3) + "rotuation_u");
		for_Int(i, 0, solver.uormat.size()) biwrite(ofs, CharP(solver.uormat[i].p()), solver.uormat[i].szof());
	}
}

void DMFT::do_you_also_have_these_backup(Bath& bth, NORG& solver, bool if_backuped){
	if(if_backuped) {
		{
			Str na_ose(prefill0(0, 3) + "ose");
			IFS ifs(na_ose);
			for_Int(i, 0, p.norbs) biread(ifs, CharP(bth.vec_ose[i].p()), bth.vec_ose[i].szof());
		}{
			Str na_hop(prefill0(0, 3) + "hop");
			IFS ifs(na_hop);
			for_Int(i, 0, p.norbs) biread(ifs, CharP(bth.vec_hop[i].p()), bth.vec_hop[i].szof());
		}{
			Str na_rotuation_u(prefill0(0, 3) + "rotuation_u");
			IFS ifs(na_rotuation_u);
			for_Int(i, 0, solver.uormat.size()) biread(ifs, CharP(solver.uormat[i].p()), solver.uormat[i].szof());
		}
		for_Int(i, 0, p.norbs){
			slctsort(bth.vec_ose[i], bth.vec_hop[i]);
			bth.vec_hop[i] = ABS(bth.vec_hop[i]);
		}
	}
}

// ------------------------------- since the 2023.04.15 -------------------------------

// bethe lattice part:
ImGreen DMFT::g0_loc() const {
	ImGreen g0_loc(p.nband, p);
	for_Int(n, 0, p.num_omg) {
		for_Int(i, 0, p.norbs) {
			g0_loc[n][i][i] = (g0_loc.z(n) + bethe_mu - SQRT(SQR(g0_loc.z(n) + bethe_mu) - 4 * SQR(bethe_t[i]))) / (2 * bethe_t[i] * bethe_t[i]);
		}
	}
	return g0_loc;
}

void DMFT::set_parameter() {
	bethe_mu = 0.15;
	bethe_t.reset(Vec{0.5, 0.25});
	p.U = bethe_u = 3.0;
	p.Uprm = bethe_u12 = bethe_u - 0.3;
	p.jz = 0.;
}
