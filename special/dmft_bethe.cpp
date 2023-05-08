/*
coded by Jia-Ming Wang (jmw@ruc.edu.cn, RUC, China) date 2022
*/
//! not finished yet 


#include "dmft_bethe.h"


void DMFT::set_parameter() {
	bethe_mu = p.mu;
	bethe_t.reset(Vec{0.5, 0.25});
	bethe_u = p.U;
	bethe_u12 = p.Uprm;
	p.jz = 0.;
	p.bandw = 2 * SQRT(SQR(bethe_u) + SQR(bethe_u12) + SUM(bethe_t * bethe_t));
	p.derive_ImGreen();
}

DMFT::DMFT(const MyMpi& mm_i, Prmtr& prmtr_i, const Int mode) :
	mm(mm_i), p(prmtr_i), iter_cnt(0), g_loc(p.nband, p), n_eles(p.norbs, 0.),
	gloc_err({ 1.e99,1.e98,1.e97 ,1.e96 ,1.e95}), se(p.nband, p), imp_backup(p.imp_backup)
{
	// make random seed output together
	{ mm.barrier(); SLEEP(1); }
	IFS fitdata(prefill0(p.nI2B[0], 2) + ".ose_hop"), mbgfdata("mb.gfimp");
	log("initial");	set_parameter();
	Bath bth(mm, p); 
	Impurity imp(mm, p, bth);
	if(fitdata) {bth.read_ose_hop();		imp.update();	if (mm) imp.write_H0info(bth, -1, iter_cnt);}
	if(mbgfdata) {g_loc = ImGreen(p.nband, p, "mb.gfimp");	if (mm) g_loc.write("g_loc", iter_cnt);}
	else {g_loc = g0_loc();									if (mm) g_loc.write("g_0loc", iter_cnt);}
	
    VEC<MatReal> norg_tempU;
	while (iter_cnt < p.iter_max && !converged()) 
	{
		++iter_cnt;ImGreen hb(p.nband, p);
		if(!(fitdata && iter_cnt == 1)){
			if(mode == 1) {hb = find_hb_by_se(se);					if (mm) hb.write("hb", iter_cnt);}
			if(mode == 0) {hb = find_hb(g_loc);						if (mm) hb.write("hb", iter_cnt);}
			bth.bath_fit(hb, VecInt{1,2});							if (mm) bth.write_ose_hop(iter_cnt);
		}
		ImGreen hb_imp(p.nband, p);   	imp.find_hb(hb_imp); 					if (mm) hb_imp.write("hb-fit", iter_cnt);
		imp.update();															if (mm) imp.write_H0info(bth, -1, iter_cnt);
		auto_nooc("ful_pcl_sch", imp);	NORG norg(mm, p);
		norg.up_date_h0_to_solve(imp.impH, 1);	n_eles = norg.write_impurtiy_occupation(iter_cnt);
		ImGreen g0imp(p.nband, p);	imp.find_g0(g0imp);							if (mm)	g0imp.write("g0imp", iter_cnt);
		ImGreen gfimp(p.nband, p);	norg.get_gimp(gfimp);						if (mm) gfimp.write("gfimp", iter_cnt);
		ImGreen seimp(p.nband, p);	seimp = g0imp.inverse() - gfimp.inverse();	if (mm) seimp.write("seimp", iter_cnt);

		if (mode == 0) {
			append_gloc_err(gfimp.error(g_loc));				log("gerr_update");
			g_loc = 0.3 * g_loc + 0.7 * gfimp;
		}

		// obtain se from impurity model.
		if (mode == 1) {
			append_gloc_err(se.error(seimp));					log("gerr_update");
			se = seimp;	g_loc = find_gloc_by_se(se);
		}
		norg_tempU = norg.uormat;
	}

	// ! auto_nooc("for_green_calc", imp); // in here we can change the nooc space once, last chance.
	NORG finalrg(mm, p); 		finalrg.uormat = norg_tempU;	finalrg.up_date_h0_to_solve(imp.impH, 1);
	ImGreen gfimp(p.nband, p);	finalrg.get_gimp(gfimp);		if (mm) gfimp.write("gfimp", iter_cnt);
	ReGreen g0_imp_re(p.nband, p);imp.find_g0(g0_imp_re);		if (mm) g0_imp_re.write("g0fimp");
	ReGreen gfimp_re(p.nband, p);	finalrg.get_gimp(gfimp_re);	if (mm) gfimp_re.write("gfimp");
	ReGreen se = g0_imp_re.inverse() - gfimp_re.inverse();		if (mm) se.write("se_loc");

		
	// if (mm) bth.write_ose_hop();
	// ReGreen g_loc(p.nband, p);find_gloc_by_se(g_loc,se);		if (mm) g_loc.write("g_loc_re");
}


// -------------------------------------------------- bethe lattice --------------------------------------------------

ImGreen DMFT::find_gloc_by_se(const ImGreen& se_i) const
{
	ImGreen gloc_temp(g_loc.norbs, p);
	for_Int(n, 0, gloc_temp.nomgs) {
		for_Int(i,0,p.nband){
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
	ImGreen hb(p.nband, p);
	for_Int(i, 0, p.nband) {
		for_Int(n, 0, g_loc.nomgs) {
			hb[n][i][i] = -bethe_t[i] * bethe_t[i] * g_loc[n][i][i];
			// hb[n][i][i] -= real(hb[n][i][i]);
		}
	}
	return hb;
}

ImGreen DMFT::find_hb_by_se(const ImGreen& se_i) const
{
	ImGreen hb(p.nband, p);
	// for_Int(i,0,p.nband){
	// 	for_Int(n, 0, g_loc.nomgs){
	// 		Cmplx z = g_loc.z(n) + bethe_mu - se_i[n][i][i];
	// 		hb[n][i][i]= -z + g_loc.inverse()[n][i][i];
	// 	}
	// }
	for_Int(n, 0, g_loc.nomgs){
		MatCmplx z(p.nband, p.nband, 0.);
		for_Int(i, 0, p.nband) { z[i][i] = g_loc.z(n) + bethe_mu - se_i[n][i][i]; }
		hb[n] = g_loc.inverse()[n] - z;
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
	ImGreen se(p.nband, p);
	// ImGreen g_loc_inv = g_loc.inverse();
	for_Int(n, 0, se.nomgs) {
		for_Int(i, 0, p.nband){
			se[n][i][i] = g_imp.z(n) + bethe_mu + g_imp.inverse()[n][i][i] - hb_imp[n][i][i];
		}
	}
	return se;
}

// ------------------------------- since the 2023.04.15 -------------------------------

// bethe lattice part:
ImGreen DMFT::g0_loc() const {
	ImGreen g0_loc(p.nband, p);
	for_Int(n, 0, p.num_omg) {
		for_Int(i, 0, p.nband) {
			g0_loc[n][i][i] = (g0_loc.z(n) + bethe_mu - SQRT(SQR(g0_loc.z(n) + bethe_mu) - 4 * SQR(bethe_t[i]))) / (2 * bethe_t[i] * bethe_t[i]);
		}
	}
	return g0_loc;
}

void DMFT::auto_nooc(Str mode, const Impurity& imp) {
	if(mode == "ful_pcl_sch"){
		Occler opcler(mm, p);
		VEC<MatReal> uormat;
		VecInt ordeg = VecInt{ 1, 1, 2, 2 }, nppso;
		Vec<VecInt> controler(MAX(ordeg) + 1, VecInt(p.ndiv, 0));
		MatReal occnum, occweight;
		controler[0] = p.control_divs[0];
		{
			NORG norg(opcler.find_ground_state_partical(imp.impH, VecInt{1, 2}));
			uormat = norg.uormat;
			occnum = norg.occnum.mat(p.norg_sets, p.n_rot_orb / p.norg_sets);occweight = occnum;
			nppso = norg.scsp.nppso;
		}
		for_Int(i, 0, p.norg_sets) for_Int(j, 0, p.n_rot_orb/p.norg_sets) occweight[i][j] = occnum[i][j] > 0.5 ? (1 - occnum[i][j]) : occnum[i][j];

		for_Int(i, 0, MAX(VecInt{1, 2})){
			Int o(0), freze_o(0), e(0), freze_e(0), orb_rep(0), nooc_o(0), nooc_e(0);
			for_Int(j, 0, p.norg_sets) {orb_rep = j; if(ordeg[j] == i + 1) break;}
			o = nppso[orb_rep] - 1; e = p.nI2B[orb_rep] - nppso[orb_rep];
			for_Int(j, 0, o) 							if(occweight[orb_rep][j] < 1e-7) freze_o++;
			for_Int(j, nppso[orb_rep], p.nI2B[orb_rep])	if(occweight[orb_rep][j] < 1e-7) freze_e++;
			nooc_o = o - freze_o; nooc_e = e - freze_e;
			controler[i+1] = p.if_norg_imp ?  VecInt{freze_o, nooc_o, 1, 1, nooc_e, freze_e } : VecInt{1, freze_o, nooc_o, 1, nooc_e, freze_e };
		}

		// p.if_norg_imp = true; p.after_modify_prmtr(); 
		p.according_controler(controler, ordeg);
	}
}