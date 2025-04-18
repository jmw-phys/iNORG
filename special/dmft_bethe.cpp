/*
coded by Jia-Ming Wang (jmw@ruc.edu.cn, RUC, China) date 2022
*/
//! not finished yet 


#include "dmft_bethe.h"


void DMFT::set_parameter() {
	bethe_mu = p.mu;
	bethe_t = p.bethe_t;
	bethe_u = p.U;
	// p.jz = 0.25 * p.U;
	p.Uprm = bethe_u12 = p.U - p.delta;
	p.bandw = 2 * (2 * SQRT(SQR(bethe_u) + SQR(bethe_u12) + SQR(bethe_u12) + SUM(bethe_t * bethe_t)));
	p.derive_ImGreen();
	if (mm) { OFS ofss("log.parameters.txt", std::ios::app);  p.print(ofss); }
}

DMFT::DMFT(const MyMpi& mm_i, Prmtr& prmtr_i, const Int mode) :
	mm(mm_i), p(prmtr_i), iter_cnt(0), g_loc(p.nband, p), n_eles(p.norbs, 0.),
	gloc_err({ 1.e99,1.e98,1.e97 ,1.e96 ,1.e95}), se(p.nband, p), imp_backup(p.imp_backup)
{
	// make random seed output together
	{ mm.barrier(); SLEEP(1); }
	IFS fitdata("ose_hop"), mbgfdata("mb.gfloc"), mbsedata("mb.seimp");
	log("initial");	set_parameter();
	Bath bth(mm, p); 	Impurity imp(mm, p, bth);
	NORG norg(mm, p);
	if (fitdata) bth.read_ose_hop();
	if (mbsedata) { se	  = ImGreen(p.nband, p, "mb.seimp");								if (mm) se.write("seimp", iter_cnt); }
	if (mbgfdata) { g_loc = ImGreen(p.nband, p, "mb.gfloc");								if (mm) g_loc.write("gfloc", iter_cnt);}
	else {g_loc = g0_loc();																	if (mm) g_loc.write("g0loc", iter_cnt);}

    VEC<MatReal> norg_tempU;
	se_input.push_back(se);
	Int Flag_semix(0);
	// while (iter_cnt < p.iter_max && !converged()) 
	while (iter_cnt < p.iter_max) 
	{
		++iter_cnt;ImGreen hb(p.nband, p);
		if (!(fitdata && iter_cnt == 1)) {
			hb = (mode == 1) ? find_hb_by_se(se):find_hb(g_loc);							if (mm) hb.write("hb", iter_cnt);
			bth.bath_fit(hb,iter_cnt);														if (mm) bth.write_ose_hop(iter_cnt);
		}
		imp.update("behte_alpha");																if (mm) imp.write_H0info(bth, -1, iter_cnt);
		ImGreen hb_imp(p.nband, p);   	imp.find_hb(hb_imp); 								if (mm) hb_imp.write("hb-fit", iter_cnt);
		// auto_nooc("ful_pcl_sch", imp);	NORG norg(mm, p);
		if (iter_cnt > 1) norg.uormat = norg_tempU;	norg.up_date_h0_to_solve(imp.impH, 1);	n_eles = norg.write_impurtiy_occupation(iter_cnt);
		ImGreen g0imp(p.nband, p);	imp.find_g0(g0imp);										if (mm)	g0imp.write("g0imp", iter_cnt);
		ImGreen gfimp(p.nband, p);	norg.get_gimp_eigpairs(gfimp);							if (mm) gfimp.write("gfimp", iter_cnt);
		ImGreen seimp(p.nband, p);	seimp = g0imp.inverse() - gfimp.inverse();				if (mm) seimp.write("seimp", iter_cnt);

		if (mode == 0) {
			append_gloc_err(gfimp.error(g_loc));											log("gerr_update");
			// g_loc = 0.3 * g_loc + 0.7 * gfimp;
			g_loc = gfimp;
		}

		// obtain se from impurity model.
		if (mode == 1) {
			Real err_temp(se.error(seimp));
			append_gloc_err(err_temp);														log("sigerr_update");
			// se = (iter_cnt == 1) ? seimp : 0.5 * se + 0.5 * seimp;	
			// se = seimp;
			if (iter_cnt > 1 && err_temp < 0.01) Flag_semix = 1;
			pulay_mixing(seimp);
			if (!Flag_semix) {
				se = seimp;
				res_past.clear();
				se_input.clear();
				se_input.push_back(se);
			}
			g_loc = find_gloc_by_se(se);													if (mm) g_loc.write("gfloc", iter_cnt);
		}
		norg_tempU = norg.uormat;

		if (converged()) {
			Real& var_a(p.U);
			if(mm) {
				g0imp.write("U" + STR(var_a) + "mb.g0imp");
				gfimp.write("U" + STR(var_a) + "mb.gfimp");
				seimp.write("U" + STR(var_a) + "mb.seimp");
				g_loc.write("U" + STR(var_a) + "mb.gfloc");
				bth.write_ose_hop(-1, "U" + STR(var_a));
			}
			norg.write_impurtiy_occupation(-1, "U" + STR(var_a));
			ReGreen g0_imp_re(p.nband, p);	imp.find_g0(g0_imp_re);							if (mm) g0_imp_re.write("U" + STR(var_a) + "Re-g0fimp");
			ReGreen gfimp_re(p.nband, p);	norg.get_gimp_eigpairs(gfimp_re);				if (mm) gfimp_re.write("U" + STR(var_a) + "Re-gfimp");
			ReGreen se_re = g0_imp_re.inverse() - gfimp_re.inverse();						if (mm) se_re.write("U" + STR(var_a) + "Re-seimp");
			ReGreen g_loc_re(p.nband, p); g_loc_re = find_gloc_by_se(se_re);				if (mm) g_loc_re.write("U" + STR(var_a) + "Re-gfloc");


			// excitation spectrum
			// ReGreen hd_exsp(p.nband, p);	norg.get_gimp_hdQPs(hd_exsp);					if (mm)	hd_exsp.write("U" + STR(var_a) + "Re-hdex");
			
			var_a += 0.5;
			set_parameter();
			if (var_a > 4.0)
				break;
			res_past.clear();
			se_input.clear();
			Flag_semix = iter_cnt = 0;
			se_input.push_back(se);
			gloc_err = { 1.e99,1.e98,1.e97 ,1.e96 ,1.e95};
			if(mm) norg.scsp.print();
		}

		{ mm.barrier(); SLEEP(1); }
		if(IFS("norg.stop"))  break;
	}

/*
	// ! auto_nooc("for_green_calc", imp); // in here we can change the nooc space once, last chance.
	NORG finalrg(mm, p); 		finalrg.uormat = norg_tempU;								finalrg.up_date_h0_to_solve(imp.impH, 1);
	// ImGreen gfimp(p.nband, p);	finalrg.get_gimp_eigpairs(gfimp);							if (mm) gfimp.write("gfimp", iter_cnt);
	ReGreen g0_imp_re(p.nband, p);imp.find_g0(g0_imp_re);									if (mm) g0_imp_re.write("g0fimp");
	ReGreen gfimp_re(p.nband, p);	finalrg.get_gimp_eigpairs(gfimp_re);					if (mm) gfimp_re.write("gfimp");
	ReGreen se_re = g0_imp_re.inverse() - gfimp_re.inverse();								if (mm) se_re.write("se_loc");

	ReGreen g_loc_re(p.nband, p); g_loc_re = find_gloc_by_se(se_re);						if (mm) g_loc_re.write("g_loc_re");
*/
	
	// if (mm) bth.write_ose_hop();
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

ReGreen DMFT::find_gloc_by_se(const ReGreen& se_i) const
{
	ReGreen gloc_temp(g_loc.norbs, p);
	for_Int(n, 0, gloc_temp.nomgs) {
		for_Int(i,0,p.nband){
			Cmplx z = gloc_temp.z(n) + bethe_mu - se_i[n][i][i];
			Cmplx temp = SQRT(SQR(z) - 4 * SQR(bethe_t[i]));
			if(imag(temp) < 0) temp *= -1;
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
		VecInt ordeg(p.norbs, 0), nppso;
		for_Int(i, 0, p.nband) for_Int(j, 0, 2) ordeg[i * 2 + j] = i + 1;
		Vec<VecInt> controler(MAX(ordeg) + 1, VecInt(p.ndiv, 0));
		MatReal occnum, occweight;
		controler[0] = p.control_divs[0];
		{
			NORG norg(opcler.find_ground_state_partical(imp.impH, VecInt{1,1,2,2}));
			uormat = norg.uormat;
			occnum = norg.occnum.mat(p.norg_sets, p.n_rot_orb / p.norg_sets);occweight = occnum;
			nppso = norg.scsp.nppso;
			p.npartical = norg.scsp.nppso;
		}
		for_Int(i, 0, p.norg_sets) for_Int(j, 0, p.n_rot_orb/p.norg_sets) occweight[i][j] = occnum[i][j] > 0.5 ? (1 - occnum[i][j]) : occnum[i][j];

		for_Int(i, 0, MAX(ordeg)){
			Int o(0), freze_o(0), e(0), freze_e(0), orb_rep(0), nooc_o(0), nooc_e(0);
			for_Int(j, 0, p.norg_sets) {orb_rep = j; if(ordeg[j] == i + 1) break;}
			o = nppso[orb_rep] - 1; e = p.nI2B[orb_rep] - nppso[orb_rep];
			for_Int(j, 0, o) 							if(occweight[orb_rep][j] < 1e-8) freze_o++;
			for_Int(j, nppso[orb_rep], p.nI2B[orb_rep])	if(occweight[orb_rep][j] < 1e-8) freze_e++;
			nooc_o = o - freze_o; nooc_e = e - freze_e;
			controler[i+1] = p.if_norg_imp ?  VecInt{freze_o, nooc_o, 1, 1, nooc_e, freze_e } : VecInt{1, freze_o, nooc_o, 1, nooc_e, freze_e };
		}
		// if(mm) WRN(NAV(controler));
		// p.if_norg_imp = true; p.after_modify_prmtr(); 
		p.according_controler(controler, ordeg);
	}
}

void DMFT::pulay_mixing(const ImGreen& seimp){
    Int len_std = 5;
    ImGreen res = seimp - se;
    res_past.push_back(res);
    if (res_past.size() > len_std) {
        se_input.erase(se_input.begin());
        res_past.erase(res_past.begin());
    }
    Int len = res_past.size();

    VecReal beta(len);
    for_Int(i,0,len){
        // beta[len - 1 - i] = exp(-(iter_cnt - i - 1));
		beta[len - 1 - i] = 1;
    }

    MatReal temp(len + 1, len + 1, 0.);
    MatReal alpha(len + 1, 1, 0.);
    for_Int(i, 0, len) {
        for_Int(j, i, len){
            temp[i][j] = dot(res_past[i], res_past[j]);
            temp[j][i] = temp[i][j];
        }
    }
    for_Int(i, 0, len) {
        temp[len][i] = 1.;
        temp[i][len] = 1.;
    }
    alpha[len] = 1.;
    if (mm) WRN(NAV(temp));
    Int info = gaussj(temp, alpha);
	if(mm)	WRN(NAV(alpha));
    se -= se;
    for_Int(i, 0, len) {
        se += alpha[i][0] * se_input[i] + alpha[i][0] * beta[i] * res_past[i];
    }
    se_input.push_back(se);
}

// void DMFT::get_gimpim_evenodd(NORG& norg,Green& gfimp) const {
void DMFT::get_gimp_evenodd(NORG& norg,Green& gfimp) const {

	if(gfimp.type_info() == STR("ImGreen"))
	{
		ImGreen gfimp1(p.norbs, p);	norg.get_gimp_all(gfimp1);
		if (mm)	gfimp1.write("gfimp1", iter_cnt);
		ImGreen gfimp2(p.norbs, p); gfimp2 = gfimp1;
		for_Int(i, 0, p.norbs) {
			for_Int(n, 0, gfimp1.nomgs) {
				if (i < p.nband) {
					gfimp1[n][i][i] += gfimp2[n][i + p.nband][i + p.nband];
				}
				else {
					gfimp1[n][i][i] += gfimp2[n][i - p.nband][i - p.nband];
				}
				gfimp1[n][i][i] /= 2.;
			}
		}
		if (mm)	gfimp1.write("gfimp11", iter_cnt);
		for_Int(i, 0, p.nband) {
			for_Int(n, 0, gfimp.nomgs) {
				gfimp[n][i][i] = gfimp1[n][i * 2][i * 2];
			}
		}
	}

	if(gfimp.type_info() == STR("ReGreen"))
	{
		ReGreen gfimp1(p.norbs, p);	norg.get_gimp_all(gfimp1);
		//if (mm)	gfimp1.write("gfimp1", iter_cnt);
		ReGreen gfimp2(p.norbs, p); gfimp2 = gfimp1;
		for_Int(i, 0, p.norbs) {
			for_Int(n, 0, gfimp1.nomgs) {
				if (i < p.nband) {
					gfimp1[n][i][i] += gfimp2[n][i + p.nband][i + p.nband];
				}
				else {
					gfimp1[n][i][i] += gfimp2[n][i - p.nband][i - p.nband];
				}
				gfimp1[n][i][i] /= 2.;
			}
		}
		//if (mm)	gfimp1.write("gfimp11", iter_cnt);
		for_Int(i, 0, p.nband) {
			for_Int(n, 0, gfimp.nomgs) {
				gfimp[n][i][i] = gfimp1[n][i * 2][i * 2];
			}
		}
	}
}

// -------------------------------------------------- private --------------------------------------------------

	bool DMFT::converged() const {
		const Real dev = DEV(gloc_err);
		if(gloc_err[gloc_err.size()-2]<1.E-4 && gloc_err[gloc_err.size()-2]<gloc_err[gloc_err.size()-1]) return true;
		if(gloc_err[gloc_err.size()-1]<1.E-5) return true;
		if (dev > 1.E-4) { return false; }
		else if (dev > 1.E-10) {
			for_Int(i, 1, gloc_err.size()) {
				if (gloc_err[0] > gloc_err[i]) { return false; }
			}
			return true;
		}
		else{ return true; }
	}
