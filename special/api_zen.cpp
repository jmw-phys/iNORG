/*
coded by Jia-Ming Wang (jmw@ruc.edu.cn, RUC, China) date 2022 - 2023
*/

#include "api_zen.h"

APIzen::APIzen(const MyMpi& mm_i, Prmtr& prmtr_i, const Str& file) :
	mm(mm_i), p(prmtr_i),num_omg(prmtr_i.num_omg),
	num_nondegenerate(-1), dmft_cnt(0)
{
	update(file);
	Bath bth(mm, p);
	Impurity imp(mm, p, bth, or_deg_idx);
	// NORG norg(mm, prmtr_i);		// for test, not search the whole particle space.
	// OFS ofs("norg.lock");
	

	// while(true) {
	// 	if(IFS("norg.stop"))  break;
	// 	if(if_lock("norg")) {
	// if(mode == "table_mode") while(true) if(if_lock("norg")) {

	dmft_cnt++; update(file);
	ImGreen hb(nband, p);	
	for_Int(j, 0, hb.nomgs) for_Int(i, 0, nband)	hb.g[j][i][i] = -imfrq_hybrid_function[i][j];	if (mm) hb.write_zen("hb_zen", "Read");
	bth.read_ose_hop();	bth.bath_fit(hb, or_deg_idx);												if (mm) bth.write_ose_hop();
	imp.update();																					if (mm) imp.write_H0info(bth, MAX(or_deg_idx));
	ImGreen hb_imp(p.nband, p);		imp.find_hb(hb_imp); 											if (mm) hb_imp.write_zen("hb_imp", "Fit");
	auto_nooc("ful_pcl_sch", imp);
	NORG norg(mm, p);
	if (!norg.check_NTR()) norg.uormat = p.rotationU;
	else MatReal tmp_b = norg.read_NTR();
	norg.up_date_h0_to_solve(imp.impH, 1);															norg.write_impurtiy_occupation();
	MatReal tmp_e = norg.save_NTR();
	// MatReal local_multiplets_state = norg.oneedm.local_multiplets_state(norg.oneedm.ground_state);	if (mm)WRN(NAV(local_multiplets_state));
	ImGreen g0imp(p.nband, p);	imp.find_g0(g0imp);													if (mm)	g0imp.write_zen("g0imp");
	ImGreen gfimp(p.nband, p);	norg.get_gimp_eigpairs(gfimp, or_deg_idx);							if (mm) gfimp.write_zen("gfimp");
	ImGreen seimp(p.nband, p);	seimp = g0imp.inverse() - gfimp.inverse();	seimp_fixer(seimp);		if (mm) seimp.write_zen("seimp");

	

	// 	if(mm)	std::remove("norg.lock");
	// 	{ mm.barrier(); SLEEP(1); }
	// 	}
	// }
/*
	// if(mode == "realf_mode"){
		ReGreen g0imp_re(p.nband, p);	imp.find_g0(g0imp_re);													if (mm)	g0imp_re.write_zen("Re_g0imp");
		ReGreen gfimp_re(p.nband, p);	norg.get_gimp_eigpairs(gfimp_re, or_deg_idx);							if (mm) gfimp_re.write_zen("Re_gfimp");
		ReGreen seimp_re(p.nband, p);	seimp_re = g0imp_re.inverse() - gfimp_re.inverse();						if (mm) seimp_re.write_zen("Re_seimp");
		ReGreen hd_exsp(p.nband, p);	norg.get_gimp_hdQPs(hd_exsp);											if (mm)	 hd_exsp.write_zen("Re-hdex");
	// }
*/

 /*
	{
		ImGreen g_asnci(p.nband, p);
		VecInt or_deg = or_deg_idx.truncate(0,nband);
		VecInt idx(MAX(or_deg),0); Int cter(0);
		for_Int(i, 0, or_deg.size()) if(cter < or_deg[i]) idx[cter++] = i; 
		for (int &i : idx) {
			StdVecInt difference = {(i+1), -(i+1)};
			for(const auto ii: difference)	{
				// Asnci nci(norg, Int(619), ii);
				Asnci nci(norg, Int(100), ii);
				if(mm) WRN("Finish finding the space")
				// nci.asnci_gimp(g_asnci, ii);
			}
		}
		for_Int(i, 0, or_deg.size()) for_Int(n, 0, g_asnci.nomgs) g_asnci[n][i][i] = g_asnci[n][idx[or_deg[i] - 1]][idx[or_deg[i] - 1]];
		if (mm) g_asnci.write_zen("g_asnci");
	}
*/

	// if(mm) gfimp.write_occupation_info();
	// if(mm) WRN(NAV(gfimp.particle_number().diagonal()));
	// ImGreen seimp(p.nband, p);	seimp = g0imp.inverse() - gfimp.inverse();	
	// if (mm) seimp.write_zen("before_fix_seimp");						seimp = fix_se(seimp);
	// if (mm) seimp.write_zen("seimp");
}

void APIzen::test_for_fitting(const Bath& bth, const ImGreen& hby_i, Int num)
{ // test for fitting(fit_hb)
	ImGreen fit_hb(1, p);
	const VecCmplx E = cmplx(bth.ose);
	const VecCmplx V = cmplx(bth.hop);
	VecCmplx Z(p.nbath);
	for_Int(n, 0, fit_hb.nomgs)
	{
		const Cmplx z = Cmplx(0., hby_i.omg(n));
		Z = z;
		VecCmplx S = INV(E - Z);
		Cmplx hyb = SUM(V * S * V.co());
		fit_hb[n][0][0] = hyb;
	}
	fit_hb.write("hb_fitting", num + 1);
}

void APIzen::read_ZEN(const Str& file)
{
	{// norg.in
		// std::vector<double> Ed;
		// std::vector<int> Deg;

		std::string CoulombF;
		double U;
		double J;
		std::vector<int> restrain_t;
		std::vector<int> distribute_t;
		read_norg_setting("solver.norg.in", CoulombF, U, J,  restrain_t, distribute_t);
		Uc = U;	Jz = J; restrain = VecInt(restrain_t); distribute = VecInt(distribute_t);
	}

	imfrq_hybrid_function.reset(norbs,num_omg,0.);
	solver_eimp_data.reset(norbs,0.);
	or_deg_idx.reset(nband, 0);

	{// hyb.in
		Str hybdata(file + ".hyb.in");
		IFS ifs(hybdata);
		if (!ifs) {
			ERR(STR("file opening failed with ") + NAV(hybdata))
		}
		else {
			Int i(0);
			bool swicher(true);
			while (true) {
				VecReal omg(num_omg, 0.), re(num_omg, 0.), im(num_omg, 0.);
				Int drop_Int(0);
				Real drop_omg(0.), drop_re(0.), drop_im(0.), drop_err1(0.), drop_err2(0.);
				if (i >= norbs) break;
				for_Int(j, 0, num_omg) {
					if (swicher) {
						ifs >> drop_Int;
					}
					ifs >> omg[j]; ifs >> re[j]; ifs >> im[j];
					ifs >> drop_err1; ifs >> drop_err2;
					//if(i==1)DBG(NAV2(re[j], im[j]));
					if (!ifs) ERR(STR("read_ZEN-in error with ") + NAV(hybdata));
					swicher = true;
				}
				imfrq_hybrid_function[i] = cmplx(re, im);
				while (true) {
					ifs >> drop_Int;
					if (drop_Int - 1 == i + 1) {
						swicher = false;
						break;
					}
					ifs >> drop_omg; ifs >> drop_re; ifs >> drop_im; ifs >> drop_err1; ifs >> drop_err2;
					if (!ifs) break;
				}
				++i;
			}
		}
		ifs.close();
	}

	{// eimp.in
		Str eimpdata(file + ".eimp.in");
		p.eimp.reset(norbs);
		IFS ifs(eimpdata);
		if (!ifs) {
			ERR(STR("file opening failed with ") + NAV(eimpdata))
		}
		else {
			Int drop_Int(0), c(0);
			VecReal temp_emps(nband, 0.);
			for_Int(i, 0, nband) {
				ifs >> drop_Int;
				ifs >> temp_emps[i]; ifs >> or_deg_idx[i];
				if (!ifs) ERR(STR("read_ZEN-in error with ") + NAV(eimpdata));
			}
			// for_Int(j, 0, nband) for_Int(i, 0, 2) p.eimp[c++] = temp_emps[j+nband*i];
			or_deg_idx.reset(concat(or_deg_idx,or_deg_idx).mat(2,or_deg_idx.size()).tr().vec());
			p.eimp.reset(concat(temp_emps,temp_emps).mat(2,temp_emps.size()).tr().vec());
			// if (test_mode) num_nondegenerate = 1;
			num_nondegenerate = MAX(or_deg_idx);
		}
		if (num_nondegenerate <= 0)ERR(STR("read_ZEN-in error with ") + NAV2(eimpdata, num_nondegenerate));
		ifs.close();
	}


	if (mm) WRN(NAV5(restrain, distribute, Uc, Jz, p.beta))

	if (mm) WRN(NAV2(p.eimp, or_deg_idx))
}



// NORG APIzen::choose_cauculation_style(Str mode, Impurity &imp){
// 	if(mode == "ful_pcl_sch"){
// 		Occler opcler(mm, p);
// 		VEC<MatReal> uormat;
// 		VecInt ordeg(p.norbs, 0), nppso;
// 		for_Int(i, 0, p.nband) for_Int(j, 0, 2) ordeg[i * 2 + j] = i + 1;
// 		Vec<VecInt> controler(MAX(ordeg) + 1, VecInt(p.ndiv, 0));
// 		MatReal occnum, occweight;
// 		controler[0] = p.control_divs[0];
// 		{
// 			NORG norg(opcler.find_ground_state_partical(imp.impH, VecInt{1,1,2,2}));
// 			uormat = norg.uormat;
// 			occnum = norg.occnum.mat(p.norg_sets, p.n_rot_orb / p.norg_sets);occweight = occnum;
// 			nppso = norg.scsp.nppso;
// 		}
// 		for_Int(i, 0, p.norg_sets) for_Int(j, 0, p.n_rot_orb/p.norg_sets) occweight[i][j] = occnum[i][j] > 0.5 ? (1 - occnum[i][j]) : occnum[i][j];

// 		for_Int(i, 0, MAX(ordeg)){
// 			Int o(0), freze_o(0), e(0), freze_e(0), orb_rep(0), nooc_o(0), nooc_e(0);
// 			for_Int(j, 0, p.norg_sets) {orb_rep = j; if(ordeg[j] == i + 1) break;}
// 			o = nppso[orb_rep] - 1; e = p.nI2B[orb_rep] - nppso[orb_rep];
// 			for_Int(j, 0, o) 							if(occweight[orb_rep][j] < 1e-7) freze_o++;
// 			for_Int(j, nppso[orb_rep], p.nI2B[orb_rep])	if(occweight[orb_rep][j] < 1e-7) freze_e++;
// 			nooc_o = o - freze_o; nooc_e = e - freze_e;
// 			controler[i+1] = p.if_norg_imp ?  VecInt{freze_o, nooc_o, 1, 1, nooc_e, freze_e } : VecInt{1, freze_o, nooc_o, 1, nooc_e, freze_e };
// 		}
// 		// if(mm) WRN(NAV(controler));
// 		p.according_controler(controler, ordeg);
// 		// {// WRN
// 		// 	MatInt m_controler(MAX(or_deg_idx) + 1, p.ndiv);
// 		// 	for_Int(i, 0, controler.size()) m_controler[i] = controler[i];
// 		// 	if(mm) WRN(NAV3(ordeg,m_controler, p.control_divs));
// 		// }
// 		NORG frezeorb(mm, p);
// 		// IFS ifs_a("ru" + frezeorb.scsp.nppso_str() + ".bi"); frezeorb.uormat = uormat;
// 		// if (ifs_a) for_Int(i, 0, frezeorb.uormat.size()) biread(ifs_a, CharP(frezeorb.uormat[i].p()), frezeorb.uormat[i].szof());
// 		frezeorb.up_date_h0_to_solve(imp.impH, 1);
// 		if (mm)	{
// 			// OFS ofs_a;
// 			// ofs_a.open("ru" + frezeorb.scsp.nppso_str() + ".bi");
// 			// for_Int(i, 0, frezeorb.uormat.size()) biwrite(ofs_a, CharP(frezeorb.uormat[i].p()), frezeorb.uormat[i].szof());
// 		}
// 		return frezeorb;
// 	}
// 	if(mode == "one_pcl_test"){
// 		Occler opcler(mm,p);
// 		VEC<MatReal> uormat;
// 		VecInt ordeg(concat(or_deg_idx.truncate(0,nband),or_deg_idx.truncate(0,nband)).mat(2,nband).tr().vec()), nppso;
// 		Vec<VecInt> controler(MAX(or_deg_idx) + 1, VecInt(p.ndiv, 0));
// 		MatReal occnum, occweight;
// 		controler[0] = {0, -1, p.control_divs[0][2], 0, p.control_divs[0][4], 1};
// 		{
// 			Int band1(p.npartical[0]), band2(p.npartical[0]);
// 			if (nband == 5) p.npartical = { band1, band1, band1, band1, band2, band2, band1, band1, band2, band2 };
// 			if (nband == 3) { p.npartical = { band1, band1, band1, band1, band1, band1 }; /*p.npartical += 1;*/ }
// 			p.according_nppso(p.npartical);
// 			NORG norg(mm, p);
// 			IFS ifs_a("ru" + norg.scsp.nppso_str() + ".bi");
// 			if (ifs_a) for_Int(i, 0, norg.uormat.size()) biread(ifs_a, CharP(norg.uormat[i].p()), norg.uormat[i].szof());
// 			norg.up_date_h0_to_solve(imp.impH, 1);

// 			uormat = norg.uormat;
// 			occnum = norg.occnum.mat(p.norg_sets, p.n_rot_orb / p.norg_sets); occweight = occnum;
// 			nppso = norg.scsp.nppso;
// 			// if(mm) norg.write_state_info(0);
// 		}
// 		for_Int(i, 0, p.norg_sets) for_Int(j, 0, p.n_rot_orb/p.norg_sets) if(occnum[i][j] > 0.5) occweight[i][j] = 1 - occnum[i][j];

// 		for_Int(i, 0, MAX(or_deg_idx)){ //? may not suit for the "if_norg_imp = true" case.
// 			Int o(0), freze_o(0), e(0), freze_e(0), orb_rep(0), nooc_o(0), nooc_e(0);
// 			for_Int(j, 0, p.norg_sets) {orb_rep = j; if(ordeg[j] == i + 1) break;}
// 			o = nppso[orb_rep] - 1; e = p.nI2B[orb_rep] - nppso[orb_rep];
// 			for_Int(j, 0, o) 							if(occweight[orb_rep][j] < 1e-6) freze_o++;
// 			for_Int(j, nppso[orb_rep], p.nI2B[orb_rep])	if(occweight[orb_rep][j] < 1e-6) freze_e++;
// 			nooc_o = o - freze_o; nooc_e = e - freze_e;
// 			controler[i+1] = VecInt{1, freze_o, nooc_o, 1, nooc_e, freze_e };
// 		}

// 		p.according_controler(controler, ordeg);
// 		// {// WRN
// 		// 	MatInt m_controler(MAX(or_deg_idx) + 1, p.ndiv);
// 		// 	for_Int(i, 0, controler.size()) m_controler[i] = controler[i];
// 		// 	if(mm) WRN(NAV3(ordeg,m_controler, p.control_divs));
// 		// }
// 		NORG frezeorb(mm, p); frezeorb.uormat = uormat;
// 		IFS ifs_a("ru" + frezeorb.scsp.nppso_str() + ".bi");
// 		if (ifs_a) for_Int(i, 0, frezeorb.uormat.size()) biread(ifs_a, CharP(frezeorb.uormat[i].p()), frezeorb.uormat[i].szof());
// 		frezeorb.up_date_h0_to_solve(imp.impH, 1);
// 		if (mm)	{
// 			OFS ofs_a;
// 			ofs_a.open("ru" + frezeorb.scsp.nppso_str() + ".bi");
// 			for_Int(i, 0, frezeorb.uormat.size()) biwrite(ofs_a, CharP(frezeorb.uormat[i].p()), frezeorb.uormat[i].szof());
// 			// frezeorb.write_state_info(1);
// 		}
// 		return frezeorb;
// 	}
// }


/*
ImGreen APIzen::fix_se(const ImGreen& se) const{
  using namespace std;
    ImGreen se_fix(se);
    MatReal mat_nfit(se.norbs, se.norbs);
    MatReal mat_slope(se.norbs, se.norbs);
    // MatReal mat_term_2(se.norbs, se.norbs);
    // MatReal mat_term_0(se.norbs, se.norbs);

	//A*e^(\alpha(iw))
    MatReal mat_A(se.norbs, se.norbs);
    MatReal mat_alpha(se.norbs, se.norbs);
    for_Int(row, 0, se.norbs) {
        // for_Int(col, 0, se.norbs) {
			Int col = row; Int abandon_n = se.nomgs-1;
            while (abandon_n > 0) {
                abandon_n -= 1;
                if(imag(se[abandon_n][row][col]) > 0.) break;
            }
            
            Int n_fit = 0; Real max_slope(0.);
            //if(abandon_n < p.fit_max_omg){
                for_Int(n, abandon_n, p.fit_max_omg) {
                    Real slope = imag(se[n][row][col]) / se.omg(n);
                    if(slope < max_slope){
                        max_slope = slope;
                        n_fit = n;
                    }
                }
            //}
            // Real term_2 = (real(se[n_fit + 1][row][col]) - real(se[n_fit][row][col])) / (2 * se.omg(n_fit) * 2 * se.unit_omg);
            // Real term_0 = real(se[n_fit][row][col]) - term_2 * SQR(se.omg(n_fit));
            // for_Int(n, 0, n_fit) {
            //     se_fix[n][row][col] = term_0 + term_2 * SQR(se.omg(n)) + I * max_slope * se.omg(n);
            // }

            Real term_alpha = atan(imag(se[n_fit][row][col])/real(se[n_fit][row][col]))/se.omg(n_fit);
            Real term_A = imag(se[n_fit][row][col])/(sin(term_alpha*se.omg(n_fit)));
            for_Int(n, 0, n_fit) se_fix[n][row][col] = term_A * std::exp(term_alpha * se.z_omg[n_fit]);

            mat_nfit[row][col] = n_fit;
            mat_slope[row][col] = max_slope;
            mat_A[row][col] = term_A;
            mat_alpha[row][col] = term_alpha;
            // mat_term_2[row][col] = term_2;
            // mat_term_0[row][col] = term_0;
        // }
    }

	if(mm) PIO(NAVC4(mat_nfit.diagonal().mat(1,p.nband),mat_slope.diagonal().mat(1,p.nband),mat_A.diagonal().mat(1,p.nband),mat_alpha.diagonal().mat(1,p.nband)))
	
    // if(mm){
    //     OFS ofs(iox + "zic" + prefill0(iter_cnt, 3) + ".mat_fix" + ".log" + ".txt");
    //     ofs<< "mat_nfit" << "\n" << mat_nfit << endl;
    //   	ofs << iofmt("sci");
    //     ofs<< "mat_slope" << "\n" << mat_slope << endl;
    //     ofs<< "mat_term_2" << "\n" << mat_term_2 << endl;
    //     ofs<< "mat_term_0" << "\n" << mat_term_0 << endl;
    //     ofs.close();
    // }
    return se_fix;
}
*/


//---------------------------------------------- private ----------------------------------------------

// to up date the whole date for the impurity system
void APIzen::update(const Str& file) {
	{// modify the parameters from ZEN.in
		read_ZEN(file);	p.U = Uc; p.mu = mu; p.jz = Jz; p.nband = nband; p.norg_sets = p.norbs = norbs;
		p.templet_restrain = restrain; p.templet_control = distribute; p.project = NAV(nband) + "KVSb";
		p.after_modify_prmtr(); p.recalc_partical_number();
		p.Uprm = p.U - 2 * p.jz;
		p.degel = 0;
		// if (mm) p.print();
	}

	// if (mm) std::cout << std::endl;						// blank line
}

bool APIzen::if_lock(const Str file) const {
	Str lock_file = file + ".lock";

	IFS ifs(lock_file);
	{ mm.barrier(); SLEEP(1); }
	if(ifs) {
		return true;
	} else {
		{ mm.barrier(); SLEEP(5000); }
		return false;
	}
}


void APIzen::auto_nooc(Str mode, const Impurity& imp) {
	if(mode == "ful_pcl_sch"){
		Occler opcler(mm, p);
		VEC<MatReal> uormat;
		VecInt ordeg(p.norbs, 0), nppso;
		for_Int(i, 0, p.nband) for_Int(j, 0, 2) ordeg[i * 2 + j] = i + 1;
		Vec<VecInt> controler(MAX(ordeg) + 1, VecInt(p.ndiv, 0));
		MatReal occnum, occweight;
		controler[0] = p.control_divs[0];
		{
			NORG norg(opcler.find_ground_state_partical(imp.impH, or_deg_idx));
			uormat = norg.uormat;
			occnum = norg.occnum.mat(p.norg_sets, p.n_rot_orb / p.norg_sets);occweight = occnum;
			nppso = norg.scsp.nppso;
			p.npartical = norg.scsp.nppso;
			p.rotationU = uormat;
			norg.PIO_occweight(norg.occnum);
		}
		for_Int(i, 0, p.norg_sets) for_Int(j, 0, p.n_rot_orb/p.norg_sets) occweight[i][j] = occnum[i][j] > 0.5 ? (1 - occnum[i][j]) : occnum[i][j];

		for_Int(i, 0, MAX(ordeg)){
			Int o(0), freze_o(0), e(0), freze_e(0), orb_rep(0), nooc_o(0), nooc_e(0);
			Int keep_o(0), keep_e(0);
			for_Int(j, 0, p.norg_sets) {orb_rep = j; if(ordeg[j] == i + 1) break;}
			o = nppso[orb_rep] - 1; e = p.nI2B[orb_rep] - nppso[orb_rep];
			for_Int(j, 0, o) {
				if(occweight[orb_rep][j] < weight_freze) freze_o++;
				else if(occweight[orb_rep][j] < weight_nooc) nooc_o++;
			}
			for_Int(j, nppso[orb_rep], p.nI2B[orb_rep]) {
				if(occweight[orb_rep][j] < weight_freze) freze_e++;
				else if(occweight[orb_rep][j] < weight_nooc) nooc_e++;
			}
			keep_o = o - nooc_o - freze_o; keep_e = e - nooc_e - freze_e;
			controler[i+1] = p.if_norg_imp ?  VecInt{freze_o, nooc_o, 1, 1, nooc_e, freze_e } : VecInt{1, freze_o, nooc_o, keep_o, 1, keep_e, nooc_e, freze_e };
		}
		if(mm) WRN(NAV4(p.if_norg_degenerate, p.nooc_mode, weight_nooc, weight_freze));
		// if(mm) WRN(NAV(controler));
		// p.nooc_mode = STR("cpnooc");
		controler[0][1] = -0; controler[0][p.ndiv-1] = 0;
		p.according_controler(controler, ordeg);
	}
}

void APIzen::seimp_fixer(ImGreen& seimp_in) {
	for_Int(i, 0, seimp_in.nomgs) {
		for_Int(m, 0, seimp_in.norbs) {
			// if (mm) WRN(NAV3(i, real(seimp_in.g[i][m][m]), imag(seimp_in.g[i][m][m])))	
			if (imag(seimp_in.g[i][m][m]) > 0) {
				for_Int(n, 0, seimp_in.norbs) {
					seimp_in.g[i][n][m] -= I * imag(seimp_in.g[i][n][m]);
					seimp_in.g[i][m][n] -= I * imag(seimp_in.g[i][m][n]);
				}
			}
		}
	}
}


//change for the stand C++----------------------------------------------------------------------------------------------------------------------------

void APIzen::read_norg_setting(
	const std::string& filename,
	std::string& CoulombF,
	double& U,
	double& J,
	std::vector<int>& restrain,
	std::vector<int>& distribute
) {
	std::ifstream file(filename);
	std::string line;
	while (std::getline(file, line)) {
		std::istringstream iss(line);
		std::string key, eq;
		iss >> key >> eq;
		// if (eq != "=") {
		// 	ERR("Some problem here.")
		//     continue;
		// }
		if (key == "NOOC") {
			iss >> p.nooc_mode;
		}
		else if (key == "CoulombF") {
			iss >> CoulombF;
		}
		else if (key == "Minsulator") {
			iss >> p.if_norg_degenerate;
		}
		else if (key == "nband") {
			iss >> nband;
		}
		else if (key == "norbs") {
			iss >> norbs;
		}
		else if (key == "mune") {
			iss >> mu;
		}
		else if (key == "beta") {
			iss >> p.beta;
		}
		else if (key == "Uc") { 
			iss >> U;
		}
		else if (key == "Jz") { 
			iss >> J;
		}
		else if (key == "weight_nooc") { 
			iss >> weight_nooc;
		}
		else if (key == "weight_freze") { 
			iss >> weight_freze;
		}
		else if (key == "restrain") {
            std::string value;
            getline(iss, value); 
            std::istringstream value_stream(value);
            char ch;
            while (value_stream >> ch) {
                if (ch == '*') {
                    restrain.push_back(0); 
                } else if (ch == '-' || isdigit(ch)) {
                    value_stream.unget(); 
                    int num;
                    value_stream >> num;
                    restrain.push_back(num);
                }
            }
        } else if (key == "distribute") {
            int num;
            while (iss >> num) {
                distribute.push_back(num);
            }
        }
	}
}
