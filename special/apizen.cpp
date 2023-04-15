/*
coded by Jia-Ming Wang (jmw@ruc.edu.cn, RUC, China) date 2022 - 2023
*/

#include "apizen.h"

APIzen::APIzen(const MyMpi& mm_i, Prmtr& prmtr_i, const Str& file, const Int test_mode_i) :
	mm(mm_i), p(prmtr_i),num_omg(prmtr_i.num_omg),
	num_nondegenerate(-1), test_mode(test_mode_i), dmft_cnt(0)
{
	read_ZEN(file);	p.U = Uc; p.mu = mu; p.jz = Jz; p.nband = nband; p.norg_sets = p.norbs = norbs;
	p.templet_restrain = restrain; p.templet_control = distribute;
	p.after_modify_prmtr(); p.recalc_partical_number();
	if(mm) p.print();
	ImGreen hb(nband, p);
	for_Int(j, 0, hb.nomgs) for_Int(i, 0, nband)  hb.g[j][i][i] = - imfrq_hybrid_function[i][j];
	hb.write_zen("hb_zen", "Read");
	Bath bth(mm, p);
	bth.read_ose_hop();	bth.bath_fit(hb, or_deg_idx.truncate(0,nband));			if (mm)	bth.write_ose_hop(dmft_cnt);


	if (mm) std::cout << std::endl;						// blank line


	Impurity imp(mm, p, bth, or_deg_idx.truncate(0,nband));
	ImGreen hb_imp(p.nband, p);   	imp.find_hb(hb_imp); 	if (mm) hb_imp.write_zen("hb_imp", "Fit");
	// imp.update();
	if (mm) imp.write_H0info(bth, MAX(or_deg_idx));


	// NORG norg(mm, prmtr_i);		// for test
	{// for test:
		// IFS ifs_a("ru" + norg.scsp.nppso_str() + ".bi");
		// if (ifs_a) for_Int(i, 0, norg.uormat.size()) biread(ifs_a, CharP(norg.uormat[i].p()), norg.uormat[i].szof());
		// norg.up_date_h0_to_solve(imp.impH, 1);
		// if (mm) {
		// 	OFS ofs_a;
		// 	ofs_a.open("ru" + norg.scsp.nppso_str() + ".bi");
		// 	for_Int(i, 0, norg.uormat.size()) biwrite(ofs_a, CharP(norg.uormat[i].p()), norg.uormat[i].szof());
		// }		
	}



	NORG norg(choose_cauculation_style("ful_pcl_sch", imp));
	ImGreen g0imp(p.nband, p);	imp.find_g0(g0imp);									if (mm)	g0imp.write_zen("g0imp");
	ImGreen gfimp(p.nband, p);	norg.get_gimp(gfimp, or_deg_idx.truncate(0,nband));	if (mm) gfimp.write_zen("gfimp");
	ImGreen seimp(p.nband, p);	seimp = g0imp.inverse() - gfimp.inverse();			if (mm) seimp.write_zen("seimp");


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
	// ImGreen gfimp(p.nband, p);	norg.get_gimp(gfimp);				if (mm) gfimp.write_zen("gfimp");


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
		Str norgdata(file+".norg.in");
		IFS ifs(norgdata);
		if (!ifs) {
			ERR(STR("file opening failed with ") + NAV(norgdata))
		}
		else {
			Str strr;
			while (test_mode == 0){
				ifs >> strr;
				if(strr == "nband") { ifs >> strr; ifs >> nband; }
				if(strr == "norbs") { ifs >> strr; ifs >> norbs; }
				if(strr == "Uc") { ifs >> strr; ifs >> Uc; }
				if(strr == "Jz") { ifs >> strr; ifs >> Jz; }
				if(strr == "mune") { ifs >> strr; ifs >> mu; }
				if(strr == "beta") { ifs >> strr; ifs >> p.beta; }
				if(strr == "nmesh") { ifs >> strr; ifs >> p.nmesh; }
				if (strr == "restrain") {
					Int l2, l1, r1, r2;
					ifs >> strr; ifs >> strr;
					ifs >> l2;	ifs >> l1; ifs >> strr; ifs >> r1; ifs >> r2;
					restrain = { 0, l2, l1, 0, r1, r2 };
					if (mm) PIO("Finish the restrain input:" + NAV(restrain.mat(1,restrain.size())));
				}
				if (strr == "distribute") {
					Int l2, l1, m0, r1, r2;
					ifs >> strr; ifs >> strr;
					ifs >> l2;	ifs >> l1; ifs >> m0; ifs >> r1; ifs >> r2;
					distribute = { 1, l2, l1, m0, r1, r2 };
					if (mm) PIO("Finish the division distribute input:" + NAV(distribute.mat(1,distribute.size())));
				}
				if (strr == "norm_mode") {
					ifs >> strr; ifs >> norm_mode;
					if (mm) if (norm_mode) WRN("Dear cutomer,you are now in the test mode PLEASE modify the console file manually.");
				}
				if (strr == "test_mode") {
					ifs >> strr; ifs >> test_mode;
					if (mm) if (test_mode) WRN("Dear cutomer,you are now in the test mode PLEASE modify the console file manually.");
				}
				if (strr == "fast_mode") {
					ifs >> strr; ifs >> fast_mode;
					if (fast_mode)ERR("Sorry, the fast mode will coming Soon...");
				}				
				if (test_mode == 1 && strr == "bathose") {
					bathose.reset(SUM(distribute) - 1);
					Real input;
					ifs >> strr;
					for_Int(i, 0, bathose.size()) {
						ifs >> input;
						bathose[i] = input;
					}
				}
				if (test_mode == 1 && strr == "bathhop") {
					bathhop.reset(SUM(distribute) - 1);
					Real input;
					ifs >> strr;
					for_Int(i, 0, bathhop.size()) {
						ifs >> input;
						bathhop[i] = input;
					}
					if (mm) WRN("Finish the test_mode input: nimp = 2" + NAV4(Uc, Jz, bathose, bathhop));
				}
				if (!ifs) break;
			}
		}
		ifs.close();
	}

	imfrq_hybrid_function.reset(norbs,num_omg,0.);
	solver_eimp_data.reset(norbs,0.);
	or_deg_idx.reset(norbs, 0);

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
			VecReal temp_emps(norbs, 0.);
			for_Int(i, 0, norbs) {
				ifs >> drop_Int;
				ifs >> temp_emps[i]; ifs >> or_deg_idx[i];
				if (!ifs) ERR(STR("read_ZEN-in error with ") + NAV(eimpdata));
			}
			for_Int(j, 0, nband) for_Int(i, 0, 2) p.eimp[c++] = temp_emps[j+nband*i];
			// if (test_mode) num_nondegenerate = 1;
			num_nondegenerate = MAX(or_deg_idx);
		}
		if (num_nondegenerate <= 0)ERR(STR("read_ZEN-in error with ") + NAV2(eimpdata, num_nondegenerate));
		ifs.close();
	}
}



NORG APIzen::choose_cauculation_style(Str mode, Impurity &imp){
	if(mode == "fast"){
		Int band1(p.npartical[0]), band2(p.npartical[0]-2);
		p.according_nppso(p.npartical);
		if(mm) WRN(NAV(SUM(p.npartical)));
		NORG norg(mm, p);
		IFS ifs_a("ru" + norg.scsp.nppso_str() + ".bi");
		if (ifs_a) for_Int(i, 0, norg.uormat.size()) biread(ifs_a, CharP(norg.uormat[i].p()), norg.uormat[i].szof());
		norg.up_date_h0_to_solve(imp.impH, 1);
		if (mm)	{
			OFS ofs_a;
			ofs_a.open("ru" + norg.scsp.nppso_str() + ".bi");
			for_Int(i, 0, norg.uormat.size()) biwrite(ofs_a, CharP(norg.uormat[i].p()), norg.uormat[i].szof());
		}
		return norg;
	}
	if(mode == "stable"){
		Occler opcler(mm,p);
		NORG norg(opcler.find_ground_state_partical(imp.impH, or_deg_idx.truncate(0,nband)));
		return norg;
	}	
	if(mode == "ful_pcl_sch"){
		Occler opcler(mm,p);
		VEC<MatReal> uormat;
		VecInt ordeg(concat(or_deg_idx.truncate(0,nband),or_deg_idx.truncate(0,nband)).mat(2,nband).tr().vec()), nppso;
		Vec<VecInt> controler(MAX(or_deg_idx) + 1, VecInt(p.ndiv, 0));
		MatReal occnum, occweight;
		controler[0] = {0, -1, p.control_divs[0][2], 0, p.control_divs[0][4], 1};
		// if(mm) WRN(NAV(controler[0]));
		{
			NORG norg(opcler.find_ground_state_partical(imp.impH, or_deg_idx.truncate(0,nband)));
			uormat = norg.uormat;
			occnum = norg.occnum.mat(p.norg_sets, p.n_rot_orb/p.norg_sets); occweight = occnum;
			nppso = norg.scsp.nppso;
		}
		for_Int(i, 0, p.norg_sets) for_Int(j, 0, p.ndiv) if(occnum[i][j] > 0.5) occweight[i][j] = 1 - occnum[i][j];

		for_Int(i, 0, MAX(or_deg_idx)){
			Int o(0), freze_o(0), e(0), freze_e(0), orb_rep(0), nooc_o(0), nooc_e(0);
			for_Int(j, 0, p.norg_sets) {orb_rep = j; if(ordeg[j] == i + 1) break;}
			o = nppso[orb_rep] - 1; e = p.nI2B[orb_rep] - nppso[orb_rep];
			for_Int(j, 0, o) 							if(occweight[orb_rep][j] < 1e-6) freze_o++;
			for_Int(j, nppso[orb_rep], p.nI2B[orb_rep])	if(occweight[orb_rep][j] < 1e-6) freze_e++;
			nooc_o = o - freze_o; nooc_e = e - freze_e;
			controler[i+1] = VecInt{1, freze_o, nooc_o, 1, nooc_e, freze_e };
		}

		p.according_controler(controler, ordeg);
		// {// WRN
		// 	MatInt m_controler(MAX(or_deg_idx) + 1, p.ndiv);
		// 	for_Int(i, 0, controler.size()) m_controler[i] = controler[i];
		// 	if(mm) WRN(NAV3(ordeg,m_controler, p.control_divs));
		// }
		NORG frezeorb(mm, p);
		IFS ifs_a("ru" + frezeorb.scsp.nppso_str() + ".bi");
		if (ifs_a) for_Int(i, 0, frezeorb.uormat.size()) biread(ifs_a, CharP(frezeorb.uormat[i].p()), frezeorb.uormat[i].szof());
		frezeorb.up_date_h0_to_solve(imp.impH, 1);
		if (mm)	{
			OFS ofs_a;
			ofs_a.open("ru" + frezeorb.scsp.nppso_str() + ".bi");
			for_Int(i, 0, frezeorb.uormat.size()) biwrite(ofs_a, CharP(frezeorb.uormat[i].p()), frezeorb.uormat[i].szof());
		}
		return frezeorb;
	}
	if(mode == "one_pcl_test"){
		Occler opcler(mm,p);
		VEC<MatReal> uormat;
		VecInt ordeg(concat(or_deg_idx.truncate(0,nband),or_deg_idx.truncate(0,nband)).mat(2,nband).tr().vec()), nppso;
		Vec<VecInt> controler(MAX(or_deg_idx) + 1, VecInt(p.ndiv, 0));
		MatReal occnum, occweight;
		controler[0] = {0, -1, p.control_divs[0][2], 0, p.control_divs[0][4], 1};
		{
			Int band1(p.npartical[0]), band2(p.npartical[0]);
			if (nband == 5) p.npartical = { band1, band1, band1, band1, band2, band2, band1, band1, band2, band2 };
			if (nband == 3) { p.npartical = { band1, band1, band1, band1, band1, band1 }; /*p.npartical += 1;*/ }
			p.according_nppso(p.npartical);
			NORG norg(mm, p);
			IFS ifs_a("ru" + norg.scsp.nppso_str() + ".bi");
			if (ifs_a) for_Int(i, 0, norg.uormat.size()) biread(ifs_a, CharP(norg.uormat[i].p()), norg.uormat[i].szof());
			norg.up_date_h0_to_solve(imp.impH, 1);

			uormat = norg.uormat;
			occnum = norg.occnum.mat(p.norg_sets, p.n_rot_orb/p.norg_sets); occweight = occnum;
			nppso = norg.scsp.nppso;
			// if(mm) norg.write_state_info(0);
		}
		for_Int(i, 0, p.norg_sets) for_Int(j, 0, p.ndiv) if(occnum[i][j] > 0.5) occweight[i][j] = 1 - occnum[i][j];

		for_Int(i, 0, MAX(or_deg_idx)){
			Int o(0), freze_o(0), e(0), freze_e(0), orb_rep(0), nooc_o(0), nooc_e(0);
			for_Int(j, 0, p.norg_sets) {orb_rep = j; if(ordeg[j] == i + 1) break;}
			o = nppso[orb_rep] - 1; e = p.nI2B[orb_rep] - nppso[orb_rep];
			for_Int(j, 0, o) 							if(occweight[orb_rep][j] < 1e-6) freze_o++;
			for_Int(j, nppso[orb_rep], p.nI2B[orb_rep])	if(occweight[orb_rep][j] < 1e-6) freze_e++;
			nooc_o = o - freze_o; nooc_e = e - freze_e;
			controler[i+1] = VecInt{1, freze_o, nooc_o, 1, nooc_e, freze_e };
		}

		p.according_controler(controler, ordeg);
		// {// WRN
		// 	MatInt m_controler(MAX(or_deg_idx) + 1, p.ndiv);
		// 	for_Int(i, 0, controler.size()) m_controler[i] = controler[i];
		// 	if(mm) WRN(NAV3(ordeg,m_controler, p.control_divs));
		// }
		NORG frezeorb(mm, p); frezeorb.uormat = uormat;
		IFS ifs_a("ru" + frezeorb.scsp.nppso_str() + ".bi");
		if (ifs_a) for_Int(i, 0, frezeorb.uormat.size()) biread(ifs_a, CharP(frezeorb.uormat[i].p()), frezeorb.uormat[i].szof());
		frezeorb.up_date_h0_to_solve(imp.impH, 1);
		if (mm)	{
			OFS ofs_a;
			ofs_a.open("ru" + frezeorb.scsp.nppso_str() + ".bi");
			for_Int(i, 0, frezeorb.uormat.size()) biwrite(ofs_a, CharP(frezeorb.uormat[i].p()), frezeorb.uormat[i].szof());
			// frezeorb.write_state_info(1);
		}
		return frezeorb;
	}
}

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