/*
coded by Jia-Ming Wang (jmw@ruc.edu.cn, RUC, China) date 2022 - 2023
*/

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>

//for the stand C++ package ------------------------------------------------------------------------------------------------------------------------
#include "api_edmft.h"


APIedmft::APIedmft(const MyMpi& mm_i, Prmtr& prmtr_i, const Str& file) :
	mm(mm_i), p(prmtr_i),num_omg(prmtr_i.num_omg),
	num_nondegenerate(-1), dmft_cnt(0), weight_nooc(1E-4), weight_freze(1E-9)
{
	update(file);
	Bath bth(mm, p);
	Impurity imp(mm, p, bth, or_deg_idx);

	dmft_cnt++; update(file);
	ImGreen hb(nband, p);	
	for_Int(j, 0, hb.nomgs) 
		for_Int(i, 0, nband) hb.g[j][i][i] = -imfrq_hybrid_function[j][or_deg_idx[i * 2] - 1];		if (mm) hb.write_edmft("hb_read.txt", or_deg_idx);
	bth.read_ose_hop();	bth.bath_fit(hb, or_deg_idx);												if (mm) bth.write_ose_hop();
	imp.update();																					if (mm) imp.write_H0info(bth, MAX(or_deg_idx));
	ImGreen hb_imp(p.nband, p);		imp.find_hb(hb_imp); 											if (mm) hb_imp.write_edmft("hb_fit.txt", or_deg_idx);

	auto_nooc("ful_pcl_sch", imp);	NORG norg(mm, p);
	if (!norg.check_NTR()) norg.uormat = p.rotationU;
	else MatReal tmp_b = norg.read_NTR();
	norg.up_date_h0_to_solve(imp.impH, 1);															norg.write_impurtiy_occupation();
	MatReal tmp_e = norg.save_NTR();
	// MatReal local_multiplets_state = norg.oneedm.local_multiplets_state(norg.oneedm.ground_state);	if (mm)WRN(NAV(local_multiplets_state));
	ImGreen g0imp(p.nband, p);	imp.find_g0(g0imp);													if (mm)	g0imp.write_edmft("g0imp.txt", or_deg_idx);
	ImGreen gfimp(p.nband, p);	norg.get_gimp_eigpairs(gfimp, or_deg_idx);							if (mm) gfimp.write_edmft("Gf.out", or_deg_idx);
	ImGreen seimp(p.nband, p);	seimp = g0imp.inverse() - gfimp.inverse();							if (mm) seimp.write_edmft("Sig.out", or_deg_idx);

//hold for checking-----------------------------------------------------------------------------------------------------------------------------------
/*
	// if(mode == "realf_mode"){
		ReGreen g0imp_re(p.nband, p);	imp.find_g0(g0imp_re);													if (mm)	g0imp_re.write_edmft("Re_g0imp");
		ReGreen gfimp_re(p.nband, p);	norg.get_gimp_eigpairs(gfimp_re, or_deg_idx);							if (mm) gfimp_re.write_edmft("Re_gfimp");
		ReGreen seimp_re(p.nband, p);	seimp_re = g0imp_re.inverse() - gfimp_re.inverse();						if (mm) seimp_re.write_edmft("Re_seimp");
		ReGreen hd_exsp(p.nband, p);	norg.get_gimp_hdQPs(hd_exsp);											if (mm)	 hd_exsp.write_edmft("Re-hdex");
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
		if (mm) g_asnci.write_edmft("g_asnci");
	}
*/

	// if(mm) gfimp.write_occupation_info();
	// if(mm) WRN(NAV(gfimp.particle_number().diagonal()));
	// ImGreen seimp(p.nband, p);	seimp = g0imp.inverse() - gfimp.inverse();	
	// if (mm) seimp.write_edmft("before_fix_seimp");						seimp = fix_se(seimp);
	// if (mm) seimp.write_edmft("seimp");
}

void APIedmft::test_for_fitting(const Bath& bth, const ImGreen& hby_i, Int num)
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

// For this function only suit for the t2g mode.
void APIedmft::read_eDMFT(const Str& file)
{
	{
		std::vector<double> Ed;
		std::vector<int> Deg;
		double J;
		std::string CoulombF;
		double beta;
		double U;
		std::vector<int> restrain_t;
		std::vector<int> distribute_t;
		read_norg_setting("PARAMS.norg", Ed, Deg, J, CoulombF, beta, U, restrain_t, distribute_t);
		//--------------------------------------------------
		// ! Here only suit for the t2g orbital.
		nband = 3;	norbs = 6;	mu = 0;	
		if (CoulombF != STR("'Ising'")) WRN("Now we only support for the Ising type interaction." + NAV(CoulombF))
		Uc = U + (8.0/7.0) * J;	Jz= (632.0/819.0) * J;
		p.beta = beta;

		p.eimp.reset(concat(VecReal(Ed), VecReal(Ed)).mat(2, Ed.size()).tr().vec());
		or_deg_idx.reset(concat(VecInt(Deg), VecInt(Deg)).mat(2, Deg.size()).tr().vec());
		num_nondegenerate = MAX(or_deg_idx);
		restrain = VecInt(restrain_t);
		distribute = VecInt(distribute_t);
		// if (mm) WRN(NAV2(VecReal(Ed), VecInt(Deg)));
		//--------------------------------------------------

	}

	if (mm) WRN(NAV2(p.eimp, or_deg_idx))
	if (mm) WRN(NAV7(restrain, distribute, Uc, Jz, p.beta, weight_nooc, weight_freze))


	imfrq_hybrid_function.reset(num_omg,num_nondegenerate,0.);

	{// Delta.inp: to get the hyb function.
		Str hybdata("Delta.inp");
		IFS ifs(hybdata);
		if (!ifs) {
			ERR(STR("file opening failed with ") + NAV(hybdata))
		}
		else {
			for_Int(i, 0, p.nmesh) {
				Real omg(0.);
				ifs >> omg;
				for_Int(m, 0, num_nondegenerate) {
				Real re(0.), im(0.);
					ifs >> re;
					ifs >> im;
					imfrq_hybrid_function[i][m] = cmplx(re, im);
					// if (mm) WRN(NAV3(omg, real(imfrq_hybrid_function[i][m]), imag(imfrq_hybrid_function[i][m])))	
				}
			}
			if (!ifs) {
				ERR(STR("read-in error with ") + NAV(file))
			}
		}
		ifs.close();
	}
}


/*
ImGreen APIedmft::fix_se(const ImGreen& se) const{
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
void APIedmft::update(const Str& file) {
	{// modify the parameters from edmft.in
		read_eDMFT(file);	p.U = Uc; p.mu = mu; p.jz = Jz; p.nband = nband; p.norg_sets = p.norbs = norbs;
		p.templet_restrain = restrain; p.templet_control = distribute; p.project = NAV(nband) + "SrVO3";
		p.after_modify_prmtr(); p.recalc_partical_number();
		p.Uprm = p.U - 2 * p.jz;
		p.degel = 0;
		// if (mm) p.print();
	}

	// if (mm) std::cout << std::endl;						// blank line
}

bool APIedmft::if_lock(const Str file) const {
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


void APIedmft::auto_nooc(Str mode, const Impurity& imp) {
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
		// if(mm) WRN(NAV(controler));
		// p.nooc_mode = STR("cpnooc");
		controler[0][1] = -1; controler[0][p.ndiv-1] = 1;
		p.according_controler(controler, ordeg);
	}
}


//change for the stand C++----------------------------------------------------------------------------------------------------------------------------

void APIedmft::read_norg_setting(
    const std::string& filename,
    std::vector<double>& Ed,
    std::vector<int>& Deg,
    double& J,
    std::string& CoulombF,
    double& beta,
    double& U,
    std::vector<int>& restrain,
    std::vector<int>& distribute
) {
    std::ifstream file(filename);
    std::string line;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string key;
        iss >> key;
        if (key == "Ed") {
            char ch;
            while (iss >> ch && ch != ']') {
                if (ch != ',' && ch != '[') {
                    double value;
                    iss.unget();
                    iss >> value;
                    Ed.push_back(value);
                }
            }
        } else if (key == "Deg") {
            char ch;
            while (iss >> ch && ch != ']') {
                if (ch != ',' && ch != '[') {
                    int value;
                    iss.unget();
                    iss >> value;
                    Deg.push_back(value);
                }
            }
        } else if (key == "J") {
            iss >> J;
        } else if (key == "CoulombF") {
            iss >> CoulombF;
        } else if (key == "NOOC") {
            iss >> p.nooc_mode;
        } else if (key == "beta") {
            iss >> beta;
        } else if (key == "U") {
            iss >> U;
		}
		else if (key == "weight_nooc") { 
			iss >> weight_nooc;
		}
		else if (key == "weight_freze") { 
			iss >> weight_freze;
		}
		else if (key == "restrain") {
			char ch;
			while (iss >> ch && ch != ']') {
				if (ch == '0') {
					restrain.push_back(0);
				}
				else if (ch == '-' || isdigit(ch)) {
					iss.unget();
					int value;
					iss >> value;
					restrain.push_back(value);
				}
			}
		}
		else if (key == "distribute") {
			char ch;
			while (iss >> ch && ch != ']') {
				if (ch == '0') {
					distribute.push_back(0);
				}
				else if (ch == '-' || isdigit(ch)) {
					iss.unget();
					int value;
					iss >> value;
					distribute.push_back(value);
				}
			}
		}
	}
}

/*
int main() {
    std::vector<double> Ed;
    std::vector<int> Deg;
    double J;
    std::string CoulombF;
    double beta;
    double U;
    std::vector<int> restrain;
    std::vector<int> distribute;

    read_data("data.txt", Ed, Deg, J, CoulombF, beta, U, restrain, distribute);


	// Print the values
    std::cout << "Ed: ";
    for (double e : Ed) std::cout << e << " ";
    std::cout << "\nDeg: ";
    for (int d : Deg) std::cout << d << " ";
    std::cout << "\nJ: " << J << "\nCoulombF: " << CoulombF << "\nbeta: " << beta << "\nU: " << U << std::endl;
    std::cout << "restrain: ";
    for (int r : restrain) std::cout << r << " ";
    std::cout << "\ndistribute: ";
    for (int d : distribute) std::cout << d << " ";
    std::cout << std::endl;

    return 0;
}
*/