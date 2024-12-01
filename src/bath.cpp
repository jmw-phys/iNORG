/*
code developed and maintained by (jmw@ruc.edu.cn, RUC, China) date 2022 - 2024
*/

#include "bath.h"

Bath::Bath(const MyMpi& mm_i, const Prmtr& prmtr_i) :mm(mm_i), p(prmtr_i), uur(mm.id()),
	npart(p.nband), nb(npart), ni(npart), info(p.nband, 6, 0.)
{
	{ SLEEP(1); mm.barrier(); }		// make random seed output together

	//set the number of imp sites in per independent parts
	for_Int(i, 0, npart) {
		ni[i] = 1;
	}

	//set the number of bath sites in per independent parts
	for_Int(i, 0, npart) {
		nb[i] = p.nI2B[i * 2];
	}

	//initialize the fitting parameter for each independent parts
	fvb_init();

	// fix_fvb_symmetry();
	if (mm) bth_write_fvb(0);
}

void Bath::number_bath_fit(const ImGreen& hb_i, const VecInt or_deg)
{
	VEC<ImGreen> v_hb = generate_hb(hb_i, or_deg);      // vec of hyb function
	Vec<MatReal> v_bs = generate_bs(); 					// bath sum rule //! Here set for empty
	Vec<VecReal> v_a(npart), v_a_temp(MAX(or_deg));
	for_Int(i, 0, MAX(or_deg)) {
		// if (mm) WRN("HERE is fine ")
		v_a_temp[i].reset(number_bath_fit_part(v_hb[i], v_bs[i], i));
	}
	for_Int(i, 0, npart) v_a[i] = v_a_temp[or_deg[i * 2] - 1];
	fvb = arr2fvb(v_a);
	//fix_fvb_symmetry();
	// if (mm) bth_write_fvb(iter);
}

VecReal Bath::number_bath_fit_part(const ImGreen& vhb_i, const MatReal& vbs_i, Int idx)
{
	Real err;
	VecReal a;
	Int nmin;
	VecReal a_mrq = fvb2arr(idx);
	std::tie(err, a, nmin) = bath_fit_number_contest(a_mrq, vhb_i, vbs_i, idx);
	// if (mm) WRN("HERE is fine ")
	if (mm) {
		HybErr hyberr(p, vhb_i, vbs_i, nb[idx]);
		Real err = hyberr(a);
		Real err_crv = hyberr.err_curve(a);
		Real err_bsr = hyberr.err_bsr(a);
		Real err_reg = hyberr.err_reg(a);
		Real a_norm = a.norm();
		using namespace std;
		Int i = idx;
		cout << setw(4) << "  " << NAV7(nb[i],nmin, err, err_crv, err_bsr, err_reg, a_norm) << "  " << present() << endl;
		NAV6(Int(info[idx][0]=Real(nmin)), info[idx][1]=err, info[idx][2]=err_crv, info[idx][3]=err_bsr, info[idx][4]=err_reg, info[idx][5]=a_norm);
	}

	return a;
}

std::tuple<Real, VecReal, Int> Bath::bath_fit_number_contest(const VecReal& a0, const ImGreen& vhb_i, const MatReal& vbs_i, Int idx)
{
	HybErr hyberr(p, vhb_i, vbs_i, nb[idx]);

	// if (mm) WRN("HERE is fine ")
	const Int np = a0.size();
	// const Int ntry_fine = ntry / 4;
	// Int ntry = 1 * (mm.np() - 1);
	Int ntry_fine = MAX(16, mm.np() - 1);
	Int ntry = MAX(12 * ntry_fine, 516);
	//Int ntry = 2 * (mm.np() - 1);
	//const Int ntry_fine = 24;
	//const Int ntry = MAX(12 * ntry_fine, 512);
	//const Int ntry = 3*(mm.np()-1);
	//const Int ntry = mm.np()-1;
	//const Int ntry = 1;
	const Real tol = SQR(MIN(1.e-8, std::pow(10, -5 - ((nb[idx] / 4 + 1) / 2)))); //add 6 bath,add 10^-1
	// const Real tol = SQR(MIN(1.e-10, std::pow(10, -5 - ((nb[idx] / 4 + 1) / 2)))); //add 6 bath,add 10^-1
	// const Real tol = SQR(MIN(1.e-6, std::pow(10, -5 - (nb[idx] / 4) / 4)));
	Int nmin = 0;		// number of fittings reaching the minimum
	MPI_Status status;
	VecReal a(np);
	VecReal a_optm = a0;
	Real err_optm = hyberr(a_optm);
	MatReal ma(ntry, np, 0.);
	VecReal ma_err(ntry);
	// if (mm) WRN("HERE is fine "+NAV(err_optm))

	Int ntimes = 2;
	Int ntime = 0;
	Int itmax = 0;
	Int ndone = 0;

	while (ntime < ntimes) {
		if (ntime == 0) {
			ntry = ntry;
			// itmax = 50;
			// ndone = 4;
			itmax = 500;
			ndone = 16;
		}
		else if (ntime == 1)
		{
			ntry = (mm.np() - 1);
			// itmax = 150;
			// ndone = 4;
			itmax = 1500;
			ndone = 16;
		}

		if (mm) {
			Int ntot = ntry;
			//		Int ntot = 2;
			Int nsnd = 0;
			Int nrcv = 0;
			Int nfst = 1;
			Int itry = 0;
			while (nrcv < ntot) {
				if (nfst < mm.np()) {		// send the first ntot-1 sets of parameters
					if (nsnd < ntot) {
						if (ntime == 0) {
							a = next_initial_fitting_parameters(a0, ntry_fine, itry);
						}
						else if (ntime == 1) {
							a = ma[itry];
							++itry;
						}
						mm.Send(a, nfst, 1);
						++nsnd;
					}
					else {
						mm.Send(a, nfst, 0);
					}
					++nfst;
				}
				else {
					mm.Recv(a, status);

					//				WRN(NAV(nrcv));
					Int sndr = status.MPI_SOURCE;
					Real err = hyberr(a);
					if (ntime == 0) {
						ma[nrcv] = a;
						ma_err[nrcv] = err;
					}
					++nrcv;
					if (true) {
						Real err_crv = hyberr.err_curve(a);
						//Real err_crv = 1.;
						Real err_bsr = hyberr.err_bsr(a);
						Real err_reg = hyberr.err_reg(a);
						//Real err_bsr = 1.;
						Real a_norm = a.norm();
						// WRN(NAV5(sndr, ntot, nsnd, nrcv, itry) + ", " + NAV6(err_optm, err, err_crv, err_reg,err_bsr, a_norm));
					}
					if (err_optm - sqrt(tol) > err) { nmin = 0; }
					if (err_optm > err) { a_optm = a; err_optm = err; }
					nmin += err - err_optm < sqrt(tol);

					if (nsnd < ntot) {
						if (ntime == 0) {
							a = next_initial_fitting_parameters(a0, ntry_fine, itry);
						}
						else if (ntime == 1) {
							a = ma[itry];
							++itry;
						}
						mm.Send(a, sndr, 1);
						++nsnd;
					}
					else {
						mm.Send(a, sndr, 0);
					}
				}
			}
			if (ntime == 0) {
				slctsort(ma_err, ma);
			}

			++ntime;
		if(mm) WRN(NAV(itry));
		}
		else {
			while (true) {
				mm.Recv(a, status, mm.ms());
				if (status.MPI_TAG == 0) break;
				FitMrq<HybErr> mrq(hyberr.x, hyberr.y, hyberr.sig, a, hyberr, tol, itmax, ndone);
				Int mrq_fit_info = mrq.fit();
				// WRN(NAV2(itmax,ndone));
				mm.Send(mrq.a, mm.ms(), 1);
			}
			++ntime;
		}
	}


	mm.Bcast(err_optm);
	mm.Bcast(a_optm);
	mm.Bcast(nmin);

	return std::make_tuple(err_optm, a_optm, nmin);
}


VecReal Bath::next_initial_fitting_parameters(const VecReal& a0, const Int& ntry_fine, Int& itry)
{
	VecReal a = a0;
	++itry;
	if (itry == 1) {
	}
	else if (itry <= ntry_fine * 1) {
		a = perturb(a, 0.1);
	}
	else if (itry <= ntry_fine * 2) {
		a = perturb(a, 0.2);
	}
	else if (itry <= ntry_fine * 3) {
		a = perturb(a, 0.4);
	}
	else if (itry <= ntry_fine * 4) {
		a = perturb(a, 0.8);
	}
	else {
		a = perturb(a, 1.6);
	}
	return a;
}


VEC<ImGreen> Bath::generate_hb(const ImGreen& hb) const {
	VEC<ImGreen> vhb;
	for_Int(i, 0, hb.norbs) {
		ImGreen vhb_i(ni[i], p);
		vhb_i = ImGreen(hb, i, i);
		vhb.push_back(vhb_i);
	}
	return vhb;
}

VEC<ImGreen> Bath::generate_hb(const ImGreen& hb, const VecInt or_deg) const {
	VEC<ImGreen> vhb;
	Idx deg_idx(0);
	for_Int(i, 0, hb.norbs) {
		if (deg_idx < or_deg[i*2]) {
			ImGreen vhb_i(ni[i], p);
			vhb_i = ImGreen(hb, i, i);
			vhb.push_back(vhb_i);
			++deg_idx;
		}
	}
	return vhb;
}

// here export the empty bath sum rule, and we replace it with the hight term of fitting.
Vec<MatReal> Bath::generate_bs() {
	Vec<MatReal> vbs(npart);
	for_Int(i, 0, npart) {
		vbs[i].reset(ni[i], ni[i], 0.);
	}
	return vbs;
}



void Bath::fvb_init()
{
	fvb.reset(npart);
	for_Int(i, 0, npart) {
		fvb[i].reset(2);
		fvb[i][0].reset(ni[i], nb[i]);
		fvb[i][1].reset(1, nb[i]);

		for_Int(j, 0, fvb[i].size()) {
			for_Int(n, 0, fvb[i][j].nrows()) {
				uur(fvb[i][j][n]);
				fvb[i][j][n] -= 0.5;
				fvb[i][j][n] *= 4;
			}
		}
	}
}

void Bath::fix_fvb_symmetry() {

	//if (mm) WRN("!!!Be careful about fix_fvb_symmetry()")
		//for_Int(i, 0, fvb.size()) {
		//	fvb[i][1] = -fvb[i][0];
		//	fvb[i][3] = fvb[i][2];
		//	fvb[i][5] = -fvb[i][4];

		//}
}

void Bath::bth_write_fvb(Int iter_cnt, const Str bath_name) const {
	using namespace std;
	OFS ofs;
	if (iter_cnt < 0) {
		if (bath_name.empty()) ofs.open("ose_hop");
		else ofs.open(bath_name + "ose_hop.out");
	}
	if (iter_cnt > 0) ofs.open(iox + "zic" + prefill0(iter_cnt, 3) + ".ose_hop.txt");
	
	ofs << iofmt("sci");
	for_Int(i, 0, npart) {
		ofs << "set[" << i << "]" << endl;
		for_Int(j, 0, fvb[i].size()) {
			switch (j) {
			case 0:
				ofs << setw(w_Real) << "the" << endl;
				break;
			case 1:
				ofs << setw(w_Real) << "epi" << endl;
				break;
			}
			for_Int(r, 0, fvb[i][j].nrows()) {
				for_Int(c, 0, fvb[i][j].ncols()) {
					ofs << setw(w_Real) << fvb[i][j][r][c];
				}
				ofs << endl;
			}
		}
	}
	ofs.close();
}

void Bath::bth_read_fvb(const Str& file) {
	IFS ifs(file);
	if (!ifs) {
		if(mm) WRN(STR("file opening failed with ") + NAV(file));
	}
	else {
		Str t;
		for_Int(i, 0, npart) {
			ifs >> t;
			for_Int(j, 0, fvb[i].size()) {
				ifs >> t;
				for_Int(r, 0, fvb[i][j].nrows()) {
					for_Int(c, 0, fvb[i][j].ncols()) {
						ifs >> fvb[i][j][r][c];
					}
				}
			}
		}
	}
	// fix_fvb_symmetry();
	ifs.close();
}

VecReal Bath::fvb2arr(const Int idx)
{
	VEC<Real> a_i;
	for_Int(j, 0, fvb[idx].size()) {
		for_Int(r, 0, fvb[idx][j].nrows()) {
			for_Int(c, 0, fvb[idx][j].ncols()) {
				a_i.push_back(fvb[idx][j][r][c]);
			}
		}
	}
	return VecReal(a_i);
}

Vec<Vec<MatReal>> Bath::arr2fvb(const Vec<VecReal>& arr)
{
	Vec<Vec<MatReal>> fvb_t(fvb);
	for_Int(i, 0, fvb.size()) {
		Int n(0);
		for_Int(j, 0, fvb[i].size()) {
			for_Int(r, 0, fvb[i][j].nrows()) {
				for_Int(c, 0, fvb[i][j].ncols()) {
					fvb_t[i][j][r][c] = arr[i][n];
					n++;
				}
			}
		}
	}
	return fvb_t;
}

// rely on imp model's frame
//! Here only suit for the 1 impurity for 1 set now.
MatReal Bath::find_hop() const
{
    MatReal h0(0, 0, 0.);
    for_Int(i, 0, p.nband) {
        const Int nb_i = p.nI2B[i * 2];
        MatReal h0_i(1 + nb_i, 1 + nb_i, 0.);
		VecReal hop(fvb[i][0][0]), ose(fvb[i][1][0]);
		slctsort(ose, hop);
        for_Int(j, 0, nb_i) {
            h0_i[0][j + 1] 		= hop[j];
            h0_i[j + 1][0] 		= hop[j];
            h0_i[j + 1][j + 1] 	= ose[j];
        }
        h0_i.reset(direct_sum(h0_i, h0_i));
        h0.reset(direct_sum(h0, h0_i));
    }
    return h0;
}