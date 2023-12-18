#include "bath.h"

Bath::Bath(const MyMpi& mm_i, const Prmtr& prmtr_i) :
	mm(mm_i), p(prmtr_i), nb(p.nI2B[0]), hb(1, p), uur(mm.id()), ose(p.nI2B[0]), hop(p.nI2B[0]), info(p.nband, 7, 0.),
	vec_ose(p.nband), vec_hop(p.nband), osea(p.nI2B[0]), oseA(p.nI2B[0]), hopb(p.nI2B[0]), hopB(p.nI2B[0])
{
	// make random seed output together
	{ SLEEP(1); mm.barrier(); }
	// init_ose_hop();
	
}

void Bath::bath_fit(const ImGreen& hb_i, Int iter)
{
	if(iter == 1) read_ose_hop(); IFS ifs("ose_hop");

	// Int band_i = 0;

	for_Int(band_i, 0, p.nband)
	{
		if(p.nband != hb_i[0].nrows()) ERR("some thing wrong with the hybrid function.")
		for_Int(i, 0, hb_i.nomgs) hb[i] = hb_i[i][band_i][band_i];
		ose.reset(p.nI2B[band_i]);   hop.reset(p.nI2B[band_i]); 
		osea.reset(p.nI2B[band_i]); hopb.reset(p.nI2B[band_i]);
		oseA.reset(p.nI2B[band_i]); hopB.reset(p.nI2B[band_i]);
		nb = p.nI2B[band_i];

		//adopt even or odd bath
		for_Int(i, 0, nb) {
			hopB[i] = 1.;
			if (i < nb / 2) oseA[i] = -1.;
			else oseA[i] = (nb % 2 == 1 && i == nb / 2) ? 0. : 1.;
		}

		/*
		if(ifs || iter > 1) {ose = vec_ose[band_i]; hop = vec_hop[band_i];} 
		else init_ose_hop();
		regularize_ose_hop();
		*/
		if(ifs || iter > 1) {
			ose = vec_ose[band_i]; 
			hop = vec_hop[band_i];
			for_Int(i, 0, nb) {
				hopb[i] = std::log(hop[i] / hopB[i]);
				osea[i] = (nb % 2 == 1 && i == nb / 2) ? 0. : std::log(ose[i] / oseA[i]);
			}
		} else init_osea_hopb();
		regularize_osea_hopb();
		//init_osea_hopb();
		//regularize_osea_hopb();

		const VecReal a0 = concat(osea, hopb);
		// if(mm) WRN(NAV(a0));
		Real err;
		VecReal a;
		Int nmin;
		// std::tie(err, a, nmin) = bath_fit_contest(a0);
		//if(mm) WRN(NAV2(a0,band_i));
		std::tie(err, a, nmin) = bath_fit_bsr(a0, band_i);
		osea = a.truncate(0, nb);
		hopb = a.truncate(nb, nb + nb);
		regularize_osea_hopb();
		vec_ose[band_i] = ose; vec_hop[band_i] = hop;
		if (mm) {
			const HybErr hyberr(p, hb, nb, oseA, hopB, band_i);
			const VecReal a = concat(osea, hopb);
			Real err = hyberr(a);
			Real err_crv = hyberr.err_curve(a);
			Real err_regE = hyberr.err_regE(a);
			//Real err_regE = 0.;
			Real err_regV = hyberr.err_regV(a);
			//Real err_regV = 0.;
			Real err_bsr = hyberr.err_bsr(a);
			//Real err_bsr = 0.; 
			//Real a_norm = a.norm();
			Real a_norm=std::sqrt(DOT(ose,ose)+DOT(hop,hop));
			using namespace std;
			cout << setw(4) << band_i+1 << "  " << NAV7(nmin, err, err_crv, err_regE, err_regV, err_bsr, a_norm) << "  " << present() << endl;
			NAV7(Int(info[band_i][0]=Real(nmin)), info[band_i][1]=err, info[band_i][2]=err_crv, info[band_i][3]=err_regE, info[band_i][4]=err_regV, info[band_i][5]=err_bsr, info[band_i][6]=a_norm);
		}
	}
	// for_Int(band_j, 1, p.nband) {vec_ose[band_j] = ose; vec_hop[band_j] = hop;} // add for same as band 1.
}

void Bath::bath_fit(const ImGreen& hb_i, VecInt or_deg)// for Zen mode
{
	read_ose_hop();IFS ifs("ose_hop");
	for_Int(degi, 0, MAX(or_deg)) {
		Int count(0), orb_rep(-1);
		VecCmplx hb_fit(hb.nomgs); 
		for_Int(i, 0, hb_i.norbs) {
			if (or_deg[i * 2] - 1 == degi) {
				for_Int(n, 0, hb.nomgs) hb_fit[n] += hb_i.g[n][i][i];
				count++;
			}
		}
		if(p.nband != hb_i[0].nrows()) ERR("some thing wrong with the hybrid function.")
		for_Int(i, 0, hb.nomgs) hb[i] = hb_fit[i] / Real(count);
		for_Int(j, 0, p.nband) { orb_rep = j; if (or_deg[j * 2] == degi + 1) break; }

		ose.reset(p.nI2B[orb_rep]);   hop.reset(p.nI2B[orb_rep]); 
		osea.reset(p.nI2B[orb_rep]); hopb.reset(p.nI2B[orb_rep]);
		oseA.reset(p.nI2B[orb_rep]); hopB.reset(p.nI2B[orb_rep]);
		nb = p.nI2B[orb_rep];

		//adopt even or odd bath
		for_Int(i, 0, nb) {
			hopB[i] = 1.;
			if (i < nb / 2) oseA[i] = -1.;
			else oseA[i] = (nb % 2 == 1 && i == nb / 2) ? 0. : 1.;
		}

		if(ifs) {
			ose = vec_ose[orb_rep]; 
			hop = vec_hop[orb_rep];
			for_Int(i, 0, nb) {
				hopb[i] = std::log(hop[i] / hopB[i]);
				osea[i] = (nb % 2 == 1 && i == nb / 2) ? 0. : std::log(ose[i] / oseA[i]);
			}
		} else init_osea_hopb();
		regularize_osea_hopb();

		const VecReal a0 = concat(osea, hopb);
		Real err;
		VecReal a;
		Int nmin;
		std::tie(err, a, nmin) = bath_fit_contest(a0);
		osea = a.truncate(0, nb);
		hopb = a.truncate(nb, nb + nb);
		regularize_osea_hopb();
		for_Int(i, 0, p.nband) if (or_deg[i * 2] - 1 == degi) vec_ose[i] = ose;
		for_Int(i, 0, p.nband) if (or_deg[i * 2] - 1 == degi) vec_hop[i] = hop;
		// if(mm) WRN(NAV2(vec_ose.size(),vec_hop.size()));
		if (mm) {
			const HybErr hyberr(p, hb, nb, oseA, hopB);
			const VecReal a = concat(osea, hopb);
			Real err = hyberr(a);
			Real err_crv = hyberr.err_curve(a);
			Real err_reg = hyberr.err_regE(a);
			//Real err_bsr = hyberr.err_bsr(a);
			Real a_norm = a.norm();
			using namespace std;
			cout << setw(4) << degi << "  " << NAV5(nmin, err, err_crv, err_reg,/* err_bsr,*/ a_norm) << "  " << present() << endl;
			NAV5(Int(info[degi][0]=Real(nmin)), info[degi][1]=err, info[degi][2]=err_crv, info[degi][3]=err_reg, /*err_bsr,*/ info[degi][6]=a_norm);
		}
	}
}


VecReal Bath::next_initial_fitting_parameters(const VecReal& a0, const Int& ntry_fine, Int& itry)
{
	VecReal a = a0;
	++itry;
	if (itry == 1) {
	}
	else if (itry <= ntry_fine) {
		perturb(a, 0.2, 12);
	}
	else if (itry <= ntry_fine * 2) {
		perturb(a, 0.4, 8);
	}
	else if (itry <= ntry_fine * 4) {
		perturb(a, 1.0, 16);
	}
	else {
		Real exponent = 4.;
		for_Int(i, 0, a.size()) {
			Real re = uur() - 0.5;
			a[i] = SIGN(std::pow(10., (8 * ABS(re) - exponent)), a[i]);
		}
	}
	return a;
}

std::tuple<Real, VecReal, Int> Bath::bath_fit_contest(const VecReal& a0)
{
	const HybErr hyberr(p, hb, nb, oseA, hopB);
	const Int np = a0.size();
	const Int ntry_fine = MAX(16, mm.np() - 1);
	// const Int ntry = MAX(64 * ntry_fine, 200);
	const Int ntry = MAX(1024 * ntry_fine, 2000);
	const Real tol = 1.e-12;
	Int nmin = 0;		// number of fittings reaching the minimum
	MPI_Status status;
	VecReal a(np);
	VecReal a_optm = a0;
	Real err_optm = hyberr(a_optm);

	if (mm) {
		Int ntot = ntry;
//		Int ntot = 2;
		Int nsnd = 0;
		Int nrcv = 0;
		Int nfst = 1;
		Int itry = 0;
		while (nrcv < ntot) {
			if (nfst < mm.np()) {
				if (nsnd < ntot) {
					a = next_initial_fitting_parameters(a0, ntry_fine, itry);
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
				++nrcv;
				Int sndr = status.MPI_SOURCE;
				Real err = hyberr(a);
				if (false) {
					Real err_crv = hyberr.err_curve(a);
					Real err_reg = hyberr.err_regE(a);
					Real a_norm = a.norm();
					WRN(NAV5(sndr, ntot, nsnd, nrcv, itry) + ", " + NAV5(err_optm, err, err_crv, err_reg, a_norm));
				}
				if (err_optm - tol > err) { nmin = 0; }
				if (err_optm > err) { a_optm = a; err_optm = err; }
				nmin += err - err_optm < tol;

				if (nsnd < ntot) {
					a = next_initial_fitting_parameters(a0, ntry_fine, itry);
					mm.Send(a, sndr, 1);
					++nsnd;
				}
				else {
					mm.Send(a, sndr, 0);
				}
			}
		}
	}
	else {
		while (true) {
			mm.Recv(a, status, mm.ms());
			if (status.MPI_TAG == 0) break;
			FitMrq<HybErr> mrq(hyberr.x, hyberr.y, hyberr.sig, a, hyberr, tol);
			//for_Int(i, 0, a.size()/2) mrq.hold(i, a[i]);
			if (nb % 2 == 1) mrq.hold(nb / 2, 0.);
			Int mrq_fit_info = mrq.fit();
			mm.Send(mrq.a, mm.ms(), 1);
		}
	}

	mm.Bcast(err_optm);
	mm.Bcast(a_optm);
	mm.Bcast(nmin);
	return std::make_tuple(err_optm, a_optm, nmin);
}

std::tuple<Real, VecReal, Int> Bath::bath_fit_bsr(const VecReal& a0, const Int& orb_i)
{
	const HybErr hyberr(p, hb, nb, oseA, hopB, orb_i);
	const Int np = a0.size();
	const Int ntry_fine = MAX(16, mm.np() - 1);
	// const Int ntry = MAX(128 * ntry_fine, 10);
	const Int ntry = MAX(1024 * ntry_fine, 2000);
	const Real tol = 1.e-12;
	Int nmin = 0;		// number of fittings reaching the minimum
	MPI_Status status;
	VecReal a(np);
	VecReal a_optm = a0;
	Real err_optm = hyberr(a_optm);

	if (mm) {
		Int ntot = ntry;
//		Int ntot = 2;
		Int nsnd = 0;
		Int nrcv = 0;
		Int nfst = 1;
		Int itry = 0;
		while (nrcv < ntot) {
			if (nfst < mm.np()) {
				if (nsnd < ntot) {
					a = next_initial_fitting_parameters(a0, ntry_fine, itry);
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
				++nrcv;
				Int sndr = status.MPI_SOURCE;
				Real err = hyberr(a);
				if (false) {
					Real err_crv = hyberr.err_curve(a);
					Real err_regE = hyberr.err_regE(a);
					//Real err_regE = 0.;
					Real err_regV = hyberr.err_regV(a);
					//Real err_regV = 0.; 
					Real err_bsr = hyberr.err_bsr(a);
					//Real err_bsr = 0.;
					Real a_norm = a.norm();
					WRN(NAV5(sndr, ntot, nsnd, nrcv, itry) + ", " + NAV7(err_optm, err, err_crv, err_regE, err_regV, err_bsr, a_norm));
				}
				if (err_optm - tol > err) { nmin = 0; }
				if (err_optm > err) { a_optm = a; err_optm = err; }
				nmin += err - err_optm < tol;

				if (nsnd < ntot) {
					a = next_initial_fitting_parameters(a0, ntry_fine, itry);
					mm.Send(a, sndr, 1);
					++nsnd;
				}
				else {
					mm.Send(a, sndr, 0);
				}
			}
		}
	}
	else {
		while (true) {
			mm.Recv(a, status, mm.ms());
			if (status.MPI_TAG == 0) break;
			FitMrq<HybErr> mrq(hyberr.x, hyberr.y, hyberr.sig, a, hyberr, tol);
			//if ((a.size() / 2) % 2 != 0) mrq.hold(Int(a.size() / 4), 0.);
			//for_Int(i, 0, a.size()/2) mrq.hold(i, a[i]);
			/*
			mrq.hold(0, -1.);
			mrq.hold(1, -0.66666666666666);
			mrq.hold(2, -0.44444444444444);
			mrq.hold(3, 0.44444444444444);
			mrq.hold(4, 0.66666666666666);
			mrq.hold(5, 1.);
			*/
			if(nb%2==1){
				mrq.hold(nb/2, 0.);
			}
			Int mrq_fit_info = mrq.fit();
			//WRN(NAV(mrq_fit_info));
			mm.Send(mrq.a, mm.ms(), 1);
		}
	}

	mm.Bcast(err_optm);
	mm.Bcast(a_optm);
	mm.Bcast(nmin);
	return std::make_tuple(err_optm, a_optm, nmin);
}


//rely on imp model's frame
MatReal Bath::find_hop() const
{
    MatReal h0(0, 0, 0.);
    for_Int(i, 0, p.nband) {
        const Int nb_i = p.nI2B[i * 2];
        MatReal h0_i(1 + nb_i, 1 + nb_i, 0.);
        for_Int(j, 0, nb_i) {
            h0_i[0][j + 1] = vec_hop[i][j];
            h0_i[j + 1][0] = vec_hop[i][j];
            h0_i[j + 1][j + 1] = vec_ose[i][j];
        }
        h0_i.reset(direct_sum(h0_i, h0_i));
        h0.reset(direct_sum(h0, h0_i));
    }
    return h0;
}

void Bath::write_ose_hop(Int iter_cnt, const Str& bath_name) const {
	using namespace std;
	OFS ofs;
	if (iter_cnt < 0) {
		if (bath_name.empty()) ofs.open("ose_hop");
		else ofs.open(bath_name + "ose_hop.out");
	}
	if (iter_cnt > 0) ofs.open(iox + "zic" + prefill0(iter_cnt, 3) + ".ose_hop.txt");
	for_Int(band_i, 0, p.nband)	{
		ofs << iofmt("sci");
		ofs << setw(4) << band_i;
		for_Int(i, 0, p.nI2B[2 * band_i]) { ofs << "\t" << setw(w_Real) << vec_ose[band_i][i]; }
		ofs << endl;
	}
	for_Int(band_i, 0, p.nband)	{
		ofs << iofmt("sci");
		ofs << setw(4) << band_i;
		for_Int(i, 0, p.nI2B[2 * band_i]) { ofs << "\t" << setw(w_Real) << vec_hop[band_i][i]; }
		ofs << endl;
	}
	ofs << endl;ofs << endl;
}

void Bath::read_ose_hop() {
	using namespace std;
	// IFS ifs_ose;ifs_ose.open("ose.txt");
	IFS ifs("ose_hop");
	Str strr;
	if(ifs)for_Int(band_i, 0, p.nband)	{
		ifs >> strr;
		VecReal ose_t(p.nI2B[2 * band_i], 0);
		for_Int(i, 0, p.nI2B[2 * band_i]) { ifs >> ose_t[i]; }
		// vec_ose.push_back(ose_t);
		vec_ose[band_i] = ose_t;
		// if(mm) WRN(NAV(ose_t));
	}		
	if(ifs)for_Int(band_i, 0, p.nband)	{
		ifs >> strr;
		VecReal hop_t(p.nI2B[2 * band_i], 0);
		for_Int(i, 0, p.nI2B[2 * band_i]) { ifs >> hop_t[i]; }
		// vec_hop.push_back(hop_t);
		vec_hop[band_i] = hop_t;
		// if(mm) WRN(NAV(hop_t));
	}
}

//------------------------------------------------------------------ print out ------------------------------------------------------------------
void Bath::regularize_ose_hop() {
	slctsort(ose, hop);
	// after a unitary transformation, hop can always be
	// non-negative when there is only one impurity site
	hop = ABS(hop);
}

void Bath::regularize_osea_hopb() {
	for_Int(i,0,nb){
		ose[i]=oseA[i]*std::exp(osea[i]);
		hop[i]=hopB[i]*std::exp(hopb[i]);
	}
	//if(mm) WRN(NAV6(oseA,osea,ose,hopB,hopb,hop));
	VecReal E_tmp1=ose;
	VecReal E_tmp2=ose;
	VecReal E_tmp3=ose;
	VecReal E_tmp4=ose;
	slctsort(E_tmp1, osea);
	slctsort(E_tmp2, hopb);
	slctsort(E_tmp3, oseA);
	slctsort(E_tmp4, hopB);
	for_Int(i,0,nb){
		ose[i]=oseA[i]*std::exp(osea[i]);
		hop[i]=hopB[i]*std::exp(hopb[i]);
	}
	//if(mm) WRN(NAV6(oseA,osea,ose,hopB,hopb,hop));
}