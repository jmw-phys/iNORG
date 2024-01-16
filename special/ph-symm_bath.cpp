#include "ph-symm_bath.h"
/*
Bath::Bath(const MyMpi& mm_i, const Prmtr& prmtr_i) :
	mm(mm_i), p(prmtr_i), nb(p.nbath), hb(1, p), uur(mm.id()), ose(p.nbath), hop(p.nbath)
{
	// make random seed output together
	{ SLEEP(1); mm.barrier(); }
	init_ose_hop();
}
*/



Bath::Bath(const MyMpi& mm_i, const Prmtr& prmtr_i) :mm(mm_i), p(prmtr_i), uur(mm.id()),
	vec_ose(p.nband), vec_hop(p.nband), info(p.nband, 7, 0.)
{
	// make random seed output together
	{ SLEEP(1); mm.barrier(); }

	for_Int(i, 0, p.nband)
	{
		vec_ose[i].reset(p.nI2B[2*i]);
		vec_hop[i].reset(p.nI2B[2*i]);
	}
	init_vec_ose_hop();
}


void Bath::init_vec_ose_hop()
{
	for_Int(i, 0, vec_ose.size())
	{
		uur(vec_ose[i]);
		vec_ose[i] -= 0.5;
		vec_ose[i] *= 4.;
		uur(vec_hop[i]);
		vec_hop[i] -= 0.5;
		vec_hop[i] *= 4.;
		slctsort(vec_ose[i], vec_hop[i]);
		vec_hop[i] = ABS(vec_hop[i]);
	}
}

void Bath::number_bath_fit(const ImGreen& hb_i, Int iter,Int mode)
{
	// for_Int(band_i, 0, p.nband)
	for_Int(band_i, 1, 3) //! set band 0 same as band 1.
	{
		ImGreen hb(hb_i,band_i,band_i);

		VecReal a0;
		if (mode == 1)
		{
			if (vec_ose[band_i].size() % 2 == 0)
			{
				a0 = concat(vec_ose[band_i].truncate(0, vec_ose[band_i].size() / 2), vec_hop[band_i].truncate(0, vec_hop[band_i].size() / 2));
			}
			else
			{
				if (vec_ose[band_i].size() == 1)
				{
					a0 = vec_hop[band_i];
				}
				else if (vec_ose[band_i].size() > 1)
				{
					a0 = concat(vec_ose[band_i].truncate(0, vec_ose[band_i].size() / 2), vec_hop[band_i].truncate(0, vec_hop[band_i].size() / 2 + 1));
				}
				
			}
		}
		else if (mode == 0)
		{
			a0 = concat(vec_ose[band_i], vec_hop[band_i]);
		}
		
		Real err;
		VecReal a;
		Int nmin;


		std::tie(err, a, nmin) = bath_fit_number_contest(a0, vec_ose[band_i].size(), hb, band_i, mode);
		


		nb = vec_ose[band_i].size();

		if (mode == 1)
		{
			if (vec_ose[band_i].size() % 2 == 0)
			{
				vec_ose[band_i] = concat(a.truncate(0, nb / 2), -a.truncate(0, nb / 2));
				vec_hop[band_i] = concat(a.truncate(nb / 2, (nb + nb) / 2), a.truncate(nb / 2, (nb + nb) / 2));
			}
			else
			{
				VecReal tmp0(1, 0.);
				if (vec_ose[band_i].size() == 1)
				{
					vec_ose[band_i] = tmp0;
					vec_hop[band_i] = a;
				}
				else if (vec_ose[band_i].size() > 1)
				{
					vec_ose[band_i] = concat(concat(a.truncate(0, nb / 2), tmp0), -a.truncate(0, nb / 2));
					vec_hop[band_i] = concat(a.truncate(nb / 2, (nb + nb) / 2), a.truncate((nb / 2), (nb + nb) / 2 - 1));
				}
			}
		}
		else if (mode == 0)
		{
			vec_ose[band_i] = a.truncate(0, nb);
			vec_hop[band_i] = a.truncate(nb, nb + nb);
		}

		
		slctsort(vec_ose[band_i], vec_hop[band_i]);
		vec_hop[band_i] = ABS(vec_hop[band_i]);
		
		if (mm)
		{
			HybErr hyberr(p, hb, vec_ose[band_i].size(), band_i, mode);
			VecReal a;
			if (mode == 1)
			{
				if (vec_ose[band_i].size() % 2 == 0)
				{
					a = concat(vec_ose[band_i].truncate(0, vec_ose[band_i].size() / 2), vec_hop[band_i].truncate(0, vec_hop[band_i].size() / 2));
				}
				else
				{
					if (vec_ose[band_i].size() == 1)
					{
						a = vec_hop[band_i];
					}
					else if (vec_ose[band_i].size() > 1)
					{
						a = concat(vec_ose[band_i].truncate(0, vec_ose[band_i].size() / 2), vec_hop[band_i].truncate(0, vec_hop[band_i].size() / 2 + 1));
					}
				}
			}
			else if (mode == 0)
			{
				a = concat(vec_ose[band_i], vec_hop[band_i]);
			}
			
			Real err = hyberr(a);
			Real err_crv = hyberr.err_curve(a);
			Real err_reg = hyberr.err_reg(a);
			Real err_bsr = hyberr.err_bsr(a);
			Real a_norm = a.norm();
			using namespace std;
			cout << setw(4) << iter << "  " << NAV7(band_i, nmin, err, err_crv, err_reg, err_bsr, a_norm) << "  " << present() << endl;
			NAV7(Int(info[band_i][0]=Real(nmin)), info[band_i][1]=err, info[band_i][2]=err_crv, info[band_i][3]=err_reg, info[band_i][4]=0, info[band_i][5]=err_bsr, info[band_i][6]=a_norm);
		}
	{vec_ose[0] = vec_ose[1]; vec_hop[0] = vec_hop[1];} //! set band 0 same as band 1.
	}
}

std::tuple<Real, VecReal, Int> Bath::bath_fit_number_contest(const VecReal& a0, Int nb, const ImGreen&hb_i, Int orb_i,Int mode)
{

	const HybErr hyberr(p, hb_i, nb, orb_i,mode);

	const Int np = a0.size();
	const Int ntry_fine = MAX(16, mm.np() - 1);
	const Int ntry = MAX(64 * ntry_fine, 1024);
	const Real tol = 1.e-8;
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
					Real err_reg = hyberr.err_reg(a);
					Real err_bsr = hyberr.err_bsr(a);
					Real a_norm = a.norm();
					WRN(NAV5(sndr, ntot, nsnd, nrcv, itry) + ", " + NAV5(err_optm, err, err_crv, err_reg,/* err_bsr,*/ a_norm));
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
			if((a.size() % 2 == 1) && a[a.size() - 1] < 5E-5) mrq.hold(a.size() - 1, 0.0); //! TESTING!!
			Int mrq_fit_info = mrq.fit();
			mm.Send(mrq.a, mm.ms(), 1);
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