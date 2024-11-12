#pragma once

/*
code developed and maintained by (jmw@ruc.edu.cn, RUC, China) date 2022 - 2024
*/

#include "specs.h"
#include "prmtr.h"
#include "green.h"

class HybErr {
public:
	const Prmtr& p;			// parameters
	//const ImGreen& hb;	// hybridization function
	const Int nw;			// nw = p.fit_num_omg
	const Int nw_hterm;		// nw_hterm = 10
	VecInt data_n;			//number of data of hyb,bsr,reg


	VecReal vy;
	Vec<VecReal> vdyda;
	VecReal a_old;
    Int ni;
    Int nb; 		
    VecInt x;
	VecReal y;
	VecReal sig;
	Real bw;

private:

	Vec<MatReal> get_dBda() const; //idx of a, and Mat idx of B

	void get_value_mrq(const VecReal& a);

	MatReal a2B(const VecReal& a);

	void fill2array(const Vec<MatCmplx>& hb, const Vec<Mat<MatCmplx>>& dydB);

	void exchange_idx(Mat<MatCmplx>& mat) const;

	VecCmplx get_dyda(const MatCmplx& dydB);

public:
	HybErr(const Prmtr& p_i, const ImGreen& vhb_i, const MatReal& vbs_i, const Int nb_i);

	// return y = f(x; a);
	Real operator()(const Int x, const VecReal& a);

	// return y = f(x; a) and dy/da
	void operator()(const Int x, const VecReal& a, Real& y_i, VecReal& dyda_i);


// ------------------------------------------------------------------ err return function ------------------------------------------------------------------

	// return relative_err = SQRT(chi_sqr)
	Real operator()(const VecReal& a) {
		// Int rank;
		// MPI_Comm_rank(MPI_COMM_WORLD, &rank);

		VecReal fx(x.size());
		for_Int(i, 0, x.size()) {
			fx[i] = (*this)(x[i], a);
		}
		VecReal relative_dev = (fx - y) * INV(sig);
		
		// if (rank == 0) WRN(NAV3(fx, y, sig));
		return SQRT(DOT(relative_dev, relative_dev));
	}


	// return relative_err = SQRT(the curve part of chi_sqr)
	Real err_curve(const VecReal& a) {
		VecReal fx(data_n[0]);
		for_Int(i, 0, data_n[0]) {
			fx[i] = (*this)(x[i], a);
		}
		VecReal y_tmp(y.truncate(0, data_n[0]));
		VecReal sig_tmp(sig.truncate(0, data_n[0]));
		VecReal relative_dev = (fx - y_tmp) * INV(sig_tmp);
		return SQRT(DOT(relative_dev, relative_dev));
	}
	
	// return relative_err = SQRT(the part of bath sum rule)
	Real err_bsr(const VecReal& a) {
		VecReal fx(data_n[1], 0.);
		Int bgn(data_n[0]);
		for_Int(i, 0, data_n[1]) {
			fx[i] = (*this)(x[bgn + i], a);
		}
		VecReal y_tmp(y.truncate(bgn, bgn + data_n[1]));
		VecReal sig_tmp(sig.truncate(bgn, bgn + data_n[1]));
		VecReal relative_dev = (fx - y_tmp) * INV(sig_tmp);
		return SQRT(DOT(relative_dev, relative_dev));
	}

	// return relative_err = SQRT(the part of ose regularization)
	Real err_reg(const VecReal& a) {
		VecReal fx(data_n[2], 0.);
		Int bgn(data_n[0] + data_n[1]);
		for_Idx(i, 0, data_n[2]){
			fx[i] = (*this)(x[bgn + i], a);
		}
		VecReal y_tmp(y.truncate(bgn, bgn + data_n[2]));
		VecReal sig_tmp(sig.truncate(bgn, bgn + data_n[2]));
		VecReal relative_dev = (fx - y_tmp) * INV(sig_tmp);
		return SQRT(DOT(relative_dev , relative_dev));
	}

	void write_xysig(Int iter_cnt) const {
		OFS ofs(iox + "zic" + prefill0(iter_cnt, 3) + ".xysig.txt");
		using namespace std;
		ofs << setw(6) << "x" << "  " << setw(w_Real) << "y" << "  " << setw(w_Real) << "sig" << endl;
		ofs << iofmt("sci");
		for_Int(i, 0, x.size()) {
			ofs << setw(6) << x[i] << "  " << setw(w_Real) << y[i] << "  " << setw(w_Real) << sig[i] << endl;
		}
	}

};

// ------------------------------------------------------------------🪦 code grave 🪦------------------------------------------------------------------
/*



*/