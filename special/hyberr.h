#pragma once

/*
code	by	Rong-Qiang He (rqhe@ruc.edu.cn, RUC, China) date 2013 - 2017
modify	by Jia-Ming Wang (jmw@ruc.edu.cn, RUC, China) date 2022 - 2023
*/

#include "specs.h"
#include "prmtr.h"
#include "green.h"

class HybErr {
public:
	const Prmtr& p;			// parameters
	const ImGreen& hb;		// hybridization function
	const Int nw;			// nw = p.fit_num_omg
	Int nb;					// nb = p.nbath
	VecInt x;
	VecReal y;
	VecReal sig;
	Real expect_err;
	Real regV_b;

	VecReal oseA;
	VecReal hopB;
public:
	HybErr(const Prmtr& p_i, const ImGreen& hb_i, const Int nb_i, VecReal oseA_i, VecReal hopB_i);
	HybErr(const Prmtr& p_i, const ImGreen& hb_i, const Int nb_i, VecReal oseA_i, VecReal hopB_i, Int orb_i);
	
	// return y = f(x; a);
	Real operator()(const Int x, const VecReal& a) const {
		VecReal temp;
		if (x < 2 * nw) {
			const Int n = x < nw ? x : x - nw;
			const Real omgn = p.Imomg(n);
			//const VecCmplx E = cmplx(temp.sm(nb, a.p()));
			VecCmplx E(oseA.size());
			for_Int(i,0,oseA.size()){
				E[i]=oseA[i]*std::exp(a[i]);
			}
			const VecCmplx S = INV(E - VecCmplx(nb, I * omgn));
			//const VecCmplx V = cmplx(temp.sm(nb, a.p() + nb));
			VecCmplx V(hopB.size());
			for_Int(i,0,hopB.size()){
				V[i]=hopB[i]*std::exp(a[oseA.size()+i]);
			}
			const VecCmplx Vco = V.co();
			Cmplx hyb = DOT(Vco, S * Vco);
			return x < nw ? real(hyb) : imag(hyb);
		}
		else if (x == 2 * nw + 0) {
			//const VecReal E = temp.sm(nb, a.p());
			VecReal E(oseA.size());
			for_Int(i,0,oseA.size()){
				E[i]=oseA[i]*std::exp(a[i]);
			}
			const VecReal E3 = E * E * E;
			const Real coefficient = expect_err * std::pow(p.bandw, -6);
			return coefficient * DOT(E3, E3);
		}
		else if (x == 2 * nw + 1) {
			//const VecReal V = temp.sm(nb, a.p() + nb);
			VecReal V(hopB.size());
			for_Int(i,0,hopB.size()){
				V[i]=hopB[i]*std::exp(a[oseA.size()+i]);
			}
			const VecReal V2 = V * V;
			const Real coefficient = expect_err * std::pow(regV_b, 3);
			return coefficient * DOT(INV(V2), INV(V));
		}
		else if (x == 2 * nw + 2) {
			//const VecReal V = temp.sm(nb, a.p() + nb);
			VecReal V(hopB.size());
			for_Int(i,0,hopB.size()){
				V[i]=hopB[i]*std::exp(a[oseA.size()+i]);
			}
			const VecReal Vco = V.co();
			Real hyb = DOT(Vco, Vco);
			return hyb;
		}
	}
	// return y = f(x; a) and dy/da
	void operator()(const Int x, const VecReal& a, Real& y, VecReal& dyda) const;
	// return relative_err = SQRT(chi_sqr)
	Real operator()(const VecReal& a) const {
		VecReal fx(x.size());
		for_Int(i, 0, x.size()) {
			fx[i] = (*this)(x[i], a);
		}
		VecReal relative_dev = (fx - y) * INV(sig);
		return SQRT(DOT(relative_dev, relative_dev));
	}
	// return relative_err = SQRT(the curve part of chi_sqr)
	Real err_curve(const VecReal& a) const {
		VecReal fx(x.size());
		for_Int(i, 0, x.size()) {
			fx[i] = (*this)(x[i], a);
		}
		fx[2 * nw] = y[2 * nw];
		VecReal relative_dev = (fx - y) * INV(sig);
		return SQRT(DOT(relative_dev, relative_dev));
	}
	// return relative_err = SQRT(the part of ose regularization)
	Real err_regE(const VecReal& a) const {
		const Int i = 2 * nw + 0;
		Real fx = (*this)(x[i], a);
		Real relative_dev = (fx - y[i]) * INV(sig[i]);
		return SQRT(relative_dev * relative_dev);
	}
	// return relative_err = SQRT(the part of hop regularization)
	Real err_regV(const VecReal& a) const {
		const Int i = 2 * nw + 1;
		Real fx = (*this)(x[i], a);
		Real relative_dev = (fx - y[i]) * INV(sig[i]);
		return SQRT(relative_dev * relative_dev);
	}
	// return relative_err = SQRT(the part of bath sum rule)
	Real err_bsr(const VecReal& a) const {
		const Int i = 2 * nw + 2;
		Real fx = (*this)(x[i], a);
		Real relative_dev = (fx - y[i]) * INV(sig[i]);
		return SQRT(relative_dev * relative_dev);
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