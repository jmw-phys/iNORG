/*
code	by	Rong-Qiang He (rqhe@ruc.edu.cn, RUC, China) date 2013 - 2017
modify	by Jia-Ming Wang (jmw@ruc.edu.cn, RUC, China) date 2022 - 2023
*/
#include "hyberr.h"

HybErr::HybErr(const Prmtr& p_i, const ImGreen& hb_i, const Int nb_i, VecReal oseA_i, VecReal hopB_i) :
	p(p_i), hb(hb_i), nw(p.fit_num_omg), nb(nb_i), regV_b(1E-10),
	x(2 * nw + 1), y(2 * nw + 1), sig(2 * nw + 1), oseA(oseA_i), hopB(hopB_i)
{
	// set the curve expected:
	expect_err = std::pow(10, -(nb - 1) / 2);

	// curve
	{
		Real mag_real = 0.;
		Real mag_imag = 0.;
		VecReal wght(2 * nw);
		for_Int(n, 0, nw) {
			x[n] = n;
			y[n] = real(hb[n][0][0]);
			mag_real += SQR(y[n]);
			x[nw + n] = nw + n;
			y[nw + n] = imag(hb[n][0][0]);
			mag_imag += SQR(y[nw + n]);
			wght[n] = wght[nw + n] = std::pow(p.Imomg(n) + p.fit_rsd, -p.fit_pow);
			// wght[n] = wght[nw + n] = 1;
		}
		mag_real = SQRT(mag_real / nw);
		mag_imag = SQRT(mag_imag / nw);
		mag_real = MAX(mag_real, 1.e-3 * mag_imag);
		wght *= INV(SUM(wght));
		for_Int(n, 0, nw) {
			sig[nw + n] = mag_imag / SQRT(wght[nw + n]);
			sig[n] = mag_real / SQRT(wght[n]);
		}
	}

	// the part of ose regularization
	if(x.size()>= 2 * nw + 1){
		x[2 * nw + 0]	= 2 * nw;
		y[2 * nw + 0]	= 0.;
		// // old
		// sig[2 * nw + 0] = nb * std::pow(8 * p.fit_max_omg, 4);
		// x*e^(0.5*(x/bw)^2)		// old
		// x*e^(0.5*(x/(bw/8.))^2) 	// new
		sig[2 * nw + 0] = 1;
	}
	// // the part of bath sum rule
	// if (x.size() >= 2 * nw + 2) {
	// 	x[2 * nw + 1]	= 2 * nw + 1;
	// 	y[2 * nw + 1]	= 0.;
	// 	sig[2 * nw] = 0.1 * (1 + ABS(p.bsr));  
	// }
}


HybErr::HybErr(const Prmtr& p_i, const ImGreen& hb_i, const Int nb_i, VecReal oseA_i, VecReal hopB_i, Int orb_i) :
	p(p_i), hb(hb_i), nw(p.fit_num_omg), nb(nb_i), regV_b(1E-10),
	x(2 * nw + 3), y(2 * nw + 3), sig(2 * nw + 3), oseA(oseA_i), hopB(hopB_i)
{
	// set the curve expected:
	expect_err = std::pow(10, -(nb - 1) / 2);

	// curve
	{	
		Real mag_real = 0.;
		Real mag_imag = 0.;
		VecReal wght(2 * nw);
		for_Int(n, 0, nw) {
			x[n] = n;
			y[n] = real(hb[n][0][0]);
			mag_real += SQR(y[n]);
			x[nw + n] = nw + n;
			y[nw + n] = imag(hb[n][0][0]);
			mag_imag += SQR(y[nw + n]);
			wght[n] = wght[nw + n] = std::pow(p.Imomg(n) + p.fit_rsd, -p.fit_pow);
			// wght[n] = wght[nw + n] = 1;
		}
		mag_real = SQRT(mag_real / nw);
		mag_imag = SQRT(mag_imag / nw);
		mag_real = MAX(mag_real, 1.e-3 * mag_imag);
		wght *= INV(SUM(wght));
		for_Int(n, 0, nw) {
			sig[nw + n] = mag_imag / SQRT(wght[nw + n]);
			sig[n] = mag_real / SQRT(wght[n]);
		}
	}
	
	// the part of ose regularization
	if(x.size()>= 2 * nw + 1){
		x[2 * nw + 0]	= 2 * nw;
		y[2 * nw + 0]	= 0.;
		// // old
		// sig[2 * nw + 0] = nb * std::pow(8 * p.fit_max_omg, 4);
		sig[2 * nw + 0] = 1;
	}
	// the part of hop regularization
	if (x.size() >= 2 * nw + 2) {
		x[2 * nw + 1]	= 2 * nw + 1;
		y[2 * nw + 1]	= 0.;
		// sig[2 * nw + 1] = 1E9 /(std::pow(10, -(nb+1)/2) * 64.);
		sig[2 * nw + 1] = 1;
		// WRN(NAV(sig[2 * nw + 1]))
	}
	// the part of bath sum rule
	if (x.size() >= 2 * nw + 3) {
		x[2 * nw + 2]	= 2 * nw + 2;
		y[2 * nw + 2]	= p.bsr[orb_i];
		sig[2 * nw + 2] = 0.1 * (1 + p.bsr[orb_i]);
		// WRN(NAV5(2 * nw + 2,x.size(),y.size(),x[2 * nw + 2],y[2 * nw + 2]));
	}
    
}

/*
chi_sqr = sum_i [(y_i - hb(x_i)) / sig_i]^2 + [(VV^+ - bsr) / sig_bsr]^2 + ose_regularization
chi_sqr = sum_i weight_i * [(y_i - hb(x_i)) / mag] ^ 2
        + weight_bsr * [(VV^+ - bsr) / mag_bsr] ^ 2
		+ ose_regularization
E = ose;			// on-site energies for bath sites
V = hop;			// hopping coefficients between impurity site and bath sites
S = INV(diag_matrix(E - i * omgn));
hyb = hb(omgn) = V * S * V^+ = SUM_i V_i * S_{i,i} * V_i^*
D = d(hb) / d(a) = concat(d(hb) / d(E), d(hb) / d(V))
d(S_{i,i}) / d(E_i) = d(1/(E_i - i * omgn)) / d(E_i) = -(E_i - i * omgn)^{-2} = -S{i,i}^2
(D_E)_i = d(hyb)/d(E_i) = V_i * d(S_{i,i}) / d(E_i) * V_i^* = V_i * -S{i,i}^2 * V_i^*
*/
void HybErr::operator()(const Int x, const VecReal& a, Real& y, VecReal& dyda) const
{
	VecReal temp;
	if (x < 2 * nw) {
		const Int n = x < nw ? x : x - nw;
		const Cmplx iomgn = p.Imz(n);
		//const VecCmplx E = cmplx(temp.sm(nb, a.p()));
		VecCmplx E(oseA.size());
		for_Int(i,0,oseA.size()){
			E[i] = oseA[i] * std::exp(a[i]);
		}
		const VecCmplx S = INV(E - VecCmplx(nb, iomgn));
		//const VecCmplx V = cmplx(temp.sm(nb, a.p() + nb));
		VecCmplx V(hopB.size());
		for_Int(i,0,hopB.size()){
			V[i] = hopB[i] * std::exp(a[oseA.size() + i]);
		}
		const VecCmplx Vco = V.co();
		Cmplx hyb = DOT(Vco, S * Vco);
		y = x < nw ? real(hyb) : imag(hyb);
		VecCmplx D_E = (V * Cmplx(-1.) * S * S * Vco) * E;
		VecCmplx D_V = (V * S + S * Vco) * V;
		VecCmplx D = concat(D_E, D_V);
		dyda = x < nw ? real(D) : imag(D);
	}
	// x*e^(0.5*(x/bw)^2)
	// (e^((0.5 x^2)/bw^2) * (bw^2 + x^2))/bw^2
	else if (x == 2 * nw + 0) { // the part of ose regularization
		//const VecReal E = temp.sm(nb, a.p());
		VecReal E(oseA.size());
		for_Int(i, 0, oseA.size()) {
			E[i] = oseA[i] * std::exp(a[i]);
		}
		const VecReal E3 = E * E * E;
		const Real coefficient = expect_err * std::pow(p.bandw, -6);
		y = coefficient * DOT(E3, E3);
		VecReal D_E = 6. * E3 * E3;
		VecReal D_V = VecReal(nb, 0.);
		VecReal D = concat(D_E, D_V);
		dyda = coefficient * D;
	}
	else if (x == 2 * nw + 1) { // the part of hop regularization
		//const VecReal V = temp.sm(nb, a.p() + nb);
		VecReal V(hopB.size());
		for_Int(i,0,hopB.size()){
			V[i]=hopB[i]*std::exp(a[oseA.size()+i]);
		}
		const VecReal V2 = V * V;
		const Real coefficient = expect_err * std::pow(regV_b, 3);
		y = coefficient * DOT(INV(V2), INV(V));
		VecReal D_E = VecReal(nb, 0.);
		VecReal D_V = -3. * INV(V2 * V);
		VecReal D = concat(D_E, D_V);
		dyda = coefficient * D;
	}
	else if (x == 2 * nw + 2) { // the part of bath sum rule
		//const VecReal V = temp.sm(nb, a.p() + nb);
		VecReal V(hopB.size());
		for_Int(i,0,hopB.size()){
			V[i]=hopB[i]*std::exp(a[oseA.size()+i]);
		}
		const VecReal Vco = V.co();
		Real hyb = DOT(Vco, Vco);
		y = hyb;
		VecReal D_E = VecReal(nb, 0.);
		VecReal D_V = 2.*(V * Vco);
		VecReal D = concat(D_E, D_V);
		dyda = D;
	}
	else {
		ERR(NAV(x));
	}
}