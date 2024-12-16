/*
code developed and maintained by (jmw@ruc.edu.cn, RUC, China) date 2022 - 2024
*/

#include "hyberr.h"


HybErr::HybErr(const Prmtr& p_i, const ImGreen& hb, const MatReal& bs, const Int nb_i) :
    p(p_i), nw(p.fit_num_omg), nb(nb_i), bw(p.fit_max_omg), data_n(3), nw_hterm(1)
{
    Int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);


    ni = hb[0].nrows();
    data_n[0] = nw * hb[0].nrows() * (hb[0].nrows() + 1);       //only pick the upper triangular part of the matrix
    // data_n[0] = 2 * nw * hb[0].nrows() * hb[0].nrows();
    // data_n[1] = bs.nrows() * (bs.nrows() + 1) / 2;
    data_n[1] = nw_hterm * hb[0].nrows() * (hb[0].nrows() + 1); // set the last nw_max omega data (real + imag part)
    data_n[2] = nb;                                             // spin up and down could nonequivalent
    Int ndata = SUM(data_n);
    x.reset(ndata);
    y.reset(ndata);
    sig.reset(ndata);
    vy.reset(ndata, 0.);
    vdyda.reset(ndata);
    a_old.reset((ni + 1) * nb, 0.);

    Int curve_shift = 0;
	{// curve
		//only pick the upper triangular part of the matrix
		VecReal mag_real(data_n[0] / (2 * nw), 0.);
		VecReal mag_imag(mag_real);
		VecReal wght(data_n[0], 0.);
		Int idx(0);
		for_Int(i, 0, ni) {
			for_Int(j, i, ni) {
				for_Int(n, 0, nw) {
					Int re_i(idx * 2 * nw + n);
					Int im_i(idx * 2 * nw + n + nw);
					x[re_i] = re_i;
                    y[re_i] = p.fit_points.empty() ? real(hb[n][i][j]):real(hb[p.fit_points[n]][i][j]);
					mag_real[idx] += SQR(y[re_i]);
					x[im_i] = im_i;
                    y[im_i] = p.fit_points.empty() ? imag(hb[n][i][j]):imag(hb[p.fit_points[n]][i][j]);
					mag_imag[idx] += SQR(y[im_i]);
                    if(n < curve_shift) wght[re_i] = wght[im_i] = (i == j) ? 1E-5 : 2E-5;
                    else wght[re_i] = wght[im_i] = (i == j) ? 1 : 2;
					//wght[re_i] = wght[im_i] = 1.;
				}
				mag_real[idx] = SQRT(1 + mag_real[idx] / nw);
				mag_imag[idx] = SQRT(1 + mag_imag[idx] / nw);
				++idx;
			}
		}
		wght *= INV(SUM(wght));
		idx = 0;
		for_Int(i, 0, ni) {
			for_Int(j, i, ni) {
				for_Int(n, 0, nw) {
					sig[idx * 2 * nw + n] = mag_real[idx] / SQRT(wght[idx * 2 * nw + n]);
					sig[idx * 2 * nw + nw + n] = mag_imag[idx] / SQRT(wght[idx * 2 * nw + nw + n]);
				}
				++idx;
			}
		}
        // if (rank == 0) WRN(NAV3(mag_real, mag_imag, wght));
	}

	/*
    // the part of bath sum rule
    //only pick the upper triangular part of the matrix
    Int bgn = data_n[0];
    idx = 0;
    for_Int(i, 0, ni) {
        for_Int(j, i, ni) {
            x[bgn + idx] = bgn + idx;
            y[bgn + idx] = bs[i][j];
            Real sig_val = std::pow(10, 6 + (nb / 2 - 4) / 11) * sqrt(16) * (1 + abs(bs[i][j]));
            sig[bgn + idx] = (i == j) ? sig_val : sig_val / SQRT(2);
            //sig[bgn + idx] = sig_val;
            ++idx;
        }
    }
	*/
	{// the part of hight ev term of hybber curve
        Int bgn = data_n[0];
		VecReal mag_real(data_n[1] / (2 * nw_hterm), 0.);
		VecReal mag_imag(mag_real);
		VecReal wght(data_n[1], 0.);
		Int idx(0);
		for_Int(i, 0, ni) {
			for_Int(j, i, ni) {
                for_Int(n, 0, nw_hterm) {
					Int re_i(bgn + idx * 2 * nw_hterm + n);
					x[re_i] = re_i;
                    y[re_i] = real(hb[p.nmesh - 1 - nw_hterm + n][i][j]);
					mag_real[idx] += SQR(y[re_i]);

					Int im_i(bgn + idx * 2 * nw_hterm + n + nw_hterm);
					x[im_i] = im_i;
                    y[im_i] = imag(hb[p.nmesh - 1 - nw_hterm + n][i][j]);
					mag_imag[idx] += SQR(y[im_i]);

					wght[idx * 2 * nw_hterm + n] = wght[idx * 2 * nw_hterm + n + nw_hterm] = (i == j) ? 1 : 2;
					//wght[re_i] = wght[im_i] = 1.;
				}
				mag_real[idx] = SQRT(1 + mag_real[idx] / nw_hterm);
				mag_imag[idx] = SQRT(1 + mag_imag[idx] / nw_hterm);
				++idx;
			}
		}
		wght *= INV(SUM(wght));
		idx = 0;
		for_Int(i, 0, ni) {
			for_Int(j, i, ni) {
				for_Int(n, 0, nw_hterm) {
					sig[bgn + idx * 2 * nw_hterm + n] = 2E-3 * mag_real[idx] / SQRT(wght[idx * 2 * nw_hterm + n]);
					sig[bgn + idx * 2 * nw_hterm + nw_hterm + n] = 2E-3 * mag_imag[idx] / SQRT(wght[idx * 2 * nw_hterm + nw_hterm + n]);
                    // sig[bgn + idx * 2 * nw_hterm + n] *= std::pow(10, 2 + (nb / 2 - 4) / 11) * sqrt(nb) * 4 * (bw);
                    // sig[bgn + idx * 2 * nw_hterm + nw_hterm + n] *= std::pow(10, 2 + (nb / 2 - 4) / 11) * sqrt(nb) * 4 * (bw);
				}
				++idx;
			}
		}

        // if (rank == 0) WRN(NAV3(mag_real, mag_imag, wght));
	}

	{ // the part of ose regularization
		Int bgn = data_n[0] + data_n[1];
		for_Idx(i, 0, nb) {
			x[i + bgn] = i + bgn;
			y[i + bgn] = 0.;
			//sig[i + bgn] = std::pow(4 * 256 * p.fit_max_omg, 2);  //Change to the part of ose regularization
			sig[i + bgn] = 1E2 * std::pow(10, 2 + (nb / 2 - 4) / 11) * sqrt(nb) * 4 * (bw);
			// if(rank == 0 && i == 0)  WRN(NAV(sig[i+bgn])) 
			// sig[i + bgn] =0.0000001;
		}
	}

    // if (rank == 0) WRN(NAV4(nw,x, y, sig));
    
}

MatReal HybErr::a2B(const VecReal& a)
{
    // Int rank;
    // MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MatReal B(ni + 1, nb);
	MatReal the = a.truncate(0, ni * nb).mat(ni, nb);
    MatReal eps = a.truncate(ni * nb, ni * nb + nb).mat(1, nb);
	for_Int(c, 0, nb) {
		for_Int(r, 0, ni) {
			B[r][c] = the[r][c];
		}
        B[ni][c] = eps[0][c];
	}
    return B;
}

void HybErr::exchange_idx(Mat<MatCmplx>& mat) const {
    Mat<MatCmplx> mat_t(mat[0][0].nrows(), mat[0][0].ncols(), MatCmplx(mat.nrows(), mat.ncols()));
    for_Int(i, 0, mat.nrows()) {
        for_Int(j, 0, mat.ncols()) {
            for_Int(m, 0, mat[i][j].nrows()) {
                for_Int(n, 0, mat[i][j].ncols()) {
                    mat_t[m][n][i][j] = mat[i][j][m][n];
                }
            }
        }
    }
    mat.reset(mat_t);
}


VecCmplx HybErr::get_dyda(const MatCmplx& dydB) {
    Vec<MatReal> dBda = get_dBda();
    VecCmplx dyda(dBda.size());
    for_Int(m, 0, dBda.size()) {
        dyda[m] = (dydB * cmplx(dBda[m].ct())).trace();
    }
    return dyda;
}

Vec<MatReal> HybErr::get_dBda() const
{
    Int i = 0;
    Vec<MatReal> dBda(a_old.size(), MatReal(ni + 1, nb, 0.));
	for_Int(r, 0, ni) {
        for_Int(c, 0, nb) {
            dBda[i][r][c] = 1.;
			i++;
		}
	}
    for_Int(c, 0, nb) {
        dBda[i][ni][c] = 1.;
        i++;
    }
    return dBda;
}

Real HybErr::operator()(const Int x, const VecReal& a) {
    VecReal dyda(a.size(), 0.);

    if (a == a_old) {
    }
    else {
        if (x != 0) WRN("aaaa")
            a_old = a;
        get_value_mrq(a);
    }
    dyda = vdyda[x];

    return vy[x];
}

void HybErr::operator()(const Int x, const VecReal& a, Real& y_i, VecReal& dyda_i)
{
    // Int rank;
    // MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (a == a_old) {
    }
    else {
        a_old = a;
        get_value_mrq(a);
    }
    y_i = vy[x];
    dyda_i = vdyda[x];
}

void HybErr::get_value_mrq(const VecReal& a)
{
    // Int rank;
    // MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MatReal B = a2B(a);
    MatCmplx V = cmplx(B.truncate(0, 0, ni, nb));
    MatCmplx E = cmplx(dmat(B[ni]));
    //[nw:nw+nw_hterm) is bsr, [nw+nw_hterm] is reg
    Vec<MatCmplx> hb(nw + nw_hterm + 1);
    // idx of z; idx of mat hb; idx of mat B
    Vec<Mat<MatCmplx>> dydB(nw + nw_hterm + 1);
    for_Int(n, 0, nw + nw_hterm + 1) {
        Mat<MatCmplx> dydV;
        Mat<MatCmplx> dydE;
        if (n < nw) {                   // the part of lower ev term of hybber curve
            Cmplx tem_z = p.fit_points.empty() ? p.Imz(n):p.Imz(p.fit_points[n]);
            MatCmplx Z = dmat(E.nrows(), tem_z);
            MatCmplx S = matinvlu(Z - E);
            MatCmplx hb_n = V * S * V.ct();
            dydV = Mat<MatCmplx>(V.nrows(), V.ncols(), MatCmplx(ni, ni, 0.));
            for_Int(r, 0, V.nrows()) {
                for_Int(c, 0, V.ncols()) {
                    MatCmplx dVdv(V.nrows(), V.ncols(), 0.);
                    dVdv[r][c] = 1.;
                    dydV[r][c] = dVdv * S * V.ct() + V * S * dVdv.ct();
                }
            }
            exchange_idx(dydV);

            dydE = Mat<MatCmplx>(1, E.ncols(), MatCmplx(ni, ni, 0.));
            for_Int(r, 0, E.nrows()) {
                MatCmplx dEde(E.nrows(), E.ncols(), 0.);
                dEde[r][r] = 1.;
                dydE[0][r] = V * S * dEde * S * V.ct();
            }
            exchange_idx(dydE);

            hb[n].reset(hb_n);
        }
        else if (n < nw + nw_hterm) {	// the part of hight ev term of hybber curve
            MatCmplx Z = dmat(E.nrows(), p.Imz(p.nmesh - 1 - nw_hterm - nw + n));
            MatCmplx S = matinvlu(Z - E);
            MatCmplx hb_n = V * S * V.ct();
            dydV = Mat<MatCmplx>(V.nrows(), V.ncols(), MatCmplx(ni, ni, 0.));
            for_Int(r, 0, V.nrows()) {
                for_Int(c, 0, V.ncols()) {
                    MatCmplx dVdv(V.nrows(), V.ncols(), 0.);
                    dVdv[r][c] = 1.;
                    dydV[r][c] = dVdv * S * V.ct() + V * S * dVdv.ct();
                }
            }
            exchange_idx(dydV);

            dydE = Mat<MatCmplx>(1, E.ncols(), MatCmplx(ni, ni, 0.));
            for_Int(r, 0, E.nrows()) {
                MatCmplx dEde(E.nrows(), E.ncols(), 0.);
                dEde[r][r] = 1.;
                dydE[0][r] = V * S * dEde * S * V.ct();
            }
            exchange_idx(dydE);

            hb[n].reset(hb_n);
        }
		/*
        else if (n < nw + 1) {
            // bath sum rule
            MatCmplx bsr = V * V.ct();
            dydV = Mat<MatCmplx>(V.nrows(), V.ncols(), MatCmplx(ni, ni, 0.));
            for_Int(r, 0, V.nrows()) {
                for_Int(c, 0, V.ncols()) {
                    MatCmplx dVdv(V.nrows(), V.ncols(), 0.);
                    dVdv[r][c] = 1.;
                    dydV[r][c] = dVdv * V.ct() + V * dVdv.ct();
                }
            }
            exchange_idx(dydV);
            dydE = Mat<MatCmplx>(ni, ni, MatCmplx(1, E.ncols(), 0.));
            hb[n].reset(bsr);
        }
		*/
        else {            // regularization,y=E*exp(0.5*(E/bw)^2)
            MatCmplx reg(1, E.ncols());
            reg[0] = E.diagonal();
            dydE = Mat<MatCmplx>(1, E.ncols(), MatCmplx(1, E.ncols(), 0.));
            Vec<Cmplx> reg_reg = Cmplx(INV(bw)) * E.diagonal();                 //regularized onsite energy
            for_Int(r, 0, E.ncols()) {
                reg[0][r] = reg[0][r] * EXP(0.5 * SQR(reg_reg[r]));
            }
            Vec<Cmplx> temp(E.ncols(), 0.);
            for_Int(r, 0, E.nrows()) {
                MatCmplx dEde(1, E.ncols(), 0.);
                dEde[0][r] = 1.;
                dydE[0][r] = (1. + SQR(reg_reg[r])) * EXP(0.5 * SQR(reg_reg[r])) * dEde;
            }
            exchange_idx(dydE);
            dydV = Mat<MatCmplx>(1, E.ncols(), MatCmplx(V.nrows(), V.ncols(), 0.));
            hb[n].reset(reg);
        }

        dydB[n].reset(hb[n].nrows(), hb[n].ncols());
        for_Int(r, 0, hb[n].nrows()) {
            for_Int(c, 0, hb[n].ncols()) {
                dydB[n][r][c].reset(concat(dydV[r][c].vec(), dydE[r][c].vec()).mat(ni + 1, nb));
            }
        }
    }

    fill2array(hb, dydB);
    // if (rank == 0) WRN("HERE is fine " + NAV3(B, V, E))
}


void HybErr::fill2array(const Vec<MatCmplx>& hb, const Vec<Mat<MatCmplx>>& dydB)
{
    // Int rank;
    // MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // vy
    Int idx(0);
    for_Int(r, 0, ni) {
        for_Int(c, r, ni) {
            for_Int(n, 0, nw) {
                vy[idx] = real(hb[n][r][c]);
                vy[idx + nw] = imag(hb[n][r][c]);
                VecCmplx dyda = get_dyda(dydB[n][r][c]);
                vdyda[idx] = real(dyda);
                vdyda[idx + nw] = imag(dyda);
                idx++;
            }
            idx += nw;
        }
    }
    for_Int(r, 0, ni) {
        for_Int(c, r, ni) {
            for_Int(n, nw, nw_hterm+nw) {
                vy[idx] = real(hb[n][r][c]);
                vy[idx + nw_hterm] = imag(hb[n][r][c]);
                VecCmplx dyda = get_dyda(dydB[n][r][c]);
                vdyda[idx] = real(dyda);
                vdyda[idx + nw_hterm] = imag(dyda);
                idx++;
            }
            idx += nw_hterm;
        }
    }
    for_Int(c, 0, nb) {
        vy[idx] = real(hb[nw + nw_hterm][0][c]);
        VecCmplx dyda = get_dyda(dydB[nw + nw_hterm][0][c]);
        vdyda[idx] = real(dyda);
        idx++;
    }
}


// ------------------------------------------------------------------🪦 code grave 🪦------------------------------------------------------------------
        /* old regularization before 20241216
        else {                          // regularization
            MatCmplx reg(1, E.ncols());
            reg[0] = E.diagonal();
            dydE = Mat<MatCmplx>(1, E.ncols(), MatCmplx(1, E.ncols(), 0.));
            for_Int(r, 0, E.nrows()) {
                MatCmplx dEde(1, E.ncols(), 0.);
                dEde[0][r] = 1.;
                dydE[0][r] = dEde;
            }
            exchange_idx(dydE);
            dydV = Mat<MatCmplx>(1, E.ncols(), MatCmplx(V.nrows(), V.ncols(), 0.));
            hb[n].reset(reg);
        }
        */