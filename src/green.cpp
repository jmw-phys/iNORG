#include "green.h"

/*
code developed and maintained by (jmw@ruc.edu.cn, RUC, China) date 2022 - 2024
*/

using namespace std;

void Green::read(const Str& file) {
	IFS ifs(file);
	if (!ifs) {
		WRN(STR("file opening failed with ") + NAV(file))
	}
	else {
		// for_Int(i, 0, (norbs * norbs + 1) * 2) {
		for_Int(i, 0, (norbs + 1) * 2) {
			Str strr;
			ifs >> strr;
		}
		Real omg;
		MatReal re(norbs, norbs, 0.);
		MatReal im(norbs, norbs, 0.);
		for_Int(i, 0, nomgs) {
			ifs >> omg;
			for_Int(m, 0, norbs) {
				// for_Int(n, 0, norbs) {
				Int n = m;
				ifs >> re[m][n];
				// }
			}
			ifs >> omg;
			for_Int(m, 0, norbs) {
				// for_Int(n, 0, norbs) {
				Int n = m;
				ifs >> im[m][n];
				// }
			}
			g[i] = cmplx(re, im);
		}
		if (!ifs) {
			ERR(STR("read-in error with ") + NAV(file))
		}
	}
	ifs.close();
}


Green::Green(const Prmtr& p, Int nomgs_i, const VecCmplx& z_omg_i, const Str& file) :
	nomgs(nomgs_i), norbs(p.nband), z_omg(z_omg_i), g(nomgs, MatCmplx(norbs, norbs, 0.)) {
	if (!file.empty()) {
		read(file);
	}
}

Green::Green(Int norbs_i, const Prmtr& p, Int nomgs_i, const VecCmplx& z_omg_i, const Str& file) :
	nomgs(nomgs_i), norbs(norbs_i), z_omg(z_omg_i), g(nomgs, MatCmplx(norbs, norbs, 0.)) {
	if (!file.empty()) {
		read(file);
	}
}

Green::Green(const Green& a, Int orb0, Int orb1) :
	nomgs(a.nomgs), norbs(1), z_omg(a.z_omg), g(nomgs, MatCmplx(norbs, norbs, 0.)) {
	for_Int(i, 0, nomgs) {
		g[i] = MatCmplx(1, 1, a.g[i][orb0][orb1]);
	}
}


Green::Green(const Prmtr& p, Int orb0, Int orb1, const Green& a) :
	nomgs(a.nomgs), norbs(p.nband), z_omg(a.z_omg), g(nomgs, MatCmplx(norbs, norbs, 0.)) {
	for_Int(i, 0, nomgs) {
		g[i][orb0][orb1] = a.g[i][0][0];
	}
}

void ImGreen::write_occupation_info() const {
	using namespace std;
	OFS ofs; ofs.open("nmat.txt");
	VecReal counter(3);
	VecReal particals(particle_number().diagonal());

	ofs << "#   < n_i >   data:" << endl;
	for_Int(orb_i, 0, norbs * 2) {
		ofs << iofmt();
		std::string temp = (orb_i % 2) == 0 ? STR(Int(orb_i / 2) + 1) + "up" : STR(Int(orb_i / 2) + 1) + "dn";
		ofs << setw(6) << temp << setw(p_Real) << particals[Int(orb_i / 2)] << endl;
	}
	counter[0] = SUM(particals);
	counter[1] = SUM(particals);

	counter[2] = counter[0] + counter[1];
	ofs << setw(6) << "sup" << setw(p_Real) << counter[0] << endl;
	ofs << setw(6) << "sdn" << setw(p_Real) << counter[1] << endl;
	ofs << setw(6) << "sum" << setw(p_Real) << counter[2];
	ofs.close();
}

void Green::write(const Str& green_name, Int iter_cnt) const {
	OFS ofs;
	// if (iter_cnt == 999) { ofs.open(iox + green_name + ".txt"); }
	if (iter_cnt == 999) { ofs.open(green_name + ".out"); }
	else { ofs.open(iox + "zic" + prefill0(iter_cnt, 3) + ".mb." + green_name + ".txt"); }
	Str iter_str = iter_cnt == 999 ? "" : STR(green_name, iter_cnt) + "_";
	// OFS ofs(iox + "zic" + prefill0(iter_cnt, 3) + ".mb." + green_name + ".txt");
	ofs << iofmt("sci");
	// real part of green
	ofs << setw(w_Real) << "w";
	for_Int(m, 0, norbs) {
		// for_Int(n, 0, norbs) {
		Int n = m;
		ofs << "  " << setw(w_Real) << iter_str + "Re" + STR(m + 1) + STR(n + 1);
		// }
	}
	// imag part of green
	ofs << "  " << setw(w_Real) << "w";
	for_Int(m, 0, norbs) {
		// for_Int(n, 0, norbs) {
		Int n = m;
		ofs << "  " << setw(w_Real) << iter_str + "Im" + STR(m + 1) + STR(n + 1);
		// }
	}
	ofs << endl;
	for_Int(i, 0, nomgs) {
		// real part of green
		ofs << setw(w_Real) << omg(i);
		for_Int(m, 0, norbs) {
			// for_Int(n, 0, norbs) {
			Int n = m;
			ofs << "  " << setw(w_Real) << real(g[i][m][n]);
			// }
		}
		// imag part of green
		ofs << "  " << setw(w_Real) << omg(i);
		for_Int(m, 0, norbs) {
			// for_Int(n, 0, norbs) {
			Int n = m;
			ofs << "  " << setw(w_Real) << imag(g[i][m][n]);
			// }
		}
		ofs << endl;
	}
}


void Green::write(const Str& green_name, const Str& rowname, Int iter_cnt) const {
	OFS ofs;
	// if (iter_cnt == 999) { ofs.open(iox + green_name + rowname + ".txt"); }
	if (iter_cnt == 999) { ofs.open(green_name + rowname + ".txt"); }
	else { ofs.open(iox + "zic" + prefill0(iter_cnt, 3) + ".mb." + green_name + rowname + ".txt"); }
	Str iter_str = iter_cnt == 999 ? STR(rowname) + "_" : STR(rowname) + STR(iter_cnt) + "_";
	// OFS ofs(iox + "zic" + prefill0(iter_cnt, 3) + ".mb." + green_name + ".txt");
	ofs << iofmt("sci");
	// real part of green
	ofs << setw(w_Real) << "w";
	for_Int(m, 0, norbs) {
		for_Int(n, 0, norbs) {
			ofs << "  " << setw(w_Real) << iter_str + "Re" + STR(m + 1) + STR(n + 1);
		}
	}
	// imag part of green
	ofs << "  " << setw(w_Real) << "w";
	for_Int(m, 0, norbs) {
		for_Int(n, 0, norbs) {
			ofs << "  " << setw(w_Real) << iter_str + "Im" + STR(m + 1) + STR(n + 1);
		}
	}
	ofs << endl;
	for_Int(i, 0, nomgs) {
		// real part of green
		ofs << setw(w_Real) << omg(i);
		for_Int(m, 0, norbs) {
			for_Int(n, 0, norbs) {
				ofs << "  " << setw(w_Real) << real(g[i][m][n]);
			}
		}
		// imag part of green
		ofs << "  " << setw(w_Real) << omg(i);
		for_Int(m, 0, norbs) {
			for_Int(n, 0, norbs) {
				ofs << "  " << setw(w_Real) << imag(g[i][m][n]);
			}
		}
		ofs << endl;
	}
}

void Green::write_edmft(const Str& green_name, const VecInt& or_deg) const {
	OFS ofs;
	ofs.open(green_name);
	ofs << iofmt("sci");
	for_Int(i, 0, nomgs) {
		Int deg_idx(0);
		// real part of green
		ofs << setw(w_Real) << omg(i);
		for_Int(m, 0, norbs) {
			if (deg_idx < or_deg[m * 2]) {
				ofs << "  " << setw(w_Real) << real(g[i][m][m]);
				ofs << "  " << setw(w_Real) << imag(g[i][m][m]);
				deg_idx++;
			}
		}
		ofs << endl;
	}
}

void Green::read_edmft(const Str& file, const VecInt& or_deg) {
    IFS ifs(file);
    if (!ifs) {
        WRN(STR("file opening failed with ") + NAV(file));
    } else {
        Str line;
        // skip the "#" line
        while (std::getline(ifs, line)) {
            if (line.empty() || line[0] != '#') {
                std::istringstream iss(line);
                Real omg;
				Int n = 0;
                MatReal re(norbs, norbs, 0.);
                MatReal im(norbs, norbs, 0.);
                VecReal re_temp(MAX(or_deg)), im_temp(MAX(or_deg));
                iss >> omg;
                for_Int(m, 0, MAX(or_deg)) {
                    iss >> re_temp[m];
                    iss >> im_temp[m];
                }
                for_Int(m, 0, norbs) {
                    re[m][m] = re_temp[or_deg[m * 2] - 1];
                    im[m][m] = im_temp[or_deg[m * 2] - 1];
                }
                g[n] = cmplx(re, im);
                break;
            }
        }

        for_Int(n, 1, nomgs) {
            Real omg;
            MatReal re(norbs, norbs, 0.);
            MatReal im(norbs, norbs, 0.);
            ifs >> omg;
            VecReal re_temp(MAX(or_deg)), im_temp(MAX(or_deg));
            for_Int(m, 0, MAX(or_deg)) {
                ifs >> re_temp[m];
                ifs >> im_temp[m];
            }
            for_Int(m, 0, norbs) {
                re[m][m] = re_temp[or_deg[m * 2] - 1];
                im[m][m] = im_temp[or_deg[m * 2] - 1];
            }
            g[n] = cmplx(re, im);
        }

        if (!ifs) {
            ERR(STR("read-in error with ") + NAV(file));
        }
    }
    ifs.close();
}

void Green::write_zen(const Str& green_name, Int nspin, Int iter_cnt) const {
	OFS ofs;
	// if (iter_cnt == 999) { ofs.open(iox + green_name + ".txt"); }
	if (iter_cnt == 999) { ofs.open(green_name + ".txt"); }
	else { ofs.open(iox + "zic" + prefill0(iter_cnt, 3) + ".mb." + green_name + ".txt"); }
	Str iter_str = iter_cnt == 999 ? "" : STR(iter_cnt) + "_";
	// OFS ofs(iox + "zic" + prefill0(iter_cnt, 3) + ".mb." + green_name + ".txt");
	ofs << iofmt("sci");
	for_Int(nspin_i, 0, nspin) {
		// real part of green
		ofs << setw(w_Real) << "w";
		for_Int(m, 0, norbs) {
			for_Int(n, 0, norbs) {
				ofs << "  " << setw(w_Real) << iter_str + "Re" + STR(m + 1) + STR(n + 1);
			}
		}
		// imag part of green
		ofs << "  " << setw(w_Real) << "w";
		for_Int(m, 0, norbs) {
			for_Int(n, 0, norbs) {
				ofs << "  " << setw(w_Real) << iter_str + "Im" + STR(m + 1) + STR(n + 1);
			}
		}
		ofs << endl;
		for_Int(i, 0, nomgs) {
			// real part of green
			ofs << setw(w_Real) << omg(i);
			for_Int(m, 0, norbs) {
				for_Int(n, 0, norbs) {
					ofs << "  " << setw(w_Real) << real(g[i][m][n]);
				}
			}
			// imag part of green
			ofs << "  " << setw(w_Real) << omg(i);
			for_Int(m, 0, norbs) {
				for_Int(n, 0, norbs) {
					ofs << "  " << setw(w_Real) << imag(g[i][m][n]);
				}
			}
			ofs << endl;
		}
		ofs << endl;						// blank line
		ofs << endl;						// blank line
	}
}


void Green::write_zen(const Str& green_name, const Str& rowname, Int nspin, Int iter_cnt) const {
	OFS ofs;
	// if (iter_cnt == 999) { ofs.open(iox + green_name + rowname + ".txt"); }
	if (iter_cnt == 999) { ofs.open(green_name + rowname + ".txt"); }
	else { ofs.open(iox + "zic" + prefill0(iter_cnt, 3) + ".mb." + green_name + rowname + ".txt"); }
	Str iter_str = iter_cnt == 999 ? STR(rowname) + "_" : STR(rowname) + STR(iter_cnt) + "_";
	// OFS ofs(iox + "zic" + prefill0(iter_cnt, 3) + ".mb." + green_name + ".txt");
	ofs << iofmt("sci");
	for_Int(nspin_i, 0, nspin) {
		// real part of green
		ofs << setw(w_Real) << "w";
		for_Int(m, 0, norbs) {
			for_Int(n, 0, norbs) {
				ofs << "  " << setw(w_Real) << iter_str + "Re" + STR(m + 1) + STR(n + 1);
			}
		}
		// imag part of green
		ofs << "  " << setw(w_Real) << "w";
		for_Int(m, 0, norbs) {
			for_Int(n, 0, norbs) {
				ofs << "  " << setw(w_Real) << iter_str + "Im" + STR(m + 1) + STR(n + 1);
			}
		}
		ofs << endl;
		for_Int(i, 0, nomgs) {
			// real part of green
			ofs << setw(w_Real) << omg(i);
			for_Int(m, 0, norbs) {
				for_Int(n, 0, norbs) {
					ofs << "  " << setw(w_Real) << real(g[i][m][n]);
				}
			}
			// imag part of green
			ofs << "  " << setw(w_Real) << omg(i);
			for_Int(m, 0, norbs) {
				for_Int(n, 0, norbs) {
					ofs << "  " << setw(w_Real) << imag(g[i][m][n]);
				}
			}
			ofs << endl;
		}
		ofs << endl;						// blank line
		ofs << endl;						// blank line
	}
}

// /*
Real ImGreen::error(const ImGreen& b, Real omg_rsd) const {
	VecReal wght(nomgs);
	for_Int(n, 0, nomgs) {
		wght[n] = INV(SQR(omg(n) + omg_rsd));
	}
	wght *= INV(SUM(wght) * (2 * norbs * norbs));

	MatReal mag_real(norbs, norbs, 0.);
	MatReal mag_imag(norbs, norbs, 0.);
	for_Int(n, 0, nomgs) {
		mag_real += SQR(real(g[n])) + SQR(real(b.g[n]));
		mag_imag += SQR(imag(g[n])) + SQR(imag(b.g[n]));
	}
	mag_real = SQRT(mag_real * INV(nomgs * 2));
	mag_imag = SQRT(mag_imag * INV(nomgs * 2));
	const Real mag = MAX(1.e-2, MAX(MAX(mag_real), MAX(mag_imag)));
	// const Real mag = MAX(1.e-2, MAX(mag_imag));

	Real sum = 0.;
	for_Int(n, 0, nomgs) {
		for_Int(i, 0, norbs) for_Int(j, 0, norbs) {
			Cmplx diff = b.g[n][i][j] - g[n][i][j];
			// Real diff = imag(b.g[n][i][j]) - imag(g[n][i][j]);
			sum += wght[n] * (norm(diff) / SQR(mag));
		}
	}
	return SQRT(sum);
}
// */

/*
Real ImGreen::error(const ImGreen& b, Real omg_rsd) const // only consider the imag part.
{
	VecReal wght(nomgs);
	for_Int(n, 0, nomgs) {
		wght[n] = INV(SQR(omg(n) + omg_rsd));
	}
	wght *= INV(SUM(wght) * (2 * norbs * norbs));

	// MatReal mag_real(norbs, norbs, 0.);
	MatReal mag_imag(norbs, norbs, 0.);
	for_Int(n, 0, nomgs) {
		// mag_real += SQR(real(g[n])) + SQR(real(b.g[n]));
		mag_imag += SQR(imag(g[n])) + SQR(imag(b.g[n]));
	}
	// // mag_real = SQRT(mag_real * INV(nomgs * 2));
	mag_imag = SQRT(mag_imag * INV(nomgs * 2));
	// const Real mag = MAX(1.e-2, MAX(MAX(mag_real), MAX(mag_imag)));
	const Real mag = MAX(1.e-2, MAX(mag_imag));

	Real sum = 0.;
	for_Int(n, 0, nomgs) {
		for_Int(i, 0, norbs) for_Int(j, 0, norbs) {
			// Cmplx diff = b.g[n][i][j] - g[n][i][j];
			Real diff = imag(b.g[n][i][j]) - imag(g[n][i][j]);
			sum += wght[n] * (norm(diff) / SQR(mag));
		}
	}
	return SQRT(sum);
}
*/

//w_n=(2n+1)*unit_omg 
//sum_n :e/(w_n^2+e^2)   =pi*tanh(e*pi/(2*unit_omg))/4*unit_omg
//(2/beta) sum_{n=0 inf}<<c_j c_i^+>>_(iw_n)

//if Re(gf) = b1/(w^2+a1^2)+b2/(w^2+a2^2) => Re(gf)=h=(xw^2+y)/(w^4+mw^2+n) =>
//w^2*x + y - hw^2*m - h*n =h*w^4,solve these equations with four different frequency   
Real ImGreen::sum(const VecCmplx& gf) const {
	Real electron_density = real(SUM(gf));
	VecReal w = imag(z_omg.truncate(nomgs - 4, nomgs));
	VecReal g = real(gf.truncate(nomgs - 4, nomgs));
	if (MAX(ABS(real(gf))) < 1.E-6) 	return 0.;					//if gf is 0, insensitive with omg, return 0
	MatReal A(4, 4);
	MatReal B(4, 1);
	for_Int(i, 0, 4) {
		w[i] = imag(z_omg[nomgs - 1 - 10 * i]);
		g[i] = real(gf[nomgs - 1 - 10 * i]);
		A[i][0] = w[i] * w[i];
		A[i][1] = 1.;
		A[i][2] = -w[i] * w[i] * g[i];
		A[i][3] = -g[i];
		B[i][0] = g[i] * w[i] * w[i] * w[i] * w[i];
	}
	Int info = gaussj(A, B);
	if (info > 0)	ERR("gaussj false" + NAV(info))
		Real x = B[0][0];
	Real y = B[1][0];
	Real m = B[2][0];
	Real n = B[3][0];
	VecReal a(2, 0.);
	VecReal b(2, 0.);
	a[0] = 0.5 * (m + SQRT(m * m - 4 * n));
	a[1] = 0.5 * (m - SQRT(m * m - 4 * n));
	b[0] = (x * a[0] - y) / (a[0] - a[1]);
	b[1] = (y - x * a[1]) / (a[0] - a[1]);
	for_Int(i, 0, 2) {
		if (a[i] < 0) { electron_density += b[i] * pi_Real * tan(SQRT(-a[i]) * pi_Real / (2. * unit_omg)) / (4 * unit_omg * SQRT(-a[i])); }
		else { electron_density += b[i] * pi_Real * tanh(SQRT(a[i]) * pi_Real / (2. * unit_omg)) / (4 * unit_omg * SQRT(a[i])); }
	}
	for_Int(n, 0, nomgs) {
		for_Int(i, 0, 2) {
			electron_density -= b[i] / (SQR(omg(n)) + a[i]);
		}
	}
	const Real temp = unit_omg / pi_Real;
	return 2 * temp * electron_density;
}

//<c_i^+ c_j> - <c_j c_i^+> = (2/beta) sum_{n=-inf inf}<<c_j c_i^+>>_(iw_n)
//mat_ij=<c_i^+ c_j>
MatReal ImGreen::particle_number() const {
	MatReal electron_density(norbs);
	for_Int(i, 0, norbs) {
		for_Int(j, 0, norbs) {
			VecCmplx gf(nomgs);
			for_Int(n, 0, nomgs) {
				gf[n] = g[n][i][j];
			}
			electron_density[j][i] = sum(gf);
		}
	}
	electron_density += dmat(norbs, 0.5);
	return electron_density;
}

/*

Green green_inverse(const Green& g)
{
	Green ginv(g.size());
	for_Int(i, 0, g.size()) {
		ginv[i] = matinvlu(g[i]);
	}
	return ginv;
}

void write_matsubara_green(const Prmtr& p, const Green& green, Str green_name, Int iter_cnt)
{
	OFS ofs(iox + "zic" + prefill0(iter_cnt, 3) + ".mb." + green_name + ".txt");
	ofs << iofmt("sci");
	// real part of green
	ofs << setw(w_Real) << "w";
	for_Int(m, 0, p.nband) {
		for_Int(n, 0, p.nband) {
			ofs << "  " << setw(w_Real) << "Re" + STR(m + 1) + STR(n + 1);
		}
	}
	// imag part of green
	ofs << "  " << setw(w_Real) << "w";
	for_Int(m, 0, p.nband) {
		for_Int(n, 0, p.nband) {
			ofs << "  " << setw(w_Real) << "Im" + STR(m + 1) + STR(n + 1);
		}
	}
	ofs << endl;
	for_Int(i, 0, p.num_omg) {
		Real omg = (2 * i + 1) * p.unit_omg;
		// real part of green
		ofs << setw(w_Real) << omg;
		for_Int(m, 0, p.nband) {
			for_Int(n, 0, p.nband) {
				ofs << "  " << setw(w_Real) << real(green[i][m][n]);
			}
		}
		// imag part of green
		ofs << "  " << setw(w_Real) << omg;
		for_Int(m, 0, p.nband) {
			for_Int(n, 0, p.nband) {
				ofs << "  " << setw(w_Real) << imag(green[i][m][n]);
			}
		}
		ofs << endl;
	}
}

*/

void write_retarded_green(const Prmtr& p, const Green& green, Str green_name, Int iter_cnt) {
	OFS ofs(iox + "zic" + prefill0(iter_cnt, 3) + ".rt." + green_name + ".txt");
	ofs << iofmt("sci");
	for_Int(i, 0, p.nfreq) {
		Real omg = p.freq_low + i * p.dlt_freq;
		// real part of green
		ofs << setw(w_Real) << omg;
		for_Int(m, 0, p.nband) {
			for_Int(n, 0, p.nband) {
				ofs << "  " << setw(w_Real) << real(green[i][m][n]);
			}
		}
		// imag part of green
		ofs << "  " << setw(w_Real) << omg;
		for_Int(m, 0, p.nband) {
			for_Int(n, 0, p.nband) {
				ofs << "  " << setw(w_Real) << imag(green[i][m][n]);
			}
		}
		ofs << endl;
	}
}
