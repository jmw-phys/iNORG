#pragma once

/*
code developed by
    Jia-Ming Wang (jmw@ruc.edu.cn, RUC, China) date 2020 - 2024
*/


#include "miblas.h"
#include "milapacke.h"		// to get the trd_heevr_vd(diagonalize Tridiagonal matrix function)
#include "stdpfx.h"
#include "vec.h"
#include "mat.h"
#include "random.h"
#include "mympi.h"

template<typename T>
bool residual_method_for_b(VEC<T> ltd, VEC<T> lt_sd, Real b, Int k, bool fast_mode = false)
{
    using namespace std;
    Mat<T> ev(ltd.size(), ltd.size(), 0.);
    Vec<T> va(ltd);
    Vec<T> vb(lt_sd);
    Int info = trd_heevr_qr(va, vb, ev); if (info != 0)ERR("the i-th parameter had an illegal value." + NAV3(info, va, vb));
    if ((fast_mode) && ABS(b * ev[k][0]) < 1.E-8) return true;
    if (!(fast_mode) && ABS(b * ev[k][0]) < 1.E-9) return true;
    return false;
}

template<typename T>
inline Real compare_error(T& pre, T& lst) { return (2 * ABS(lst - pre)) / (ABS(pre) + ABS(lst) + 1.); }

template<typename T>
void orthogonalization_one_to_multiple(const MyMpi& mm, Vec<T>& v, const VEC<Vec<T>>& m, Int e)
//	orthogonalization v to m[0] through m[e - 1], which are assumed to be orthonormal
{
    Real mag0, mag1;
    mag1 = sqrt(mm.Allreduce(real(DOT(v, v))));		// mag1 != 0 is supposed
    mag0 = mag1 * 1.E+16;
    while (mag1 / mag0 < 0.98) {
        for_Int(i, 0, e) v -= mm.Allreduce(DOT(m[i], v)) * m[i];
        mag0 = mag1;
        mag1 = sqrt(mm.Allreduce(real(DOT(v, v))));		// mag1 != 0 is supposed
    }
}

template<typename T, typename F>
VecInt lanczos(VecReal& evals, Mat<T>& evecs, Int& gs_dgcy, Idx n, Int evals_size, F& H,
    VecReal& inital_state, const MyMpi& mm, bool fast_mode = false, Int max_krylov = 9999)
    //	Lanczos algorithm to find several lowest eigenpairs of a Hermitian matrix
    //  only typename F(matrix) H needed with well define operator*(which to returns H * |p>)
    //	typename T can only be Real or Cmplx(NOT USED/TESTED)
    //Output:
    //  eval                is the eigenvalues                              (Dynamic size)
    //  evec                is the eigenvectors, one row contains one vector(Dynamic nrow)
    //  gs_dgcy             ground states's degeneracy
    //Input:
    //  n                   is the dimension of the H's space.
    //  evals_size          is the fix eigenstates to be found.
    //  H                   typename F(matrix) H needed with well define operator*(which to returns H * |p>)
    //  inital_state        (not necessary) if set it as blank, Then we will generate a random vector.
    //  max_krylov          (not necessary) is the maximum size of the Krylov space

{
    // if (fast_mode) PIO("Dear customer, you are now in the fast mode, in this mode you only can find the ONE eigen pair.");
#ifdef _CHECK_DIMENSION_MATCH_
    //ASSERT_EQ(n, evec[].());
#endif
    using namespace std;
    Int e(0);
    gs_dgcy = 1;
    VEC<VecReal> evec;
    VEC<Real> eval;
    VEC<Int> Krylovsize;
    VecPartition vp(mm.np(), mm.id(), inital_state.size());
    while (e < evals_size)
    {
        // for tridiagonal matrix
        VEC<Real> ltd;	    // diagonal elements / eigenvalues
        VEC<Real> lt_sd;	// sub-diagonal elements
        // prepare an initial Lanczos vector
        VecReal inital_state_p(vp.len());
        for_Int(i, 0, vp.len()) inital_state_p[i] = inital_state[i + vp.bgn()];
        Real norm_inv = INV(SQRT(mm.Allreduce(DOT(inital_state_p, inital_state_p))));
        VecReal temp = norm_inv * inital_state_p;
        if (e > 0) {
            orthogonalization_one_to_multiple(mm, temp, evec, e);
            temp *= INV(SQRT(mm.Allreduce(DOT(temp, temp))));
        }
        VecReal v0 = temp;
        VecReal v1(H * v0);
        Idx k = 0;
        Real a(0.), b(0.);
        while (true)
        {
            a = mm.Allreduce(DOT(v1, v0));
            ltd.push_back(a);
            if (!(fast_mode)) if ( k >= 60 ) if (residual_method_for_b(ltd, lt_sd, b, k, fast_mode)) break;
            if ((fast_mode))  if ( k >= 60 ) if (residual_method_for_b(ltd, lt_sd, b, k, fast_mode)) break;
            if (k >= max_krylov) {
                if(mm) WRN("Getting Wrong with someting? krylov was too large!" + NAV2(e, k));
                if (k >= 20 + max_krylov) break;
            }
            k++;
            v1 -= a * v0;
            b = SQRT(mm.Allreduce(DOT(v1, v1)));
            lt_sd.push_back(b);
            if (e > 0) orthogonalization_one_to_multiple(mm, v1, evec, e);
            v1 *= INV(SQRT(mm.Allreduce(DOT(v1, v1)))); //v1 *= (1 / b);
            SWAP(v0, v1);
            v1 *= -b;
            v1 += H * v0;
        }
        // WRN(NAV(k) + "The iteration" + NAV(e));
        VecReal va(ltd), vb(lt_sd);
#ifdef _ASSERTION_
        if (va.size() - 1 != vb.size())ERR("Wrong with lanczos" + NAV3(e, va.size(), vb.size()));
#endif
        MatReal ev(ltd.size(), ltd.size(), 0.);
        trd_heevr_qr(va, vb, ev);
        VecReal gs_kvsp(ev.tr()[0]);                //ground state at Krylov space.
        VecReal gs(v0.size(), 0.);
        v0 = temp;
        gs += gs_kvsp[0] * v0;
        v1 = H * v0;
        k = 1;
        while (true)
        {
            v1 -= ltd[k - 1] * v0;
            if (e > 0) orthogonalization_one_to_multiple(mm, v1, evec, e);
            // for_Idx(i, 0, e) v1 -= DOT(evec[i], v1) * evec[i];
            if (k == gs_kvsp.size() - 1) break;
            v1 *= INV(SQRT(mm.Allreduce(DOT(v1, v1))));
            gs += gs_kvsp[k] * v1;
            SWAP(v0, v1);
            v1 *= -lt_sd[k - 1];
            v1 += H * v0;  k++;
        }
        Krylovsize.push_back(k);
        evec.push_back(gs);
        eval.push_back(va[0]);
        e++;
        if (fast_mode) break;
        if (e == evals_size) {// verify if finding all the degeneracy
            Int ndeg = 1;
            for_Int(i, 1, eval.size()) if (compare_error(eval[0], eval[i]) < 2E-6) ndeg++;
            if (evals_size == ndeg) evals_size++;
        }
    }
    evecs.reset(evec.size(),n);
    for_Int(i, 0, evec.size()) evecs[i] = mm.Allgatherv(evec[i], vp);
    evals.reset(eval);

    slctsort(evals, evecs);
    for_Int(i, 1, evals.size()) if(compare_error(evals[0], evals[i]) < 1E-6) gs_dgcy++;
    // if(mm) WRN(NAV(evals))    
    VecInt krylov_size(Krylovsize);
    return krylov_size;
}