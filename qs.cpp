#include <climits>
#include "ClpSimplex.hpp"
#include "mdkp.h"
#include "srt.h"

bool MdKPQ::prep(MdKP &src)
/* Uses Coin-OR CLP to solve the linear relaxation, uses the dual solution (shadow prices) to calculate the utility
 * value of each item and finally sorts the items by decreasing utility value:
 * n = number of elements in LP matrix,
 * elts = element values in LP matrix,
 * obj = objective values in LP (item values in MdKP),
 * col_lwr = column lower values,
 * col_upr = column upper values,
 * row_upr = row upper values,
 * srg = surrogate weights (from dual solution),
 * utl = utility values of each item,
 * dl_row = dual (solution) row values,
 * sum = a summation for calculating utility values,
 * row_idc = row indices for LP matrix,
 * col_idc = column indices for LP matrix,
 * i, j, k = index variables.
*/
{
	ClpSimplex  mdl;
	int n;
	dims = src.dims;
	itms = src.itms;
	opt = src.opt;
	n = dims*itms;
	double elts[n], obj[itms], col_lwr[itms], col_upr[itms], row_upr[dims], utl[itms], *dl_row, sum;
	int row_idc[n], col_idc[n], idx[itms], i, j, k;
	for (i = 0; i < dims; i++)
		cap[i] = src.cap[i];
	k = 0;
	for (i = 0; i < dims; i++)
	{
		row_upr[i] = cap[i];
		for (j = 0; j < itms; j++)
		{
			elts[k] = src.wei[i][j];
			row_idc[k] = i;
			col_idc[k++] = j;
		}
	}
	for (j = 0; j < itms; j++)
	{
		obj[j] = src.val[j];
		col_lwr[j] = 0;
		col_upr[j] = 1;
	}
	CoinPackedMatrix mat(false, row_idc, col_idc, elts, n);
	mdl.setLogLevel(0);
	mdl.loadProblem(mat, col_lwr, col_upr, obj, NULL, row_upr);
	mdl.setOptimizationDirection(-1);
	mdl.primal();
	if (mdl.isProvenOptimal())
	{
        dl_row = mdl.dualRowSolution();
        for (j = 0; j < itms; j++)
        {
            sum = 0;
            for (i = 0; i < dims; i++)
                sum += dl_row[i]*src.wei[i][j];
            idx[j] = j;
            utl[j] = src.val[j]/sum;
        }
        qSrt(idx, utl, 0, itms - 1);
        for (j = 0; j < itms; j++)
        {
            k = idx[j];
            val[j] = src.val[k];
            for (i = 0; i < dims; i++)
                wei[i][j] = src.wei[i][k];
            inv[k] = j;
        }
        return true;
    }
    else
    {
        cout << "error LP: " << endl;
        if (mdl.isProvenPrimalInfeasible())
            cout << "proven primal infeasible" << endl;
        if (mdl.isProvenDualInfeasible())
            cout << "proven dual infeasible" << endl;
        if (mdl.isIterationLimitReached())
            cout << "iteration limit reached" << endl;
        return false;
    }
}

int MdKPQ::calV(int x[])
/* Calculates value of a solution. */
{
    int j, v;
    v = 0;
    for (j = 0; j < itms; j++)
        if (x[j] == 1)
            v += val[j];
    return v;
}

void MdKPQ::cpy(int x[], int y[])
/* Copies a solution. */
{
    int j;
    for (j = 0; j < itms; j++)
        y[j] = x[j];
}

void MdKPQ::rmc(int x[], int r[])
/* Calculates the remaining capacity. */
{
    int i, j;
    for (i = 0; i < dims; i++)
        r[i] = cap[i];
    for (j = 0; j < itms; j++)
        if (x[j] == 1)
            for (i = 0; i < dims; i++)
                r[i] -= wei[i][j];
}

void MdKPQ::cpyR(int s[], int t[])
/* Copies the remaining capacity. */
{
    int i;
    for (i = 0; i < dims; i++)
        t[i] = s[i];
}

void MdKPQ::cpyD(double s[], double t[])
/* Copies a solution. */
{
    int j;
    for (j = 0; j < itms; j++)
        t[j] = s[j];
}

bool MdKPQ::inf(int r[])
/* Returns true if infeasible e.g. one or more resource is overused. */
{
    int i;
    bool flg;
    flg = false;
    for (i = 0; i < dims; i++)
        if (r[i] < 0)
        {
            flg = true;
            break;
        }
    return flg;
}

bool MdKPQ::drp(int k, int x[], int r[], int &v)
/* Drop subroutine. */
{
    int i;
    bool flg;
    if (x[k] == 1)
    {
        x[k] = 0;
        for (i = 0; i < dims; i++)
            r[i] += wei[i][k];
        v -= val[k];
        flg = true;
        for (i = 0; i < dims; i++)
            if (r[i] < 0)
            {
                flg = false;
                break;
            }
        return flg;
    }
    else
        return false;
}

void MdKPQ::add(int k, int x[], int r[], int &v)
/* Add subroutine. */
{
    int i;
    bool flg;
    if (x[k] == 0)
    {
        flg = true;
        for (i = 0; i < dims; i++)
            if (wei[i][k] > r[i])
            {
                flg = false;
                break;
            }
        if (flg)
        {
            x[k] = 1;
            for (i = 0; i < dims; i++)
                r[i] -= wei[i][k];
            v += val[k];
        }
    }
}

void MdKPQ::lcl(int x[], int &v)
/* Local search. */
{
    int xl[N_MX], xd[N_MX], r[M_MX], rl[M_MX], rd[M_MX], vl, vd, i, j, k;
    bool imp;
    cpy(x, xl);
    rmc(x, r);
    vl = v;
    imp = true;
    while (imp)
    {
        imp = false;
        for (j = 0; j < itms; j++)
        {
            // TODO - could factor some of this:
            cpy(x, xd);
            cpyR(r, rd);
            vd = v;
            if (xd[j] == 1)
            {
                xd[j] = 0;
                for (i = 0; i < dims; i++)
                    rd[i] += wei[i][j];
                vd -= val[j];
                // add phase:
                for (k = 0; k < itms; k++)
                    if (k != j)
                        add(k, xd, rd, vd);
            }
            else
            {
                xd[j] = 1;
                for (i = 0; i < dims; i++)
                    rd[i] -= wei[i][j];
                vd += val[j];
                // drop phase:
                if (inf(rd))
                    for (k = itms - 1; k >= 0; k--)
                        if (k != j && drp(k, xd, rd, vd))
                            break;
            }
            if (vd > vl)
            {
                cpy(xd, xl);
                cpyR(rd, rl);
                vl = vd;
                imp = true;
            }
        }
        if (imp)
        {
            cpy(xl, x);
            cpyR(rl, r);
            v = vl;
        }
    }
}

void MdKPQ::algS(int s[], int t[], int &vi)
/* Algorithm subroutine. */
{
    int r[M_MX], k, v;
    v = calV(s);
    rmc(s, r);
    if (inf(r))
        for (k = itms - 1; k >= 0; k--)
            if (drp(k, s, r, v))
                break;
    for (k = 0; k < itms; k++)
        add(k, s, r, v);
    if (v > vi)
    {
        lcl(s, v);
        cpy(s, t);
        vi = v;
        if (v > v_sol)
        {
            cpy(s, sol);
            v_sol = v;
        }
    }
}

void MdKPQ::alg()
/* Main algorithm. */
{
    tm_pt t0;
    double t, ys[S_MX][N_MX], yt[S_MX][N_MX], yh[N_MX];
    int xs[S_MX][N_MX], xt[S_MX][N_MX], vs[S_MX], i, j, it;
    bet = 1 - alp;
    ep3 = 1 - ep1 - ep2;
    v_sol = -INT_MAX;
    for (i = 0; i < sz; i++)
    {
        for (j = 0; j < itms; j++)
        {
            t = rndDbl();
            ys[i][j] = t;
            xs[i][j] = t < rndDbl() ? 1 : 0;
        }
        vs[i] = -INT_MAX;
        algS(xs[i], xt[i], vs[i]);
    }
    it = 0;
    while (it < mx)
    {
        for (j = 0; j < itms; j++)
            yh[j] = alp*sol[j] + bet*(1 - sol[j]);
        for (i = 0; i < sz; i++)
        {
            for (j = 0; j < itms; j++)
            {
                yt[i][j] = alp*xt[i][j] + bet*(1 - xt[i][j]);
                t = ep1*ys[i][j] + ep2*yt[i][j] + ep3*yh[j];
                ys[i][j] = t;
                xs[i][j] = t < rndDbl() ? 1 : 0;
            }
            algS(xs[i], xt[i], vs[i]);
        }
        it++;
    }
}

void MdKPQ::info()
{
    cout << name() << endl;
}

void MdKPQ::ini(int sz0, int mx0, double alp0, double ep10, double ep20, bool out, int sty0, bool nwl0)
/* Initialise. */
{
    sz = sz0;
    mx = mx0;
    alp = alp0;
    ep1 = ep10;
    ep2 = ep20;
    sty = sty0;
    nwl = nwl0;
    if (out)
        cout << "sz = " << sz << ", mx = " << mx << ", alp = " << alp << ", ep1 = " << ep1 << ", ep2 = " << ep2 << endl;
}

int MdKPQ::slv(MdKP &src, bool out)
/* Initialises and solves a problem. */
{
    tm_pt t0, t1;
    t0 = timer::now();
    prep(src);
    alg();
    t1 = timer::now();
    tm = calcMs(duration(t1 - t0));
    if (out)
    {
        cout << "z = " << v_sol << (v_sol == opt ? " (opt)" : "") << ", tm = " << tm << " ms" << endl;
        outSol();
    }
    return v_sol;
}
