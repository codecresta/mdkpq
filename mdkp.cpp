#include <cmath>
#include <climits>
#include <stdlib.h>
#include "mdkp.h"
#include "srt.h"

int calcMs(millisecs dur)
{
    return dur.count();
}

int rndU(int lo, int up)
/* Returns a random integer:
 * lo = lower value,
 * up = upper value,
 * approximately uniform for up - lo << RAND_MAX.
 */
{
    return lo + (int)((double)rand()/(double)RAND_MAX*(up - lo + 1));
}

double rndDbl()
/* Returns a random number between 0 and 1. */
{
    return (double)rand()/RAND_MAX;
}

int MdKP::invIdx(int j)
/* Returns the inverse for the base class. */
{
    return j;
}

bool MdKP::rdInt(string txt, int &i)
/* Reads an integer from the file stream. */
{
    if (pfs.eof())
    {
        cout << "Error reading " << pfn << ": " << txt << "." << endl;
        return false;
    }
    pfs >> i;
    return true;
}

bool MdKP::clsFl(bool res)
/* Closes a file and returns the value passed in. */
{
    pfs.close();
    return res;
}

bool MdKP::rdPrb(string fn0)
/* Reads a problem from a data file. */
{
    int i, j;
    pfn = fn0;
    pfs.open(pfn);
    if (pfs)
    {
        if (rdInt("number of items", itms) && rdInt("number of dimensions", dims) && rdInt("optimal value", opt))
        {
            for (j = 0; j < itms; j++)
                if (!rdInt("value", val[j]))
                    return clsFl(false);
            for (i = 0; i < dims; i++)
                for (j = 0; j < itms; j++)
                    if (!rdInt("weight", wei[i][j]))
                        return clsFl(false);
            for (i = 0; i < dims; i++)
                if (!rdInt("capacity", cap[i]))
                    return clsFl(false);
            return clsFl(true);
        }
        else
            return clsFl(false);
    }
    else
        cout << "Error opening file " << pfn << "." << endl;
    return false;
}

void MdKP::uncrlPrb(int dims0, int itms0, int cap_lo, int cap_up, int val_lo, int val_up, int wei_lo, int wei_up)
/* Creates a random uncorrelated (weights to profits) problem from the parameters passed. */
{
    int i, j;
    dims = dims0;
    itms = itms0;
    for (i = 0; i < dims; i++)
        cap[i] = rndU(cap_lo, cap_up);
    for (j = 0; j < itms; j++)
    {
        val[j] = rndU(val_lo, val_up);
        for (i = 0; i < dims; i++)
            wei[i][j] = rndU(wei_lo, wei_up);
    }
}

void MdKP::crlPrb(int dims0, int itms0, int cap_lo, int cap_up, int wei0, double rng)
/* Creates a random correlated problem from the parameters passed. */
{
    int i, j, v_rng;
    dims = dims0;
    itms = itms0;
    for (i = 0; i < dims; i++)
        cap[i] = rndU(cap_lo, cap_up);
    v_rng = (int)(wei0*rng);
    for (j = 0; j < itms; j++)
    {
        val[j] = rndU(wei0 - v_rng, wei0 + v_rng);
        for (i = 0; i < dims; i++)
            wei[i][j] = rndU(1, wei0);
    }
}

void MdKP::outPrbSmry()
/* Outputs a problem summary. */
{
    cout << "m = " << dims << ", n = " << itms << endl;
}

bool MdKP::feasSol()
/* Check solution is feasible. */
{
    int i, j, c;
    for (i = 0; i < dims; i++)
    {
        c = cap[i];
        for (j = 0; j < itms; j++)
            if (sol[j] == 1)
                c -= wei[i][j];
        if (c < 0)
            return false;
    }
    return true;
}

int MdKP::vSol()
/* Recalculates solution value. */
{
    int j, r;
    r = 0;
    for (j = 0; j < itms; j++)
        if (sol[j] == 1)
            r += val[j];
    return r;
}

void MdKP::chk()
/* Checks that the solution is feasible and that the value is correct. */
{
    if (!feasSol())
        cout << "infeasible" << endl;
    if (v_sol != vSol())
        cout << "mismatch" << endl;
}

int MdKPB::varIdx(int j)
/* Variable index for base class. */
{
    return j;
}

void MdKPB::outPrb()
/* Outputs a problem. */
{
    int i, j;
    cout << "m = " << dims << ", n = " << itms << endl;
    for (j = 0; j < itms; j++)
    {
        cout << "p[" << j << "] = " << val[varIdx(j)];
        if (j < itms - 1)
            cout << ", ";
    }
    cout << endl;
    for (i = 0; i < dims; i++)
    {
        for (j = 0; j < itms; j++)
            cout << "w[" << i << ", " << j << "] = " << wei[i][varIdx(j)] << ", ";
        cout << "c[" << i << "] = " << cap[i] << endl;
    }
}

void MdKPB::outSol()
/* Outputs solution:
 * i = item index,
 * aft = flag set to true after first output.
*/
{
    int i;
    bool aft;
    if (sty == STY_NONE)
        return;
    cout << "solution:" << endl;
    aft = false;
    for (i = 0; i < itms; i++)
    {
        if (sty == STY_BAS || sty == STY_BAS_X)
        {
            if (sol[invIdx(i)] == 1)
            {
                if (aft & !nwl)
                    cout << ", ";
                if (sty == STY_BAS_X)
                    cout << "x";
                cout << i + 1;
                if (nwl)
                    cout << endl;
                aft = true;
            }
        }
        else
        {
            if (aft & !nwl)
                cout << ", ";
            if (sty == STY_VAL_X)
                cout << "x" << i + 1 << " = ";
            cout << sol[invIdx(i)];
            if (nwl)
                cout << endl;
            aft = true;
        }
    }
    if (!nwl)
        cout << endl;
}
