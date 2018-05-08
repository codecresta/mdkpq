/* Multidimensional Knapsack Problems by Daniel W. Grace
 *
 * Hybrid Quantum Particle Swarm Optimization Algorithm, see reference below!
 *
 * Notes:
 * - Outputs "z" as the objective value (profit)
 * - Uses the Coin-OR library
 * - Uses a command line interface
 *
 * TODO:
 * Could factor some code, see comments!
 *
 * Daniel W. Grace, email: danwgrace@gmail.com
 *
 * Reference:
 * Boukthir Haddar, Mahdi Khemakhem, Saïd Hanafi, Christophe Wilbaut
 * A hybrid quantum particle swarm optimization for the Multidimensional Knapsack Problem
 * 2016
*/

#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <dirent.h>
#include "sys_stk.h"
#include "mdkp.h"

#define STK_SZ 256*1024*1024
#define MX_RNS 30
#define MN_PREC 2
#define MX_PRBS 1000
#define STA_NON 0
#define STA_DIR 1
#define STA_INP 2
#define STA_SE 3
#define STA_RNS 4
#define STA_SZ 5
#define STA_MX 6
#define STA_ALP 7
#define STA_EP1 8
#define STA_EP2 9
#define STA_STY 10
#define STA_NWL 11

#define INC_STR "incorrect argument(s)"
#define ARGS_STR "arguments error: "

#define DEFA_SE 198671
#define DEFA_RNS 30
#define DEFA_SZ 20
#define DEFA_MX 500
#define DEFA_ALP 0.1
#define DEFA_EP1 0.4
#define DEFA_EP2 0.2
#define DEFA_STY 0
#define DEFA_NWL false

using namespace std;

MdKP mdkp;
MdKPQ mdkpq;
double alp, ep1, ep2;
int se, rns, sz, mx, sty, nwl;

void slv()
{
    int ses[MX_RNS], i, sum, bst, tm;
    sum = 0;
    bst = 0;
    tm = 0;
    srand(se);
	cout << fixed << setprecision(MN_PREC);
    for (i = 0; i < rns; i++)
        ses[i] = rand();
    for (i = 0; i < rns; i++)
    {
        srand(ses[i]);
        cout << "run " << i + 1 << endl;
        mdkpq.slv(mdkp, true);
        sum += mdkpq.v_sol;
        tm += mdkpq.tm;
        if (mdkpq.v_sol > bst)
            bst = mdkpq.v_sol;
    }
    cout << "bst = " << bst << ", avg = " << (double)sum/rns << ", tm = " << tm << ", atm = " << (double)tm/rns << " ms" << endl;
}

void ini()
{
    mdkpq.info();
    mdkpq.ini(sz, mx, alp, ep1, ep2, true, sty, nwl);
}

void slvFl(string fn)
/* Reads a problem from file, runs an algorithm and outputs results. */
{
    cout << "solving " << fn << ":" << endl;
    if (mdkp.rdPrb(fn))
    {
        ini();
        slv();
    }
}

void slvDir(string dn)
/* Solves all instances in a directory. */
{
    DIR* dir;
    dirent* pdir;
    vector<string> fls;
    vector<string>::size_type v, np;
    string str, dstr;
    cout << "solving instances in " << dn << ":" << endl;
    dstr = dn.c_str();
    dir = opendir(dn.c_str());
    if (dir == NULL)
        cout << "bad directory" << endl;
    else
    {
        while ((pdir = readdir(dir)) != NULL)
        {
            str = pdir->d_name;
            if (str != "." && str != ".." && str != "README.txt")
                fls.push_back(str);
        }
        sort(fls.begin(), fls.end());
        np = fls.size();
        if (np > MX_PRBS)
        {
            cout << "too many problem files, limit is " << MX_PRBS << endl;
            return;
        }
        ini();
        for(v = 0; v < np; v++)
        {
            cout << endl;
            cout << "solving file " << fls[v] << endl;
            if (mdkp.rdPrb(dstr + "/" + fls[v]))
            {
                if (mdkp.opt > 0)
                    cout << "opt = " << mdkp.opt << ", se = " << se << endl;
                slv();
            }
        }
    }
}

int main(int argc, char **argv)
{
    vector<std::string> arg;
    string err, dn, fn;
    int i, sta;
    cout << "MdKPQ program" << endl;
    if (!setStkSz(STK_SZ))
        return 0;
    err = "";
    dn = "";
    fn = "";
    se = DEFA_SE;
    rns = DEFA_RNS;
    sz = DEFA_SZ;
    mx = DEFA_MX;
    alp = DEFA_ALP;
    ep1 = DEFA_EP1;
    ep2 = DEFA_EP2;
    sty = DEFA_STY;
    nwl = DEFA_NWL;
    if (argc > 1)
    {
        arg.assign(argv + 1, argv + argc);
        sta = STA_NON;
        try
        {
            for (i = 0; i < argc - 1; i++)
            {
                if (sta == STA_DIR)
                {
                    dn = arg[i];
                    sta = STA_NON;
                }
                else if (sta == STA_INP)
                {
                    fn = arg[i];
                    sta = STA_NON;
                }
                else if (sta == STA_SE)
                {
                    se = atoi(arg[i].c_str());
                    sta = STA_NON;
                }
                else if (sta == STA_RNS)
                {
                    rns = atoi(arg[i].c_str());
                    sta = STA_NON;
                }
                else if (sta == STA_SZ)
                {
                    sz = atoi(arg[i].c_str());
                    sta = STA_NON;
                }
                else if (sta == STA_MX)
                {
                    mx = atoi(arg[i].c_str());
                    sta = STA_NON;
                }
                else if (sta == STA_ALP)
                {
                    alp = atof(arg[i].c_str());
                    sta = STA_NON;
                }
                else if (sta == STA_EP1)
                {
                    ep1 = atof(arg[i].c_str());
                    sta = STA_NON;
                }
                else if (sta == STA_EP2)
                {
                    ep2 = atof(arg[i].c_str());
                    sta = STA_NON;
                }
                else if (sta == STA_STY)
                {
                    sty = atoi(arg[i].c_str());
                    sta = STA_NON;
                }
                else if (sta == STA_NWL)
                {
                    nwl = arg[i] == "1";
                    sta = STA_NON;
                }
                else if (arg[i] == "-dir")
                    sta = STA_DIR;
                else if (arg[i] == "-inp")
                    sta = STA_INP;
                else if (arg[i] == "-se")
                    sta = STA_SE;
                else if (arg[i] == "-rns")
                    sta = STA_RNS;
                else if (arg[i] == "-sz")
                    sta = STA_SZ;
                else if (arg[i] == "-mx")
                    sta = STA_MX;
                else if (arg[i] == "-alp")
                    sta = STA_ALP;
                else if (arg[i] == "-ep1")
                    sta = STA_EP1;
                else if (arg[i] == "-ep2")
                    sta = STA_EP2;
                else if (arg[i] == "-sty")
                    sta = STA_STY;
                else if (arg[i] == "-n")
                    sta = STA_NWL;
                else
                    err = string(INC_STR) + string(": ") + string(arg[i]);
            }
            if (sta != STA_NON)
                err = INC_STR;
        }
        catch(std::exception& ex)
        {
            err = ex.what();
        }
    }
    else
        err = "no arguments given";
    if (err == "")
    {
        if (dn != "")
            slvDir(dn);
        else if (fn != "")
            slvFl(fn);
        else
            err = "directory or file name not specified";
    }
    if (err != "")
    {
        cout << ARGS_STR << err << endl;
        cout << "compulsory arguments:" << endl;
        cout << "dir - instance directory" << endl;
        cout << "(or...)" << endl;
        cout << "inp - input filename" << endl;
        cout << "optional arguments:" << endl;
        cout << "se - random number seed" << endl;
        cout << "rns - number of runs" << endl;
        cout << "sz - population size" << endl;
        cout << "mx - maximum iterations" << endl;
        // TODO document (see reference paper):
        cout << "alp - (alpha) parameter" << endl;
        cout << "ep1 - (epsilon 1) parameter" << endl;
        cout << "ep2 - (epsilon 2) parameter" << endl;
        cout << "sty - output style (0 - 4)" << endl;
        cout << "nwl - switch that turns on new lines between outputs (0 or 1)" << endl;
    }
    return 0;
}
