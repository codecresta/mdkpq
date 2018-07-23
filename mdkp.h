#include <string>
#include <fstream>
#include <chrono>
#include <iostream>

using namespace std;

typedef chrono::high_resolution_clock timer;
typedef chrono::milliseconds millisecs;
typedef chrono::high_resolution_clock::time_point tm_pt;

#define duration chrono::duration_cast<millisecs>

#define M_MX 100
#define N_MX 2500
#define S_MX 20
#define STY_NONE 0
#define STY_BAS 1
#define STY_BAS_X 2
#define STY_VAL 3
#define STY_VAL_X 4

int calcMs(millisecs dur);
int rndU(int lo, int up);
double rndDbl();

class MdKP
/* Class for Multidimensional Knapsack Problems:
 * pfn = problem filename (if needed),
 * pfs = file stream (if needed),
 * itms = number of items,
 * dims = number of dimensions (constraints),
 * val = values of items,
 * wei = weights of items (for constraints),
 * cap = capacity of the knapsack,
 * sol = incumbent solution (variable values),
 * val_sol = value of incumbent solution,
 * opt = optimum value (read from file, zero if not known).
*/
{
public:
    string pfn;
    ifstream pfs;
    int itms, dims, val[N_MX], wei[M_MX][N_MX], cap[M_MX], sol[N_MX], v_sol, opt;
    virtual int invIdx(int i);
    bool rdPrb(string fn0);
    void uncrlPrb(int dims0, int itms0, int cap_lo, int cap_up, int val_lo, int val_up, int wei_lo, int wei_up);
    void crlPrb(int dims0, int itms0, int cap_lo, int cap_up, int wei0, double rng);
    void outPrbSmry();
    bool feasSol();
    int vSol();
    void chk();
private:
    bool rdInt(string txt, int &i);
    bool clsFl(bool res);
};

class MdKPB: public MdKP
/* Abstract class for Multidimensional Knapsack Algorithm variants. */
{
public:
    int sty;
    bool nwl;
    virtual std::string name() = 0;
    virtual int varIdx(int j);
    void outPrb();
    void outSol();
};

class MdKPQ: public MdKPB
/* Hybrid Quantum Particle Swarm Optimisation Algorithm. */
{
private:
    int inv[N_MX], sz, mx;
    double alp, bet, ep1, ep2, ep3;
    bool prep(MdKP &src);
    int calV(int x[]);
    void cpy(int x[], int y[]);
    void rmc(int x[], int r[]);
    void cpyR(int s[], int t[]);
    void cpyD(double s[], double t[]);
    bool inf(int r[]);
    bool drp(int k, int x[], int r[], int &v);
    void add(int k, int x[], int r[], int &v);
    void lcl(int x[], int &v);
    void algS(int s[], int t[], int &vi);
    void alg();
public:
    int tm;
    std::string name()
    {
        return "hybrid quantum particle swarm optimisation algorithm";
    }
    void info();
    void ini(int sz0, int mx0, double alp0, double ep10, double ep20, bool out, int sty0, bool nwl0);
    int slv(MdKP &src, bool out);
};
