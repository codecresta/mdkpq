#include <iostream>
#include "sys_stk.h"

using namespace std;

bool setStkSz(rlim_t sz)
{
	struct rlimit rl;
	int res;
	res = getrlimit(RLIMIT_STACK, &rl);
	if (res == 0)
	{
		if (rl.rlim_cur < sz)
		{
			rl.rlim_cur = sz;
			res = setrlimit(RLIMIT_STACK, &rl);
			if (res != 0)
			{
				cout << "setrlimit returned result: " << res << endl;
				return false;
			}
		}
	}
	else
	{
		cout << "getrlimit returned result: " << res << endl;
		return false;
	}
	return true;
}
