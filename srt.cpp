#include "srt.h"

int prt(int idx[], double val[], int lo, int up)
/* Partition for quick sort:
 * idx = index array,
 * val = value array,
 * lo = lower index,
 * up = upper index.
*/
{
	int i, k, t;
	double a;
	i = lo;
	a = val[idx[up]];
	for (k = lo; k < up; k++)
	{
		t = idx[k];
		if (val[t] > a)
		{
			idx[k] = idx[i];
			idx[i++] = t;
		}
	}
	t = idx[i];
	idx[i] = idx[up];
	idx[up] = t;
	return i;
}

void qSrt(int idx[], double val[], int lo, int up)
/* Quick sort algorithm. Sorts by value in reverse order.
 * idx = index array,
 * val = value array,
 * lo = lower index,
 * up = upper index.
*/
{
	int p;
	if (lo < up)
	{
		p = prt(idx, val, lo, up);
		qSrt(idx, val, lo, p - 1);
		qSrt(idx, val, p + 1, up);
	}
}
