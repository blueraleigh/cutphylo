/* 
** Generate all pairwise combinations of the integer sequences 
** 0:(m-1) and 0:(n-1)
**
** Here is an example:
**
** int m = 2, n = 3, b[2];
** for (i = combn_begin(m, n, b); i; i = combn_next(m, n, b))
** {
**      // do something with b
** }
**
** where the loop sets b, in turn, to the arrays
**
** [0,0]
** [0,1]
** [0,2]
** [1,0]
** [1,1]
** [1,2]
*/
static int combn_begin(int m, int n, int *p)
{
    if (m <= 0 || n <= 0)
        return 0;
    p[0] = p[1] = 0;
    return 1;
}

static int combn_next(int m, int n, int *p)
{
    if (p[0] == m-1 && p[1] == n-1)
        return 0;
    if (p[1] == n-1)
    {
        p[0] += 1;
        p[1] = 0;
    }
    else
        p[1] += 1;
    return 1;
}
