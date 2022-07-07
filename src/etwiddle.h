/* https://stackoverflow.com/questions/15577651/generate-all-compositions-of-an-integer-into-k-parts
**
** Generate all k-part compositions of in integer n.
** Here is an example:
**
** int n = 5, k = 3, b[k];
** for (i = etwiddle_begin(k, n, b); i; i = etwiddle_next(k, n, b))
** {
**      // do something with b
** }
**
** where the loop sets b, in turn, to the arrays
**
** [1,1,3]
** [1,2,2]
** [1,3,1]
** [2,1,2]
** [2,2,1]
** [3,1,1]
*/
static int etwiddle_begin(int k, int n, int *p)
{
    if (k > n || k <= 0)
        return 0;
    int i;
    for (i = 0; i < k-1; ++i)
        p[i] = 1;
    p[k-1] = n-k+1;
    return 1;
}

static int etwiddle_next(int k, int n, int *p)
{
    if (p[0] == n-k+1)
        return 0;
    
    int last = k - 1;

    while (p[last] == 1)
        --last;
    
    // turn    a b ...   y   z 1 1 ...   1
    //                       ^ last
    // into    a b ... (y+1) 1 1 1 ... (z-1)

    // be careful, there may be no 1's at the end

    int z = p[last];
    p[last - 1] += 1;
    p[last] = 1;
    p[k - 1] = z - 1;
    
    return 1;
}
