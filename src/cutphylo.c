#include <R.h>
#include <Rinternals.h>

#include "set.h"

// compute the total residual sum of squares on the partition
// given by the union of the current cutset and the candidate
// branch that is being cut.
static double SSE(
    int v      // candidate cut
    , int ntip
    , int d
    , int nedge
    , const int *parent
    , const int *child
    , double *w
    , double *x
    , double *s
    , iset *cutset  // current cutset
)
{
    int i;
    int j;
    int u;
    const int *kid;
    double w1, w2, s1, s2;
    // loop over all internal nodes in postorder sequence
    for (i = 0, u = parent[i]; i < nedge; i += 2, u = parent[i])
    {
        kid = child + i;
        s[u] = s[kid[0]] + s[kid[1]];
        if (iset_contains(kid[0], cutset) || kid[0] == v)
        {
            w[u] = w[kid[1]];
            for (j = 0; j < d; ++j)
                x[j+u*d] = x[j+kid[1]*d];
        }
        else if (iset_contains(kid[1], cutset) || kid[1] == v)
        {
            w[u] = w[kid[0]];
            for (j = 0; j < d; ++j)
                x[j+u*d] = x[j+kid[0]*d];
        }
        else
        {
            s1 = 0;
            s2 = 0;
            w1 = w[kid[0]];
            w2 = w[kid[1]];
            w[u] = w1 + w2;
            for (j = 0; j < d; ++j)
                x[j+u*d] = (w1*x[j+kid[0]*d] + w2*x[j+kid[1]*d]) / (w1 + w2);
            for (j = 0; j < d; ++j)
            {
                s1 += (x[j+kid[0]*d] - x[j+u*d]) * (x[j+kid[0]*d] - x[j+u*d]);
                s2 += (x[j+kid[1]*d] - x[j+u*d]) * (x[j+kid[1]*d] - x[j+u*d]);
            }
            s[u] += w1*s1 + w2*s2;
        }
    }
    return s[ntip];
}

// greedy search for the best cutset
static void do_cuts(
    int ntip
    , int nnode
    , int d
    , int K
    , int nedge
    , const int *parent
    , const int *child
    , const int *sib
    , double *w
    , double *x
    , double *s
    , double *RSS
    , iset *cutset
)
{
    int i;
    int v;
    int k = 1;
    int best_cut;
    double min_s;
    double score;
    iset *p = malloc(sizeof(*p));
    if (!p || !iset_malloc(ntip+nnode, p))
        error("cutset allocation failure");
    for (i = 0; i < (ntip+nnode); ++i)
        iset_add(i, p);
    iset_remove(ntip, p);
    RSS[0] = SSE(-1,ntip,d,nedge,parent,child,w,x,s,cutset);
    while (k < K)
    {
        best_cut = -1;
        min_s = R_PosInf;
        for (i = 0; i < p->n; ++i)
        {
            v = p->dense[i];
            score = SSE(v,ntip,d,nedge,parent,child,w,x,s,cutset);
            if (score < min_s)
            {
                min_s = score;
                best_cut = v;
            }
        }
        iset_remove(best_cut, p);
        iset_remove(sib[best_cut], p);
        iset_add(best_cut, cutset);
        RSS[k] = min_s;
        Rprintf("\rProgress: %d / %d", k++, K-1);
    }
    Rprintf("\n");
    iset_free(p);
    free(p);
}


static void phyloweights(
    int ntip
    , int nnode
    , int nedge
    , const int *parent
    , const int *child
    , const double *brlen
    , double *w)
{
    int i;
    int p;
    const int *kid;
    double *v = malloc((ntip+nnode)*sizeof(*v));
    if (!v)
        error("memory allocation failure");
    memcpy(v, brlen, (ntip+nnode)*sizeof(*v));
    for (i = 0; i < nedge; i += 2)
    {
        p = parent[i];
        kid = child + i;
        v[p] += (v[kid[0]]*v[kid[1]]) / (v[kid[0]]+v[kid[1]]);
    }
    w[ntip] = (double)ntip;
    for (i = nedge-1; i >=0; i -= 2)
    {
        p = parent[i];
        kid = child + i - 2 + 1;
        w[kid[0]] = w[p] * v[kid[1]] / (v[kid[0]]+v[kid[1]]);
        w[kid[1]] = w[p] * v[kid[0]] / (v[kid[0]]+v[kid[1]]);
    }
    free(v);
}


static SEXP dp(
    int K
    , int d
    , int ntip
    , int nnode
    , int nedge
    , const int *parent
    , const int *child
    , const double *brlen
    , double *w
    , double *x)
{
    int i;
    int u;
    int p;
    int n = ntip + nnode;

    double *s = calloc(n, sizeof(*s));

    int *sib = malloc(n * sizeof(*sib));

    if (!s || !sib)
        error("memory allocation failure");

    for (i = 0; i < nedge; i += 2)
    {
        sib[child[i]] = child[i+1];
        sib[child[i+1]] = child[i];
    }

    iset *cutset = malloc(sizeof(*cutset));

    if (!cutset || !iset_malloc(n, cutset))
        error("cutset allocation failure");
    
    if (*w < 0)
        phyloweights(ntip,nnode,nedge,parent,child,brlen,w);


    SEXP cuts = PROTECT(allocVector(INTSXP, K));
    SEXP score = PROTECT(allocVector(REALSXP,K));
    SEXP subtrees = PROTECT(allocMatrix(INTSXP, n, K));

    do_cuts(ntip,nnode,d,K,nedge,parent,child,sib,w,x,s,REAL(score),cutset);

    INTEGER(cuts)[0] = ntip+1;
    for (i = 1; i < K; ++i)
        INTEGER(cuts)[i] = cutset->dense[i-1] + 1;

    for (i = 0; i < K; ++i)
    {
        p = 1;
        INTEGER(subtrees)[ntip+i*n] = p;
        for (u = nedge-1; u >= 0; u -= 2)
        {
            if (iset_contains(child[u], cutset) && 
                cutset->sparse[child[u]] < i)
            {
                INTEGER(subtrees)[child[u]+i*n] = ++p;
            }
            else
            {
                INTEGER(subtrees)[child[u]+i*n] = INTEGER(subtrees)[parent[u]+i*n];
            }
            if (iset_contains(child[u-1], cutset) && 
                cutset->sparse[child[u-1]] < i)
            {
                INTEGER(subtrees)[child[u-1]+i*n] = ++p;
            }
            else
            {
                INTEGER(subtrees)[child[u-1]+i*n] = INTEGER(subtrees)[parent[u-1]+i*n];
            }
        }
    }

    iset_free(cutset);
    free(cutset);
    free(s);
    free(sib);
    UNPROTECT(3);
    return list3(score, cuts, subtrees);
}


SEXP C_cutphylo(
    SEXP K          // max number of clusters to consider
    , SEXP ntip     // number of terminals
    , SEXP nnode    // number internal nodes
    , SEXP nedge    // number of edges
    , SEXP parent   // parent array in postorder sequence
    , SEXP child    // child array in postorder sequence
    , SEXP brlen    // branch lengths
    , SEXP w        // tip weights
    , SEXP x)       // trait data for terminals
{    
    return dp(
        *INTEGER(K)
        , *INTEGER(getAttrib(x, R_DimSymbol))
        , *INTEGER(ntip)
        , *INTEGER(nnode)
        , *INTEGER(nedge)
        , INTEGER(parent)
        , INTEGER(child)
        , REAL(brlen)
        , REAL(w)
        , REAL(x)
    );
}
