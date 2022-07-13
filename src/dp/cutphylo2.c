#include <R.h>
#include <Rinternals.h>

#include "etwiddle.h"
#include "combn.h"

typedef struct iset {
    int n;
    int size;
    int *dense;
} iset;

static int iset_malloc(int size, iset *s)
{
    if (s)
    {
        s->n = 0;
        s->size = size;
        s->dense = malloc(size * sizeof(int));
        if (!s->dense)
            return 0;
        return 1;
    }
    return 0;
}

static void iset_free(iset *s)
{
    if (s)
    {
        free(s->dense);
        s->dense = NULL;
    }
}

static int iset_add(int i, iset *s)
{
    if (s->n < s->size)
    {
        s->dense[(s->n)++] = i;
        return 1;
    }
    return 0;
}

static void iset_union(iset *a, iset *b, iset *c)
{
    int i;
    for (i = 0; i < a->n; ++i)
        iset_add(a->dense[i], c);
    for (i = 0; i < b->n; ++i)
        iset_add(b->dense[i], c);
}

static void iset_clear(iset *s)
{
    s->n = 0;
}

typedef struct cutset {
    // size of cutset
    int k;

    // length of x
    int d;

    // root of subtree
    int v;

    // set of branches to cut that partition
    // subtree into k subtrees
    iset *set;

    // this is the mean value of the
    // rootmost partition subtree
    double *x;

    // weight of the rootmost
    // partition subtree
    double w;

    // total residual sum of squares of
    // partition subtrees induced by cutset
    double score;

    // next cutset
    struct cutset *next;

} cutset;

static int cutset_malloc(int k, int d, int v, cutset *s)
{
    if (s)
    {
        s->k = k;
        s->d = d;
        s->v = v;
        s->w = 0;
        s->score = R_PosInf;
        s->next = NULL;
        s->set = malloc(sizeof(iset));
        if (!s->set)
            return 0;
        s->x = calloc(d, sizeof(double));
        if (!s->x)
        {
            free(s->set);
            s->set = NULL;
            return 0;
        }
        if (!iset_malloc(k, s->set))
        {
            free(s->set);
            free(s->x);
            s->set = NULL;
            s->x = NULL;
            return 0;
        }
        return 1;
    }
    return 0;
}

static void cutsetcpy(cutset *dest, const cutset *src)
{
    memcpy(dest->x, src->x, src->d * sizeof(double));
    iset_clear(dest->set);
    for (int i = 0; i < src->set->n; ++i)
        iset_add(src->set->dense[i], dest->set);
    dest->v = src->v;
    dest->w = src->w;
    dest->score = src->score;
    dest->next = NULL;
}

static cutset *cutset_copy(const cutset *a)
{
    cutset *b = malloc(sizeof(*b));
    if (!cutset_malloc(a->k, a->d, a->v, b))
        error("cutset allocation failure");
    memcpy(b->x, a->x, a->d * sizeof(double));
    for (int i = 0; i < a->set->n; ++i)
        iset_add(a->set->dense[i], b->set);
    b->w = a->w;
    b->score = a->score;
    b->next = NULL;
    return b;
}

static void cutset_free(cutset *s)
{
    if (s)
    {
        iset_free(s->set);
        free(s->set);
        free(s->x);
        s->next = NULL;
    }
}

static void cutset_clear(cutset *s)
{
    cutset *tmp;
    while (s->next)
    {
        tmp = s->next;
        s->next = tmp->next;
        cutset_free(tmp);
        free(tmp);
    }
    s->w = 0;
    s->score = R_PosInf;
    s->next = NULL;
    iset_clear(s->set);
}

static void cutset_add(cutset *a, const cutset *b)
{
    cutset *bb = cutset_copy(b);
    bb->next = a->next;
    a->next = bb;
}

static void cutset_merge(const cutset *a, const cutset *b, cutset *c, int mtype)
{
    int i;
    double s1 = 0;
    double s2 = 0;
    double w1 = a->w;
    double w2 = b->w;
    double w = w1 + w2;

    if (mtype == 3) // merge both
    {
        c->w = w;
        for (i = 0; i < a->d; ++i)
            c->x[i] = (w1 * a->x[i] + w2 * b->x[i]) / w;
        for (i = 0; i < a->d; ++i)
        {
            s1 += (a->x[i] - c->x[i]) * (a->x[i] - c->x[i]);
            s2 += (b->x[i] - c->x[i]) * (b->x[i] - c->x[i]);
        }
        iset_union(a->set, b->set, c->set);
        c->score = a->score + b->score + w1*s1 + w2*s2;
    }
    else if (mtype == 2) // merge right
    {
        c->w = w2;
        for (i = 0; i < b->d; ++i)
            c->x[i] = b->x[i];
        iset_union(a->set, b->set, c->set);
        iset_add(a->v, c->set);
        c->score = a->score + b->score;
    }
    else if (mtype == 1) // merge left
    {
        c->w = w1;
        for (i = 0; i < a->d; ++i)
            c->x[i] = a->x[i];
        iset_union(a->set, b->set, c->set);
        iset_add(b->v, c->set);
        c->score = a->score + b->score;
    }
}

// c is the parent cutset, a and b are the child cutsets
static void foo(cutset *a, cutset *b, cutset *c, int mtype)
{
    double min_s = c->score;

    cutset *aa = a;
    cutset *bb = b;

    cutset *tmp = malloc(sizeof(*tmp));
    if (!cutset_malloc(a->k, a->d, c->v, tmp))
        error("cutset allocation failure");

    while (aa)
    {
        bb = b;
        while (bb)
        {
            cutset_merge(aa, bb, tmp, mtype);
            if (tmp->score < min_s)
            {
                min_s = tmp->score;
                cutset_clear(c);
                cutsetcpy(c, tmp);
            }
            else if (min_s < R_PosInf && tmp->score == min_s)
            {
                cutset_add(c, tmp);
            }
            cutset_clear(tmp);
            bb = bb->next;
        }
        aa = aa->next;
    }
    cutset_free(tmp);
    free(tmp);
}


static void downpass(
    int ntip
    , int nnode
    , int K
    , int nedge
    , const int *parent
    , const int *child
    , const double *brlen
    , cutset *s)
{

    int v, i, k, l, i1, i2, i3, m, ks[2], r[2];

    const int *kid;

    cutset *a, *b, *c;

    // loop over all internal nodes in postorder sequence
    for (i = 0, v = parent[i], l = 0; i < nedge; i += 2, v = parent[i])
    {
        kid = child + i;
        for (k = 1; k <= K; ++k)
        {
            for (i1 = 1; i1 <= 3; ++i1)
            {
                c = s + (i1-1) + k*3 + v*3*(K+1);
                m = (i1 == 3) ? 2 : 1;
                // loop over all compositions of the integer k+m-1 into 2 
                // parts. subtree v will be partitioned into k subtrees. each
                // component of the composition will be the number of subtrees 
                // that v's child subtrees are partitioned into.
                for (i2 = etwiddle_begin(2,k+m-1,ks);
                     i2;
                     i2 = etwiddle_next(2,k+m-1,ks))
                {
                    for (i3 = combn_begin(3,3,r); i3; i3 = combn_next(3,3,r))
                    {
                        a = s + r[0] + ks[0]*3 + kid[0]*3*(K+1);
                        b = s + r[1] + ks[1]*3 + kid[1]*3*(K+1);
                        foo(a, b, c, i1);
                    }
                }
            }
        }
        Rprintf("\rProgress: %d / %d", ++l, nnode);
    }
    Rprintf("\n");
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
    , const double *x)
{
    int u;
    int v;
    int k;
    int n = ntip + nnode;

    cutset *a;
    cutset *s = malloc(3*(K+1)*n*sizeof(*s));

    if (!s)
        error("cutset allocation failure");
    
    if (*w < 0)
        phyloweights(ntip,nnode,nedge,parent,child,brlen,w);      

    // initialize the cutsets
    for (v = 0; v < n; ++v)
    {
        if (v < ntip)
        {
            for (k = 1; k <= K; ++k)
            {
                for (u = 0; u < 3; ++u)
                {
                    if (!cutset_malloc(k, d, v, s + u + k*3 + v*3*(K+1)))
                        error("cutset allocation failure");
                    if (u == 2 && k == 1)
                    {
                        a = s + u + k*3 + v*3*(K+1);
                        a->w = w[v];
                        a->score = 0;
                        for (u = 0; u < d; ++u)
                            a->x[u] = x[u+v*d];
                    }
                }
            }
        }
        else
        {
            for (k = 1; k <= K; ++k)
            {
                for (u = 0; u < 3; ++u)
                {
                    if (!cutset_malloc(k, d, v, s + u + k*3 + v*3*(K+1)))
                        error("cutset allocation failure");
                }
            }
        }
    }

    SEXP dim = PROTECT(allocVector(INTSXP, 3));
    INTEGER(dim)[0] = n;
    INTEGER(dim)[1] = 3;
    INTEGER(dim)[2] = K;

    SEXP cutset = PROTECT(allocArray(INTSXP, dim));
    SEXP score = PROTECT(allocMatrix(REALSXP, 3, K));

    memset(INTEGER(cutset), 0, n*3*K*sizeof(int));

    downpass(ntip, nnode, K, nedge, parent, child, brlen, s);

    for (k = 1; k <= K; ++k)
    {
        for (u = 0; u < 3; ++u)
        {
            a = s + u + k*3 + ntip*3*(K+1);

            REAL(score)[u+(k-1)*3] = a->score;
            //REAL(score)[u+(k-1)*3] = 0;
            //while (a)
            //{
            //    REAL(score)[u+(k-1)*3] += 1;
            //    a = a->next;
            //}
        }
    }

    for (v = 0; v < n; ++v)
    {
        for (k = 1; k <= K; ++k)
        {
            for (u = 0; u < 3; ++u)
            {
                a = s + u + k*3 + v*3*(K+1);
                cutset_clear(a);
                cutset_free(a);
            }
        }
    }

    free(s);

    UNPROTECT(3);
    //return list2(score, cutset);
    return score;
}

SEXP C_cutphylo2(
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