#include <R.h>
#include <Rinternals.h>

#include "etwiddle.h"
#include "combn.h"

/* 
** For each k from 1...K, we use a dynamic programming algorithm to find the
** set of branches in the tree that, when cut, yield the optimal 
** partitioning of the data into k subtrees.
**
** The following reference served as inspiration:
**
** Ying. Xu, Victor Olman, Dong Xu. Clustering gene expression data using
** a graph-theoretic approach: an application of minimum spanning trees. 
** Bioinformatics 18(4), 536-545, 2002.
*/
typedef struct scorearr {

    // number of internal and terminal nodes
    int nnode;
    
    // max number of clusters plus 1
    int K;
    
    // data dimension
    int d;
    
    // each s(,,v) is a 3-by-K matrix
    //
    //    0  1  2  3 ... (K-1)
    // 10 -  -  -  - ... -
    // 01 -  -  -  - ... -
    // 11 -  -  -  - ... -
    //
    // that records the minimum RSS possible
    // on the subtree rooted at v under the constraints
    // specified by the row and column combination.
    // the columns indicate the number of subtrees
    // that v's subtree is partitioned into. the
    // rows indicate which of v's child branches belong
    // to v's partition subtree. 11 indicates both branches, 01
    // indicates the right branch, 10 the left branch.
    // dim: (3,K,nnode)
    double *s;

    // weighted means of each partition subtree
    // dim: (d,3,K,nnode)
    double *x;

    // weights of each partition subtree (measured as sum of tip weights)
    // dim: (3,K,nnode)
    double *w;

    // bookkeeping array for backtracking to actually cut
    // the tree into optimal subtrees once the algorith completes.
    // for each cell in the s(,,v) table, we need to record the
    // cell identities in v's child tables that were used in its
    // formation. so, each cell needs 4 ints: the first two will
    // be the (row,col) of the left child and the last two will be
    // the (row,col) indices of the right child.
    // dim: (4,3,K,nnode)
    int *b;
} scorearr;


#define SCORE(i,k,v) s->s[(i) + (k)*3 + (v)*3*s->K]
#define WEIGHT(i,k,v) s->w[(i) + (k)*3 + (v)*3*s->K]
#define MEAN(i,j,k,v) s->x[(i) + (j)*s->d + (k)*s->d*3 + (v)*s->d*3*s->K]
#define BOOK(i,j,k,v) s->b[(i) + (j)*4 + (k)*4*3 + (v)*4*3*s->K]

static void downpass(
  int ntip
  , int nnode
  , int K
  , int nedge
  , const int *parent
  , const int *child
  , const double *brlen
  , scorearr *s)
{
  int i;
  int j;
  int k;
  int d;
  int v;
  int n;
  int u;
  int l;
  const int *kid;

  int i1, i2, i3, m;
  int b[2];
  int r[2];
  int ks[2];

  double w, w1, w2;
  double tmp;
  double bs[2];
  double score;
  double min_s;
  double X[s->d];

  // loop over all internal nodes in postorder sequence
  for (i = 0, v = parent[i], l = 0; i < nedge; i += 2, v = parent[i])
  {
    kid = child + i;
    for (k = 1; k <= K; ++k)
    {
      for (i1 = 1; i1 <= 3; ++i1)
      {
        switch (i1)
        {
          case 1:
            b[0] = 1;
            b[1] = 0;
            m = 1;
            break;
          case 2:
            b[0] = 0;
            b[1] = 1;
            m = 1;
            break;
          case 3:
            b[0] = 1;
            b[1] = 1;
            m = 2;
            break;
        }
        min_s = SCORE(i1-1, k, v);
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
            bs[0] = bs[1] = 0;
            if (b[0] && b[1])
            {
              // both children are being combined into v's partition subtree
              w1 = WEIGHT(r[0],ks[0],kid[0]);
              w2 = WEIGHT(r[1],ks[1],kid[1]);
              w = w1 + w2;
              for (d = 0; d < s->d; ++d)
              {
                X[d] = (w1*MEAN(d,r[0],ks[0],kid[0]) +
                        w2*MEAN(d,r[1],ks[1],kid[1])) / w;
                
                bs[0] += (MEAN(d,r[0],ks[0],kid[0]) - X[d]) *
                         (MEAN(d,r[0],ks[0],kid[0]) - X[d]);
                bs[1] += (MEAN(d,r[1],ks[1],kid[1]) - X[d]) *
                         (MEAN(d,r[1],ks[1],kid[1]) - X[d]);
              }
              tmp = w1*bs[0] + w2*bs[1];
              score = SCORE(r[0],ks[0],kid[0]) + SCORE(r[1],ks[1],kid[1]) + tmp;
            }
            else
            {
              u = b[0] ? 0 : 1;
              for (d = 0; d < s->d; ++d)
                X[d] = MEAN(d,r[u],ks[u],kid[u]);
              w = WEIGHT(r[u],ks[u],kid[u]);
              score = SCORE(r[0],ks[0],kid[0]) + SCORE(r[1],ks[1],kid[1]);
            }
            if (score < min_s)
            {
              min_s = score;
              SCORE(i1-1, k, v) = min_s;
              WEIGHT(i1-1, k, v) = w;
              BOOK(0, i1-1, k, v) = r[0];
              BOOK(1, i1-1, k, v) = ks[0];
              BOOK(2, i1-1, k, v) = r[1];
              BOOK(3, i1-1, k, v) = ks[1];
              for (d = 0; d < s->d; ++d)
                MEAN(d, i1-1, k, v) = X[d];
            }
          }
        }
      }
    }
    Rprintf("\rProgress: %d / %d", ++l, nnode);
  }
  Rprintf("\n");
}

static void backtrack(
    int ntip
    , int nnode
    , int K
    , int nedge
    , const int *parent
    , const int *child
    , int *cutset
    , double *score
    , scorearr *s)
{
    int i;
    int j;
    int k;
    int u;
    int v;
    int n = ntip + nnode;
    const int *kid;

    int p;
    int root = parent[nedge-1];

    int *ks = malloc(n * sizeof(int));
    int *rs = malloc(n * sizeof(int));

    for (k = 1; k <= K; ++k)
    {
        for (u = 0; u < 3; ++u)
        {
            if (SCORE(u, k, root) < R_PosInf)
            {
                p = 1;
                cutset[root+u*n+(k-1)*n*3] = p;
                ks[root] = k;
                rs[root] = u;
                score[u+3*(k-1)] = SCORE(u, k, root);
                for (i = nedge-1; i >= 0; i -= 2)
                {
                    v = parent[i];
                    kid = child+i-2+1;
                    switch (rs[v])
                    {
                        case 0:
                            cutset[kid[0]+u*n+(k-1)*n*3] = cutset[v+u*n+(k-1)*n*3];
                            cutset[kid[1]+u*n+(k-1)*n*3] = ++p;
                            break;
                        case 1:
                            cutset[kid[1]+u*n+(k-1)*n*3] = cutset[v+u*n+(k-1)*n*3];
                            cutset[kid[0]+u*n+(k-1)*n*3] = ++p;
                            break;
                        case 2:
                            cutset[kid[0]+u*n+(k-1)*n*3] = cutset[v+u*n+(k-1)*n*3];
                            cutset[kid[1]+u*n+(k-1)*n*3] = cutset[v+u*n+(k-1)*n*3];
                            break;
                    }
                    rs[kid[0]] = BOOK(0,rs[v],ks[v],v);
                    ks[kid[0]] = BOOK(1,rs[v],ks[v],v);
                    rs[kid[1]] = BOOK(2,rs[v],ks[v],v);
                    ks[kid[1]] = BOOK(3,rs[v],ks[v],v);
                }
            }
            else
            {
                score[u+3*(k-1)] = R_PosInf;
            }
        }
    }
    free(ks);
    free(rs);
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
    scorearr a;
    scorearr *s = &a;
    a.nnode = n;
    a.K = K+1;
    a.d = d;
    a.s = malloc(3*a.K*n*sizeof(double));
    a.w = calloc(3*a.K*n, sizeof(double));
    a.x = calloc(a.d*3*a.K*n, sizeof(double));
    a.b = malloc(4*3*a.K*n*sizeof(int));
    
    if (*w < 0)
        phyloweights(ntip,nnode,nedge,parent,child,brlen,w);      


    // initialize the scorearr
    for (v = 0; v < n; ++v)
    {
        if (v < ntip)
        {
            for (k = 0; k <= K; ++k)
            {
                if (k == 1)
                {
                    SCORE(0, k, v) = R_PosInf;
                    SCORE(1, k, v) = R_PosInf;
                    SCORE(2, k, v) = 0;
                    WEIGHT(0, k, v) = 0;
                    WEIGHT(1, k, v) = 0;
                    WEIGHT(2, k, v) = w[v];
                    for (u = 0; u < d; ++u)
                    {
                        MEAN(u, 0, k, v) = 0;
                        MEAN(u, 1, k, v) = 0;
                        MEAN(u, 2, k, v) = x[u+v*d];
                    }
                }
                else
                {
                    SCORE(0, k, v) = R_PosInf;
                    SCORE(1, k, v) = R_PosInf;
                    SCORE(2, k, v) = R_PosInf;
                }
            }
        }
        else
        {
            for (k = 0; k <= K; ++k)
            {
                SCORE(0, k, v) = R_PosInf;
                SCORE(1, k, v) = R_PosInf;
                SCORE(2, k, v) = R_PosInf;
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
    backtrack(ntip, nnode, K, nedge, parent, child, 
              INTEGER(cutset), REAL(score), s);

    free(a.s);
    free(a.w);
    free(a.x);
    free(a.b);

    UNPROTECT(3);
    return list2(score, cutset);
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
