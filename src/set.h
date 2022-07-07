/*
  Reference:

  P. Briggs and L. Torczon. An efficient representation for sparse sets.
  ACM Letters on Programming Languages and Systems 2(1-4), 59-69, 1993.
*/

typedef struct iset {
    int n;
    int size;
    int *sparse;
    int *dense;
} iset;

static int iset_malloc(int size, iset *s)
{
    if (s)
    {
        s->n = 0;
        s->size = size;
        s->sparse = NULL;
        s->dense = NULL;
        s->sparse = malloc(size * sizeof(int));
        if (!s->sparse)
            return 0;
        s->dense = malloc(size * sizeof(int));   
        if (!s->dense)
        {
            free(s->sparse);
            s->sparse = NULL;
            return 0;
        }
        while (--size >= 0)
        {
            s->sparse[size] = -1;
            s->dense[size] = -1;
        }
        return 1;
    }
    return 0;
}

static void iset_free(iset *s)
{
    if (s)
    {
        free(s->sparse);
        free(s->dense);
        s->sparse = NULL;
        s->dense = NULL;
    }
}

static int iset_contains(int i, iset *s)
{
    int a = s->sparse[i];
    if (a >= 0 && a < s->n && s->dense[a] == i)
        return 1;
    return 0;
}

static int iset_get(int i, iset *s)
{
    if (i > 0 && i <= s->n)
        return s->dense[--i];
    return -1;
}

static int iset_size(iset *s)
{
    return s->n;
}

static void iset_clear(iset *s)
{
    s->n = 0;
}

static int iset_add(int i, iset *s)
{
    if (!iset_contains(i, s))
    {
        s->sparse[i] = s->n;
        s->dense[(s->n)++] = i;
        return 1;
    }
    return 0;
}

static int iset_remove(int i, iset *s)
{
    if (iset_contains(i, s))
    {
        int a = s->sparse[i];
        int e = s->dense[--(s->n)];
        s->dense[a] = e;
        s->sparse[e] = a;
        return 1;
    }
    return 0;
}

static int iset_pop(iset *s)
{
    if (s->n)
        return s->dense[--(s->n)];
    return -1;
}

static void iset_union(iset *a, iset *b, iset *c)
{
    int i;
    for (i = 0; i < a->n; ++i)
        iset_add(a->dense[i], c);
    for (i = 0; i < b->n; ++i)
        iset_add(b->dense[i], c);
}

static void iset_intersect(iset *a, iset *b, iset *c)
{
    int i;

    if (a->n < b->n)
    {
        for (i = 0; i < a->n; ++i)
        {
            if (iset_contains(a->dense[i], b))
                iset_add(a->dense[i], c);
        }
    }
    else
    {
        for (i = 0; i < b->n; ++i)
        {
            if (iset_contains(b->dense[i], a))
                iset_add(b->dense[i], c);
        }
    }
}
