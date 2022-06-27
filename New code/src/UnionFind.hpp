#ifndef UNIONFIND_HPP
#define UNIONFIND_HPP

class UnionFind
{
    size_t n;
    int* parent;
    int* rank;

    public:

    UnionFind(size_t n) : n{n}, parent{new int[n]}, rank{new int[n]()}
    {
        for (int i = 0; i < n; i++)
        {
            parent[i] = i;
        }
    }

    ~UnionFind()
    {
        delete[] parent;
        delete[] rank;
    }

    inline int size()
    {
        return n;
    }

    int find(int x)
    {
        while (x != parent[x])
        {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }

        return x;
    }

    void merge(const int x, const int y)
    {
        int p = find(x);
        int q = find(y);

        if (p == q)
        {
            return;
        }

        if (rank[p] > rank[q])
        {
            parent[q] = p;
        }
        else
        {
            parent[p] = q;

            if (rank[p] == rank[q])
            {
                rank[q]++;
            }
        }
    }
};

#endif
