#include <queue>
#include <algorithm>
#include <vector>
#include <cmath>
#include <cstdio>
#include <glpk.h>

using namespace std;

const double infty   = 1000000000.0;
const double epsilon = 0.001;

class node
{
public:
    vector<int>  tour;
    vector<int> visited;
    double       length;
    double       bound;
    int          cVisited;
};

bool operator<(const node& lhs, const node& rhs)
{
    // We want min heap behaviour!
    return lhs.bound > rhs.bound;
}

// {{{ TCPTour class
class TCPTour
{
public:
    TCPTour(const vector<pair<double, double> > &nodes, double dist)
        : m_dist(dist), m_best(infty)
    {
        printf("Creating TCPTour\n");
        for (vector<pair<double, double> >::const_iterator it = nodes.begin();
                it != nodes.end();
                ++it)
        {
            m_nodes.push_back(*it);
            m_nearby.push_back(vector<int>());
        }
        printf("Stored %d vertices\n", m_nodes.size());
        fflush(stdout);

        // Initialize nearby vertices
        for (int i = 0; i < m_nodes.size(); ++i)
        {
            m_nearby[i].push_back(i);
            for (int j = 0; j < i; ++j)
                if (getDist(i, j) <= m_dist + epsilon)
                {
                    m_nearby[i].push_back(j);
                    m_nearby[j].push_back(i);
                }
        }
//        for (int i=0;i<m_nodes.size();++i)
//        {
//            printf("m_nearby[%d]; ", i);
//            for (vector<int>::iterator it = m_nearby[i].begin();
//                    it != m_nearby[i].end();
//                    ++it)
//                printf("%d ", *it);
//            printf("\n");
//        }
//        fflush(stdout);
        printf("Done initializing\n");
        fflush(stdout);
    }

    void run()
    {
        TCP();
    }

    double getLength()
    {
        return m_best;
    }

    vector<int> getTour()
    {
        return m_bestTour;
    }

    int getNodeCount()
    {
        return m_cNodes;
    }

private:
    // List of all the vertices. Note that all vertices are connected
    vector<pair<double, double> > m_nodes;
    // m_nearby[i] is all vertices at most m_dist distance away from vertex i
    vector<vector<int> > m_nearby;
    // Distance to view a monument
    double      m_dist;
    double      m_best;
    vector<int> m_bestTour;
    int         m_cNodes;

    double getDist(int i, int j)
    {
        // Euclidian distance
        double dist = sqrt((m_nodes[i].first  - m_nodes[j].first)  *
                    (m_nodes[i].first  - m_nodes[j].first)  +
                    (m_nodes[i].second - m_nodes[j].second) *
                    (m_nodes[i].second - m_nodes[j].second));
//        printf("dist(%d, %d) = %.6lf\n", i, j, dist);
        return dist;
    }

    double getBound(const node& n)
    {
        glp_prob *pLP = glp_create_prob();
        glp_set_obj_dir(pLP, GLP_MIN);
        glp_add_cols(pLP, m_nodes.size() * m_nodes.size());
        for (int i = 0; i < m_nodes.size(); ++i)
        {
            // Decision variables. Relaxed to 0 <= x_ij <= 1
            for (int j = 0; j < i; ++j)
            {
                glp_set_obj_coef(pLP, i*m_nodes.size() + j+1, getDist(i, j));
                glp_set_col_bnds(pLP, i*m_nodes.size() + j+1, GLP_DB, 0.0, 1.0);
            }
        }
        int ind[m_nodes.size()][m_nodes.size()];
        double coef[m_nodes.size()][m_nodes.size()];
        glp_add_rows(pLP, m_nodes.size());
        for (int i = 0; i < m_nodes.size(); ++i)
        {
            // All nodes should have exactly 0 or 2 edges incident. Relaxed
            // to 0 <= deg <= 2
            glp_set_row_bnds(pLP, i+1, GLP_DB, 0.0, 2.0);
            for (int j = 0; j < i; ++j)
            {
                ind[i][j+1]  = i*m_nodes.size() + j;
                coef[i][j+1] = 1.0;
            }
            for (int j = i+1; j < m_nodes.size(); ++j)
            {
                ind[i][j]  = j*m_nodes.size() + i;
                coef[i][j] = 1.0;
            }
            glp_set_mat_row(pLP, i+1, m_nodes.size()-1, ind[i], coef[i]);
        }
        for (int i = 0; i < m_nodes.size(); ++i)
        {
            // Add vicinity constraints for each node
//            int k = glp_add_rows(pLP, m_nearby[i].size());
//            int ind2[(m_nodes.size()-1) * m_nearby
        }
        glp_smcp params;
        glp_init_smcp(&params);
        params.msg_lev = GLP_MSG_ERR;
        glp_simplex(pLP, &params);
        double dRet = glp_get_obj_val(pLP);
        glp_delete_prob(pLP);
//        return n.length;
        return dRet;
    }

    void TCP()
    {
        m_cNodes = 0;
        m_best   = infty;
        priority_queue<node> pq;
        // Create the initial nodes
        printf("Creating initial nodes\n");
        fflush(stdout);
        for (int i = 0; i < m_nodes.size(); ++i)
        {
            node n;
            for (int j = 0; j < m_nodes.size(); ++j)
                n.visited.push_back(0);
            n.length = 0.0;
//            printf("vis[%d] = ", i);
            n.cVisited = 0;
            for (vector<int>::const_iterator it = m_nearby[i].begin();
                    it != m_nearby[i].end();
                    ++it)
            {
//                printf("%d ", *it);
                n.visited[*it] = 1;
                n.cVisited++;
            }
//            printf("\n");
            n.tour.push_back(i);
            n.bound = getBound(n);
//            printf("%d\n", n.visited.size());
            pq.push(n);
//            printf("%d\n", pq.top().visited.size());
        }
        printf("Created initial nodes.\n");
        fflush(stdout);
        while (!pq.empty())
        {
            node n = pq.top();
            pq.pop();
            ++m_cNodes;
            printf("last: %d, length %.6lf. nodes: %d, vis: %d, bound: %.6lf\n",
                    n.tour.back(), n.length, n.tour.size(), n.cVisited,
                    n.bound);
            if (n.bound > m_best)
                continue;
            if (n.cVisited == m_nodes.size())
            {
//                printf("Has visited all, heading home\n");
                // Full tour
                double length = n.length +
                    getDist(n.tour.front(), n.tour.back());
                if (length < m_best)
                {
                    m_best     = length;
                    m_bestTour = n.tour;
                }
                continue;
            }
//            printf("%d\n", n.visited.size());
//            for (int i = 0; i < m_nodes.size(); ++i)
//                printf("%d ", n.visited[i]);
//            printf("\n");
//            fflush(stdout);
            for (int i = 0; i < m_nodes.size(); ++i)
            {
//                printf("trying %d: %d\n", i, n.visited[i]);
                if (!n.visited[i])
                {
                    node child;
                    child.cVisited = n.cVisited;
                    child.visited  = n.visited;
                    for (vector<int>::const_iterator it = m_nearby[i].begin();
                            it != m_nearby[i].end();
                            ++it)
                    {
                        if (!child.visited[*it])
                        {
                            child.visited[*it] = 1;
                            child.cVisited++;
                        }
                    }
                    child.length   = n.length + getDist(i, n.tour.back());
                    child.tour     = n.tour;
                    child.tour.push_back(i);
                    child.bound    = getBound(child);
                    pq.push(child);
                }
            }
        }
    }

};
// }}}

int main()
{
    vector<pair<double, double> > nodes;
    FILE *pIn = fopen("Graph3.txt", "r");
    int n;
    double d;
    fscanf(pIn, "%d %lf", &n, &d);
    for (int i = 0; i < n; ++i)
    {
        double a,b;
        fscanf(pIn, "%lf %lf", &a, &b);
        nodes.push_back(pair<double, double>(a, b));
    }
    printf("Read input:\n");
    for (vector<pair<double, double> >::iterator it = nodes.begin();
            it != nodes.end();
            ++it)
        printf("%.6lf, %.6lf\n", it->first, it->second);

    printf("Creating solver\n");
    TCPTour solver(nodes, d);

    printf("Running solver\n");
    solver.run();

    printf("Solution found after evaluating %d nodes\n", solver.getNodeCount());
    printf("Length: %.6f\n", solver.getLength());
    vector<int> tour = solver.getTour();
    for (vector<int>::const_iterator it = tour.begin();
            it != tour.end();
            ++it)
        printf("%d => ", *it);
    printf("%d\n", tour.front());
}
