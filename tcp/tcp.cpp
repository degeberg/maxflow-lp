#include <queue>
#include <algorithm>
#include <vector>
#include <cmath>
#include <cstdio>
#include <glpk.h>
#include <set>

#define MAX_VERTICES 100

using namespace std;

const double infty   = 1000000000.0;
const double epsilon = 0.0001;

class node
{
public:
    vector<int>  tour;
    set<int>     visited;
    double       length;
    double       bound;
};

bool operator<(const node& lhs, const node& rhs)
{
    // We want min. heap behaviour!
    double diff = lhs.bound - rhs.bound;
    if (diff <= epsilon && diff >= -epsilon)
        return lhs.visited.size() < rhs.visited.size();
    return diff > 0.0;
}

// {{{ TCPTour class
class TCPTour
{
public:
    TCPTour(const vector<pair<double, double> > &nodes, double dist)
        : m_dist(dist), m_best(infty)
    {
        // Store information about the problem
        m_cNodes = nodes.size();
        m_nodes  = nodes;
        printf("Stored %d vertices\n", m_nodes.size());

        // Initialize nearby vertices
        for (int i = 0; i < m_cNodes; ++i)
        {
            m_distances[i][i] = 0.0;
            m_nearby[i].push_back(i);
            for (int j = 0; j < i; ++j)
            {
                m_distances[i][j] = m_distances[j][i] = getDist(i,j);
                // No need to calculate this more than once
                if (m_distances[i][j] <= m_dist + epsilon)
                {
                    m_nearby[i].push_back(j);
                    m_nearby[j].push_back(i);
                }
            }
        }
        // Find potential start vertex if we can define one.
        m_start = -1;
        for (int i = 0; i < m_cNodes; ++i)
            if (m_nearby[i].size() == 1)
                m_start = i;
        // Used by lp solver.
        for (int i = 0; i < m_cNodes * m_cNodes; ++i)
        {
            m_coef[i] = 1.0;
            if (i <= m_cNodes - 1 || i > 3 * m_cNodes - 3)
                m_coef2[i] = 1.0;
            else
                m_coef2[i] = -1.0;
        }
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
        return m_cCheck;
    }

    /* Used to create data for tcp.mod */
    void printDistances()
    {
        for (int i = 0; i < m_cNodes; ++i)
        {
            for (int j = 0; j < m_cNodes; ++j)
            {
                int a = i == m_start ? 0 : i == 0 ? m_start : i;
                int b = j == m_start ? 0 : j == 0 ? m_start : j;
                if (i == j)
                    printf("%d\t%d\t%d\n", a + 1, a + 1, 1000000);
                else
                    printf("%d\t%d\t%d\n", a + 1, b + 1, (int)(m_distances[i][j] * 1000));
            }
        }
    }

private:
    // List of all the vertices. Note that all vertices are connected
    // so we don't need an adjacency list/matrix
    vector<pair<double, double> > m_nodes;
    // m_nearby[i] is all vertices at most m_dist distance away from vertex i
    vector<int> m_nearby[MAX_VERTICES];
    // Distance to view a monument
    double      m_dist;     // d
    double      m_best;     // Length of the best tour
    vector<int> m_bestTour; // best tour
    int         m_cCheck;   // Nodes evaluated
    int         m_cNodes;   // n
    // m_distances[i][j] = getDist(i,j) = The euclidian distance between
    // nodes i and j.
    double      m_distances[MAX_VERTICES][MAX_VERTICES];
    // A huge array of just 1.0's. Just so we don't have to create it each
    // time we create a linear problem.
    double      m_coef[MAX_VERTICES * MAX_VERTICES];
    double      m_coef2[MAX_VERTICES * MAX_VERTICES];
    // If any vertex has only itself as nearby monuments it MUST be in the
    // path. m_start is this vertex, as we can pick it as the start of the
    // path. If m_start == -1 no vertex is by itself.
    int         m_start;

    double getDist(int i, int j)
    {
        double dist = sqrt((m_nodes[i].first  - m_nodes[j].first)  *
                    (m_nodes[i].first  - m_nodes[j].first)  +
                    (m_nodes[i].second - m_nodes[j].second) *
                    (m_nodes[i].second - m_nodes[j].second));
        return dist;
    }

    // {{{ getBound_err
    // This method was slow and produced some weird output.
    // It tries to add the subtour constraint to the bound as well.
    double getBound_err(const node& n)
    {
        glp_prob *pLP = glp_create_prob();
        glp_set_obj_dir(pLP, GLP_MIN);
        // There is n*(n-1) / 2 pairs (i,j), i > j. These are the variables
        int pairs = m_cNodes * m_cNodes;
#define v(i, j) (i * m_cNodes + j)
        glp_add_cols(pLP, pairs);
        for (int i = 0; i < m_cNodes; ++i)
        {
            // Decision variables. Relaxed to 0 <= x_ij <= 1
            for (int j = 0; j < m_cNodes; ++j)
            {
                glp_set_obj_coef(pLP, v(i,j) + 1, m_distances[i][j]);
                glp_set_col_bnds(pLP, v(i,j) + 1, GLP_DB, 0.0, 1.0);
            }
        }
        // Edges already in the tour _must_ be included.
        for (int i = 1; i < n.tour.size(); ++i)
            glp_set_col_bnds(pLP, v(n.tour[i-1], n.tour[i]) + 1, GLP_FX, 1.0, 1.0);

        int ind[m_cNodes][m_cNodes * 2 + 1];
        glp_add_rows(pLP, m_cNodes * 2);
        for (int i = 0; i < m_cNodes; ++i)
        {
            // All nodes should have exactly 0 or 2 edges incident. Relaxed
            // to 0 <= deg <= 2
            bool inPath = false;
            for (vector<int>::const_iterator it = n.tour.begin();
                    it != n.tour.end(); ++it)
                if (*it  == i)
                    inPath = true;
            if (inPath)
                glp_set_row_bnds(pLP, 2*i+1, GLP_FX, 2.0, 2.0);
            else
                glp_set_row_bnds(pLP, 2*i+1, GLP_DB, 0.0, 2.0);
            glp_set_row_bnds(pLP, 2*i+2, GLP_FX, 0.0, 0.0);

            for (int j = 0; j < i; ++j)
            {
                ind[i][j+1]          = v(i,j) + 1;
                ind[i][j + m_cNodes] = v(j,i) + 1;
            }
            for (int j = i+1; j < m_cNodes; ++j)
            {
                ind[i][j]                = v(i,j) + 1;
                ind[i][j + m_cNodes - 1] = v(j,i) + 1;
            }
            glp_set_mat_row(pLP, 2*i + 1, 2*m_cNodes - 2, ind[i], m_coef);
            glp_set_mat_row(pLP, 2*i + 2, 2*m_cNodes - 2, ind[i], m_coef2);
        }
        for (int i = 0; i < m_cNodes; ++i)
        {
            // Add vicinity constraints for each node. We do this another
            // way than in the ILP because we don't have vertex variables.
            // Instead we say the nearby vertices must total have 2 edges
            // incident.
            // We have already calculated all the indices for these vertices,
            // so we just filter duplicates and store them in a bigger array.
            int    x = glp_add_rows(pLP, 1);
            set<int> s;
            for (int j = 0; j < m_nearby[i].size(); ++j)
                for (int k = 1; k < 2 * m_cNodes - 1; ++k)
                    s.insert(ind[m_nearby[i][j]][k]);
            int    ind2[s.size() + 1];
            int j = 1;
            for (set<int>::const_iterator it = s.begin();
                    it != s.end();
                    ++it, ++j)
                ind2[j] = *it;
            glp_set_mat_row(pLP, x, s.size(), ind2, m_coef);
            glp_set_row_bnds(pLP, x, GLP_LO, 2.0, 2.0);
        }
        // Subtour constraint. this should be polynomial in size.
        if (m_start != -1)
        {
            int k = glp_add_cols(pLP, pairs);
            for (int i = 0; i < m_cNodes; ++i)
            {
                for (int j = 0; j < m_cNodes; ++j)
                {
                    int y = v(i,j);
                    glp_set_col_bnds(pLP, k + y, GLP_LO, 0.0, 0.0);
                    int x = glp_add_rows(pLP, 1);
                    glp_set_row_bnds(pLP, x, GLP_UP, 0.0, 0.0);
                    int ind2[]  = {0, k + y, y + 1};
                    double coef2[] = {0, 1.0, (double)(-m_cNodes + 1)};
                    glp_set_mat_row(pLP, x, 2, ind2, coef2);
                }
            }
            int x = glp_add_rows(pLP, m_cNodes);
            int ind2[MAX_VERTICES * MAX_VERTICES + 3 * MAX_VERTICES];
            for (int i = 0; i < m_cNodes; ++i)
            {
                glp_set_row_bnds(pLP, x + i, GLP_FX, 0.0, 0.0);
                for (int j = 0; j < i; ++j)
                {
                    ind2[j + 1]   = k + v(j,i);
                    ind2[j + m_cNodes] = k + v(i,j);
                    ind2[j + 2*m_cNodes - 1] = v(i,j) + 1;
                }
                for (int j = i+1; j < m_cNodes; ++j)
                {
                    ind2[j]   = k + v(j,i);
                    ind2[j + m_cNodes - 1] = k + v(i,j);
                    ind2[j + 2*m_cNodes - 2] = v(i,j) + 1;
                }
                if (i != m_start)
                    glp_set_mat_row(pLP, x + i, 3 * m_cNodes - 3, ind2, m_coef2);
                else
                {
                    int m = 3 * m_cNodes - 3 + 1;
                    for (int j = 0; j < m_cNodes; ++j)
                    {
                        for (int l = 0; l < m_cNodes; ++l)
                        {
                            if (j == i || l == i)
                                continue;
                            ind2[m++] = v(j,l) + 1;
                        }
                    }

                    glp_set_mat_row(pLP, x + i, m - 3 * m_cNodes + 3,
                            ind2, m_coef2);
                    glp_set_row_bnds(pLP, x + i, GLP_FX, 1.0, 1.0);
                }
            }
        }
#undef v
        // We don't want all the output from glpk.
        glp_smcp params;
        glp_init_smcp(&params);
        params.msg_lev = GLP_MSG_ERR;
        glp_simplex(pLP, &params);
        double dRet = glp_get_obj_val(pLP);
        glp_delete_prob(pLP);
        return dRet;
    }
    // }}}

    // {{{ getBound
    // Get bound using the method described in the report
    double getBound(const node& n)
    {
        glp_prob *pLP = glp_create_prob();
        glp_set_obj_dir(pLP, GLP_MIN);
        // There is n*(n-1) / 2 pairs (i,j), i > j. These are the variables
        glp_add_cols(pLP, (m_cNodes * (m_cNodes - 1))/2);
        for (int i = 0; i < m_cNodes; ++i)
        {
            // Decision variables. Relaxed to 0 <= x_ij <= 1
            for (int j = 0; j < i; ++j)
            {
                glp_set_obj_coef(pLP, (i * (i-1))/2 + j + 1,
                        m_distances[i][j]);
                glp_set_col_bnds(pLP, (i * (i-1))/2 + j + 1, GLP_DB, 0.0, 1.0);
            }
        }
        // Edges already in the tour _must_ be included.
        for (int i = 1; i < n.tour.size(); ++i)
        {
            int a = max(n.tour[i], n.tour[i-1]);
            int b = min(n.tour[i], n.tour[i-1]);
            glp_set_col_bnds(pLP, (a * (a-1))/2 + b + 1, GLP_FX, 1.0, 1.0);
        }

        int ind[m_cNodes][m_cNodes];
        glp_add_rows(pLP, m_cNodes);
        for (int i = 0; i < m_cNodes; ++i)
        {
            // All nodes should have exactly 0 or 2 edges incident. Relaxed
            // to 0 <= deg <= 2
            glp_set_row_bnds(pLP, i+1, GLP_DB, 0.0, 2.0);
            for (int j = 0; j < i; ++j)
                ind[i][j+1]  = (i * (i-1))/2 + j + 1;
            for (int j = i+1; j < m_cNodes; ++j)
                ind[i][j]  = (j * (j-1))/2 + i + 1;
            glp_set_mat_row(pLP, i + 1, m_cNodes - 1, ind[i], m_coef);
        }
        for (int i = 0; i < m_cNodes; ++i)
        {
            // Add vicinity constraints for each node. We do this another
            // way than in the ILP because we don't have vertex variables.
            // Instead we say the nearby vertices must total have 2 edges
            // incident.
            // We have already calculated all the indices for these vertices,
            // so we just filter duplicates and store them in a bigger array.
            int    x = glp_add_rows(pLP, 1);
            set<int> s;
            for (int j = 0; j < m_nearby[i].size(); ++j)
                for (int k = 1; k < m_cNodes; ++k)
                    s.insert(ind[m_nearby[i][j]][k]);
            int    ind2[s.size() + 1];
            int j = 1;
            for (set<int>::const_iterator it = s.begin();
                    it != s.end();
                    ++it, ++j)
                ind2[j] = *it;
            glp_set_mat_row(pLP, x, s.size(), ind2, m_coef);
            glp_set_row_bnds(pLP, x, GLP_LO, 2.0, 2.0);
        }
        // We don't want all the output from glpk.
        glp_smcp params;
        glp_init_smcp(&params);
        params.msg_lev = GLP_MSG_ERR;
        glp_simplex(pLP, &params);
        double dRet = glp_get_obj_val(pLP);
        glp_delete_prob(pLP);
        return dRet;
    }
    // }}}

    // {{{ BnB procedure
    void TCP()
    {
        m_cCheck = 0;
        m_best   = infty;
        priority_queue<node> pq;
        // Create the initial nodes
        printf("Creating initial nodes\n");
        fflush(stdout);
        for (int i = 0; i < m_nodes.size(); ++i)
        {
            node n;
            n.length = 0.0;
            for (vector<int>::const_iterator it = m_nearby[i].begin();
                    it != m_nearby[i].end();
                    ++it)
                n.visited.insert(*it);
            n.tour.push_back(i);
            n.bound = getBound(n);
            pq.push(n);
        }
        // Run the algorithm with the initial nodes in the pq.
        printf("Created initial nodes.\n");
        fflush(stdout);
        while (!pq.empty())
        {
            node n = pq.top();
            pq.pop();
            // Status info :)
            printf("b: %.3lf, l: %.3lf, n: %d, v: %d - best: %.3lf\n",
                    n.bound, n.length, n.tour.size(), n.visited.size(),
                    m_best);
            ++m_cCheck;
            if (n.bound > m_best - epsilon)
                // If this bound is worse than the best, all the rest will be
                // A node will never have a better lower bound than its parent
                break;
            for (int i = 0; i < m_cNodes; ++i)
            {
                // Is this node in the tour?
                bool found = false;
                for (vector<int>::const_iterator it = n.tour.begin();
                        it != n.tour.end(); ++it)
                    if (*it == i)
                        found = true;
                if (!found)
                {
                    // It wasn't. Create the child node
                    node child;
                    child.tour    = n.tour;
                    child.tour.push_back(i);
                    child.length  = n.length + m_distances[i][n.tour.back()];
                    child.visited = n.visited;
                    for (vector<int>::const_iterator it = m_nearby[i].begin();
                            it != m_nearby[i].end(); ++it)
                        child.visited.insert(*it);
                    if (child.visited.size() == m_cNodes)
                    {
                        // Full tour. Evaluating this here rather than later
                        // can save us a small amount of space and time.
                        double length = child.length +
                            m_distances[child.tour.front()][i];
                        if (length < m_best)
                        {
                            m_best     = length;
                            m_bestTour = child.tour;
                        }
                        continue;
                    }
                    child.bound   = getBound(child);
                    pq.push(child);
                }
            }
        }
    }
    // }}}

};
// }}}

int main(int argc, char **argv)
{
    vector<pair<double, double> > nodes;
    FILE *pIn = fopen(argv[1], "r");
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

    printf("Creating solver\n");
    TCPTour solver(nodes, d);

//    solver.printDistances();

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
