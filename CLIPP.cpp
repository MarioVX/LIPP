#include <iostream>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphml.hpp>
#include <boost/graph/graph_traits.hpp>
#include <fstream>
#include <unordered_map>
#include <stack>
#include <algorithm>
#include <queue>
#include <chrono>

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> Graph;
typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
typedef boost::graph_traits<Graph>::adjacency_iterator AdjacencyIterator;

static Graph readGraphML(const std::string& filename) {
    Graph g;
    std::ifstream file(filename.c_str());
    if (!file.is_open()) {
        std::cerr << "Error: Could not open the file." << std::endl;
        return g;
    }

    // Reading the GraphML file into the graph g
    boost::dynamic_properties dp;
    dp.property("id", boost::get(boost::vertex_index, g)); // to read vertex ids
    boost::read_graphml(file, g, dp);

    return g;
}

static std::vector<Vertex> getNeighbors(const Graph& g, Vertex v) {
    std::vector<Vertex> neighbors;
    std::pair<AdjacencyIterator, AdjacencyIterator> adjItr = boost::adjacent_vertices(v, g);
    for (AdjacencyIterator itr = adjItr.first; itr != adjItr.second; ++itr) {
        neighbors.push_back(*itr);
    }
    return neighbors;
}

struct State;
static std::vector<Vertex> getStateNeighbors(const Graph& g, State s, Vertex v);
static std::vector<bool> getConnected(const Graph& g, size_t num_vertices, State s, Vertex root);

// The core part of a graph state, used to identify it
struct State {
    std::vector<bool> F;
    Vertex v;
    State(size_t num_vertices, const Graph& g, Vertex present) : F(num_vertices,true), v(present) {
        F = ::getConnected(g, num_vertices, *this, present);
    }
    State(State other, size_t num_vertices, const Graph& g, Vertex new_present) : F(other.F), v(new_present) {
        for (Vertex neighbor : ::getStateNeighbors(g, other, other.v)) {
            F[neighbor] = false;
        }
        F[new_present] = true;
        F[other.v] = false;
        F = ::getConnected(g, num_vertices, *this, new_present);
    }
    std::vector<Vertex> vertices(const Graph& g) {
        std::vector<Vertex> result;
        std::pair<Graph::vertex_iterator, Graph::vertex_iterator> gverts = boost::vertices(g);
        for (Graph::vertex_iterator vi = gverts.first; vi != gverts.second; ++vi) {
            if (this->F[*vi]) {
                result.push_back(*vi);
            }
        }
        return result;
    }
};

bool operator== (const State& left, const State& right) {
    return (left.F == right.F && left.v == right.v);
}

struct StateHash {
    std::size_t operator()(const State& s) const {
        std::size_t hash = 0;
        for (bool bit : s.F) {
            hash = (hash << 1) ^ std::hash<bool>()(bit);
        }
        return hash / (1 + s.v);
    }
};

static std::vector<Vertex> getStateNeighbors(const Graph& g, State s, Vertex v) {
    std::vector<Vertex> sneighbors;
    std::vector<Vertex> neighbors = getNeighbors(g, v);
    for (Vertex neighbor : neighbors) {
        if (s.F[neighbor]) {
            sneighbors.push_back(neighbor);
        }
    }
    return sneighbors;
}

// get all vertices in the same connected component in the graph state as the root vertex. used to prune unreachable vertices from states.
static std::vector<bool> getConnected(const Graph& g, size_t num_vertices, State s, Vertex root) {
    std::vector<bool> seen(num_vertices, false);
    std::stack<Vertex> discovered;
    discovered.push(root);
    seen[root] = true;
    while (!discovered.empty()) {
        Vertex v = discovered.top();
        discovered.pop();
        std::vector<Vertex> neighbors = getStateNeighbors(g, s, v);
        for (Vertex v2 : neighbors) {
            if (!seen[v2]) {
                discovered.push(v2);
                seen[v2] = true;
            }
        }
    }
    return seen;
}

// compute an upper bound on the final maximum path length based on the graph state's degree sequence
static unsigned int degreeSequenceUB(const Graph& g, State s) {
    std::vector<Vertex> V = s.vertices(g);
    size_t len = V.size();
    if (len == 0) {
        return 0;
    }
    std::vector<unsigned int> ds(len, 0);
    for (int i = 0; i < len; ++i) {
        ds[i] = getStateNeighbors(g, s, V[i]).size();
    }
    std::sort(ds.begin(), ds.end());
    if (ds.back() == 0) {
        return 1;
    }
    while (!ds.empty() && ds.front() == 0) {
        ds.erase(ds.begin());
    }
    unsigned int used = 0;
    unsigned int reserved_ones = 0;
    while (!ds.empty() && ds.front() == 1) {
        if (used < 2) {
            used++;
        }
        else {
            reserved_ones++;
        }
        ds.erase(ds.begin());
    }
    while (!ds.empty() && ds.front() == 2) {
        used++;
        ds.erase(ds.begin());
    }
    if (ds.empty()) {
        return used;
    }
    for (int i = ds.size() - 1; i >= 0; --i) {
        std::vector<unsigned int> Using(ds.begin(), ds.begin() + i);
        std::vector<unsigned int> Reserved(ds.begin() + i, ds.end());
        for (unsigned int j = 0; j < reserved_ones; ++j) {
            Reserved.emplace(Reserved.begin(), 1);
        }
        while (!Using.empty()) {
            int y = Using.back() - 2;
            if (Reserved.size() < y) {
                break;
            }
            Using.pop_back();
            for (int k = 0; k < y; ++k) {
                Reserved[Reserved.size() - 1 - k] -= 1;
            }
            std::sort(Reserved.begin(), Reserved.end());
            while (!Reserved.empty() && Reserved.front() == 0) {
                Reserved.erase(Reserved.begin());
            }
        }
        if (Using.empty()) {
            return used + i;
        }
    }
}

// generate all the successor states of a given state
std::vector<std::shared_ptr<State>> successors(const Graph& g, State s, size_t num_vertices) {
    std::vector<std::shared_ptr<State>> out;
    for (Vertex newStart : getStateNeighbors(g, s, s.v)) {
        out.emplace_back(std::make_shared<State>(s, num_vertices, g, newStart));
    }
    return out;
}


// graph state augmented with quality info
struct EvaluatedState {
    std::shared_ptr<State> s;
    unsigned int past;
    unsigned int ub;

    EvaluatedState(const Graph& g, std::shared_ptr<State> fixed_state, unsigned int past) : past(past), s(fixed_state), ub(degreeSequenceUB(g, *fixed_state)) {}
};

struct StateCompare {
    bool operator()(EvaluatedState left, EvaluatedState right) {
        if (left.past + left.ub > right.past + right.ub) {
            return false;
        }
        else if (left.past + left.ub < right.past + right.ub) {
            return true;
        }
        else {
            if (left.past >= right.past) {
                return false;
            }
            else if (left.past < right.past) {
                return true;
            }
        }
    }
};

int main() {
    std::string filename = "C:/Users/pc/Desktop/Mario_LIPP/lip_crn/usair.graphml";
    Graph g = readGraphML(filename);

    auto started = std::chrono::high_resolution_clock::now();

    std::pair<Graph::vertex_iterator, Graph::vertex_iterator> vertices = boost::vertices(g);
    
    size_t num_vertices = 0;
    for (Graph::vertex_iterator vi = vertices.first; vi != vertices.second; ++vi) {
        num_vertices++;
    }
    
    std::priority_queue <EvaluatedState, std::vector<EvaluatedState>, StateCompare> queue;
    for (Graph::vertex_iterator vi = vertices.first; vi != vertices.second; ++vi) {
        queue.emplace(EvaluatedState(g, std::make_shared<State>(num_vertices, g, *vi), 0));
    }

    std::unordered_map<State, EvaluatedState, StateHash> idstruct;

    while (!queue.empty()) {
        EvaluatedState es = queue.top();
        queue.pop();
        if (idstruct.find(*es.s) != idstruct.end()) {
            // v,F already exists
            if (es.past <= idstruct.find(*es.s)->second.past) {
                continue;
            }
            //else {
                //std::cout << "warning: re-found same v,F with better past at later point." << std::endl;
            //}
        }
        idstruct.insert(std::make_pair(*es.s, es));
        unsigned int nF = es.s->vertices(g).size();
        if (nF == 1) {
            // state is terminal. successful completion, return maximum path length
            auto ended = std::chrono::high_resolution_clock::now();
            std::cout << es.past + 1 << " " << std::chrono::duration_cast<std::chrono::milliseconds>(ended-started).count() << std::endl;
            break;
        }
        else if (nF < 1) {
            std::cout << "warning: got to a state with less than one vertices. this should not happen." << std::endl;
            return 1;
        } else {
            // state is nonterminal. expand and enqueue its successor states.
            for (auto s2 : successors(g, *es.s, num_vertices)) {
                queue.emplace(EvaluatedState(g, s2, es.past + 1));
            }
        }
    }

    return 0;
}
