import networkx as nx
from heapq import heappush, heappop
from time import monotonic
from scipy.optimize import milp, LinearConstraint, Bounds
from itertools import permutations
from scipy.sparse import dok_array
import numpy as np
import pynauty

def int_from_set(inputset : set[int]) -> int:
    return sum(2**x for x in inputset)

def set_from_int(inputint : int) -> set[int]:
    return set(i for i in range(len(bin(inputint)[:1:-1])) if bin(inputint)[:1:-1][i] == '1')

class SetTrie:

    def __init__(self):
        self.d = [dict(), False, None]

    def insert(self, set_of_ints, val):
        a = sorted(set_of_ints, reverse = True)
        t = self.d
        if self.d[2] is None:
            self.d[2] = val
        else:
            self.d[2] = max(self.d[2], val)
        while a:
            n = a.pop()
            if n not in t[0]:
                t[0][n] = [dict(), False, val]
            else:
                t[0][n][2] = max(t[0][n][2], val)
            t = t[0][n]
        t[1] = True
        return None

    def beatingSuperset(self, set_of_ints, val) -> int:
        """
        0 - no superset contained or all are worse than val
        1 - a superset with val, no better superset
        2 - a superset better than val
        """
        if not set_of_ints:
            if self.d[2] > val:
                return 2
            elif self.d[2] == val:
                return 1
            else:
                return 0
        a = sorted(set_of_ints, reverse = True)
        q = [[self.d, 0, len(a)-1],]
        b = -1
        while q and b<=val:
            t, n_min, idx = q.pop()
            if idx == 0 and a[0] in t[0]:
                b = max(b, t[0][a[0]][2])
                continue
            for n in range(n_min, a[idx]+1):
                if n in t[0]:
                    q.append([t[0][n], n+1, idx-int(n==a[idx])])
        if b > val:
            return 2
        if b == val:
            return 1
        return 0

def rcertificate(input_graph: nx.Graph|nx.DiGraph, root: int) -> tuple[bytes,tuple[int]]:
    """
    For the given rooted graph, returns a rooted certificate and the corresponding isomorphism.
    """
    V = sorted(input_graph.nodes)
    ri = V.index(root)
    N = len(input_graph)
    g = nx.relabel_nodes(input_graph, dict((V[i],i) for i in range(N)))
    g = nx.to_dict_of_lists(g)
    g = pynauty.Graph(N, directed=(type(input_graph)==nx.DiGraph), adjacency_dict=g, vertex_coloring=[{ri,},])
    cl = pynauty.canon_label(g)
    assert cl[0] == ri
    return (pynauty.certificate(g),tuple(V[x] for x in cl))

def LP_master(inputgraph:nx.Graph|nx.DiGraph, directed:int, bisect:bool, flow:bool, clint:int, cutsets:list[set[tuple[int,int]]]=[], start=None, target=None, doc_time=False) -> float:
    """
    Model the LIP/DMIS as a linear program using various settings.

    Parameters
    ----------
    inputgraph : nx.Graph|nx.DiGraph
        the input graph to operate on.
    directed : int
        0 - undirected model. input graph must be undirected or all arcs reciprocal.
        1 - directed model, no back chords allowed.
        2 - directed model, back chords allowed. This requires flow to be turned on.
    bisect : bool
        Whether to split edge constraints. The lower constraints are replaced by the split,
        as these are strictly tighter.
        The upper constraints are added (vertex bound), since they are neither strictly
        stronger nor weaker than the fused version.
        This is always an improvement and efficient, False is only supported for comparison.
        Must be true to admit directed=2.
    flow : bool
        Whether to augment the model with flow variables and constraints.
        These ensure global connectivity and are necessary to distinguish back chords.
    clint : int
        0 - LP relaxation, no clique constraints.
        1 - LP relaxation, added clique constraints.
        2 - Integer solution. No clique constraints necessary.
    start : TYPE, optional
        a vertex, a subset of vertices, or None to imply any.
    target : TYPE, optional
        a vertex, a subset of vertices, or None to imply any.

    Returns
    -------
    tuple[dict, float]
        nonzero variables in the solution, function value in the solution.

    """
    if doc_time:
        time0, time1 = monotonic(), monotonic()
    # unifying start and target lists and validating input
    sta : list
    if start is None:
        sta = sorted(inputgraph.nodes)
    elif start in inputgraph:
        sta = [start,]
    elif set(start) <= set(inputgraph.nodes):
        sta = sorted(start)
    else:
        assert False
    tar : list
    if target is None:
        tar = sorted(inputgraph.nodes)
    elif target in inputgraph:
        tar = [target,]
    elif set(target) <= set(inputgraph.nodes):
        tar = sorted(target)
    else:
        assert False
    connected = False
    for s in sta:
        for t in tar:
            if nx.has_path(inputgraph, s, t):
                connected = True
                break
        if connected:
            break
    assert connected
    del connected
    stas = set(sta)
    tars = set(tar)
    if doc_time:
        time0, time1 = time1, monotonic()
        print("connectivity check:", time1-time0)
    graph : nx.DiGraph | nx.Graph
    # building cliques. in the directed case, only pairwise biconnected counts.
    if clint == 1:
        ug : nx.Graph
        if type(inputgraph) == nx.DiGraph:
            ug = inputgraph.to_undirected(reciprocal=True)
        else:
            ug = inputgraph.copy()
        cliques = list(x for x in nx.find_cliques(ug) if len(x)>=3)
        #print(cliques)
        if directed:
            del ug
        if doc_time:
            time0, time1 = time1, monotonic()
            print("clique check:", time1-time0)
    if not directed:
        if type(inputgraph) == nx.DiGraph:
            if clint != 1:
                ug = inputgraph.to_undirected(reciprocal=True)
            # if the input graph was directed type, all arcs must be symmetric
            # to admit modeling it as undirected.
            assert ug == inputgraph.to_undirected(reciprocal=False)
        elif clint != 1:
            ug = inputgraph.copy()
        graph = ug
        del ug
    else:
        if type(inputgraph) == nx.Graph:
            graph = inputgraph.to_directed()
        else:
            graph = inputgraph.copy()
    assert (not directed and type(graph) == nx.Graph) or (directed and type(graph) == nx.DiGraph)
    V = sorted(graph.nodes)
    E = sorted(graph.edges)
    # print(E)
    Esym = list(set(tuple(sorted(e)) for e in E))
    if doc_time:
        time0, time1 = time1, monotonic()
        print("graph type casting & node list sorting", time1-time0)
    # now let's build the variable list alongside integrality and objective.
    vari, integ, objective, varid = list(), list(), list(), dict()
    ns, nt, nn, ne, nv, nes = len(sta), len(tar), len(graph), len(graph.edges), 0, len(Esym)
    for s in sta:
        variable = ("s", s)
        varid[variable] = nv
        vari.append(variable)
        nv += 1
    for t in tar:
        variable = ("t", t)
        varid[variable] = nv
        vari.append(variable)
        nv += 1
    for v in V:
        variable = ("y", v)
        varid[variable] = nv
        vari.append(variable)
        nv += 1
    for e in E:
        variable = ("x", e)
        varid[variable] = nv
        vari.append(variable)
        nv += 1
    integ += [0,]*(ns+nt) + [clint==2,]*nn + [0,]*ne
    objective += [0,]*(ns+nt) + [-1,]*nn + [0,]*ne
    if flow:
        for v in V:
            for e in E:
                variable = ("f", e, v)
                varid[variable] = nv
                vari.append(variable)
                nv += 1
                variable = ("r", e, v)
                varid[variable] = nv
                vari.append(variable)
                nv += 1
                integ += [0,0]
                objective += [0,0]
            for s in sta:
                variable = ("rs", s, v)
                varid[variable] = nv
                vari.append(variable)
                nv += 1
                integ.append(0)
                objective.append(0)
            for t in tar:
                variable = ("ft", t, v)
                varid[variable] = nv
                vari.append(variable)
                nv += 1
                integ.append(0)
                objective.append(0)
    if doc_time:
        time0, time1 = time1, monotonic()
        print("variable, objective, integrality initialization", time1-time0)
    # ==== now we build constraints ====
    cons = list()
    # start and end
    A = dok_array((2,nv),dtype=np.uint8)
    for v in sta:
        A[0,varid[("s", v)]] = 1
    for v in tar:
        A[1,varid[("t", v)]] = 1
    cons.append(LinearConstraint(A, lb=1, ub=1))
    if doc_time:
        time0, time1 = time1, monotonic()
        print("start+end constraints", time1-time0)
    # coupling node and edge variables
    A = dok_array(((1+bool(directed))*nn,nv), dtype=np.int8)
    for i in range(nn):
        if V[i] in tar:
            A[i, varid[("t",V[i])]] = 1
        if V[i] in sta:
            A[i+bool(directed)*nn, varid[("s", V[i])]] = 1
        for e in graph.edges(V[i]):
            A[i,varid[("x", (tuple(sorted(e)), e)[bool(directed)])]] = 1
        if directed:
            for e in graph.in_edges(V[i]):
                A[nn+i, varid[("x", e)]] = 1
        A[i, varid[("y",V[i])]] = -(1 + int(not bool(directed)))
        if directed:
            A[nn+i, varid[("y", V[i])]] = -1
    cons.append(LinearConstraint(A,lb=0,ub=0))
    if doc_time:
        time0, time1 = time1, monotonic()
        print("coupling node & edge variables", time1-time0)
    # upper edge constraint
    assert (bisect and flow) or directed<2
    if directed<2:
        A = dok_array((nes,nv), dtype=np.int8)
        for i in range(nes):
            A[i, varid[("y", Esym[i][0])]] = 1
            A[i, varid[("y", Esym[i][1])]] = 1
            if ("x", Esym[i]) in varid:
                A[i, varid[("x", Esym[i])]] = -1
            if ("x", (Esym[i][1], Esym[i][0])) in varid:
                A[i, varid[("x", (Esym[i][1], Esym[i][0]))]] = -1
        cons.append(LinearConstraint(A,lb=-np.inf, ub=1))
        if doc_time:
            time0, time1 = time1, monotonic()
            print("fused upper edge constraint", time1-time0)
    # lower edge constraint
    A = dok_array((nes*(1+int(bisect)),nv), dtype=np.int8)
    for i in range(nes):
        A[i, varid[("y", Esym[i][0])]] = 1
        A[i+int(bisect)*nes, varid[("y", Esym[i][1])]] = 1
        if ("x", Esym[i]) in varid:
            A[i, varid[("x", Esym[i])]] = -(2-int(bisect))
            if bisect:
                A[nes+i, varid[("x", Esym[i])]] = -1
        if ("x", (Esym[i][1], Esym[i][0])) in varid:
            A[i, varid[("x", (Esym[i][1], Esym[i][0]))]] = -(2-int(bisect))
            if bisect:
                A[nes+i, varid[("x", (Esym[i][1], Esym[i][0]))]] = -1
    cons.append(LinearConstraint(A,lb=0,ub=np.inf))
    if doc_time:
        time0, time1 = time1, monotonic()
        print("lower edge constraint", time1-time0)
    # clique constraints
    if clint == 1:
        nc = len(cliques)
        A = dok_array((nc*(1+int(bisect)),nv), dtype=np.uint8)
        for i in range(nc):
            for e in permutations(cliques[i], r=2):
                if ("x", e) in varid:
                    A[i, varid[("x", e)]] = 1
            if bisect:
                for v in cliques[i]:
                    A[nc+i, varid[("y", v)]] = 1
        ubu = [1,]*nc
        if bisect:
            ubu += [2,]*nc
        cons.append(LinearConstraint(A,lb=-np.inf, ub=ubu))
        if doc_time:
            time0, time1 = time1, monotonic()
            print("clique constraints", time1-time0)
    # cutset constraints
    if cutsets:
        ncs = len(cutsets)
        A = dok_array((ncs,nv), dtype=np.uint8)
        for i in range(ncs):
            for e in cutsets[i]:
                A[i, varid[("x", e)]] = 1
        cons.append(LinearConstraint(A,lb=1, ub=1))
        if doc_time:
            time0, time1 = time1, monotonic()
            print("cutset constraints", time1-time0)
    # flow constraints
    if flow:
        for v in V:
            # edges
            A = dok_array((2*ne,nv), dtype=np.int8)
            for i in range(ne):
                A[i, varid[("f", E[i], v)]] = 1
                A[ne+i, varid[("r", E[i], v)]] = 1
                A[i, varid[("x", E[i])]] = -1
                A[ne+i, varid[("x", E[i])]] = -1
            cons.append(LinearConstraint(A, lb=-np.inf, ub=0))
            # starts
            A = dok_array((ns,nv), dtype=np.int8)
            for i in range(ns):
                A[i, varid[("rs", sta[i], v)]] = 1
                A[i, varid[("s", sta[i])]] = -1
            cons.append(LinearConstraint(A, lb=-np.inf, ub=0))
            # ends
            A = dok_array((nt,nv), dtype=np.int8)
            for i in range(nt):
                A[i, varid[("ft", tar[i], v)]] = 1
                A[i, varid[("t", tar[i])]] = -1
            cons.append(LinearConstraint(A, lb=-np.inf, ub=0))
            # nodes
            A = dok_array((2*nn,nv), dtype=np.int8)
            for i in range(nn):
                if V[i] == v:
                    A[i, varid[("y", v)]] = 1
                    A[nn+i, varid[("y", v)]] = 1
                for e in graph.edges(V[i]):
                    #if ("x", e) in varid:
                    A[i, varid[("f", e, v)]] = -1
                    A[nn+i, varid[("r", e, v)]] = 1
                    # elif ("x", (e[1], e[0])) in varid:
                    #     A[i, varid[("f", (e[1], e[0]), v)]] = -1
                    #     A[nn+i, varid[("r", (e[1], e[0]), v)]] = 1
                    # else:
                    #     assert False
                for e in graph.in_edges(V[i]):
                    A[i, varid[("f", e, v)]] = 1
                    A[nn+i, varid[("r", e, v)]] = -1
                if V[i] in stas:
                    #A[i, varid[("rs", V[i], v)]] = 1
                    A[nn+i, varid[("rs", V[i], v)]] = -1
                if V[i] in tars:
                    A[i, varid[("ft", V[i], v)]] = -1
                    #A[nn+i, varid[("ft", V[i], v)]] = 1
            cons.append(LinearConstraint(A, lb=0, ub=0))
        # edges
        if flow and directed == 2:
            for ed in E:
                A = dok_array((2*(nn-1),nv), dtype=np.int8)
                i=0
                j=2*nn-3
                for v in V:
                    if v != ed[1]:
                        for e in graph.edges(ed[0]):
                            A[i, varid[("f", e, v)]] = -1
                        if ed[0] in tars:
                            A[i, varid[("ft", ed[0], v)]] = -1
                        for e in graph.edges(ed[1]):
                            A[i, varid[("f", e, v)]] = 1
                        if ed[1] in tars:
                            A[i, varid[("ft", ed[1], v)]] = 1
                        A[i, varid[("y", ed[0])]] = 1
                        i += 1
                    if v != ed[0]:
                        for e in graph.in_edges(ed[0]):
                            A[j, varid[("r", e, v)]] = 1
                        if ed[0] in stas:
                            A[j, varid[("rs", ed[0], v)]] = 1
                        for e in graph.in_edges(ed[1]):
                            A[j, varid[("r", e, v)]] = -1
                        if ed[1] in stas:
                            A[j, varid[("rs", ed[1], v)]] = -1
                        A[j, varid[("y", ed[1])]] = 1
                        j -= 1
                cons.append(LinearConstraint(A, lb=-np.inf, ub=1))
    if doc_time and flow:
        time0, time1 = time1, monotonic()
        print("flow constraints", time1-time0)
    pass
    # ===== finally solve the problem and return the result
    bnds : Bounds
    if bisect:
        bnds = Bounds(lb=0, ub=1)
    else:
        ubn = [1,]*nv
        for i in range(nv):
            if vari[i][0] == "y":
                ubn[i]=np.inf
        bnds = Bounds(lb=0, ub=ubn)
    if doc_time:
        time0, time1 = time1, monotonic()
        print("constructing variable bounds", time1-time0)
    #res = milp(objective, integrality=integ, bounds=bnds, constraints=cons)
    res = milp(objective, integrality=int(clint==2), bounds=bnds, constraints=cons)
    #print(res)
    if doc_time:
        time0, time1 = time1, monotonic()
        print("solving the LP", time1-time0)
    assert res.success
    return -res.fun

def deg_seq_ub(graph:nx.Graph) -> int:
    ds = sorted([d for n, d in graph.degree()], reverse=True)
    if not ds:
        return 0
    if not ds[0]:
        return 1
    while ds and ds[-1] == 0:
        ds.pop()
    used_nodes = 0
    reserved_ones = 0
    while ds and ds[-1] == 1:
        if used_nodes < 2:
            used_nodes += 1
        else:
            reserved_ones += 1
        ds.pop()
    while ds and ds[-1] == 2:
        used_nodes += 1
        ds.pop()
    if not ds:
        return used_nodes
    used_nodes += len(ds)
    for x in range(len(ds)+1):
        reserved, using = ds[:x]+[1,]*reserved_ones, ds[x:]
        while using:
            y = using.pop()
            if len(reserved) < y - 2:
                using.append(0)
                break
            for z in reserved[:y-2]:
                z -= 1
            reserved.sort(reverse=True)
            while not reserved[-1]:
                reserved.pop()
        if not using:
            return used_nodes - x
    print("ds failed on graph:",str(nx.to_dict_of_lists(graph)),"; remaining ds:",ds)
    assert False

def solve(inputgraph : nx.Graph, settings, timeout = 1200, inputstart = None):
    assert "sup" in settings
    use_supercomp = settings["sup"]
    assert "iso" in settings
    iso_upto = settings["iso"]
    assert "ub" in settings and (settings["ub"] in ("trivial", "ds", "fLP", "pLP") or type(settings["ub"])==int)
    # "trivial" - len(F)
    # "ds" - degree sequence
    # "fLP" - flat LP. this gets overwritten to pLP if decomposition is done
    # "pLP" - flat LP on pointed graph
    # if a decomposition is done, the min with longest path from decomp is taken
    assert "dec" in settings and settings["dec"] in range(3)
    # 0 - no decomposition at all
    # 1 - do the condensation and contraction, compute node weights with ds
    # 2 - also take its line graph, compute line weights with LP
    assert "bbt" in settings
    assert (not settings["bbt"]) or settings["dec"] > 0
    if settings["bbt"]:
        assert "bbt_with_ub" in settings and (settings["bbt_with_ub"] in ("ds", "LP") or type(settings["bbt_with_ub"])==int)
        assert "bbt_without_ub" in settings and (settings["bbt_without_ub"] in ("ds", "LP") or type(settings["bbt_without_ub"])==int)
    is_consistent = settings["ub"]!= "ds" and type(settings["ub"]==str) and settings["dec"]!=1 and (not settings["bbt"] or (settings["bbt_with_ub"]=="LP" and settings["bbt_without_ub"]=="LP"))
    inputnodenames = sorted(inputgraph.nodes)
    numinputnodes = len(inputnodenames)
    ig : nx.Graph = nx.relabel_nodes(inputgraph, dict(((inputnodenames[i], i) for i in range(numinputnodes))))
    start : set[int]
    if inputstart is None:
        start = set(range(numinputnodes))
    elif type(inputstart) == set:
        start = inputstart
    else:
        start = {inputnodenames.index(inputstart),}
    queue = []
    searchnodes = 0
    def mature(past : int, present : int, immature_future : set[int], parent : tuple[int,int]) -> None:
        nonlocal ig, queue, searchnodes
        # trim disconnected
        searchnodes += 1
        wg = nx.induced_subgraph(ig, nx.node_connected_component(nx.induced_subgraph(ig, immature_future), present)).copy()
        if len(wg) == 1:
            heappush(queue, (-((numinputnodes+1)*(past+1)+1), (past, present, int_from_set(wg.nodes), parent)))
            return None
        # point
        dg = None
        if settings["ub"]=="pLP" or settings["dec"]>0 or (type(settings["ub"])==int and len(wg) <= settings["ub"]):
            dg = wg.to_directed()
            for n in wg.neighbors(present):
                dg.remove_edge(n, present)
            changed = True
            while changed:
                changed = False
                for e in tuple(dg.edges):
                    dg2 = nx.induced_subgraph(dg, {e[0],}|(set(dg.nodes)-set(wg.neighbors(e[1]))))
                    if present not in dg2.nodes or not nx.has_path(dg2, present, e[0]):
                        dg.remove_edges_from((e,))
                        changed = True
        if settings["dec"]>0:
            con : nx.DiGraph = nx.condensation(dg)
            node_map = dict((v,set(con.nodes[v]['members'])) for v in con.nodes)
            in_node_map = con.graph['mapping'].copy()
            edge_map = dict((e,set()) for e in con.edges)
            con_edge_set = set(con.edges)
            for e in dg.edges:
                ce = (in_node_map[e[0]], in_node_map[e[1]])
                if ce in con_edge_set:
                    edge_map[ce].add(e)
            del in_node_map, con_edge_set
            internal_cuts = dict((v,list()) for v in con.nodes)
            changed = True
            while changed:
                changed = False
                for ce in tuple(con.edges):
                    if con.out_degree(ce[0]) == 1 and con.in_degree(ce[1]) == 1:
                        nx.contracted_edge(con, ce, self_loops=False, copy=False)
                        del con.nodes[ce[0]]['contraction']
                        con.nodes[ce[0]]['members'] |= node_map[ce[1]]
                        node_map[ce[0]] |= node_map[ce[1]]
                        internal_cuts[ce[0]].append(edge_map[ce])
                        internal_cuts[ce[0]] += internal_cuts[ce[1]]
                        del node_map[ce[1]]
                        del edge_map[ce]
                        del internal_cuts[ce[1]]
                        for e in tuple(edge_map.keys()):
                            if e[0] == ce[1]:
                                edge_map[(ce[0],e[1])] = edge_map[e]
                                del edge_map[e]
                        changed = True
                        assert con.nodes[ce[0]]['members'] == node_map[ce[0]]
                        break
            starting_comp : int
            for comp in con.nodes:
                if present in node_map[comp]:
                    starting_comp = comp
                    break
            if settings["dec"]>1:
                boundary_map = dict()
                for ce in edge_map:
                    boundary_map[ce] = (set(e[0] for e in edge_map[ce]), set(e[1] for e in edge_map[ce]))
                # at this point, boundary_map : (c0, c1) -> (set[exits of c0], set[entries of c1])
                # edge_map : (c0, c1) -> set[edges crossing from c0 to c1]
                # node_map : c -> contained vertices
                # internal_cuts : c -> list[set[tuple[int,int]]]
                lg = con.copy()
                ls : tuple[tuple[str,int],int]
                for c in con.nodes:
                    lg.add_node(("end",c))
                    lg.add_edge(c, ("end", c))
                    if present in node_map[c]:
                        lg.add_node(("start",c))
                        lg.add_edge(("start",c), c)
                        ls = (("start",c),c)
                lg : nx.DiGraph = nx.line_graph(lg)
                # nodes of this line graph have the form (c0, c1)
                # edges of this line graph have the form ((c0,c1),(c1,c2))
                # this allows computing an edge weight LP by restricting node set to c1,
                # start set to the boundary of c1 against c0,
                # and target set to the boundary of c1 against c2.
                # in the special cases where c0 = ("start", c), start is restricted to present,
                # and where c2 = ("end", c), target set is open.
                lines_from_comp : dict[int, list[tuple]] = dict()
                for line in lg.edges:
                    assert line[0][1] in con.nodes
                    if line[0][1] not in lines_from_comp:
                        lines_from_comp[line[0][1]] = list()
                    lines_from_comp[line[0][1]].append(line)
                linenodes_from_comp : dict[int, list[tuple[int,int]]] = dict()
                for linenode in lg.nodes:
                    if type(linenode[0]) == tuple:
                        assert linenode[0][1] in con.nodes
                        if linenode[0][1] not in linenodes_from_comp:
                            linenodes_from_comp[linenode[0][1]] = list()
                        linenodes_from_comp[linenode[0][1]].append(linenode)
                    else:
                        assert linenode[0] in con.nodes
                        if linenode[0] not in linenodes_from_comp:
                            linenodes_from_comp[linenode[0]] = list()
                        linenodes_from_comp[linenode[0]].append(linenode)
                    if type(linenode[1]) == tuple:
                        assert linenode[1][1] in con.nodes
                        if linenode[1][1] not in linenodes_from_comp:
                            linenodes_from_comp[linenode[1][1]] = list()
                        linenodes_from_comp[linenode[1][1]].append(linenode)
                    else:
                        assert linenode[1] in con.nodes
                        if linenode[1] not in linenodes_from_comp:
                            linenodes_from_comp[linenode[1]] = list()
                        linenodes_from_comp[linenode[1]].append(linenode)
                # compute the component LPs and assign weight
                for line in tuple(lg.edges):
                    comp = line[0][1]
                    assert comp == line[1][0]
                    if type(line[0][0]) == tuple:
                        comp_start = present
                    else:
                        comp_start = boundary_map[line[0]][1]
                    if type(line[1][1]) == tuple:
                        comp_target = None
                    else:
                        comp_target = boundary_map[line[1]][0]
                    line_weight = int(1e-6 + LP_master(nx.induced_subgraph(dg, node_map[comp]), 1, True, False, 1, cutsets=(internal_cuts[comp], [])[type(line[1][1])==tuple], start=comp_start, target=comp_target))
                    lg.add_edge(line[0], line[1], weight=line_weight)
            else:
                # compute node weights with ds instead
                for c in tuple(con.nodes):
                    node_weight = deg_seq_ub(nx.induced_subgraph(wg, node_map[c]))
                    con.add_node(c, weight=node_weight)
        # B&B trim
        if settings["bbt"]:
            changed = True
            while changed:
                changed = False
                for comp in con.nodes:
                    if comp == starting_comp:
                        continue
                    # compute an upper bound WITH the component
                    ub_with : int
                    if settings["bbt_with_ub"]=="LP" or type(settings["bbt_with_ub"])==int:
                        entering_edges = set()
                        for ce in con.in_edges(comp):
                            entering_edges |= edge_map[ce]
                    con_with = nx.induced_subgraph(con, {comp,}|set(nx.descendants(con, comp))|set(nx.ancestors(con, comp))).copy()
                    assert starting_comp in con_with
                    if settings["dec"]>1:
                        lg_with = lg.copy()
                    dg_with = dg.copy()
                    for comp1 in set(con.nodes)-set(con_with.nodes):
                        dg_with.remove_nodes_from(node_map[comp1])
                        if settings["dec"]>1:
                            lg_with.remove_nodes_from(linenodes_from_comp[comp1])
                    # flat computation
                    assert present in dg_with
                    if settings["bbt_with_ub"]=="LP" or (type(settings["bbt_with_ub"])==int and len(dg_with) <= settings["bbt_with_ub"]):
                        ub_with = past + int(1e-6 + LP_master(dg_with, 1, True, False, 1, cutsets=[entering_edges,], start=present))
                    elif settings["bbt_with_ub"]=="ds" or (type(settings["bbt_with_ub"])==int and len(dg_with) > settings["bbt_with_ub"]):
                        ub_with = past + deg_seq_ub(dg_with.to_undirected(reciprocal=False, as_view=False))
                    else:
                        assert False
                    # decomposed computation: longest weighted path in lg_with
                    longest_from : dict[int,int] = dict()
                    decomp_max = 0
                    if settings["dec"]>1:
                        for v in reversed(list(nx.topological_sort(lg_with))):
                            longest_from[v] = max((longest_from[s]+lg_with.edges[(v,s)]['weight'] for s in lg_with.successors(v)), default=0)
                            if type(v[0]) == tuple:
                                assert v[0][0] == "start"
                                decomp_max = max(decomp_max, longest_from[v])
                    else:
                        for v in reversed(list(nx.topological_sort(con_with))):
                            longest_from[v] = max((longest_from[s]+con_with.nodes[v]['weight'] for s in con_with.successors(v)), default=con_with.nodes[v]['weight'])
                            if v == starting_comp:
                                decomp_max = max(decomp_max, longest_from[v])
                    decomp_max += past
                    ub_with = min(ub_with, decomp_max)
                    # attempt greedy solve WITHOUT the component
                    lb_without = past + 1
                    con_without = nx.induced_subgraph(con, nx.node_connected_component(nx.induced_subgraph(con, set(con.nodes)-{comp,}).to_undirected(reciprocal=False, as_view=False), starting_comp)).copy()
                    dg_without : nx.DiGraph = dg.copy()
                    for comp1 in set(con.nodes)-set(con_without.nodes):
                        dg_without.remove_nodes_from(node_map[comp1])
                    start_without = present
                    if settings["bbt_without_ub"]=="LP" or (type(settings["bbt_without_ub"])==int and len(dg_without) <= settings["bbt_without_ub"]):
                        ub_without = past + int(1e-6 + LP_master(dg_without, 1, True, False, 1, start=start_without))
                    elif settings["bbt_without_ub"]=="ds" or (type(settings["bbt_without_ub"])==int and len(dg_without) > settings["bbt_without_ub"]):
                        ub_without = past + deg_seq_ub(dg_without.to_undirected(reciprocal=False, as_view=False))
                    else:
                        assert False
                    while lb_without <= ub_with and ub_with < ub_without:
                        ssu = set(dg_without.successors(start_without))
                        if not ssu:
                            break
                        sud : dict[int,tuple[int,nx.DiGraph]] = dict()
                        for succ in ssu:
                            dgwo2 = nx.induced_subgraph(dg_without, nx.node_connected_component(nx.induced_subgraph(dg_without.to_undirected(reciprocal=False, as_view=True), {succ,}|(set(dg_without.nodes)-(ssu|{start_without,}))), succ))
                            if settings["bbt_without_ub"]=="LP" or (type(settings["bbt_without_ub"])==int and len(dgwo2) <= settings["bbt_without_ub"]):
                                ubwo2 = lb_without + LP_master(dgwo2, 1, True, False, 1, start=succ)
                            elif settings["bbt_without_ub"]=="ds" or (type(settings["bbt_without_ub"])==int and len(dgwo2) > settings["bbt_without_ub"]):
                                ubwo2 = lb_without + deg_seq_ub(dgwo2.to_undirected(reciprocal=False, as_view=True))
                            sud[succ] = (ubwo2, dgwo2.copy())
                        start_without = max(sud, key= lambda x : sud[x][0])
                        lb_without += 1
                        dg_without = sud[start_without][1]
                        ub_without = int(1e-6 + sud[start_without][0])
                        assert type(ub_without) == int
                    if lb_without > ub_with:
                        # Mr. Trimmy McTrimTrim
                        con.remove_node(comp)
                        dg.remove_nodes_from(node_map[comp])
                        if settings["dec"]>1:
                            lg.remove_nodes_from(linenodes_from_comp[comp])
                        for c2 in set(con.nodes) - ({starting_comp,}|set(nx.descendants(con, starting_comp))):
                            con.remove_node(c2)
                            dg.remove_nodes_from(node_map[c2])
                            if settings["dec"]>1:
                                lg.remove_nodes_from(linenodes_from_comp[c2])
                        assert starting_comp in con
                        changed = True
                        break
        # upper bound
        upper_bound : int = past
        # flat LP
        match settings["ub"]:
            case "trivial":
                if settings["dec"] == 0:
                    upper_bound += len(wg)
                else:
                    upper_bound += len(dg)
            case "ds":
                if settings["dec"] == 0:
                    upper_bound += deg_seq_ub(wg)
                else:
                    upper_bound += deg_seq_ub(dg.to_undirected())
            case "fLP":
                if settings["dec"] == 0:
                    upper_bound += int(1e-6 + LP_master(wg, 1, True, False, 1, start=present))
                else:
                    upper_bound += int(1e-6 + LP_master(dg, 1, True, False, 1, start=present))
            case "pLP":
                upper_bound += int(1e-6 + LP_master(dg, 1, True, False, 1, start=present))
            case _:
                assert type(settings["ub"])==int
                if dg is None:
                    if len(wg) > settings["ub"]:
                        upper_bound += deg_seq_ub(wg)
                    else:
                        upper_bound += int(1e-6 + LP_master(wg, 1, True, False, 1, start=present))
                else:
                    if len(dg) > settings["ub"]:
                        upper_bound += deg_seq_ub(dg.to_undirected())
                    else:
                        upper_bound += int(1e-6 + LP_master(dg, 1, True, False, 1, start=present))
        # decomposed computation: longest weighted path in lg
        if settings["dec"] > 0:
            longest_from : dict[int,int] = dict()
            decomp_max = 0
            if settings["dec"]>1:
                for v in reversed(list(nx.topological_sort(lg))):
                    longest_from[v] = max((longest_from[s]+lg.edges[(v, s)]['weight'] for s in lg.successors(v)), default=0)
                    if type(v[0]) == tuple:
                        assert v[0][0] == "start"
                        decomp_max = max(decomp_max, longest_from[v])
            else:
                for v in reversed(list(nx.topological_sort(con))):
                    longest_from[v] = max((longest_from[s]+con.nodes[v]['weight'] for s in con.successors(v)), default=con.nodes[v]['weight'])
                    if v == starting_comp:
                        decomp_max = max(decomp_max, longest_from[v])
            decomp_max += past
            upper_bound = min(upper_bound, decomp_max)
        if dg is not None:
            heappush(queue, (-((numinputnodes+1)*upper_bound+len(dg)), (past, present, int_from_set(dg.nodes), parent)))
        else:
            heappush(queue, (-((numinputnodes+1)*upper_bound+len(wg)), (past, present, int_from_set(wg.nodes), parent)))
        return None
    solvetime = monotonic()
    if timeout is not None:
        stoptime = solvetime + timeout
    for s in start:
        if timeout is not None and monotonic() > stoptime:
            break
        mature(0, s, ig.nodes, None)
    opt_length = -1
    opt_seeds = set()
    if use_supercomp:
        domcomp : list[SetTrie] = list()
        for i in range(numinputnodes):
            domcomp.append(SetTrie())
    if iso_upto:
        idstruct : list[dict[int,list[int,list[tuple[int,int]],list[tuple[bytes,tuple[int]]]]]] = list()
    else:
        idstruct : list[dict[int,list[int,list[tuple[int,int]]]]] = list()
    for i in range(numinputnodes):
        idstruct.append(dict())
    if iso_upto:
        if not pynauty_imported:
            import pynauty
        isostruct : dict[bytes,tuple[int,list[tuple[int]]]] = dict()
        def iso_backprop(present : int, ifuture : int, paths : bool) -> set[tuple[int,int]]|list[tuple[int]]:
            """
            Propagates the given state back up the ancestor tree to apply isomorphisms.
            If paths = False, returns a set of isomorphic states.
            If paths = True, assumes singleton future, returns a set of paths.
            """
            nonlocal isostruct, idstruct
            ol, nl = {(present, ifuture):{(present,)+tuple(set_from_int(ifuture)-{present,}),},}, dict()
            term = False
            while True:
                for head in ol:
                    for parent in idstruct[head[0]][head[1]][1]:
                        if parent is None:
                            term = True
                            break
                        if parent not in nl:
                            nl[parent] = set()
                        F : tuple[int] = idstruct[parent[0]][parent[1]][2][0][1]
                        # the position of each element in F tells us where to find
                        # its replacement in iso_parent.
                        F = dict((F[i],i) for i in range(len(F)))
                        for iso_parent in isostruct[idstruct[parent[0]][parent[1]][2][0][0]][1]:
                            # iso_parent is a tuple of vertex ints.
                            ip = (iso_parent[0], int_from_set(iso_parent))
                            if ip not in nl:
                                nl[ip] = set()
                            # use them to apply a monomorphism
                            for x in ol[head]:
                                y = list(iso_parent[F[z]] for z in x)
                                if paths:
                                    y.append(iso_parent[0])
                                #else:
                                #    y.sort()
                                nl[ip].add(tuple(y))
                    if term:
                        break
                if term:
                    break
                ol, nl = nl, dict()
            assert bool(ol) and not bool(nl)
            # finally, pool the results
            if paths:
                out : set[tuple[int]] = set()
                for head in ol:
                    out |= ol[head]
                return sorted(out)
            else:
                out : set[tuple[int,int]] = set()
                for head in ol:
                    for tail in ol[head]:
                        out.add((tail[0], int_from_set(tail)))
                return out
    success = False
    while queue and (timeout is None or monotonic() < stoptime):
        if (-queue[0][0])//(numinputnodes+1) < opt_length:
            success = True
            break
        hp = heappop(queue)
        past,present,ifuture,parent = hp[1]
        assert (-hp[0])//(numinputnodes+1) >= past
        future = set_from_int(ifuture)
        #print((-hp[0])//(numinputnodes+1), past, present, "{", *future, "}")
        # check dominance
        if use_supercomp:
            status = domcomp[present].beatingSuperset(future, past)
            if status == 2:
                continue
            if status == 1 and ifuture in idstruct[present]:
                assert idstruct[present][ifuture][0] == past
                assert parent not in idstruct[present][ifuture][1]
                idstruct[present][ifuture][1].append(parent)
                continue
            assert (status == 1 and ifuture not in idstruct[present]) or status == 0
            if is_consistent:
                assert not (status==0 and ifuture in idstruct[present])
        else:
            if ifuture in idstruct[present]:
                assert past <= idstruct[present][ifuture][0] or not is_consistent
                if idstruct[present][ifuture][0] == past:
                    idstruct[present][ifuture][1].append(parent)
                if past <= idstruct[present][ifuture][0]:
                    continue
                else:
                    del idstruct[present][ifuture]
        assert ifuture not in idstruct[present]
        # integrate
        if iso_upto:
            cert, isom = rcertificate(nx.induced_subgraph(ig, future), present)
            if cert in isostruct:
                assert past <= isostruct[cert][0] or not is_consistent
                if past == isostruct[cert][0]:
                    isostruct[cert][1].append(isom)
                    idstruct[present][ifuture] = [past, [parent,], [(cert, isom),]]
                if len(future) == 1:
                    opt_seeds.add(present)
                if use_supercomp and past == isostruct[cert][0]:
                    domcomp[present].insert(future, past)
                    if past == isostruct[cert][0] and use_supercomp >= 2:
                        for x in iso_backprop(present, ifuture, False):
                            domcomp[x[0]].insert(set_from_int(x[1]), past)
                if past <= isostruct[cert][0]:
                    continue
                else:
                    del isostruct[cert]
            isostruct[cert] = (past, [isom,])
            idstruct[present][ifuture] = [past, [parent,], [(cert, isom),]]
            if use_supercomp:
                domcomp[present].insert(future, past)
                if past == isostruct[cert][0] and use_supercomp >= 2:
                    for x in iso_backprop(present, ifuture, False):
                        domcomp[x[0]].insert(set_from_int(x[1]), past)
        else:
            idstruct[present][ifuture] = [past, [parent,]]
        if use_supercomp and not iso_upto:
            domcomp[present].insert(future, past)
        if len(future) == 1:
            # terminal state
            assert opt_length == -1 or past + 1 == opt_length
            if opt_length == -1:
                opt_length = past + 1
                #print("opt_length =", opt_length)
            opt_seeds.add(present)
            continue
        # expand!
        neighbors = set(ig.neighbors(present)) & future
        for v in neighbors:
            mature(past+1, v, ({v,}|(future-neighbors))-{present,}, (present, ifuture))
    solvetime = monotonic() - solvetime
    if not queue:
        success = True
    # retrieve all solutions
    if success:
        if iso_upto:
            out = set()
            for x in opt_seeds:
                out |= set(iso_backprop(x, int_from_set({x,}), True))
            out = sorted(tuple(reversed(tuple(inputnodenames[i] for i in p))) for p in out)
        else:
            out = list()
            pathq = list()
            for x in opt_seeds:
                pathq.append([(x, int_from_set({x,})),])
            while pathq:
                path = pathq.pop()
                for p in idstruct[path[-1][0]][path[-1][1]][1]:
                    if p is not None:
                        pathq.append(path[:-1] + [path[-1][0],] + [p,])
                    else:
                        path = path[:-1] + [path[-1][0],]
                        assert len(path) == opt_length
                        out.append(tuple(path))
            out = list(tuple(reversed(tuple(inputnodenames[i] for i in p))) for p in out)
            out.sort()
        return out, searchnodes, solvetime
    elif queue:
        return [], searchnodes, (-queue[0][0])//(numinputnodes+1)
    else:
        return [], searchnodes, 0

#s = {"iso":1, "sup":1, "dec":1, "bbt":1, "bbt_with_ub":-1, "bbt_without_ub":-1, "ub":10}
#karate_paths = [('12', '3', '2', '28', '31', '25', '23', '29', '26'), ('16', '5', '0', '1', '30', '32', '23', '25', '24'), ('16', '5', '0', '1', '30', '32', '23', '27', '24'), ('16', '5', '0', '1', '30', '33', '23', '25', '24'), ('16', '5', '0', '1', '30', '33', '27', '24', '25'), ('16', '5', '0', '2', '28', '33', '23', '25', '24'), ('16', '5', '0', '2', '9', '33', '23', '25', '24'), ('16', '5', '0', '31', '24', '27', '23', '29', '26'), ('16', '6', '0', '1', '30', '32', '23', '25', '24'), ('16', '6', '0', '1', '30', '32', '23', '27', '24'), ('16', '6', '0', '1', '30', '33', '23', '25', '24'), ('16', '6', '0', '1', '30', '33', '27', '24', '25'), ('16', '6', '0', '2', '28', '33', '23', '25', '24'), ('16', '6', '0', '2', '9', '33', '23', '25', '24'), ('16', '6', '0', '31', '24', '27', '23', '29', '26'), ('17', '1', '2', '28', '31', '25', '23', '29', '26'), ('19', '1', '2', '28', '31', '25', '23', '29', '26'), ('21', '1', '2', '28', '31', '25', '23', '29', '26'), ('24', '25', '23', '32', '30', '1', '0', '5', '16'), ('24', '25', '23', '32', '30', '1', '0', '6', '16'), ('24', '25', '23', '33', '28', '2', '0', '5', '16'), ('24', '25', '23', '33', '28', '2', '0', '6', '16'), ('24', '25', '23', '33', '30', '1', '0', '5', '16'), ('24', '25', '23', '33', '30', '1', '0', '6', '16'), ('24', '25', '23', '33', '9', '2', '0', '5', '16'), ('24', '25', '23', '33', '9', '2', '0', '6', '16'), ('24', '27', '23', '32', '30', '1', '0', '5', '16'), ('24', '27', '23', '32', '30', '1', '0', '6', '16'), ('25', '24', '27', '33', '30', '1', '0', '5', '16'), ('25', '24', '27', '33', '30', '1', '0', '6', '16'), ('26', '29', '23', '25', '31', '28', '2', '1', '17'), ('26', '29', '23', '25', '31', '28', '2', '1', '19'), ('26', '29', '23', '25', '31', '28', '2', '1', '21'), ('26', '29', '23', '25', '31', '28', '2', '1', '30'), ('26', '29', '23', '25', '31', '28', '2', '3', '12'), ('26', '29', '23', '25', '31', '28', '2', '8', '30'), ('26', '29', '23', '27', '24', '31', '0', '1', '30'), ('26', '29', '23', '27', '24', '31', '0', '5', '16'), ('26', '29', '23', '27', '24', '31', '0', '6', '16'), ('26', '29', '23', '27', '24', '31', '0', '8', '30'), ('30', '1', '0', '31', '24', '27', '23', '29', '26'), ('30', '1', '2', '28', '31', '25', '23', '29', '26'), ('30', '8', '0', '31', '24', '27', '23', '29', '26'), ('30', '8', '2', '28', '31', '25', '23', '29', '26')]
#res = solve(nx.read_graphml("graphs/lip_crn/karate.graphml"), s, timeout=10)
#if bool(res[0]) and res[0] != karate_paths:
#    print(s)
#    assert set(res[0]) < set(karate_paths)
#    print(*sorted(set(karate_paths)-set(res[0])), sep='\n')
#    assert False
#print(s, res[1], res[2])
#gg = nx.grid_2d_graph(3,3)
# print(s, solve(gg, s)[1:])
#sb = {"iso":0, "sup":0, "dec":0, "bbt":0, "ub":"trivial"}
#pb = solve(gg, sb)[0]
#print()
#p = solve(gg, s)[0]
