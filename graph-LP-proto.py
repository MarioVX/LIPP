# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 18:01:41 2024

@author: pc
"""
import numpy as np
from scipy.optimize import linprog, milp, LinearConstraint, Bounds
from scipy.sparse import dok_array
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import networkx as nx


class Graph:
    def __init__(self, addict):
        self.V = list(addict.keys())
        self.V.sort()
        self.E = np.zeros((len(self.V),len(self.V)),dtype=np.uint8)
        for i in range(len(self.V)):
            for x in addict[self.V[i]]:
                self.E[i,self.V.index(x)] = 1
    
    def edges(self):
        for i in range(len(self.V)):
            for j in range(len(self.V)):
                if self.E[i,j] and i<j:
                    yield (self.V[i],self.V[j])

def grid(m,n):
    a = list((i,j) for j in range(n) for i in range(m))
    b = set(a)
    c = dict((x,[]) for x in a)
    for x in a:
        for y in ((1,0),(-1,0),(0,1),(0,-1)):
            x2 = list(x)
            x2[0] += y[0]
            x2[1] += y[1]
            x2 = tuple(x2)
            if x2 in b:
                c[x].append(x2)
    return c

def complete(n):
    return dict((i,list(set(range(n))-{i,})) for i in range(n))

def LP_base(g: Graph, start=None, end=None):
    assert start is None or start in g.V
    assert end is None or end in g.V
    V = g.V.copy()
    E = sorted(list(g.edges()))
    var_names = V + E
    c = np.array([0,]*len(V)+[-1,]*len(E))
    if start is None and end is None:
        A_eq = np.array([[1,]*len(V)+[0,]*len(E),])
        b_eq = np.array([2,])
    elif start is not None and end is not None:
        A_eq = np.zeros((2,len(V)+len(E)))
        A_eq[0,var_names.index(start)] = 1
        A_eq[1,var_names.index(end)] = 1
        b_eq = np.array([1,1])
    else:
        A_eq = np.array([[1,]*len(V)+[0,]*len(E), [0,]*(len(V)+len(E))])
        if start is not None:
            A_eq[1,var_names.index(start)] = 1
        if end is not None:
            A_eq[1,var_names.index(end)] = 1
        b_eq = np.array([2,1])
    b_ub = np.array([0,2]*len(E))
    A_ub = np.zeros((2*len(E),len(var_names)))
    for i in range(len(E)):
        A_ub[2*i, len(V)+i] = 2
        A_ub[2*i, V.index(E[i][0])] = -1
        A_ub[2*i+1, V.index(E[i][0])] = 1
        A_ub[2*i, V.index(E[i][1])] = -1
        A_ub[2*i+1, V.index(E[i][1])] = 1
        for f in E:
            if (f[0]==E[i][0] or f[0]==E[i][1] or f[1]==E[i][0] or f[1]==E[i][1]) and f!=E[i]:
                A_ub[2*i, var_names.index(f)] = -1
                A_ub[2*i+1, var_names.index(f)] = 1
    res = linprog(c, A_ub, b_ub, A_eq, b_eq, bounds=(0,1))
    assert res.success
    x = dict((var_names[i],res.x[i]) for i in range(len(var_names)) if res.x[i]>1e-6)
    return x, -res.fun#, A_ub, b_ub

def LP_base_vertex(g: Graph, start=None, end=None):
    assert start is None or start in g.V
    assert end is None or end in g.V
    V = g.V.copy()
    E = sorted(list(g.edges()))
    var_names = V + E
    c = np.array([0,]*len(V)+[-1,]*len(E))
    if start is None and end is None:
        A_eq = np.array([[1,]*len(V)+[0,]*len(E),])
        b_eq = np.array([2,])
    elif start is not None and end is not None:
        A_eq = np.zeros((2,len(V)+len(E)))
        A_eq[0,var_names.index(start)] = 1
        A_eq[1,var_names.index(end)] = 1
        b_eq = np.array([1,1])
    else:
        A_eq = np.array([[1,]*len(V)+[0,]*len(E), [0,]*(len(V)+len(E))])
        if start is not None:
            A_eq[1,var_names.index(start)] = 1
        if end is not None:
            A_eq[1,var_names.index(end)] = 1
        b_eq = np.array([2,1])
    b_ub = np.array([0,2]*len(E)+[2,]*len(V))
    A_ub = np.zeros((2*len(E)+len(V),len(var_names)))
    for i in range(len(E)):
        A_ub[2*i, len(V)+i] = 2
        A_ub[2*i, V.index(E[i][0])] = -1
        A_ub[2*i+1, V.index(E[i][0])] = 1
        A_ub[2*i, V.index(E[i][1])] = -1
        A_ub[2*i+1, V.index(E[i][1])] = 1
        for f in E:
            if (f[0]==E[i][0] or f[0]==E[i][1] or f[1]==E[i][0] or f[1]==E[i][1]) and f!=E[i]:
                A_ub[2*i, var_names.index(f)] = -1
                A_ub[2*i+1, var_names.index(f)] = 1
    for i in range(2*len(E),2*len(E)+len(V)):
        for e in E:
            if V[i-2*len(E)] in e:
                A_ub[i,var_names.index(e)] = 1
    res = linprog(c, A_ub, b_ub, A_eq, b_eq, bounds=(0,1))
    assert res.success
    x = dict((var_names[i],res.x[i]) for i in range(len(var_names)) if res.x[i]>1e-6)
    return x, -res.fun#, A_ub, b_ub

def LP_masterdraft(inputgraph:nx.Graph|nx.DiGraph, dirlevel:int, vertexc:bool, clint:int, start=None, target=None) -> tuple[dict,float]:
    """
    Compute an LP solution of the given instance using various model parameters.

    Parameters
    ----------
    inputgraph : TYPE
        The input graph, can be undirected or directed.
    dirlevel : int
        0 - only undirected
        1 - directed model, distinguishing in and out locally
        2 - multi-commodity flow, ensures connectedness
    vertexc : bool
        Whether to include constraints for each vertex too.
    clint : int
        0 - relaxed, no clique constraints
        1 - relaxed, with clique constraints
        2 - integer solution
    start : vertex, or container of vertices, or None
        the admitted start of the path. None means anywhere.
    target : vertex, or container of vertices, or None
        the admitted end of the path. None means anywhere.

    Returns
    -------
    tuple[dict,float]
        The dict is a mapping from variable names to values in the optimal solution.
        This admits visualizing the found solutions.
        Zero values are filtered out.
        The float denotes the objective value.

    """
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
    graph : nx.DiGraph | nx.Graph
    # building cliques. in the directed case, only pairwise biconnected counts.
    if clint == 1:
        ug : nx.Graph
        if type(inputgraph) == nx.DiGraph:
            ug = inputgraph.to_undirected(reciprocal=True)
        else:
            ug = inputgraph.copy()
        cliques = list(x for x in nx.find_cliques(ug) if len(x)>=3)
        if dirlevel:
            del ug
    if not dirlevel:
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
    assert (not dirlevel and type(graph) == nx.Graph) or (dirlevel and type(graph) == nx.DiGraph)
    # detect vertices relevant for constraints
    if vertexc:
        VC = [sorted(set(v for v in graph.nodes if graph.edges(v))|set(tar)),]
        if dirlevel:
            VC.append(sorted(set(v for v in graph.nodes if graph.in_edges(v))|set(sta)))
    # now let's build the variable list alongside integrality and objective. This only depends on dirlevel.
    vari = list()
    integ = list()
    objective = list()
    ns, nt, ne = len(sta), len(tar), len(list(graph.edges))
    merged_ends = (sta == tar) and (dirlevel == 0)
    if merged_ends:
        nt = 0
        vari += [("xs", s) for s in sta] + [("x", e) for e in graph.edges]
    else:
        vari += [("xs", s) for s in sta] + [("xt", t) for t in tar] + [("x", e) for e in graph.edges]
    integ += [int(clint==2),]*(ns+nt+ne)
    objective += [0,]*(ns+nt)+[-1,]*ne
    if dirlevel == 2:
        for v in graph.nodes:
            vari += [("fs", v, s) for s in sta] + [("ft", v, t) for t in tar] + [("f", v, e) for e in graph.edges]
            integ += [0,]*(ns+nt+ne)
            objective += [0,]*(ns+nt+ne)
    nv = len(vari)
    # ===== from here on, we build all the constraints.
    cons = list()
    # start and end
    if merged_ends:
        A = dok_array((1,nv),dtype=np.uint8)
        for i in range(nv):
            if vari[i][0] == "xs":
                A[0,i] = 1
        cons.append(LinearConstraint(A, lb=2, ub=2))
    else:
        A = dok_array((2,nv),dtype=np.uint8)
        for i in range(nv):
            if vari[i][0] == "xs":
                A[0,i] = 1
            elif vari[i][0] == "xt":
                A[1,i] = 1
        cons.append(LinearConstraint(A, lb=1, ub=1))
    
    # ===== finally solve the problem and return the result
    res = milp(objective, integrality=integ, bounds=Bounds(lb=0, ub=1), constraints=cons)
    assert res.success
    return list((x,res.x[x]) for x in vari if res.x[x]>1e-6), -res.fun

def LP_master(inputgraph:nx.Graph|nx.DiGraph, directed:int, bisect:bool, flow:bool, clint:int, start=None, target=None) -> tuple[dict, float]:
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
    # now let's build the variable list alongside integrality and objective.
    vari, integ, objective = list(), list(), list()
    ns, nt, nn, ne = len(sta), len(tar), len(graph), len(graph.edges)
    vari += [("s", s) for s in sta] + [("t", t) for t in tar] + [("y", v) for v in graph] + [("x", e) for e in graph.edges]
    integ += [0,]*(ns+nt) + [1,]*nn + [0,]*ne
    objective += [0,]*(ns+nt+nn) + [-1,]*ne
    if flow:
        pass
    nv = len(vari)
    # ==== now we build constraints ====
    cons = list()
    # start and end
    A = dok_array((2,nv),dtype=np.uint8)
    for i in range(nv):
        if vari[i][0] == "s":
            A[0,i] = 1
        elif vari[i][0] == "t":
            A[1,i] = 1
    cons.append(LinearConstraint(A, lb=1, ub=1))
    # coupling node and edge variables
    A = dok_array(((1+bool(directed))*nn,nv), dtype=np.int8)
    for i in range(nn):
        if V[i] in tar:
            A[i, vari.index(("t",V[i]))] = 1
        if V[i] in sta:
            A[i+bool(directed)*nn, vari.index(("s", V[i]))] = 1
        for e in graph.edges(V[i]):
            A[i,vari.index(("x", tuple(sorted(e))))] = 1
        if directed:
            for e in graph.in_edges(V[i]):
                A[nn+i, vari.index(("x", e))] = 1
        A[i, vari.index(("y",V[i]))] = -(1 + int(not bool(directed)))
        if directed:
            A[nn+i, vari.index(("y", V[i]))] = -1
    #print(A)
    cons.append(LinearConstraint(A,lb=0,ub=0))
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
    res = milp(objective, integrality=integ, bounds=bnds, constraints=cons)
    assert res.success
    return list((vari[i],res.x[i]) for i in range(nv) if res.x[i]>1e-6), -res.fun

def degree_sequence_ub(graph:nx.Graph):
    ds = sorted([d for n, d in graph.degree()])
    while ds[0] == 0:
        ds = ds[1:]
    used_nodes = 0
    while ds[0] == 1 and used_nodes <2:
        used_nodes += 1
        ds = ds[1:]
    while ds[0] == 1:
        ds = ds[1:]
    while ds[0] == 2:
        used_nodes += 1
        ds = ds[1:]
    while len(ds) > 1:
        ds[-1] -= ds[0] - 2
        used_nodes += 1
        ds = ds[1:]
        if ds[-1] <= 0 and len(ds) > 1:
            ds[-2] += ds[-1]
            ds = ds[:-1]
    return used_nodes -1

def plotgrid(m, n=None,start=None,end=None, solver=LP_base):
    g = Graph(grid(m,(n,m)[n is None]))
    res = solver(g,start,end)[0]
    b = list(res.keys())
    while type(b[0][0]) == int:
        b.pop(0)
    c = list(str(min(1.0,max(0.0,1.0-res[x]))) for x in b)
    print(res)
    lc = LineCollection(b,colors=c, linewidth=2)
    fig,ax = plt.subplots()
    ax.set_xlim(-1,m)
    ax.set_ylim(-1,(n,m)[n is None])
    ax.add_collection(lc)
    plt.show()
    return None

print(*LP_master(nx.gnp_random_graph(10, 0.2, seed=42), 1, True, False, 0), sep="\n")