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
from time import monotonic
from itertools import permutations


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

def LP_master(inputgraph:nx.Graph|nx.DiGraph, directed:int, bisect:bool, flow:bool, clint:int, start=None, target=None, doc_time=False) -> tuple[dict, float]:
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
    objective += [0,]*(ns+nt+nn) + [-1,]*ne
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
    return sorted((vari[i],res.x[i]) for i in range(nv) if res.x[i]>1e-6), -res.fun

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

def deg_seq_2(graph:nx.Graph) -> int:
    ds = sorted([d for n, d in graph.degree()], reverse=True)
    if not ds:
        return 0
    if not ds[0]:
        return 1
    while ds and ds[-1] == 0:
        ds.pop()
    used_nodes = 0
    while ds and ds[-1] == 1 and used_nodes < 2:
        used_nodes += 1
        ds.pop()
    while ds and ds[-1] == 2:
        used_nodes += 1
        ds.pop()
    if not ds:
        return used_nodes
    starting_un = used_nodes
    starting_ds = ds.copy()
    while len(starting_ds) > 1:
        starting_ds[0] -= starting_ds[-1] - 2
        starting_un += 1
        starting_ds.pop()
        if len(starting_ds) > 1 and starting_ds[0] <= 0:
            starting_ds[1] += starting_ds[0]
            starting_ds.pop(0)
    starting_un -= 1
    starting_un -= used_nodes
    # print("starting used nodes:", starting_un)
    del starting_ds
    used_nodes += len(ds)
    for x in range(len(ds)-starting_un, len(ds)+1):
        # print("x =",x)
        reserved, using = ds[:x], ds[x:]
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
        if using:
            continue
        return used_nodes - x
    assert False

def plotgrid(res):
    b = res.copy()
    xmax = 0
    ymax = 0
    for v in b.copy():
        if v[0][0] != 'x':
            b.remove(v)
        if v[0][0] == 'y':
            xmax = max(xmax, v[0][1][0])
            ymax = max(ymax, v[0][1][1])
    #print(b)
    c = list(str(min(1.0,max(0.0,1.0-b[i][1]))) for i in range(len(b)))
    #print(b)
    #print(c)
    lc = LineCollection(list(x[0][1] for x in b),colors=c, linewidth=2)
    fig,ax = plt.subplots()
    ax.set_xlim(0,xmax)
    ax.set_ylim(0,ymax)
    ax.add_collection(lc)
    plt.xticks(range(xmax+1))
    plt.yticks(range(ymax+1))
    plt.grid()
    plt.show()
    return None

def find_discrepancy(generator):
    for graph in generator():
        if len(graph) > 2 and nx.is_connected(graph):
            yf, yi = LP_master(graph, 0, True, False, 1), LP_master(graph, 0, True, False, 2)
            if yf[1] - yi[1] > 0:
                print(yf)
                print(yi)
                yield graph

#tgrid = nx.grid_2d_graph(3,4)
#res = LP_master(tgrid, 1, True, False, 1)
#print(res[1])
#plotgrid(res[0])
#print(LP_master(nx.grid_2d_graph(3,3), 1, True, True, 1)[1])
#disexample = nx.from_dict_of_lists({0:[1,],1:[0,2,3],2:[1,],3:[1,],4:[5,],5:[4,]})
