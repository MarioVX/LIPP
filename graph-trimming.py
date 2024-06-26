# -*- coding: utf-8 -*-
"""
Created on Sat Jun 22 10:50:17 2024

@author: pc
"""

import networkx as nx

def dead(g: nx.Graph, s, t):
    V = tuple(sorted(g.nodes))
    assert s in V and t in V
    res = set()
    stack = [[s,set(V)],]
    while stack:
        state = stack.pop()
        assert state[0] in state[1]
        if t in state[1] and nx.has_path(nx.induced_subgraph(g, state[1]), state[0], t):
            res.add(state[0])
            if state[0] != t:
                neighbors = set(nx.all_neighbors(g, state[0])) & state[1]
                for v in neighbors:
                    stack.append([v,state[1] - ({state[0],}|(neighbors-{v,}))])
                    #print(stack[-1])
    return sorted(set(V)-res)

from itertools import combinations, chain

def minsep(gr: nx.Graph, s, t):
    assert s in gr.nodes and t in gr.nodes and nx.has_path(gr,s,t)
    g = (nx.induced_subgraph(gr, nx.node_connected_component(gr,s))).copy()
    V = set(g.nodes) - {s,t}
    res = set()
    for S in chain.from_iterable(combinations(V, r) for r in range(nx.node_connectivity(g,s,t), len(V)+1)):
        if not nx.has_path(nx.induced_subgraph(g, (V|{s,t})-set(S)), s, t):
            minimal = True
            for x in res:
                if x <= frozenset(S):
                    minimal = False
                    break
            if minimal:
                res.add(frozenset(S))
                yield frozenset(S)

def recog_dead(gr: nx.Graph, s, t, draw=False):
    V = set(gr.nodes)
    assert s in V and t in V and nx.has_path(gr, s, t)
    if gr.has_edge(s,t):
        return sorted(V-{s,t})
    g = nx.induced_subgraph(gr, nx.node_connected_component(gr, s))
    res = V - set(g.nodes)
    V = set(g.nodes)
    assert s in V and t in V and nx.has_path(g, s, t)
    #S = nx.all_node_cuts(g)
    # for sep in minsep(g, s, t):
    #     if s not in sep:
    #         A = nx.node_connected_component(nx.induced_subgraph(g, V-set(sep)), s)
    #         j = V.copy()
    #         B = j - (A|sep)
    #         for v in sep:
    #             j &= {v,}|set(g.neighbors(v))
    #         j -= A
    #         for v in j&B:
    #             if v == t or not nx.has_path(nx.induced_subgraph(g, V-{v,}), s, t):
    #                 res |= j&B - {v,}
    assert s not in res
    assert t not in res
    assert nx.has_path(nx.induced_subgraph(g, V-res), s, t)
    dg = g.to_directed()
    for n in g.neighbors(s):
        dg.remove_edges_from(((n, s),))
        dg.remove_edges_from(list((n2,n) for n2 in set(g.neighbors(n))-{s,}))
    for n in g.neighbors(t):
        dg.remove_edges_from(((t, n),))
        dg.remove_edges_from(list((n,n2) for n2 in set(g.neighbors(n))-{t,}))
    changed = True
    while changed:
        changed = False
        for e in tuple(dg.edges):
            # if s in e or t in e:
            #     continue
            dg2 = nx.induced_subgraph(dg, set(dg.nodes)-(set(g.neighbors(e[0]))-{e[1],}))
            dg3 = nx.induced_subgraph(dg, set(dg.nodes)-(set(g.neighbors(e[1]))-{e[0],}))
            if t not in dg2.nodes or not nx.has_path(dg2, e[1], t):
                dg.remove_edges_from((e,))
                changed = True
                continue
            if s not in dg3.nodes or not nx.has_path(dg3, s, e[0]):
                dg.remove_edges_from((e,))
                changed = True
    if draw:
        nx.draw(dg, with_labels=True)
    dg = nx.condensation(dg)
    # if draw:
    #     nx.draw(dg)
    for c in dg.nodes:
        if not nx.has_path(dg, dg.graph['mapping'][s], c) or not nx.has_path(dg, c, dg.graph['mapping'][t]):
            res |= set(dg.nodes[c]['members'])
    return sorted(res)

def recog2(gr: nx.Graph, s, t, depth=1, report_vertex=None):
    V = set(gr.nodes)
    assert s in V and t in V and nx.has_path(gr, s, t)
    if gr.has_edge(s, t):
        return sorted(V-{s,t})
    g = nx.induced_subgraph(gr, nx.node_connected_component(gr, s))
    res = V - set(g.nodes)
    V = set(g.nodes)
    for v in V-{s,t}:
        stack = [[[v,], depth, g.copy()],]
        while stack:
            node = stack.pop()
            if (node[0][0] == s and node[0][-1] == t) or node[1]==0:
                stack.append(True)
                if report_vertex is not None and v == report_vertex:
                    print(node[0])
                break
                #if (node[0][0]==s or nx.has_path(nx.induced_subgraph(node[2], set(node[2].nodes)-set(node[2].neighbors(node[0][-1]))), s, node[0][0])) and (node[0][-1]==t or nx.has_path(nx.induced_subgraph(node[2], set(node[2].nodes)-set(node[2].neighbors(node[0][0]))), node[0][-1], t)):
            else:
                if node[0][0] == s:
                    if not node[2].neighbors(node[0][-1]):
                        continue
                    for su in node[2].neighbors(node[0][-1]):
                        newnode = [node[0]+[su,], node[1]-1, node[2].copy()]
                        newnode[2].remove_nodes_from(set(node[2].neighbors(node[0][-1]))-{su,})
                        if newnode[0][-1] == t or ({newnode[0][-1],t} <= set(nx.induced_subgraph(newnode[2], set(newnode[2].nodes)-set(newnode[2].neighbors(newnode[0][0]))).nodes) and nx.has_path(nx.induced_subgraph(newnode[2], set(newnode[2].nodes)-set(newnode[2].neighbors(newnode[0][0]))), newnode[0][-1], t)):
                            stack.append(newnode.copy())
                elif node[0][-1] == t:
                    if not node[2].neighbors(node[0][0]):
                        continue
                    for pr in node[2].neighbors(node[0][0]):
                        newnode = [[pr,]+node[0], node[1]-1, node[2].copy()]
                        newnode[2].remove_nodes_from(set(node[2].neighbors(node[0][0]))-{pr,})
                        if newnode[0][0] == s or ({s, newnode[0][0]} <= set(nx.induced_subgraph(newnode[2], set(newnode[2].nodes)-set(newnode[2].neighbors(newnode[0][-1]))).nodes) and nx.has_path(nx.induced_subgraph(newnode[2], set(newnode[2].nodes)-set(newnode[2].neighbors(newnode[0][-1]))), s, newnode[0][0])):
                            stack.append(newnode.copy()) 
                elif node[0][0] == node[0][-1]:
                    if not node[2].neighbors(node[0][0]):
                        continue
                    for pr in node[2].neighbors(node[0][0]):
                        if not set(node[2].neighbors(node[0][-1]))-({pr,}|set(node[2].neighbors(pr))):
                            continue
                        for su in set(node[2].neighbors(node[0][-1]))-({pr,}|set(node[2].neighbors(pr))):
                            newnode = [[pr,]+node[0]+[su,], node[1]-1, node[2].copy()]
                            newnode[2].remove_nodes_from((set(node[2].neighbors(node[0][0]))|set(node[2].neighbors(node[0][-1])))-{pr,su})
                            if (newnode[0][0] == s or ({s, newnode[0][0]} <= set(nx.induced_subgraph(newnode[2], set(newnode[2].nodes)-set(newnode[2].neighbors(newnode[0][-1]))).nodes) and nx.has_path(nx.induced_subgraph(newnode[2], set(newnode[2].nodes)-set(newnode[2].neighbors(newnode[0][-1]))), s, newnode[0][0]))) and (newnode[0][-1] == t or ({newnode[0][-1],t} <= set(nx.induced_subgraph(newnode[2], set(newnode[2].nodes)-set(newnode[2].neighbors(newnode[0][0]))).nodes) and nx.has_path(nx.induced_subgraph(newnode[2], set(newnode[2].nodes)-set(newnode[2].neighbors(newnode[0][0]))), newnode[0][-1], t))):
                                stack.append(newnode.copy())
                else:
                    if not set(node[2].neighbors(node[0][0]))-set(node[2].neighbors(node[0][-1])):
                        continue
                    for pr in set(node[2].neighbors(node[0][0]))-set(node[2].neighbors(node[0][-1])):
                        if not set(node[2].neighbors(node[0][-1]))-({pr,}|set(node[2].neighbors(pr))|set(node[2].neighbors(node[0][0]))):
                            continue
                        for su in set(node[2].neighbors(node[0][-1]))-({pr,}|set(node[2].neighbors(pr))|set(node[2].neighbors(node[0][0]))):
                            newnode = [[pr,]+node[0]+[su,], node[1]-1, node[2].copy()]
                            newnode[2].remove_nodes_from((set(node[2].neighbors(node[0][0]))|set(node[2].neighbors(node[0][-1])))-{pr,su})
                            if (newnode[0][0] == s or ({s, newnode[0][0]} <= set(nx.induced_subgraph(newnode[2], set(newnode[2].nodes)-set(newnode[2].neighbors(newnode[0][-1]))).nodes) and nx.has_path(nx.induced_subgraph(newnode[2], set(newnode[2].nodes)-set(newnode[2].neighbors(newnode[0][-1]))), s, newnode[0][0]))) and (newnode[0][-1] == t or ({newnode[0][-1],t} <= set(nx.induced_subgraph(newnode[2], set(newnode[2].nodes)-set(newnode[2].neighbors(newnode[0][0]))).nodes) and nx.has_path(nx.induced_subgraph(newnode[2], set(newnode[2].nodes)-set(newnode[2].neighbors(newnode[0][0]))), newnode[0][-1], t))):
                                stack.append(newnode.copy())
        if not stack:
            res.add(v)
    return sorted(res)

def unrecog(g: nx.Graph, func=recog_dead):
    res = list()
    for s in g.nodes:
        for t in set(g.nodes)-{s,}:
            if nx.has_path(g, s, t):
                d = set(dead(g, s, t))
                rd = set(func(g, s, t))
                assert d >= rd
                if d != rd:
                    res.append((s,t,tuple(sorted(d-rd))))
    return res

from itertools import product
import numpy as np

def all_undirected_simple_connected(start=1):
    l = start
    while True:
        for x in product((0, 1), repeat=(l*(l-1))//2):
            A = np.zeros((l,l))
            k = 0
            for i in range(l):
                for j in range(i+1, l):
                    A[i,j] = x[k]
                    A[j,i] = x[k]
                    k += 1
            assert k == len(x)
            #print(A, end='\n\n')
            A = nx.from_numpy_array(A)
            if nx.is_connected(A):
                yield A
        print("finished",l)
        l += 1

def random_graph(n, p):
    while True:
        x = nx.gnp_random_graph(n,p)
        x.remove_nodes_from(set(x.nodes)-set(max(nx.connected_components(x), key=len)))
        yield x

def find_unrecog(gen=all_undirected_simple_connected, args=[], func=recog_dead):
    for g in gen(*args):
        u = unrecog(g, func=func)
        if u:
            yield [g, u]
