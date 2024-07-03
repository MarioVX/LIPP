import pynauty
import networkx as nx

def pynaut_graph(input_graph : nx.Graph|nx.DiGraph) -> pynauty.Graph:
    a : dict[int,list[int]]
    if sorted(input_graph.nodes) != list(range(len(input_graph))):
        g = nx.relabel_nodes(input_graph, dict((sorted(input_graph.nodes)[i], i) for i in range(len(input_graph))))
        a = nx.to_dict_of_lists(g)
    else:
        a = nx.to_dict_of_lists(input_graph)
    return pynauty.Graph(len(input_graph), directed=(type(input_graph)==nx.DiGraph), adjacency_dict=a)


def canon_color(graph:pynauty.Graph) -> tuple[int]:
    vc = graph.vertex_coloring
    cl = pynauty.canon_label(graph)
    out = [0,] * len(cl)
    for c in range(len(vc)):
        for v in vc[c]:
            out[cl.index(v)] = c
    return tuple(out)

def actual_certificate(graph:pynauty.Graph):
    return pynauty.certificate(graph), canon_color(graph)
