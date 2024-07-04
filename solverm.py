import networkx as nx
from heapq import heappush, heappop

def int_from_set(inputset : set[int]) -> int:
    return sum(2**x for x in inputset)

def set_from_int(inputint : int) -> set[int]:
    return set(i for i in range(len(bin(inputint)[:1:-1])) if bin(inputint)[:1:-1][i] == '1')

def solve(inputgraph : nx.Graph, settings = dict(), inputstart = None):
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
    def mature(past : int, present : int, immature_future : set[int]) -> None:
        nonlocal ig, queue
        wg = nx.induced_subgraph(ig, nx.node_connected_component(nx.induced_subgraph(ig, immature_future), present)).copy()
        # point
        pass
        # B&B trim
        pass
        # upper bound
        upper_bound : float
        upper_bound = past + len(wg) - 1
        pass
        heappush(queue, (-upper_bound, (past, present, int_from_set(wg.nodes))))
        return None
    for s in start:
        mature(0, s, ig.nodes)
    opt_length = -1
    while queue:
        if -queue[0][0] < opt_length:
            break
        past,present,ifuture = heappop(queue)[1]
        future = set_from_int(ifuture)
        if len(future) == 1:
            # terminal state
            assert opt_length == -1 or past == opt_length
            if opt_length == -1:
                opt_length = past
        pass

solve(nx.read_graphml("graphs/lip_crn/karate.graphml"), settings={"iso":9e9})
