import networkx as nx
from heapq import heappush, heappop
from time import monotonic

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
        q = list()
        b = -1
        for n in range(a[-1]+1):
            if n in self.d[0]:
                q.append([self.d[0][n], n+1, len(a)-1-int(n==a[-1])])
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
    def mature(past : int, present : int, immature_future : set[int], parent : tuple[int,int]) -> None:
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
        heappush(queue, (-upper_bound, (past, present, int_from_set(wg.nodes), parent)))
        return None
    solvetime = monotonic()
    searchnodes = 0
    for s in start:
        mature(0, s, ig.nodes, None)
    opt_length = -1
    domcomp : list[SetTrie] = list()
    idstruct : list[dict[int,list[int,list[tuple[int,int]]]]] = list()
    opt_seeds = set()
    for i in range(numinputnodes):
        domcomp.append(SetTrie())
        idstruct.append(dict())
    while queue:
        if -queue[0][0] < opt_length:
            break
        past,present,ifuture,parent = heappop(queue)[1]
        future = set_from_int(ifuture)
        status = domcomp[present].beatingSuperset(future, past)
        if status == 2:
            continue
        if status == 1 and ifuture in idstruct[present]:
            assert idstruct[present][ifuture][0] == past
            assert parent not in idstruct[present][ifuture][1]
            idstruct[present][ifuture][1].append(parent)
            continue
        assert (status == 1 and ifuture not in idstruct[present]) or status == 0
        idstruct[present][ifuture] = [past, [parent,]]
        domcomp[present].insert(future, past)
        if len(future) == 1:
            # terminal state
            assert opt_length == -1 or past == opt_length
            if opt_length == -1:
                opt_length = past
            opt_seeds.add(present)
            continue
        # expand!
        searchnodes += 1
        neighbors = set(ig.neighbors(present)) & future
        for v in neighbors:
            mature(past+1, v, ({v,}|(future-neighbors))-{present,}, (present, ifuture))
    solvetime = monotonic() - solvetime
    # retrieve all solutions
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
                assert len(path) == opt_length + 1
                path.reverse()
                out.append(tuple(path))
    out.sort()
    return out, searchnodes, solvetime

res = solve(nx.read_graphml("graphs/lip_crn/karate.graphml"), settings={"iso":9e9})
karate_paths = [(4, 23, 12, 21, 25, 18, 16, 22, 19), (8, 29, 0, 1, 24, 26, 16, 20, 17), (8, 29, 0, 1, 24, 27, 20, 17, 18), (8, 29, 0, 25, 17, 20, 16, 22, 19), (8, 30, 0, 1, 24, 26, 16, 20, 17), (8, 30, 0, 1, 24, 27, 20, 17, 18), (8, 30, 0, 25, 17, 20, 16, 22, 19), (9, 1, 12, 21, 25, 18, 16, 22, 19), (11, 1, 12, 21, 25, 18, 16, 22, 19), (14, 1, 12, 21, 25, 18, 16, 22, 19), (17, 18, 16, 26, 24, 1, 0, 30, 8), (17, 18, 16, 27, 21, 12, 0, 30, 8), (17, 18, 16, 27, 24, 1, 0, 30, 8), (17, 18, 16, 27, 33, 12, 0, 30, 8), (17, 20, 16, 26, 24, 1, 0, 30, 8), (18, 17, 20, 27, 24, 1, 0, 30, 8), (19, 22, 16, 18, 25, 21, 12, 1, 9), (19, 22, 16, 18, 25, 21, 12, 1, 11), (19, 22, 16, 18, 25, 21, 12, 1, 14), (19, 22, 16, 18, 25, 21, 12, 23, 4), (19, 22, 16, 18, 25, 21, 12, 32, 24), (19, 22, 16, 20, 17, 25, 0, 30, 8), (19, 22, 16, 20, 17, 25, 0, 32, 24), (24, 1, 0, 25, 17, 20, 16, 22, 19), (24, 1, 12, 21, 25, 18, 16, 22, 19), (24, 32, 0, 25, 17, 20, 16, 22, 19), (24, 32, 12, 21, 25, 18, 16, 22, 19)]
assert res[0] == karate_paths
print(res[1], res[2])
