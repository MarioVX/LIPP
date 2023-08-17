"""
Created on Fri Nov 19 16:25:25 2021.

@author: Mario
"""
import json
import matplotlib.pyplot as plt
from time import process_time

lowerbound = 0

class InputGraph:
    
    def __init__(self, inputdict: dict):
        self.vertices = list(inputdict.keys())
        self.vertices.sort()
        self.vertices = tuple(self.vertices)
        self.goal = None
        self.succ = list(frozenset(self.vertices.index(x) for x in inputdict[self.vertices[i]]) for i in range(len(self.vertices)))
        self.pred = list(frozenset(i for i in range(len(self.vertices)) if j in self.succ[i]) for j in range(len(self.vertices)))
    
    def setgoal(self, goal):
        if goal not in self.vertices and goal is not None:
            raise KeyError("invalid goal: "+str(goal))
        if goal is None:
            self.goal = None
            return None
        self.goal = self.vertices.index(goal)
        return None


def load(filename: str) -> InputGraph:
    with open(filename, "r") as inputfile:
        G = json.load(inputfile)
    return InputGraph(G)


class SetTrie:
    
    def __init__(self):
        self.d = [dict(), False]
    
    def insert(self, bitstring):
        bs = tuple(int(x) for x in bitstring)
        if sum(bs) == 0:
            self.d[1] = True
            return None
        lastidx = len(bs) - 1
        while not bs[lastidx]:
            lastidx -= 1
        t = self.d[0]
        idx = 0
        while True:
            while not bs[idx]:
                idx += 1
            if idx == lastidx:
                if idx not in t:
                    t[idx] = [dict(), True]
                else:
                    t[idx][1] = True
                return None
            else:
                if idx not in t:
                    t[idx] = [dict(), False]
                t = t[idx][0]
                idx += 1
    
    def delete(self, bitstring):
        idxseq = list()
        for i in range(len(bitstring)):
            if int(bitstring[i]):
                idxseq.append(i)
        if not idxseq:
            self.d[1] = False
            return None
        nodes = [self.d,]
        isi = 0
        while isi < len(idxseq):
            nodes.append(nodes[-1][0][idxseq[isi]])
            isi += 1
        isi -= 1
        nodes[-1][1] = False
        while nodes and not nodes[-1][1] and not nodes[-1][0]:
            nodes.pop()
            if nodes:
                del nodes[-1][0][idxseq[isi]]
            isi -= 1
        return None
    
    def existsSuperset(self, bitstring) -> bool:
        bs = tuple(int(x) for x in bitstring)
        if sum(bs) == 0:
            if self.d[1] or self.d[0]:
                return True
        lastidx = len(bs) -1
        while not bs[lastidx]:
            lastidx -= 1
        firstidx = 0
        while not bs[firstidx]:
            firstidx += 1
        q = list()
        for i in range(firstidx+1):
            if i in self.d[0]:
                q.append([self.d[0], i])
        while q:
            t, idx = q.pop()
            if idx == lastidx:
                return True
            nextidx = idx + 1
            while not (bs[nextidx] or (nextidx == lastidx)):
                nextidx += 1
            for i in range(idx+1, nextidx+1):
                if i in t[idx][0]:
                    q.append([t[idx][0], i])
        return False
    
    def getAllSubsets(self, bitstring) -> list:
        # print("Own content:")
        # print(self.d)
        # print("Looking for subsets of:", bitstring)
        idxseq = list()
        for i in range(len(bitstring)):
            if int(bitstring[i]):
                idxseq.append(i)
        foundsubs = list()
        if self.d[1]:
            foundsubs.append(tuple())
        q = list()
        for i in range(len(idxseq)):
            if idxseq[i] in self.d[0]:
                q.append([self.d[0], (idxseq[i],)])
        # print("idxseq =", idxseq)
        # print("initial foundsubs:", foundsubs)
        # print("initial q:", q)
        while q:
            t, iseq = q.pop()
            # print("t:", t)
            # print("iseq:", iseq)
            if t[iseq[-1]][1]:
                foundsubs.append(iseq)
            for i in range(iseq[-1]+1, len(idxseq)):
                if idxseq[i] in t[iseq[-1]][0]:
                    q.append([t[iseq[-1]][0], iseq+(idxseq[i],)])
            # print("updated q:", q)
        formatted = list()
        for iseq in foundsubs:
            x = "0"*len(bitstring)
            for i in iseq:
                x = x[:i] + "1" + x[i+1:]
            formatted.append(x)
        return formatted


def prio2(bitstring, pathlength) -> int:
    space = 0
    for x in bitstring:
        if int(x):
            space += 1
    return space*len(bitstring) + pathlength

def prio3(bitstring, pathlength) -> int:
    space = 0
    for x in bitstring:
        if int(x):
            space += 1
    return space + pathlength

def transmit2(prev, oldstart: int, newstart: int, graph: InputGraph):
    goal = graph.goal
    if oldstart == goal:
        return None
    if newstart == goal:
        return "0"*goal + "1" + "0"*(len(prev)-1-goal)
    prevs = frozenset(i for i in range(len(prev)) if int(prev[i]))
    news = set()
    removed = {oldstart,} | ((prevs & graph.succ[oldstart]) - {newstart,})
    mustcheck = {newstart,}
    if goal is None:
        while mustcheck:
            v = mustcheck.pop()
            news.add(v)
            mustcheck.update((prevs & graph.succ[v]) - (news | removed))
    else:
        # directed case. forward & backward flood without bypassing
        unsure = set()
        # forward pass
        while mustcheck:
            v = mustcheck.pop()
            unsure.add(v)
            if v in graph.pred[goal]:
                if goal not in removed:
                    unsure.add(goal)
            else:
                mustcheck.update((prevs & graph.succ[v]) - (unsure | removed))
        if goal not in unsure:
            return None
        # backward pass
        mustcheck = {goal, }
        while mustcheck:
            v = mustcheck.pop()
            news.add(v)
            if v in graph.succ[newstart]:
                if newstart not in removed:
                    news.add(newstart)
            else:
                mustcheck.update((prevs & graph.pred[v]) - (news | removed))
        if newstart not in news:
            return None
        news &= unsure
    # leaf pruning
    if goal is not None:
        unsure = set()
        while news != unsure:
            unsure = news.copy()
            for v in unsure:
                if v != newstart and v != goal:
                    # now v must have a predeccsor-successor-pair that are neither
                    # identical nor neighbours, otherwise remove it from news.
                    approved = False
                    for p in (news & graph.pred[v]):
                        for s in (news & graph.succ[v]):
                            if p!=s and (s not in (unsure & graph.succ[p])) and (p not in (unsure & graph.pred[s])):
                                approved = True
                                break
                        if approved:
                            break
                    if not approved:
                        news.discard(v)
    newbs = ""
    for i in range(len(prev)):
        newbs += str(int(i in news))
    return newbs


class ParetoSet:
    priof = prio2
    transf = transmit2
    
    def __init__(self, start: int, graph: InputGraph):
        self.start = start
        self._graph = graph
        self.prio = 0
        self.inbox = dict((v, list()) for v in self._graph.pred[self.start])
        self.outbox = dict((v,list()) for v in self._graph.succ[self.start])
        self._unsent = dict() # prio -> set
        self._st = dict() # pathlength -> SetTrie
    
    def integrate(self):
        for sender in self.inbox:
            while self.inbox[sender]:
                bs, pl = self.inbox[sender].pop()
                # look for non-strict supersets of bs with length >= pl. proceed only if none are found
                dominated = False
                for l in self._st:
                    if l >= pl:
                        if self._st[l].existsSuperset(bs):
                            dominated = True
                            break
                if not dominated:
                    # add to SetTrie and unsent, delete all subsets with length <= pl.
                    pr = ParetoSet.priof(bs, pl)
                    self.prio = max(self.prio, pr)
                    if pr not in self._unsent:
                        self._unsent[pr] = set()
                    self._unsent[pr].add((bs, pl))
                    for l in self._st:
                        if l <= pl:
                            for sub in self._st[l].getAllSubsets(bs):
                                subpr = ParetoSet.priof(sub, l)
                                self._st[l].delete(sub)
                                if subpr in self._unsent:
                                    self._unsent[subpr].discard(sub)
                                    if not self._unsent[subpr]:
                                        del self._unsent[subpr]
                    if pl not in self._st:
                        self._st[pl] = SetTrie()
                    self._st[pl].insert(bs)
        return None
    
    def prepsend(self, priority: int):
        global lowerbound
        if priority > self.prio:
            return None
        elif priority < self.prio:
            raise ValueError("Miscoordination! This node still had a higher unsent priority="
                             + str(self.prio) + " than what was just called: "+str(priority))
        while self._unsent[priority]:
            bs, pl = self._unsent[priority].pop()
            for v in self.outbox:
                if int(bs[v]):
                    x = (ParetoSet.transf(bs, self.start, v, self._graph), pl+1)
                    if x[0] is not None:
                        lowerbound = max(lowerbound, x[1])
                        if x[1]+x[0].count('1') >= lowerbound:
                            self.outbox[v].append(x)
        del self._unsent[priority]
        if not self._unsent:
            self.prio = 0
        else:
            self.prio = max(self._unsent)
        return None


class OperationGraph:
    
    def __init__(self, igraph: InputGraph, start, goal):
        self.igraph = igraph
        self.igraph.setgoal(goal)
        self.start = start
        self.goal = goal
        self.psets = list(ParetoSet(i, igraph) for i in range(len(igraph.vertices)))
        if start is None:
            for i in range(len(self.psets)):
                next(iter(self.psets[i].inbox.values())).append(("1"*len(self.psets), 0))
                self.psets[i].integrate()
        else:
            si = self.igraph.vertices.index(start)
            next(iter(self.psets[si].inbox.values())).append(("1"*len(self.psets), 0))
            self.psets[si].integrate()
        self.prio = ParetoSet.priof("1"*len(self.psets), 0)
        self.maxlength = None
    
    def solve(self):
        while self.prio > 0:
            for i in range(len(self.psets)):
                self.psets[i].prepsend(self.prio)
                for j in self.igraph.succ[i]:
                    self.psets[j].inbox[i] = self.psets[i].outbox[j].copy()
                    self.psets[i].outbox[j].clear()
            for i in range(len(self.psets)):
                self.psets[i].integrate()
            self.prio = max(x.prio for x in self.psets)
            # print("priority:", self.prio//len(self.psets), self.prio%len(self.psets))
        l = max(max(x._st) for x in self.psets)
        print("maximum path length:", l)
        self.maxlength = l
        return l
    
    def extractpaths(self):
        if self.maxlength is None:
            raise RuntimeError("The graph must be solved before paths can be extracted.")
        finishedpaths = list()
        wippaths = list()
        # initialize
        if self.goal is None:
            for i in range(len(self.psets)):
                if max(self.psets[i]._st) == self.maxlength:
                    wippaths.append([i,])
        else:
            assert self.maxlength == max(self.psets[self.igraph.vertices.index(self.goal)]._st)
            wippaths.append([self.igraph.vertices.index(self.goal),])
        # extend
        while wippaths:
            path = wippaths.pop()
            l = self.maxlength - len(path)
            if l == -1:
                path.reverse()
                finishedpaths.append(tuple(path))
            else:
                exclude = set(path)
                for x in path[:-1]:
                    exclude.update(self.igraph.pred[x])
                pathbs = ""
                for i in range(len(self.psets)):
                    pathbs += str(int(i in path))
                for v in self.igraph.pred[path[-1]] - exclude:
                    if l in self.psets[v]._st:
                        if self.psets[v]._st[l].existsSuperset(pathbs):
                            wippaths.append(path + [v,])
        return finishedpaths

def showgridpath(igraph: InputGraph, path):
    p = list(igraph.vertices[i] for i in path)
    p = list(eval(x) for x in p)
    p = list(zip(*p))
    width, height = eval(igraph.vertices[-1])
    fig, ax = plt.subplots()
    ax.plot(p[0], p[1], 'ro-')
    ax.set_xticks(range(0, width+1))
    ax.set_yticks(range(0, height+1))
    ax.grid(lw=1.5)
    plt.show()
    return None

t = process_time()
ig = load("prison.json")
# # print("|V| =", len(ig.vertices))
og = OperationGraph(ig, None, None)
og.solve()
temp = og.extractpaths()
t = process_time() - t
print (t, "seconds")
# for x in temp:
#     showgridpath(ig, x)
print("Number of tied optimal paths:",len(temp))
# print("End vertices:", set(ig.vertices[x[-1]] for x in temp))