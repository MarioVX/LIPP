
def adjacent_settings(inpset : dict[str,int|str]) -> list[dict[str,int|str]]:
    out : list[dict[str,int|str]] = list()
    # mirror iso
    s = inpset.copy()
    s["iso"] = 1 - s["iso"]
    out.append(s)
    # mirror sup
    s = inpset.copy()
    s["sup"] = 1 - s["sup"]
    out.append(s)
    # relax upper bound
    if inpset["ub"] != "trivial":
        s = inpset.copy()
        match s["ub"]:
            case "ds":
                s["ub"] = "trivial"
            case 3:
                s["ub"] = "ds"
            case _:
                assert type(s["ub"]) == int and s["ub"] > 3
                s["ub"] -= 1
        out.append(s)
    # tighten upper bound
    s = inpset.copy()
    match s["ub"]:
        case "trivial":
            s["ub"] = "ds"
        case "ds":
            s["ub"] = 3
        case _:
            assert type(s["ub"]) == int and s["ub"] > 2
            s["ub"] += 1
    out.append(s)
    if inpset["dec"]:
        # deactivate decomposition
        s = inpset.copy()
        s["dec"] = 0
        s["bbt"] = 0
        out.append(s)
        # toggle decomposition
        s = inpset.copy()
        s["dec"] = 3 - s["dec"]
        out.append(s)
        if inpset["bbt"]:
            # deactivate bbt
            s = inpset.copy()
            s["bbt"] = 0
            out.append(s)
            # relax ub with
            if inpset["bbt_with_ub"] != "ds":
                s = inpset.copy()
                match s["bbt_with_ub"]:
                    case 3:
                        s["bbt_with_ub"] = "ds"
                    case _:
                        assert type(s["bbt_with_ub"]) == int and s["bbt_with_ub"] > 3
                        s["bbt_with_ub"] -= 1
                out.append(s)
            # tighten ub with
            s = inpset.copy()
            match s["bbt_with_ub"]:
                case "ds":
                    s["bbt_with_ub"] = 3
                case _:
                    assert type(s["bbt_with_ub"]) == int
                    s["bbt_with_ub"] += 1
            out.append(s)
            # relax ub without
            if inpset["bbt_without_ub"] != "ds":
                s = inpset.copy()
                match s["bbt_without_ub"]:
                    case 3:
                        s["bbt_without_ub"] = "ds"
                    case _:
                        assert type(s["bbt_without_ub"]) == int and s["bbt_without_ub"] > 3
                        s["bbt_without_ub"] -= 1
                out.append(s)
            # tighten ub without
            s = inpset.copy()
            match s["bbt_without_ub"]:
                case "ds":
                    s["bbt_without_ub"] = 3
                case _:
                    assert type(s["bbt_without_ub"]) == int
                    s["bbt_without_ub"] += 1
            out.append(s)
        else:
            # no bbtrim so far, try cheapest one
            s = inpset.copy()
            s["bbt"] = 1
            s["bbt_with_ub"] = "ds"
            s["bbt_without_ub"] = "ds"
            out.append(s)
    else:
        # no decomposition so far, try the cheapest one
        s = inpset.copy()
        s["dec"] = 1
        s["bbt"] = 0
        out.append(s)
    return out

import networkx as nx
def gen_hypercubes():
    n = 3
    while True:
        yield [nx.hypercube_graph(n),]
        n += 1

def gen_square_grids():
    n = 3
    while True:
        yield [nx.grid_2d_graph(n,n),]
        n += 1

def gen_random(exp_deg : int, sample_size : int):
    n = 10
    while True:
        yield list(nx.gnp_random_graph(n, exp_deg/(n-1)) for _ in range(sample_size))
        n += 10

def gen_barabasi(m : int, sample_size : int):
    n = 10
    while True:
        yield list(nx.barabasi_albert_graph(n, m) for _ in range(sample_size))
        n += 10


def tup_from_dict(s : dict) -> tuple:
    return tuple(s[k] for k in sorted(s.keys()))

def dict_from_tup(t : tuple) -> dict:
    h = sorted(["iso", "sup", "ub", "dec", "bbt", "bbt_with_ub", "bbt_without_ub"])
    return dict((h[i], t[i]) for i in range(len(h)))

import solverm
import csv
import os
from datetime import datetime

def experiment(name, generator, gen_args, init_settings={"iso":0, "sup":0, "ub":"ds", "dec":0, "bbt":0, "bbt_with_ub":"ds", "bbt_without_ub":"ds"}):
    print(name)
    header = sorted(init_settings.keys()) + ["V", "success", "prio/length", "searchnodes", "time"]
    with open('experiment_'+name+'_'+datetime.now().strftime('%H%M%S')+'.csv', 'w', newline='') as resultfile:
        writer = csv.writer(resultfile)
        writer.writerow(header)
        pettings = init_settings
        for graphpack in generator(*gen_args):
            print(len(graphpack[0]), pettings)
            to_test = [pettings,]
            scores = dict()
            bestscore = None
            all_timeouts = True
            while to_test:
                settings = to_test.pop()
                tettings = tup_from_dict(settings)
                if tettings in scores:
                    continue
                scores[tettings] = 0
                for graph in graphpack:
                    res = solverm.solve(graph, settings)
                    line = list(settings[k] for k in header[:7])
                    line.append(len(graph))
                    if not res[0]:
                        line += [0, res[2], res[1], 1200.0]
                    else:
                        line += [1, len(res[0][0]), res[1], res[2]]
                        all_timeouts = False
                    writer.writerow(line)
                    scores[tettings] +=  line[-1] + line[-3]
                if bestscore is None or scores[tettings] <= bestscore:
                    bestscore = scores[tettings]
                    to_test = adjacent_settings(settings)
            if all_timeouts:
                break
            pettings = dict_from_tup(min(scores.keys(), key= lambda x:scores[x]))
    return None

if __name__ == '__main__':
    experiment("squares", gen_square_grids, [])
    experiment("hypercubes", gen_hypercubes, [])
    experiment("ER_D=3", gen_random, [3, 10])
    experiment("ER_D=4", gen_random, [4, 10])
    experiment("ER_D=5", gen_random, [5, 10])
    experiment("BA_m=3", gen_barabasi, [3, 10])
    print("experiments successfully completed!")
