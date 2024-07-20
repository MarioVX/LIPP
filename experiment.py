import csv
import os
import solverm
from datetime import datetime
import networkx as nx

graph_dir_name = 'graphs/lip_bal/ba_100_3'
exp_name = 'ba_100_3'

settings = list()
for iso in (0,1):
    for sup in (0,1):
        for dec in (0,1):
            for bbt in (0,1):
                for ub in ("ds", 5, 10, 15):
                    if bbt <= dec:
                        s = {"iso":iso, "sup":sup, "dec":dec, "bbt":bbt, "ub":ub}
                        if bbt:
                            for bwu in ("ds",):
                                for bwou in ("ds",):
                                    s2 = s.copy()
                                    s2["bbt_with_ub"] = bwu
                                    s2["bbt_without_ub"] = bwou
                                    settings.append(s2)
                        else:
                            s["bbt_with_ub"] = 0
                            s["bbt_without_ub"] = 0
                            settings.append(s)

header = sorted(settings[0].keys()) + ["success", "prio/length", "nodes", "time"]

with open('experiment_'+exp_name+'_'+datetime.now().strftime('%Y%m%d%H%M%S')+'.csv', 'w', newline='') as resultfile:
    writer = csv.writer(resultfile)
    writer.writerow(header)
    for setting in settings:
        print(setting)
        for graphfile in os.listdir(graph_dir_name):
            print(graphfile)
            g = nx.read_graphml(graph_dir_name+"/"+str(graphfile))
            res = solverm.solve(g, setting)
            line = list(setting[k] for k in header[:7])
            if not res[0]:
                line += [0, res[2], res[1], 1200.0]
            else:
                line += [1, len(res[0][0]), res[1], res[2]]
            writer.writerow(line)
            input()
