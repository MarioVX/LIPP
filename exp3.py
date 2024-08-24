import solverm
import csv
import os
from datetime import datetime

settings = {"iso":0, "sup":0, "ub":"ds", "dec":0, "bbt":0, "bbt_with_ub":"ds", "bbt_without_ub":"ds"}
graph_dir_name = 'graphs/lip_crn'

with open('experiment_' + datetime.now().strftime('%H%M%S') + '.csv', 'w', newline='') as resultfile:
    writer = csv.writer(resultfile)
    writer.writerow(['name','success', 'prio/length', 'searchnodes', 'time'])
    for graphfile in os.listdir(graph_dir_name):
        g = nx.read_graphml(graph_dir_name + "/" + str(graphfile))
        res = solverm.solve(g, settings)
        line = [str(graphfile),]
        if not res[0]:
            line += [0, res[2], res[1], 1200.0]
        else:
            line += [1, len(res[0][0]), res[1], res[2]]
        writer.writerow(line)
