"""
LIPP Kit
Solver

The program responsible for actually solving one graph with one parameter setting.
It's supposed to be invoked by the controller with appropriate command line arguments.
Upon successfull execution it should print results and metrics to stdout,
which can be caught by the controller who runs this as a subprocess.

@author: Mario
"""

class InputGraph:
    
    def __init__(self, inputdict: dict):
        self.vertices = list(inputdict.keys())
        self.vertices.sort()
        self.vertices = tuple(self.vertices)
        self.goal = None
        self.succ = list(frozenset(self.vertices.index(x) for x in inputdict[self.vertices[i]]) for i in range(len(self.vertices)))
        self.pred = list(frozenset(i for i in range(len(self.vertices)) if j in self.succ[i]) for j in range(len(self.vertices)))
    
    def setgoal(self, goal) -> None:
        assert goal is None or goal in self.vertices
        if goal is None:
            self.goal = None
        else:
            self.goal = self.vertices.index(goal)
        return None


if __name__ == '__main__':
    from interface import Instructions
    
    instructions = Instructions(True)
    
    print("Solver is called on graph:", instructions.graphfile,
          "with the following list of settings:", instructions.slist, sep='\n')
    del instructions.slist
    
    # --- read in the given graph ---
    def dict_from_graphml() -> dict:
        with open('graphs/' + instructions.graphfile, "r") as f:
            for _ in range(3):
                line = f.readline()
            line = line.strip()
            line = line.split()
            undirected = (line[-1][13:-2] == "undirected")
            V = dict()
            line = f.readline()
            line = line.strip()
            while line:
                line = line.split()
                if line[0] == "<node":
                    V[line[1][4:-1]] = list()
                elif line[0] == "<edge":
                    s = line[2][8:-1]
                    t = line[3][8:-1]
                    V[s].append(t)
                    if undirected:
                        V[t].append(s)
                line = f.readline()
                line = line.strip()
        return V
    
    inputgraph = InputGraph(dict_from_graphml())
    if instructions.has_goal:
        inputgraph.setgoal(instructions.goal)
    
    print(inputgraph.succ)