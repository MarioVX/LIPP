"""
LIPP Kit
Controller

The main program through which to use the toolkit. Handles user interaction.
Experiments can be specified through sets of parameter settings and graphs,
results will be collected and saved.

@author: Mario
"""

import subprocess
from interface import Instructions

graph: str = 'lip_crn/494bus.graphml'
instructions = Instructions(False)

result = subprocess.run(["python", "solver.py", graph] +
                        instructions.slist, capture_output=True, text=True)
print(result.stdout)
print(result.stderr)
