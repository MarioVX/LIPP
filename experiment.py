import csv
import os
import subprocess
from datetime import datetime

graph_dir_name = 'graphs/lip_bal/

with open('experiment_'+datetime.now().strftime('%Y%m%d%H%M%S')+'.csv', 'w', newline='') as resultfile:
    writer = csv.writer(resultfile)
    
