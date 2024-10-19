# Longest Induced Path Problem Solver

The code submission accompanying by Master's thesis, "Exploration of Search-based Methods for the Longest Induced Path Problem".

The Python version requires Nauty and Pynauty installed to employ isomorphism reduction. It also needs networkx under any configuration.
The C++ version requires the Boost Graph Library to be installed and linked to compile.

Python version usage:
1. Construct or parse the graph in which the LIPs are to be sought using the networkx module, as a networkx.Graph object.
2. Import the solve function from solverm.py
3. Call the solve function with the graph object and a parameter settings dictionary, with the following items:
   - "iso": isomorphism checking, 1 or 0
   - "sup": superset inclusion checking, 1 or 0
   - "dec":
     - 0: no decomposition
     - 1: compute condensation and contraction, node weights from degree sequence method
     - 2: also computer line graph, edge weights from linear programs
   - "ub": upper bound method:
     - "trivial": |F|
     - "ds": degree sequence method
     - "fLP": flat LP
     - "pLP": pointed LP
     - int: use degree sequence above this number of remaining vertices, LP below
   - "bbt": B&B trimming, 1 or 0. 1 requires dec>0 and the following params:
     - "bbt_with_ub": upper bound method when enforcing the selection of a region. "ds", "LP" or int
     - "bbt_without_ub": upper bound method when prohibiting the selection of a region. "ds", "LP" or int
4. (optional) kwargs: timeout in seconds or None, inputstart to restrict eligible path starts (a specific vertex, a set, or None for any).
5. output:
   - if successful within time limit: (sorted list of longest induced paths, number of generated search nodes, time taken to compute solution)
   - if unsuccessful: (empty list, number of generated search nodes, upper bound at timeout)

C++ version usage:
1. Save the graph in which the LIPs are to be sought as a .graphml file and denote its file path.
2. Edit the CLIPP.cpp, setting the denoted file path in the first line of the main function.
3. Compile CLIPP with the necessary Boost libraries and link with BGL for parsing the .graphml file.
4. Execute the program.
5. Standard output: longest induced path length, time taken to compute solution
