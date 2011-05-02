#!/usr/bin/python2

from pulp import *

graph = [
(( 0,  8), 150),
(( 0,  11), 150),
(( 1,  8), 150),
(( 1,  9), 150),
(( 2,  3), 150),
(( 2,  9), 150),
(( 2,  15), 30),
(( 3,  16), 100),
(( 4,  16), 100),
(( 4,  18), 150),
(( 5,  17), 30),
(( 6,   7), 150),
(( 7,  10), 150),
(( 7,  25), 30),
(( 7,  28), 30),
(( 8,  12), 100),
(( 9,  14), 100),
((10,  11), 100),
((10,  21), 30),
((11,  12), 30),
((12,  13), 100),
((12,  21), 100),
((13,  20), 30),
((13,  14), 100),
((14,  15), 100),
((14,  19), 30),
((15,  16), 100),
((15,  19), 30),
((16,  19), 100),
((17,  18), 30),
((17,  24), 30),
((18,  19), 30),
((18,  23), 150),
((19,  20), 30),
((19,  22), 30),
((20,  21), 30),
((20,  27), 30),
((21,  26), 30),
((22,  23), 30),
((22,  27), 30),
((22,  28), 30),
((23,  24), 30),
((23,  29), 30),
((24,  29), 30),
((25,  26), 30),
((26,  27), 30),
((27,  28), 30),
((28,  29), 30),
]

vertexconstraints = [
    (0, 150),
    (1, 150),
    (2, 150),
    (3, 150),
    (4, 150),
    (6, 150),
    (7, 150),
    (8, 150),
    (9, 150),
    (18, 150),
    (10, 100),
    (11, 100),
    (12, 100),
    (13, 100),
    (14, 100),
    (15, 100),
    (16, 100),
    (17, 30),
    (19, 30),
    (20, 30),
    (21, 30),
    (22, 30),
    (23, 30),
    (24, 30),
    (25, 30),
    (5, 30),
]

sources = range(0, 7)
sinks = range(26, 30)

rest = []

ingoing = {}
outgoing = {}

for ((u, v), c) in graph:
    if u not in sources and u not in sinks:
        rest.append(u)
    if v not in sources and v not in sinks:
        rest.append(v)

    outvar = LpVariable("node_"+str(v)+"_"+str(u), 0, c, LpInteger)
    invar = LpVariable("node_"+str(u)+"_"+str(v), 0, c, LpInteger)

    if u not in outgoing:
        outgoing[u] = []
        ingoing[u] = []
    if v not in outgoing:
        outgoing[v] = []
        ingoing[v] = []

    ingoing[u].append(invar)
    outgoing[u].append(outvar)
    ingoing[v].append(outvar)
    outgoing[v].append(invar)

rest = set(rest)

prob = LpProblem("Copenhagen", LpMaximize)

astuff = []
astuff2 = []
for s in sources:
    for v in outgoing[s]:
        astuff.append(v)
    for v in ingoing[s]:
        astuff2.append(v)
prob += lpSum(astuff2) - lpSum(astuff), "source_stuff"

for u in rest:
    prob += lpSum(ingoing[u]) - lpSum(outgoing[u]) == 0, "node " + str(u)

for u, c in vertexconstraints:
    prob += lpSum(ingoing[u]) <= c, "v_" + str(u)

prob.writeLP("part2.lp")

prob.solve()
print "Status: ", LpStatus[prob.status]

for v in prob.variables():
    print v.name, "=", v.varValue

print "Total Cost = ", value(prob.objective)
