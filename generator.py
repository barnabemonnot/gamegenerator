from gurobipy import *
import numpy as np
import itertools
import random as rd
import gambit
import gambit.nash
from decimal import *

# Game A represented as array of payoff matrices (one per player)
# Payoff matrices are multidimensional arrays (as many dimensions as there are players)

def cartesian(arrays, out=None):
    # Returns cartesian product of arrays
    # Props to http://stackoverflow.com/questions/1208118/using-numpy-to-build-an-array-of-all-combinations-of-two-arrays/1235363#1235363
    arrays = [np.asarray(x) for x in arrays]
    dtype = arrays[0].dtype

    n = np.prod([x.size for x in arrays])
    if out is None:
        out = np.zeros([n, len(arrays)], dtype=dtype)

    m = n / arrays[0].size
    out[:,0] = np.repeat(arrays[0], m)
    if arrays[1:]:
        cartesian(arrays[1:], out=out[0:m,1:])
        for j in xrange(1, arrays[0].size):
            out[j*m:(j+1)*m,1:] = out[0:m,1:]
    return out

def parseGame(A):
    # Returns shape (number of players, number of strategies of P1, number of strategies of P2...)
    # number of players
    # pure moves as tuples of strategies
    shape = np.shape(A)
    num_players = shape[0]
    arrays = []
    for i in range(1, len(shape)):
        arrays.append(range(0, shape[i]))
    pure_moves = cartesian(arrays)
    return (shape, num_players, pure_moves)

def mixedNashEquilibria(A):
    # Returns all mixed NE of game A
    (shape, num_players, pure_moves) = parseGame(A)
    g = gambit.new_table(list(shape[1:]))
    for move in pure_moves:
        for player in range(0, num_players):
            g[tuple(move)][player] = Decimal(A[(player,)+tuple(move)])
    solver = gambit.nash.ExternalEnumMixedSolver()
    a = solver.solve(g)
    ne = []
    for eq in a:
        new_eq = []
        for player in range(0, num_players):
            eq_pl = []
            for strat in range(0, shape[player+1]):
                eq_pl.append(eq[g.players[player].strategies[strat]])
            new_eq.append(eq_pl)
        ne.append(new_eq)
    return ne
    
def pureNashEquilibria(A):
    # Needs more testing
    (shape, num_players, pure_moves) = parseGame(A)
    ne = []
    for move in pure_moves:
        is_ne = 1
        for player in range(1, num_players+1):
            payoff = A[player-1,]
            for strat in range(0, shape[player]):
                if strat != move[player-1]:
                    dev = np.copy(move)
                    dev[player-1] = strat
                    if payoff[tuple(dev)] < payoff[tuple(move)]:
                        is_ne = 0
        if is_ne == 1:
            new_ne = []
            for player in range(1, num_players+1):
                ne_str = []
                for strat in range(0, shape[player]):
                    if strat == move[player-1]:
                        ne_str.append(1)
                    else:
                        ne_str.append(0)
                new_ne.append(ne_str)
                        
            ne.append(new_ne)
    print ne

def getCorrelatedEquilibria(A1, A2):
    a = np.size(A1,0)
    b = np.size(A1,1)
    m = Model('correlated')
    profiles = tuplelist([(x, y) for x in range(1,a+1) for y in range(1,b+1)])
    cost_1 = { (x, y): A1[x-1,y-1] for x in range(1,a+1) for y in range(1,b+1) }
    cost_2 = { (x, y): A2[x-1,y-1] for x in range(1,a+1) for y in range(1,b+1) }
    
    printconstr = 0
    p = {}
    for profile in profiles:
        p[profile] = m.addVar(name=('p(%d,%d)' % profile), obj=(cost_1[profile]+cost_2[profile]))
    
    m.update()
    m.setParam('OutputFlag', False)
    m.setParam('FeasibilityTol', 1e-9)
    
    for i in range(1, a+1):
       for k in range(1, a+1):
           if i != k:
               m.addConstr(quicksum(p[profile]*cost_1[profile] for profile in profiles.select(i,'*')) <= quicksum(p[profiles.select(i,'*')[t]] * cost_1[profiles.select(k,'*')[t]] for t in range(0,b)), name='p1constr'+str(i)+'->'+str(k))
               if printconstr == 1:
                   print '--------------------------------------------------------------'
                   print 'i = %d, k = %d' % (i, k)
                   print quicksum(p[profile]*cost_1[profile] for profile in profiles.select(i,'*'))
                   print '<='
                   print quicksum(p[profiles.select(i,'*')[t]] * cost_1[profiles.select(k,'*')[t]] for t in range(0,b))
    for j in range(1, b+1):
       for l in range(1, b+1):
           if j != l:
               m.addConstr(quicksum(p[profile]*cost_2[profile] for profile in profiles.select('*',j)) <= quicksum(p[profiles.select('*',j)[t]] * cost_2[profiles.select('*',l)[t]] for t in range(0,a)), name='p2constr'+str(j)+'->'+str(l))
               if printconstr == 1:
                   print '--------------------------------------------------------------'
                   print 'j = %d, l = %d' % (j, l)
                   print quicksum(p[profile]*cost_2[profile] for profile in profiles.select('*',j))
                   print '<='
                   print quicksum(p[profiles.select('*',j)[t]] * cost_2[profiles.select('*',l)[t]] for t in range(0,a))
               
    m.addConstr(quicksum(p[profile] for profile in profiles) == 1, name='proba')
    
    m.optimize()

    resp = np.array([v.x for v in m.getVars()])
    resobj = m.objVal
    
    slack = m.getAttr(GRB.attr.Slack, m.getConstrs())
    if np.product(slack >= 0):
        return (resobj, resp)
    else:
        return (0, 0)

def getSocialCost(A, eq):
    (shape, num_players, pure_moves) = parseGame(A)
    sc = 0
    for move in pure_moves:
        mult = 1
        cost = 0
        for player in range(0, num_players):
            mult *= eq[player][move[player]]
            cost += A[(player,)+tuple(move)]
        sc = sc + mult * cost
    return sc

def twoPlayersNashEquilibria(A1, A2):
    ne = np.empty([0,4])
    if A1[0,0] <= A1[1,0] and A2[0,0] <= A2[0,1]:
        ne = np.vstack((ne, [[1,0,1,0]]))
    if A1[0,1] <= A1[1,1] and A2[0,1] <= A2[0,0]:
        ne = np.vstack((ne, [[1,0,0,1]]))
    if A1[1,0] <= A1[0,0] and A2[1,0] <= A2[1,1]:
        ne = np.vstack((ne, [[0,1,1,0]]))
    if A1[1,1] <= A1[0,1] and A2[1,1] <= A2[1,0]:
        ne = np.vstack((ne, [[0,1,0,1]]))
    x = (A2[1,1]-A2[1,0])/(A2[0,0]+A2[1,1]-A2[0,1]-A2[1,0])
    y = (A1[1,1]-A1[0,1])/(A1[0,0]+A1[1,1]-A1[0,1]-A1[1,0])
    if x <= 1 and x >= 0 and y <= 1 and y >= 0:
        ne = np.vstack((ne, [[x, 1-x, y, 1-y]]))
    return ne

def twoPlayersGetSocialCost(A1, A2, eq):
    cost = np.dot(np.dot(eq[[0,1]],A1),eq[[2,3]])+np.dot(np.dot(eq[[0,1]],A2),eq[[2,3]])
    return cost

def generateRandomGame(num_players, strategies):
    return np.random.random_sample((num_players,)+tuple(strategies))

def reversePayoff(A):
    max_payoff = np.amax(A)
    min_payoff = np.amin(A)
    (shape, num_players, pure_moves) = parseGame(A)
    new_A = np.copy(A)
    for move in pure_moves:
        for player in range(0, num_players):
            new_A[(player,)+tuple(move)] = -A[(player,)+tuple(move)]
    return new_A   

def generateRandomGame2x2():
    A1 = np.array([[rd.random(), rd.random()], [rd.random(), rd.random()]])
    A2 = np.array([[rd.random(), rd.random()], [rd.random(), rd.random()]])
    return (A1, A2)
