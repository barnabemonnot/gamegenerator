from gurobipy import *
import numpy as np
import itertools
import random as rd

def twoPlayersNashEquilibrium(A1, A2):
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

def getCorrelatedEquilibria(A1, A2):
    a = np.size(A1,0)
    b = np.size(A1,1)
    m = Model('correlated')
    profiles = tuplelist([(x, y) for x in range(1,a+1) for y in range(1,b+1)])
    cost_1 = { (x, y): A1[x-1,y-1] for x in range(1,a+1) for y in range(1,b+1) }
    cost_2 = { (x, y): A2[x-1,y-1] for x in range(1,a+1) for y in range(1,b+1) }
    
    printconstr = False
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
               if printconstr:
                   print '--------------------------------------------------------------'
                   print 'i = %d, k = %d' % (i, k)
                   print quicksum(p[profile]*cost_1[profile] for profile in profiles.select(i,'*'))
                   print '<='
                   print quicksum(p[profiles.select(i,'*')[t]] * cost_1[profiles.select(k,'*')[t]] for t in range(0,b))
    for j in range(1, b+1):
       for l in range(1, b+1):
           if j != l:
               m.addConstr(quicksum(p[profile]*cost_2[profile] for profile in profiles.select('*',j)) <= quicksum(p[profiles.select('*',j)[t]] * cost_2[profiles.select('*',l)[t]] for t in range(0,a)), name='p2constr'+str(j)+'->'+str(l))
               if printconstr:
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

def twoPlayersGetSocialCost(A1, A2, eq):
    cost = np.dot(np.dot(eq[[0,1]],A1),eq[[2,3]])+np.dot(np.dot(eq[[0,1]],A2),eq[[2,3]])
    return cost

def generateRandomGame():
    A1 = np.array([[rd.random(), rd.random()], [rd.random(), rd.random()]])
    A2 = np.array([[rd.random(), rd.random()], [rd.random(), rd.random()]])
    return (A1, A2)
