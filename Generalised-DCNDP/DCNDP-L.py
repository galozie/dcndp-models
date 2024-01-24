import networkx as nx
#import matplotlib.pyplot as plt
from gurobipy.gurobipy import Model, GRB, LinExpr
#import collections
import dist_connectivity as dc

#-------Lazy cut callback to separate path constraints-------------------
def mycut(model, where):
    global savebnd, cutcount, bndcheck 
    cost = {} #dictionary for solution of x variables
    connect={} # dictionary for solution of y variables
    
    #-----separates integer solution -----
    def depthKsp1(graph, cost, connect, L, cutlimit):# separates integer solution 
        roots=[n for (n,attr) in cost.items() if attr < 1-1e-5]
        input_graph=graph.subgraph(roots) 
        for rt in roots:
            cutcount=0   
            
            length, path = nx.single_source_dijkstra(input_graph, rt, 
                                                     cutoff=L, weight='weight')
            for v,distance in length.items():
                if rt!=v:
                     i=min([rt,v])
                     j=max([rt,v])   
                     if connect[(i, j, distance)] < 1: 
                         model.cbLazy(sum(model._x_delete[n] for n in path[v])\
                         + model._u_connect[i,j,distance] >= 1)
                         cutcount+=1
                         if cutcount ==cutlimit:
                             break
                
    #--------separates fractional solution                 
    def depthKsp2(graph, cost, connect, L, cutlimit): 
        roots=[n for (n,attr) in cost.items() if attr < 1-1e-5]
        input_graph=graph.subgraph(roots) 
        for rt in roots:
            cutcount=0
            length, path = nx.single_source_dijkstra(input_graph, rt,
                                                     cutoff=L, weight='weight')
            for v,distance in length.items():
                if rt!=v:
                     i=min([rt,v])
                     j=max([rt,v]) 
                     if sum(cost[n] for n in path[v]) + \
                     connect[(i, j, distance)] < 1:
                         model.cbLazy(sum(model._x_delete[n] for n in path[v])\
                         +  model._u_connect[i,j,distance] >= 1)
                         cutcount+=1
                         if cutcount ==cutlimit:
                             break
                         
    #----if integer solution    
    if where == GRB.Callback.MIPSOL:
         savebnd=model.cbGet(GRB.callback.MIPSOL_OBJBND)
         nodecnt = model.cbGet(GRB.Callback.MIPSOL_NODCNT)         
         for j in G.nodes():
            cost[j]=abs(model.cbGetSolution(model._x_delete[j]))
            for i in G.nodes():
                if i<j:
                    for l in range(1,L+1):
                        connect[(i,j,l)] = \
                        abs(model.cbGetSolution(model._u_connect[i,j,l]))
         depthKsp1(G, cost, connect, L, GRB.INFINITY) 

    #-----if fractional solution          
    elif where ==GRB.Callback.MIPNODE: 
         if model.cbGet(GRB.Callback.MIPNODE_STATUS) == GRB.Status.OPTIMAL:
            currentbnd=model.cbGet(GRB.callback.MIPNODE_OBJBND)
            nodecnt=int(model.cbGet(GRB.callback.MIPNODE_NODCNT))
            if savebnd==currentbnd:
                bndcheck+=1                  
                if bndcheck >=5 or nodecnt>0:
                   bndcheck=0
                else:
                    for j in G.nodes():
                        cost[j]=abs(model.cbGetNodeRel(model._x_delete[j]))
                        for i in G.nodes():
                            if i<j:
                                for l in range(1,L+1):
                                    connect[(i,j,l)] = \
                                    abs(model.cbGetNodeRel(model._u_connect[i,j,l]))
                    depthKsp2(G, cost, connect, L, 300)
                        
            else:
                savebnd=currentbnd
                for j in G.nodes():
                        cost[j]=abs(model.cbGetNodeRel(model._x_delete[j]))
                        for i in G.nodes():
                            if i<j:
                                for l in range(1,L+1):
                                    connect[(i,j,l)] = \
                                    abs(model.cbGetNodeRel(model._u_connect[i,j,l]))
                depthKsp2(G, cost, connect, L, 300)
        

#-----Minimize distance-based connectivity metric specified by f((l)---------
def main(H,k,C):
    model = Model('Minimize distance-based connectivity objective')
    
   #Decision variables
    x_delete = {} #x varaibles defined for every node i
    u_connect = {} #y varaibles defined for every node pair (i,j)
    
    for j in H.nodes():
        if H.degree[j]==1: #fix value of x-variable of degree 1 nodes to zero
            x_delete[j] = model.addVar(lb=0.0, ub=0.0, 
                    vtype=GRB.BINARY, name ="x[%d]" %(j))
        else: #define domain of x-variable for non-degree 1 nodes
            x_delete[j] = model.addVar(lb=0.0, ub=1.0, 
                    vtype=GRB.BINARY, name ="x[%d]" %(j))
            
        #----define domain of x-variable for non-degree 1 nodes    
        for i in H.nodes():
            if i <j:
                for l in range(1,L+1):
                    u_connect[i, j, l] =model.addVar(lb=0.0, ub=1.0,
                             vtype=GRB.CONTINUOUS, name="u[%d,%d,%d]" %(i,j,l))  

    #----objective function
    obj = LinExpr(0)
    for j in H.nodes():
        for i in H.nodes():
            if i<j:
                obj.add(f[0] * u_connect[i, j, 1])
                for l in range(1,L):
                    obj.add(f[l] * (u_connect[i, j, l + 1] - u_connect[i, j, l]))

    #-----constraint on number of critical nodes, that is constraint (15)
    model.addConstr(sum((x_delete[j]) for j in H.nodes())<=C) 


    #----constraints on adjacent node pairs (i,j) in E 
    for (i,j) in H.edges():
         weight_ij=H.edges[i, j]['weight']
         if i<j:
            model.addConstr(u_connect[i,j,weight_ij] + x_delete[i]+ \
                            x_delete[j] >= 1) #constraint (10) for edges
            weights_other=[H.edges[i,t]['weight'] for t in set(H[i])-set([j])]
            weights_other.append(weight_ij-1)
            weight=min(weights_other)
            for l in range(1, weight+1):
                model.addConstr(u_connect[i,j,l] == 0) #constraint (13)
            for l in range(weight_ij+1, L+1):   
                model.addConstr(u_connect[i,j,l] == u_connect[i,j,weight_ij]) #constraint (12) for edges
         else: #that is j<i
             model.addConstr(u_connect[j,i,weight_ij] + \
                             x_delete[j]+ x_delete[i] >= 1) #constraint (10) for edges
             weights_other=[H.edges[j,t]['weight'] for t in set(H[j])-set([i])]
             weights_other.append(weight_ij-1)
             weight=min(weights_other)
             for l in range(1, weight+1):
                model.addConstr(u_connect[j,i,l] == 0) #constraint (13) for edges
             for l in range(weight_ij+1, L+1):
                model.addConstr(u_connect[j,i,l] == u_connect[j,i,weight_ij]) #constraint (12)                                    
           
    #----constraints on non-adjacent node pairs (i,j) not in E
    for j in H.nodes():
        for i in H.nodes():
            if i not in H.neighbors(j) and i<j:
                for l in range(2, L): 
                    model.addConstr(u_connect[i,j,l] <= u_connect[i,j,l+1]) #constraint (11)     

    #------update model and set gurobi model parameters
    model.update()
    model.setObjective(obj, GRB.MINIMIZE)
    model._x_delete=x_delete
    model._u_connect=u_connect
    model.setParam(GRB.param.Cuts, 0)
    model.setParam('LogToConsole', 0)
    model.setParam(GRB.param.PreCrush, 1)
    model.setParam('LazyConstraints', 1)
    model.setParam('TimeLimit', 3600)
    model.optimize(mycut)
    run_time=model.Runtime
    xval=model.getAttr('x', x_delete)
   
    #-----retrieve critical nodes(that is nodes with value of x-variable equal to one)
    critical_nodes=[i for i in xval.keys() if xval[i]>=1- 1e-4]
    
    #-----retrieve value of objective function
    opt_obj = 0
    for j in H.nodes():
        for i in H.nodes():
            if i < j:
                opt_obj+=f[0] * u_connect[i, j, 1].X
                for l in range(1,L):
                    opt_obj+=f[l] * (u_connect[i, j, l + 1].X - \
                              u_connect[i, j, l].X)  
    
    return critical_nodes,opt_obj, run_time, model.Runtime

#------------------Main body------------------------------
#-----read graphs
#G=nx.read_edgelist(path="hi_tech.edgelist", nodetype=int)
#G=nx.karate_club_graph()
#G=nx.read_edgelist(path="mexican.edgelist", nodetype=int)
#G=nx.read_edgelist(path="chesapeake.edgelist", nodetype=int)
#G=nx.read_gml(path="lesmiserable.gml", label='id')
#G=nx.read_edgelist(path="Sawmill.edgelist", nodetype=int)
#G=nx.read_gml(path="dolphins.gml", label='id')
#G=nx.read_edgelist(path="santafe.edgelist", nodetype=int)
#G=nx.read_edgelist(path="Sanjuansur2.edgelist", nodetype=int)
#G=nx.read_edgelist(path="attiro.edgelist", nodetype=int)
G=nx.read_edgelist(path="SmallWorld.edgelist", nodetype=int)


#-----assign weights to edges based on edge betweenness------------
edge_btw=nx.edge_betweenness_centrality(G, normalized=True)
#-----set parameter b to 1 , 2, 3, 4, 5 and 6 respectively for Smallworld graph 
#----(Table 2 of paper) and set b=6 for the rest of the graphs (Table 1 of paper)
b=6
for e in G.edges():
    G.edges[e[0],e[1]]['weight'] = min(b,max(1,round(0.1/edge_btw[e]))) 
    
bndcheck=0
n=G.number_of_nodes()
L=dc.cal_diameter(G) #diameter of the graph
C=int(0.05*n) #budget on critical nodes B=0.05n in paper
ind=0

#----define distance connectivity function (eg f(d)=1/d)
f=[]
for l in range(L+1):
    f.append(1/float(l+1))
    
#----find the critical nodes
critical_nodes,opt_obj,run_time,cpu_time=main(G,L,C)
#-----print
print("critical nodes are :", critical_nodes)
print("Running Time = {:.2f} seconds".format(run_time))
print('Final objective = {:.2f}'.format(opt_obj))
print('Final objective percentage = {:.2f} %'.format(2*100*opt_obj/(n*(n-1))))

    