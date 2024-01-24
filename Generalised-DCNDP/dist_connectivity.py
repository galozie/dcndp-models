"""
We use Breadth-First-Search algorithm to calculate the proportion of node pairs
#connected by hop distance less than or equal to k(for example k=3) that is DCNDP-1 objective
"""

import networkx as nx
import collections

def k_distBFS(input_graph, root,k): 
    visited, queue =set([root]), collections.deque([root])  #keep track of all visited nodes and nodes to be checked
    levels = {}         # this dict keeps track of levels
    levels[root]= 0    # depth of root node is 0
   
    pred={}  #this dict keeps track of predecessors
    pred[root]=-1 #predecessor of root node is -1
    while queue: #keep looping until there are no nodes still to be checked
        vertex = queue.popleft() #pop first node from the queue
        for neighbour in input_graph[vertex]:
            if neighbour not in visited: 
                newlevel=levels[vertex]+1
                if  newlevel>k:
                    break #return pred#, levels
                else:
                   levels[neighbour]=newlevel
                   pred[neighbour]=vertex 
                   visited.add(neighbour) #mark neighbours of node as visited to avoid revisiting
                   queue.append(neighbour) #add neighbours of node to queue 
        else:
            continue
        break
    return pred.keys()#, levels

#-----------------------------------
#-----calculate DCNDP-1 objective for graph G for specified k
def k_connectivity(G,k): 
    count=0
    for root in G.nodes():
        k_conn_nodes=len(k_distBFS(G,root,k))-1
        count+=k_conn_nodes
    #print("The k-connectivity of the graph is %i" %(count*0.5))
    #print('The percentage k-connectivity of the graph is = {:.2f}%'.format(100*count/(G.number_of_nodes()*(G.number_of_nodes()-1))))
    k_conn=count*0.5
    k_conn_percent=round((100*count/(G.number_of_nodes()*(G.number_of_nodes()-1))),1)
    return k_conn,k_conn_percent

#--------------------------------------------
#--calculate DCNDP-1 objective for subgraph H of G induced by nodes in the node list
#subgraph_nodes for specified k
def subG_kconnectivity(G,k,subgraph_nodes): 
    H=G.subgraph(subgraph_nodes) 
    count=0
    for root in H.nodes():
        k_conn_nodes=len(k_distBFS(H,root,k))-1
        count+=k_conn_nodes
    k_conn=count*0.5
    #k_conn=round((100*count/(G.number_of_nodes()*(G.number_of_nodes()-1))),3)
    return k_conn

#-------------------------------
#---calculate diameter of a weighted graph
def cal_diameter(g):
    length = dict(nx.all_pairs_dijkstra_path_length(g,weight='weight'))
    eccent=[max(v.values()) for k,v in length.items()]
    #mid=[np.median(list(v.values())) for k,v in length.items()]
    return int(max(eccent))





   
      