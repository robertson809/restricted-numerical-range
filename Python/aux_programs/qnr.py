# NR: Numerical Range Plot
#
# Author: Thomas R. Cameron
# Date: 5/30/2019
import numpy as np
from matplotlib import pyplot as plt
import networkx as nx
import itertools
from math import pi as pi
import sys
from operator import itemgetter
##Global variables are realllly bad practice but otherwise
##I'd have to pass things all over the place
singleton_index_list = []
###############################################
###             Numerical Range             ###
###############################################
def nr(a, graph_num):
    """Plots numerical range of a matrix with its eigenvalues."""
    nv = 120 #some kind of precision
    m, n = a.shape
    if(m!=n):
        print('Warning: matrix is non-square')
        return
    else:
        e = np.linalg.eigvals(a)
        f = []
        for k in range(1,nv+1):
            z = np.exp(2*pi*1j*(k-1)/nv) #roots of unity
            a1 = z*a
            a2 = (a1 + np.transpose(np.conjugate(a1)))/2 #A_2 is Hermitian part, and the NR of A_2 is the Real part of NR(A)
            w, v = np.linalg.eig(a2) 
            ind = np.argsort(w)
            w = w[ind] #sorts the array
            v = v[:,ind]
            v = v[:,n-1]
            #gives a point on the nr
            f.append(np.dot(np.conjugate(v),np.dot(a,v))/np.dot(np.conjugate(v),v))
        f.append(f[0])
        f = np.array(f)
        if is_singleton(f):
            singleton_index_list.append(graph_num)
        plt.subplot(122)
        plt.plot(np.real(f),np.imag(f))
        plt.plot(np.real(e),np.imag(e),'r*')
        #plt.gca().set_aspect('equal', adjustable='box')

        
###############################################
###             Q Numerical Range           ###
###############################################
def qnr(l, graph_num, title):
    """Plots q numerical range of graph Laplacian."""
    m, n = l.shape
    if(m!=n):
        print('Warning: matrix is non-square')
        return
    q = np.zeros((n,n-1))
    for j in range(n-1):
        q[0:j+1,j] = 1
        q[j+1,j] = -(j+1)
        q[:,j] = q[:,j]/np.linalg.norm(q[:,j])
    m = np.dot(np.transpose(q), np.dot(l,q))
    print('normal size off by', np.linalg.norm(m.transpose() @ m - m
                   @ m.transpose(), 'fro'))
    if abs(np.linalg.norm(m.transpose() @ m - m @ m.transpose(), 'fro')) <= \
            np.linalg.norm(m) * sys.float_info.epsilon:
        title = title + ' [Q(N)] '
    print("******** Q Matrix {} ********".format(graph_num))
    print(m)
    nr(m, graph_num)
    return title

###############################################
###             Unique Permutation          ###
###############################################
def unique_perm(a,n):
    ind = range(n)
    perm = []
    for i in itertools.permutations(ind):
        b = a[i,:]
        b = b[:,i]
        add = not any(np.array_equal(b,x) for x in perm)
        if(add):
            perm.append(b)
    return perm


def is_singleton(f):
    """
    checks to see if the shape described by the points in the array
    are really all the same point
    """
    return np.var(f) < 0.001
    
    
###############################################
###             allQNR                      ###
###############################################
def allQNR(n, r_graph_num = -1):
    """
    r_graph_num is for requested graph number, it will only
    print the requested graph-numerical range pair
    """
    adj = [] #adjacency matrix of all i-unique graphs
    #find all isomorphically unique graphs
    for i in itertools.product([0, 1], repeat = n*n):
        a = np.reshape(np.array(i),(n,n))
        if(np.trace(a)==0):
            perm = unique_perm(a,n)
            add = not any(np.array_equal(p,x) for p in perm for x in adj)
            if(add):
                adj.append(a)
    print('There are {} isomorphically unique graphs '
        'with {} vertices'.format(len(adj), n)) 
        
    if r_graph_num != -1: #turn adj into a singleton to only return one graph
        if type(r_graph_num) == int:
            g = adj[r_graph_num]
            adj = [g]
        if type(r_graph_num) == list:
            adj = [adj[num] for num in r_graph_num]
        
    graph_num = 0  
    for a in adj:
        x = np.array([np.sum(a[i,:]) for i in range(n)])
        l = np.diag(x) - a
        print("******** Graph Laplacian {} ********".format(graph_num))
        print(l)
        title = ''
        print('normal size is off by {}'.format(np.linalg.norm(l.transpose() @ l - l
                               @ l.transpose(), 'fro')))
        if np.linalg.norm(l.transpose() @ l - l
                               @ l.transpose(), 'fro') == 0:
            title = title + ' [L(N)] '
        g = nx.DiGraph(a)
        title = qnr(l, graph_num, title)
        nx.draw(g, with_labels=True, ax = plt.subplot(121))
        if r_graph_num == -1:
            plt.title(title + 'Graph {}'.format(graph_num))
        elif type(r_graph_num) == int:
            plt.title('Graph {}'.format(r_graph_num))
        elif type(r_graph_num) == list:
            plt.title('Graph {}'.format(r_graph_num[graph_num]))

        plt.show()
        #print('showing the numerical range of just L')
        #nr(l, graph_num)
        graph_num += 1 

#The imploding stars for 3 graphs are 0, 6, 13 and 15
#for 4 they are 0, 76, 176, 213, and 217
###############################################
###             impStar                     ###
###############################################
def impStar(n):
    a = np.zeros((n,n))
    for i in range(n-1):
        a[i,n-1] = -1
    g = nx.DiGraph(a)
    nx.draw(g)
    plt.show()
    x = np.array([np.sum(a[i,:]) for i in range(n)])
    l = np.diag(x) - a
    qnr(l)
    
baseball = [15,20, 39]
pseduocycles = [27, ]
lines = [21]
three_IS = [0, 6, 13, 15]
four_IS = [0, 76, 176, 213, 217]
#two exploding star for five vertices. 
#look at k exploding stars and lines 
#what are lines? We have a theory that they are either exploding stars
# or they are k -- imploding stars unioned with disjoint isolated vertices.

#check to see if the 2 exploding star on n vertices is always the line from 0 to n
# or it's the 3 exploding star, or it's the n -2 exploding

#allQNR(4, r_graph_num= 76)
# allQNR(4, r_graph_num = four_IS)
allQNR(3)
alex_graph = [[2, -1, 0,0,-1],
              [0,1,-1,0,0],
              [0,0,1,-1,0],
              [0,0,0,1,-1],
              [-1,0,0,0,1]]

test_graph = np.array([[1,-1,0,0,0],
              [0,1,-1,0,0],
              [-1,0,1,0,0],
              [0,0,0,1,-1],
              [0,0,0,-1,1]])
#qnr(test_graph, 1)
#k = test_graph.transpose() @ test_graph - test_graph @ test_graph.transpose()
#print(k, np.linalg.norm(k, 'fro'))
print('the singletons are at', singleton_index_list)

#shouldn't be possible that graph
# 3 on 3 vertices has normal laplacian but not normal 1



