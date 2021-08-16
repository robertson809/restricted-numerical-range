# Polygonal
#
# This module determines graphs with polygonal numerical range
#
# Author: Thomas R. Cameron
# Date: 3/16/2020
import numpy as np
import networkx as nx
import sys
from multiprocessing import Process
import warnings
warnings.filterwarnings("ignore", category=UserWarning)
from math import pi as pi
from matplotlib import pyplot as plt

EPS = np.finfo(float).eps       ## 2**(-52) for double precision
TOL = EPS*2**4

###############################################
###             Cmplx Convex Hull           ###
###############################################
def cmplxConvHull(points):
    # Sort the points lexicographically (tuples are compared lexicographically).
    # Remove duplicates to detect the case we have just one unique point.
    points = sorted(set(points))
    # Boring case: no points or a single point, possibly repeated multiple times.
    if(len(points) <= 1):
        return points
    # 2D cross product of OA and OB vectors, i.e. z-component of their 3D cross product.
    # Returns a positive value, if OAB makes a counter-clockwise turn,
    # negative for clockwise turn, and zero if the points are collinear.
    def _cross(o, a, b):
        return (a.real - o.real)*(b.imag - o.imag) - (a.imag - o.imag)*(b.real-o.real)
    # Build lower hull
    lower = []
    for p in points:
        while((len(lower)>=2) and (_cross(lower[-2],lower[-1],p)<EPS)):
            lower.pop()
        lower.append(p)
    # Build upper hull
    upper = []
    for p in reversed(points):
        while((len(upper)>=2) and (_cross(upper[-2],upper[-1],p)<EPS)):
            upper.pop()
        upper.append(p)
    # Concatenation of the lower and upper hulls gives the convex hull.
    # Last point of each list is omitted because it is repeated at the beginning of the other list.
    return lower[:-1] + upper[:-1]
###############################################
###             Numerical Range Main        ###
###############################################
def nr_main(a):
    # set size and norm
    n,_ = a.shape
    norm = np.linalg.norm(a,ord='fro')
    # Hermitian test, don't need to check for normal if it's Hermitian
    err = np.linalg.norm(a - np.conj(np.transpose(a)),ord='fro')
    if(err < max(TOL*norm,TOL)):
        #print("Hermitian Case")
        eig = np.linalg.eigvalsh(a)
        ind = np.argsort(eig)
        eig = eig[ind]
        return [eig[0],eig[-1]], eig
    # Normal test
    err = np.linalg.norm(np.dot(np.conj(np.transpose(a)),a) - np.dot(a,np.conj(np.transpose(a))),ord='fro')
    if(err < max(TOL*norm,TOL)):
        print("Normal Case")
        eig = np.linalg.eigvals(a)
        convHull = cmplxConvHull(eig)
        convHull.append(convHull[0])
        return convHull, eig
    # eigenvalues and inscribed polygon vertices (for non-normal matrices)
    e = []; p = []; r = n
    for k in range(r):
        # multiply by complex exponential
        a1 = np.exp(2*pi*1j*k/r)*a
        # compute hermitian part
        a2 = (a1 + np.transpose(np.conj(a1)))/2
        # compute eigenvalues and eigenvectors of hermitian part
        eigval, eigvec = np.linalg.eigh(a2)
        # sort eigenvalues and store maximum eigenvalue and associated inscribed polygonal vertex
        ind = np.argsort(eigval)
        eigval = eigval[ind]
        e.append(eigval[n-1])
        eigvec = eigvec[:,ind]
        v = eigvec[:,n-1]/np.linalg.norm(eigvec[:,n-1])
        p.append(np.dot(np.conj(v),np.dot(a,v)))
    # complete cycle
    e.append(e[0])
    p.append(p[0])
    # compute circumscribed polygon vertices
    q = []
    for k in range(r):
        q.append(np.exp(-2*pi*1j*k/r)*(e[k]+1j*(e[k]*np.cos(2*pi/r)-e[k+1])/np.sin(2*pi/r)))
    # complete cycle
    q.append(q[0])
    # compute area
    s1 = 0; s2 = 0
    for k in range(r):
        s1 = s1 + np.conj(p[k])*p[k+1]
        s2 = s2 + np.conj(q[k])*q[k+1]
    s1 = 0.5*np.imag(s1)
    s2 = 0.5*np.imag(s2)
    # compute error
    err = (s2 - s1)/min(-1,s2)
    # refine approximation
    return _nr_sub(a,e,p,q,n,r,err), np.linalg.eigvals(a)
###############################################
###             Numerical Range Sub         ###
###############################################
def _nr_sub(a,e,p,q,n,r,err):
    while(err>1E-6 and r<2**12):
        # update r
        r = 2*r        
        for k in range(1,r,2):
            # multiply by complex exponential
            a1 = np.exp(2*pi*1j*k/r)*a
            # compute hermitian part
            a2 = (a1 + np.transpose(np.conj(a1)))/2
            # compute eigenvalues and eigenvectors of hermitian part
            eigval, eigvec = np.linalg.eigh(a2)
            # sort eigenvalues and store maximum eigenvalue and associated inscribed polygonal vertex
            ind = np.argsort(eigval)
            eigval = eigval[ind]
            e.insert(k,eigval[n-1])     ## all inserts could be made more efficient
            eigvec = eigvec[:,ind]
            v = eigvec[:,n-1]/np.linalg.norm(eigvec[:,n-1])
            p.insert(k,np.dot(np.conj(v),np.dot(a,v)))      ## all inserts could be made more efficient
        # compute circumscribed polygon vertices
        for k in range(1,r,2):
            q.insert(k,np.exp(-2*pi*1j*k/r)*(e[k]+1j*(e[k]*np.cos(2*pi/r)-e[k+1])/np.sin(2*pi/r)))      ## all inserts could be made more efficient
        # compute area
        s1 = 0; s2 = 0
        for k in range(r):
            s1 = s1 + np.conj(p[k])*p[k+1]
            s2 = s2 + np.conj(q[k])*q[k+1]
        s1 = 0.5*np.imag(s1)
        s2 = 0.5*np.imag(s2)
        # compute error
        err = (s2 - s1)/min(-1,s2)
    # return
    return p
###############################################
###         Restricted Numerical Range      ###
###############################################
def qnr(l):
    n, _ = l.shape
    q = np.zeros((n,n-1),dtype=float)
    for j in range(n-1):
        q[0:j+1,j] = 1
        q[j+1,j] = -(j+1)
        q[:,j] = q[:,j]/np.linalg.norm(q[:,j])
    a = np.dot(np.transpose(q),np.dot(l,q))
    return nr_main(a)
###############################################
###             Check Hull                  ###
###############################################
def checkHull(convHull,z):
    # set counter variables
    k = 0
    n = len(convHull)
    # set initial point and res
    r0 = convHull[-1]
    res = False
    # while there is still another point in convHull
    # and (x,y) is not between a previous pair of points
    while((k<n) and not(res)):
        # set r1
        r1 = convHull[k]
        res = res or (abs(abs(r1-r0) - (abs(r1-z) + abs(z-r0))) < max(TOL,abs(r1-r0)*TOL))
        # update r0 and k
        r0 = r1
        k = k + 1
    # return res
    return res

def is_singleton(f):
    """
    @param f: numpy complex array with points representing the numerical range
    @return: true if all the points are within (our) tolerance of each other,
    false otherwise. Tolerance is very lax so as not to lose any singletons.
    """
    return np.var(f) < 0.001


def is_line(f):
    """
    checks if the shape of the points f is a real line
    @param f: numpy complex array with points representing the numerical range
    @return: true if all the points lie on the real line, false otherwise.
    """
    return np.sum(np.absolute(np.imag(f))) < .0000001

###############################################
###             Is Polygon                  ###
###############################################
def isPolygon(f,e):
    # compute convex hull from eigenvalues
    convHull = cmplxConvHull(e)
    # set res
    res = True
    # for each point in f, check if it lies between two points in convHull
    for z in f:
        res = res and checkHull(convHull,z)
    # return res
    return res


###############################################
###             Normality Tests             ###
###############################################
def normality(l):
    # check if l is normal
    ln = False
    norm = np.linalg.norm(l,ord='fro')
    err = np.linalg.norm(np.dot(np.conj(np.transpose(l)),l) - np.dot(l,np.conj(np.transpose(l))),ord='fro')
    if(err < max(TOL,TOL*norm)):
        ln = True
    # build a = q^{T}lq
    n, _ = l.shape
    q = np.zeros((n,n-1),dtype=float)
    for j in range(n-1):
        q[0:j+1,j] = 1
        q[j+1,j] = -(j+1)
        q[:,j] = q[:,j]/np.linalg.norm(q[:,j])
    a = np.dot(np.transpose(q),np.dot(l,q))
    # check if a is normal
    an = False
    norm = np.linalg.norm(a,ord='fro')
    err = np.linalg.norm(np.dot(np.conj(np.transpose(a)),a) - np.dot(a,np.conj(np.transpose(a))),ord='fro')
    if(err < max(TOL*norm,TOL)):
        an = True
    # return
    return ln, an
###############################################
###             Poly Graphs                 ###
###############################################
def polyGraphs(n,pt):
    # open adjacency matrix file
    if pt == -1:
        mfile = open('matrices/directed_adjacency_mats/adj%d.txt' % n)
    else:
        mfile = open("matrices/adj%d_pt%d_cleaned.txt" % (n,pt))
    # read lines
    lineList = mfile.readlines()
    # paramters
    i = 0
    poly_count = 0
    a = np.zeros((n,n),dtype=float)

    singleton_adj = []
    line_adj = []
    poly_adj = []

    # read all adjacency matrices
    for row in lineList:
        row = row.split(" ")
        if(len(row)<n):
            # build laplacian
            x = np.array([np.sum(a[i,:]) for i in range(n)])
            l = np.diag(x) - a
            # restricted numerical range
            f, e = qnr(l)
            # check if graph is polygonal
            if(isPolygon(f,e)):
                poly_count += 1
                # draw graph
                g = nx.DiGraph(a)
                nx.draw_shell(g,with_labels=True,ax=plt.subplot(221))
                # draw restricted numerical range
                plt.subplot(222)
                plt.title('Restricted Numerical Range')
                plt.plot(np.real(f),np.imag(f),"b")
                plt.plot(np.real(e),np.imag(e),"r*")
                # draw Laplacian
                plt.subplot(223)
                plt.axis("off")
                plt.matshow(l,fignum=False,cmap="winter")
                for i in range(n):
                    for j in range(n):
                        plt.text(i,j,l[j,i],va="center",ha="center",color="white")
                # add text
                ln,an = normality(l)
                plt.subplot(224)
                plt.axis("off")
                if(ln):
                    plt.text(0.5,0.5,"Normal Laplacian",size=12,ha="center")
                    plt.text(0.5,0.0,"Normal Restricted Laplacian",size=12,ha="center")
                elif(an):
                    plt.text(0.5,0.5,"Non-normal Laplacian",size=12,ha="center")
                    plt.text(0.5,0.0,"Normal Restricted Laplacian",size=12,ha="center")
                else:
                    plt.text(0.5,0.5,"Non-normal Laplacian",size=12,ha="center")
                    plt.text(0.5,0.0,"Non-normal Restricted Laplacian",size=12,ha="center")

                # save figure and matrix
                fig = plt.gcf()
                if is_singleton(f):
                    fig.savefig("figures/%d_all/singletonGraph%d/singeltonGraph%d.png" % (n, n, poly_count), dpi=400)
                    singleton_adj.append(l)
                elif is_line(f):
                    fig.savefig("figures/%d_all/lineGraph%d/lineGraph%d.png" % (n, n, poly_count), dpi=400)
                    line_adj.append(l)
                # already determined the graph was polygonal, so if its not a line or singleton
                # its a polygon
                else:
                    fig.savefig("figures/%d_all/polyGraph%d/polyGraph%d.png" % (n, n, poly_count), dpi=400)
                    poly_adj.append(l)

                # clear plot
                plt.clf()
            # reset the i paramter
            i = 0
        else:
            # store row
            a[i,:] = [eval(row[j]) for j in range(n)]
            # update i
            i = i + 1
    write_out((singleton_adj, line_adj, poly_adj), n)


def write_out(adjs, n):
    """
    Write out the adjacency matrices for singleton, line, and polygon matrices
    @param adjs: tuple of three adjacency matrices, singleton, line and polygon
    @param n: rank of matrix
    """
    names = ['singleton', 'line', 'polygon']
    for name, mat in zip(names, adjs):
        # outfile = open('matrices/polygon_adjs/{}/{}_adjs.txt'.format(n, name), 'w+')
        outfile = open('matrices/polygon_adjs/{}/{}_adjs.txt'.format(n, name), 'w+')
        for g in mat:
            for i in range(n):
                for j in range(n - 1):
                    outfile.write("%d " % g[i, j])
                outfile.write("%d\n" % g[i, n - 1])
            outfile.write("\n")
        outfile.close()
    
###############################################
###             main                        ###
###############################################
def main(argv):

    n = int(argv[0])
    pt = int(argv[1])
    polyGraphs(n,pt)


    

    
if __name__ == '__main__':
    main([4,-1])
    # main(sys.argv[1:])