# Polygonal Graphs
#
# This module finds graphs with polygonal numerical range
#
# Author: Michael D. Robertson and Thomas R. Cameron
# Date: 3/16/2020

import numpy as np
import networkx as nx
import warnings
import time
from math import pi as pi
from matplotlib import pyplot as plt
import sys

warnings.filterwarnings("ignore", category=UserWarning)

EPS = np.finfo(float).eps
graph_num = 1
title = ""
verbose = False




def convex_hull(points):
    """Computes the convex hull of a set of 2D points.

    Input: an iterable sequence of (x, y) pairs representing the points.
    Output: a list of vertices of the convex hull in counter-clockwise order,
      starting from the vertex with the lexicographically smallest coordinates.
    Implements Andrew's monotone chain algorithm. O(n log n) complexity.
    """

    # Sort the points lexicographically (tuples are compared lexicographically).
    # Remove duplicates to detect the case we have just one unique point.
    points = sorted(set(points))

    # Boring case: no points or a single point, possibly repeated multiple times.
    if len(points) <= 1:
        return points

    # 2D cross product of OA and OB vectors, i.e. z-component of their 3D cross product.
    # Returns a positive value, if OAB makes a counter-clockwise turn,
    # negative for clockwise turn, and zero if the points are collinear.
    def cross(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    # Build lower hull 
    lower = []
    for p in points:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)

    # Build upper hull
    upper = []
    for p in reversed(points):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)

    # Concatenation of the lower and upper hulls gives the convex hull.
    # Last point of each list is omitted because it is repeated at the beginning of the other list. 
    return lower[:-1] + upper[:-1]


def check_hull(conv_hull, x, y):
    """
    Returns true if an (x,y), or (real, complex) point lies on a convex hull
    represented by the conv_hull variable
    @param conv_hull: array of points representing the convex hull
    @param x: real or x part of the point in question
    @param y: imaginary or y part of the point in question
    @return: true if the point lies on the convex hull
    """
    tol = 1000 * EPS
    r0 = conv_hull[-1]
    res = False
    for k in range(len(conv_hull)):
        r1 = conv_hull[k]
        a0 = [r1[0] - r0[0], r1[1] - r0[1]]
        a1 = [x - r0[0], y - r0[1]]
        if np.linalg.norm(a0) <= tol:
            res = res or (np.linalg.norm(a1) <= tol)
        else:
            a = np.zeros((2, 2))
            a[:, 0] = a0
            a[:, 1] = a1
            s = np.linalg.svd(a, compute_uv=False)
            res = res or (s[1] < tol * s[0])
        r0 = r1
    return res


def nr(mat):
    """
    Numerical Range -- returns the numerical range of a matrix mat, as an
    array of complex values.

    @todo what is being done to calculate the numerical range here
    @param mat: square matrix
    @return: numerical range of adj and the eigenvalues of adj
    """
    num_values = 360
    n, _ = mat.shape
    numerical_range = []
    for k in range(1, num_values + 1):
        z = np.exp(2 * pi * 1j * (k - 1) / num_values) # roots of unity
        a1 = z * mat
        a2 = (a1 + np.transpose(np.conjugate(a1))) / 2 # A_2 is Hermitian part, and the NR of
        # A_2 is the Real part of NR(A)
        w, v = np.linalg.eig(a2)
        ind = np.argsort(w)
        w = w[ind]
        v = v[:, ind]
        v = v[:, n - 1]
        val = np.dot(np.conjugate(v), np.dot(mat, v)) / np.dot(np.conjugate(v), v)

        # kill the roundoff
        if 0 < abs(val.real) < EPS * 10:
            val = val - val.real
        if 0 < abs(val.imag) < EPS * 10:
            val = val - val.imag * 1j

        numerical_range.append(val)
    numerical_range.append(numerical_range[0])
    numerical_range = np.array(numerical_range)
    evals = np.around(np.linalg.eigvals(mat), decimals=13)
    return numerical_range, evals


def qnr(lap):
    """
    Finds the Q numerical range of a laplacian matrix
    @param lap: laplacian matrix
    @return: the restircted numerical range of lap
    """
    n, _ = lap.shape

    # form Q matrix in order to restrict
    q = np.zeros((n, n - 1))
    for j in range(n - 1):
        q[0:j + 1, j] = 1
        q[j + 1, j] = -(j + 1)
        q[:, j] = q[:, j] / np.linalg.norm(q[:, j])

    # restrict lap to form a
    a = np.dot(np.transpose(q), np.dot(lap, q))


    # plot Q matrix, we decided not to do this
    # plt.subplot(224)
    # plt.matshow(a, fignum=False, cmap='winter')
    # for i in range(n):
    #     for j in range(n):
    #         c = lap[j, i]
    #         plt.text(i, j, str(c), va='center', ha='center', color='white')

    global title
    if abs(np.linalg.norm(a.transpose() @ a - a @ a.transpose(), 'fro')) <= np.linalg.norm(a, 'fro') * \
            sys.float_info.epsilon * 100:
        title = title + ' N_Q '
    return nr(a)


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


def is_polygon(f, e):
    """
    checks if the shape of the points in f is the convex hull of the eigenvalues e
    @param f: numpy complex array with points representing the numerical range
    @param e: array of the eigenvalues of the matrix
    @return: true if the numerical range is a polygon that is the convex hull
    of the eigenvalues.
    """
    edges = []
    for eig in e:
        edges.append((eig.real, eig.imag))
    conv_hull = convex_hull(edges)
    res = True
    for z in f:
        res = res and check_hull(conv_hull, z.real, z.imag)
    return res


def restricted_nr(lap):
    """
    Plots both the restricted and unrestricted numerical ranges
    @param lap: lapalcian matrix
    @return: the numerical range eigenvalue pair for the restricted
    and unrestricted numerical ranges
    """
    nr_restricted, e_restricted = qnr(lap)
    nr_unrestricted, e_unrestricted = nr(lap)

    # restricted numerical range
    # plt.subplot(121)
    # plt.title('NR(Q)')
    plt.plot(np.real(nr_restricted), np.imag(nr_restricted),'#03a5e0', linewidth=2.5)
    plt.plot(np.real(e_restricted), np.imag(e_restricted), '*',
             linestyle='None', marker='*', color='#0337e0', markersize = 8)
    plt.fill(np.real(nr_restricted), np.imag(nr_restricted), '#03a5e0')  # fill
    plt.show()

    # unrestricted numerical range
    # plt.subplot(122)
    # plt.title('NR(L)')
    # plt.plot(np.real(nr_unrestricted), np.imag(nr_unrestricted))
    # plt.plot(np.real(e_unrestricted), np.imag(e_unrestricted), 'g*')
    # plt.show()
    return nr_restricted, e_restricted


def determine_polygon(adj, count_vect=[0, 0, 0], singleton_adj=[], line_adj=[], poly_adj=[], pt=-1, disp=False):
    """
    Determines whether an an adjacency matrix has a singleton, line, or polygon numerical range, and saves
    four-part figures which include a visualization of the graph, the Lapalcian matrix, and the restricted
    and unrestricted numerical ranges.
    @param adj: adjacency matrix
    @param count_vect: three vector with a position to count the singletons, lines and polygons
    @param singleton_adj: list of adjacency matrices with singleton numerical ranges
    @param line_adj: list of adjacency matrices with line numerical ranges
    @param poly_adj: list of adjacency matrices with polygonal numerical range
    @param pt: used to allow for parallel computing with six or more vertices
    @param disp: if true, displays the figure and exits rather than saving it
    """
    global graph_num
    global title
    title = ""
    g = nx.DiGraph(adj)
    n = len(adj)

    # include graph in figure
    # nx.draw_shell(g, with_labels=True, ax=plt.subplot(221))

    # form Laplacian
    diag = np.array([np.sum(adj[i, :]) for i in range(n)])
    lap = np.diag(diag) - adj
    if not disp:
        print(lap)
        print(graph_num)

    # test for normality
    # if np.linalg.norm(lap.transpose() @ lap - lap @ lap.transpose(), 'fro') == 0:  # normality test
    #     title = title + ' N_L '
    # plt.subplot(221)

    # include graph number in figure
    if pt != -1:
        plt.title(title + "Graph Number %d.%d" % (pt, graph_num))
    if disp:
        print('could include a title in line 267')
        # plt.title()
    # else:
    #     plt.title(title + "Graph Number %d" % graph_num)

    # plt.subplot(222)

    # include Laplacian in figure
    # plt.matshow(lap, fignum=False, cmap='winter')
    # for k in range(n):
    #     for j in range(n):
    #         c = lap[j, k]
            # plt.text(k, j, str(c), va='center', ha='center', color='white')

    # find unrestricted numerical ranges and include both restricted and unrestricted in figure
    nr_restricted, eig_restricted = restricted_nr(lap)
    fig = plt.gcf()

    if disp:
        plt.show()
        exit()
    # in figure nad save figure
    if is_singleton(nr_restricted):
        fig.savefig("figures/%d_all/singletonGraph%d/singletonGraph%d.png" % (n, n, graph_num), dpi=400)
        count_vect[0] += 1
        singleton_adj.append(lap)
    elif is_line(nr_restricted):
        fig.savefig("figures/%d_all/lineGraph%d/lineGraph%d.png" % (n, n, graph_num), dpi=400)
        count_vect[1] += 1
        line_adj.append(lap)
    elif is_polygon(nr_restricted, eig_restricted):
        count_vect[2] += 1
        poly_adj.append(lap)
        fig.savefig("figures/%d_all/polyGraph%d/polyGraph%d.png" % (n, n, graph_num), dpi=400)
    else:
        fig.savefig("figures/%d_all/other/graph%d.png" % (n, graph_num), dpi=400)

    plt.clf()
    graph_num = graph_num + 1


def write_out(adjs, n):
    """
    Write out the adjacency matrices for singleton, line, and polygon matrices
    @param adjs: tuple of three adjacency matrices, singleton, line and polygon
    @param n: rank of matrix
    """
    names = ['singleton', 'line', 'polygon']
    for name, mat in zip(names, adjs):
        # outfile = open('matrices/polygon_adjs/{}/{}_adjs.txt'.format(n, name), 'w+')
        outfile = open('defense_figs/polygon_adjs/{}/{}_adjs.txt'.format(n, name), 'w+')
        for g in mat:
            for i in range(n):
                for j in range(n - 1):
                    outfile.write("%d " % g[i, j])
                outfile.write("%d\n" % g[i, n - 1])
            outfile.write("\n")
        outfile.close()


def find_poly_graphs(n, pt=-1):
    """
    Given appropriately formatted (see nauty_cleaner) files of adjacency matrices of directed graphs, will find the
    matrices which have singleton, line, and polygon numerical range, and save images containing the laplacian,
    visualization of the graph, and restricted and unrestricted numerical ranges.
    @param n: number of vertices in the graphs
    @param pt: arguement necessary when the text files have been split into many parts in order to keep track of the
    graph indexing
    """
    start_time = time.time()
    if pt != -1:
        infile = open('matrices/directed_adjacency_mats/adj' + str(n) + '_pt' + str(pt) + '.txt')
    else:
        infile = open('matrices/directed_adjacency_mats/adj' + str(n) + '.txt')
    # infile = open('matrices/directed_adjacency_mats/4complete.txt')


    line_list = infile.readlines()
    a = np.zeros((n, n))
    mat_line_count = 0
    count_vect = np.zeros(3)  # counts the number of singletons, lines, and polygons in the 0th, 1st, and 2nd positions
    singleton_adj = []
    line_adj = []
    poly_adj = []
    # print(line_list)
    for row in line_list:
        if graph_num % 1000 == 0:
            print('Examined graph {}.{}'.format(pt, graph_num))
            print('We have {} singletons, {} lines, and {} polygons'.format(count_vect[0], count_vect[1],
                                                                            count_vect[2]))
            print('Total time elapsed ---- %.1f hours' % ((time.time() - start_time) / 3600))
        row = row.split(" ")
        if len(row) < n:  # we've read a full adjacency matrix
            print(a)
            determine_polygon(a, count_vect, singleton_adj, line_adj, poly_adj, pt)
            a = np.zeros((n, n))  # reset adjacency matrix
            mat_line_count = 0  # reset line count
        else:  # not the end of a matrix in infile, read the next line
            a[mat_line_count, :] = [eval(row[j]) for j in range(n)]
            mat_line_count = mat_line_count + 1

    write_out((singleton_adj, line_adj, poly_adj), n)

    report_file = open('reports/%d_report.txt' % n, 'w+')
    report_file.write('%d singletons \n%d lines \n%d polygons' % (count_vect[0], count_vect[1], count_vect[2]))


def plt_nr(l, restricted = True, save_num = 0):
    """
    Plots the numerical range of a laplacian using matplotlib
    @param a: square matrix
    """
    if restricted:
        f, e = qnr(l)  # make it q instead
        # print('the eigenvalues are', e)
    else:
        f, e = nr(l)

    # color options
    # print('the eigen values are', e)
    # close blue #0337e0
    # 03a5e0 - 'calming blue
    # outside green - #8da63a
    # 03e0ad
    # contrast orange - #e03e03
    # ac1a2f davidson wildcat red
    # plt.plot(np.real(f), np.imag(f),'#03a5e0', linewidth = 2.5) #boundary
    # plt.plot(np.real(e), np.imag(e), linestyle = 'None', marker =  '*',
    #          color = '#0337e0', markersize=8) #evals
    # plt.fill(np.real(f), np.imag(f), '#03a5e0') #fill
    #
    # plt.show()
    # print(a)

    # plt.plot(np.real(f), np.imag(f), '#03a5e0', linewidth=2.5)  # boundary
    plt.plot(np.real(e), np.imag(e), linestyle='None', marker='*', color='#0337e0', markersize=8)
    plt.plot(np.real(e), np.imag(e), linestyle='None', marker='*', color='#0337e0', markersize = 8)  # evals
    plt.title('Adrian\'s Restricted Numerical Range')
    plt.fill(np.real(f), np.imag(f), '#03a5e0')# fill

    # plt.xlim(-1, 4)
    # plt.ylim(-2, 2)
    # plt.gca().set_aspect('equal', adjustable='box')

    # either save or
    # fig = plt.gcf()
#     fig.savefig("nice_poly_nr_figs/polyGraph%d.png" % save_num, dpi=400)
#     plt.clf()

    # display
    plt.show()

    


def main():
    find_poly_graphs(4)


if __name__ == '__main__':
    s61 = np.array([[1,0,0,0,0,-1],[0,1,0,0,0,-1],[0,0,1,0,0,-1],[0,0,0,1,0,-1],
                    [0,0,0,0,1,-1],[0,0,0,0,0,0]])
    s62 = np.array([[2,0,0,0,-1,-1],[0,2,0,0,-1,-1],[0,0,2,0,-1,-1],
                    [0,0,0,2,-1,-1],[0,0,0,0,1,-1],[0,0,0,0,-1,1]])
    six_star = np.array([[0,0,0,1,1,1],[0,0,0,1,1,1],[0,0,0,1,1,1],[0,0,0,0,1,1],
                         [0,0,0,1,0,1],[0,0,0,1,1,0]])
    six_star_deg = np.array([[0,0,1,0,1,1],[0,0,0,1,1,1],[0,0,0,1,1,1],[0,0,0,0,1,1],
                         [0,0,0,1,0,1],[0,0,0,1,1,0]])
    cam_ex = np.array([[0,1,0,1],[0,0,1,1],[1,0,0,1],[0,0,0,0]])
    empty = np.array([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]])
    complete = np.array([[0,1,1,1],[1,0,1,1],[1,1,0,1],[1,1,1,0]])
    four_cycle = np.array([[0,1,0,0],[0,0,1,0],[0,0,0,1],[1,0,0,0]])
    reg_torn = np.array([[0,1,1,0,0],[0,0,1,1,0],[0,0,0,1,1],[1,0,0,0,1],[1,1,0,0,0]])
    cam_ex = np.array(
        [[2, -1, 0, 0, -1, 0], [0, 1, -1, 0, 0, 0], [-1, 0, 1, 0, 0, 0], [-1, 0, 0, 1, 0, 0], [0, 0, 0, -1, 2, -1],
         [0, 0, 0, 0, -1, 1]], dtype=float)



    # determine_polygon(four_cycle, disp=True)
    # determine_polygon(x, disp=True)
    # determine_polygon(reg_torn, disp=True)

    empty_4 = np.array([[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]])
    complete_2 = np.array([[1,-1],[-1,1]])
    complete_3 = np.array([[2,-1,-1],[-1,2,-1],[-1,-1,2]])
    complete_4 = np.array([[3,-1, -1, -1], [-1, 3,-1, -1], [-1, -1, 3,-1],[-1,-1,-1,3]])
    four_cycle = np.array([[1, -1, 0, 0], [0, 1, -1, 0], [0, 0, 1, -1], [-1, 0, 0, 1]])
    five_cycle = np.array([[1, -1, 0, 0, 0], [0, 1, -1, 0, 0], [0, 0, 1, -1, 0],
                          [0, 0, 0, 1, -1], [-1, 0, 0, 0, 1]])
    six_cycle = np.array([[1, -1, 0, 0, 0, 0], [0, 1, -1, 0, 0, 0], [0, 0, 1, -1, 0, 0],
                          [0, 0, 0, 1, -1, 0], [0, 0, 0, 0, 1, -1], [-1, 0, 0, 0, 0, 1]])
    s61 = np.array([[1, 0, 0, 0, 0, -1], [0, 1, 0, 0, 0, -1], [0, 0, 1, 0, 0, -1], [0, 0, 0, 1, 0, -1],
                    [0, 0, 0, 0, 1, -1], [0, 0, 0, 0, 0, 0]])
    s62 = np.array([[2, 0, 0, 0, -1, -1], [0, 2, 0, 0, -1, -1], [0, 0, 2, 0, -1, -1],
                    [0, 0, 0, 2, -1, -1], [0, 0, 0, 0, 1, -1], [0, 0, 0, 0, -1, 1]])
    s63 = np.array([[3, 0, 0, -1, -1, -1], [0, 3, 0, -1, -1, -1], [0, 0, 3, -1, -1, -1],
                    [0, 0, 0, 2, -1, -1], [0, 0, 0, -1, 2, -1], [0, 0, 0, -1, -1, 2]])
    reg_torn_5 = np.array([[2, -1, -1, 0, 0], [0, 2, -1, -1, 0], [0, 0, 2, -1, -1],
                         [-1, 0, 0, 2, -1], [-1, -1, 0, 0, 2]])
    reg_torn_7 = np.array([[3, -1, -1, -1, 0,0,0], [0, 3, -1, -1, -1,0,0], [0, 0, 3, -1, -1,-1,0],
                           [0, 0, 0, 3, -1,-1,-1], [-1, 0, 0, 0, 3,-1,-1],[-1,-1,0,0,0,3,-1],
                           [-1,-1,-1,0,0,0,3]])
    reg_torn_9 = np.array([[4, -1, -1, -1,-1,0,0,0, 0], [0, 4, -1, -1, -1, -1,0,0,0],
                           [0, 0, 4, -1, -1,-1,-1,0,0],
                           [0, 0, 0, 4, -1,-1,-1,-1,0], [0, 0, 0, 0, 4,-1,-1,-1,-1],
                           [-1,0,0,0,0,4,-1,-1,-1],[-1,-1,0,0,0,0,4,-1,-1],[-1,-1,-1,0,0,0,0,4,-1],
                           [-1,-1,-1, -1,0,0,0,0,4]])
    restricted_nr(reg_torn_9)
    
    
    three_balanced_4 = np.array([[2,0,-1,-1],[0,2,-1,-1],[0,0,0,0],
                                 [0,0,0,0]])
    three_balanced_5 = np.array([[1,-1,0,0,0],[-1,2,-1,0,0],[0,-1,1,0,0],
                                 [-1,-1,-1,3,0],[-1,-1,-1,0,3]])

    case_1 = np.array([[2,0,0,-1,0,-1],[0,1,0,0,0,-1],[0,0,1,0,0,-1],[0,0,0,2,-1,-1],
                      [-1,0,0,0,2,-1],[-1,-1,-1,-1,-1,5]])
    case_2 = np.array([[1,0,0,-1,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,1,-1,0],
                       [-1,0,0,0,1,0],[-1,-1,-1,-1,-1,5]])

    case_3 = np.array([[2,0,-1,-1,0,0],[0,2,-1,-1,0,0],[-1,-1,4,-1,-1,0],
                       [-1,-1,-1,4,0,-1],[-1,-1,0,-1,4,-1],[-1,-1,-1,0,-1,4]])

    case_3_comp = np.array([[3, -1, 0, 0, -1, -1], [-1, 3, 0, 0, -1, -1], [0, 0, 1, 0, 0, -1],
                       [0, 0, 0, 1, -1, 0], [0, 0, -1, 0, 1, 0], [0, 0, 0, -1, 0, 1]])


    six_star = np.array([[3, 0, 0, -1, -1, -1], [0, 3, 0, -1, -1, -1], [0, 0, 3, -1, -1, -1], [0, 0, 0, 2, -1, -1],
                         [0, 0, 0, -1, 2, -1], [0, 0, 0, -1, -1, 2]])
    six_star_deg = np.array([[3, 0, -1, 0, -1, -1], [0, 3, 0, -1, -1, -1], [0, 0, 3, -1, -1, -1], [0, 0, 0, 2, -1, -1],
                             [0, 0, 0, -1, 2, -1], [0, 0, 0, -1, -1, 2]])
    # plt_nr(three_balanced_4)
    # plt_nr(six_cycle)
    # plt_nr(complete_3)
    # plt_nr(complete_4, restricted=False)
    # plt_nr(complete_4)
    # plt_nr(case_3)
    # plt_nr(case_3_comp)
    # plt_nr(six_star, restricted=False)
    # plt_nr(six_star_deg, restricted=False)
    # plt.show()
    # plt_nr(case_2)
    # plt_nr(case_3)


    # plt_nr(np.transpose(x))
    # plt_nr(astar)
    # find_poly_graphs(4)
    
    
    six_cycle = np.array([[1, -1, 0, 0, 0, 0], [0, 1, -1, 0, 0, 0], [0, 0, 1, -1, 0, 0],
                          [0, 0, 0, 1, -1, 0], [0, 0, 0, 0, 1, -1], [-1, 0, 0, 0, 0, 1]])
    # plt_nr(six_cycle)
    
                        ######  #######    #     #####  #     # 
                        #     # #         # #   #     # #     # 
                        #     # #        #   #  #       #     # 
                        ######  #####   #     # #       ####### 
                        #   #   #       ####### #       #     # 
                        #    #  #       #     # #     # #     # 
                        #     # ####### #     #  #####  #     # 
                        
    arsi = np.array([
        [2,0,-1,-1,0,0],
        [-1,4,0,-1,-1,-1],
        [-1,-1,3,0,-1,0],
        [-1,-1,-1,4,0,-1],
        [0,-1,-1,0,3,-1],
        [0,-1,0,0,-1,2],
    ])
    dj = np.array([
        [2, -1, 0, -1],
        [0, 0, 0, 0],
        [-1, -1, 3, -1],
        [-1, -1, -1, 3]
    ])
    joshua = np.array([
        [3, -1, 0, -1, -1],
        [-1, 2, 0, -1, 0],
        [-1, -1, 4, -1, -1],
        [-1, -1, -1, 3, 0],
        [0, -1, -1, -1, 3]])
    christian = np.array([[1,0,-1,0,0],
    [0,2,-1,-1,0],
    [0,-1,3,-1,-1],
    [0,0,-1,1,0],
    [-1,0,-1,0,2]])
    lucas = np.array([
        [2, 0, 0, -1, 0, 0, -1],
    [0, 2, -1, 0, 0, 0, -1],
    [-1, -1, 4, -1, 0, -1, 0],
    [-1, -1, 0, 3, -1, 0, 0],
    [0, 0, -1, -1, 2, 0, 0],
    [0, 0, 0, 0, 0, -1, -1],
    [0, 0, 0, -1, 0, -1, 2]])

    adrian = np.array([
    [2,-1,-1,0,0],
    [-1,1,0,0,0],
    [0,-1,2,-1,0],
    [-1,-1,-1,4,-1],
    [-1,-1,-1,-1,4]])
    

    print('plotting DJ')
    plt_nr(dj)
    print('plotting arsi')      
    plt_nr(arsi)
    print('plotting joshua')
    plt_nr(joshua)
    print('plotting christian')
    plt_nr(christian)
    print('plotting lucas')
    plt_nr(lucas)
    print('plotting adrian')
    plt_nr(adrian)

# main()
