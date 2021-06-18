# Polygonal Graphs
#
# This module determines graphs with polygonal numerical range
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


def plt_nr(a):
    """
    Plots the numerical range using matplotlib
    @param a: square matrix
    """
    f, e = nr(a)
    plt.plot(np.real(f), np.imag(f))
    plt.plot(np.real(e), np.imag(e), 'r*')
    plt.show()


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
        a2 = (a1 + np.transpose(np.conjugate(a1))) / 2 # A_2 is Hermitian part, and the NR of A_2 is the Real part of NR(A)
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
    plt.subplot(223)
    plt.title('NR(Q)')
    plt.plot(np.real(nr_restricted), np.imag(nr_restricted))
    plt.plot(np.real(e_restricted), np.imag(e_restricted), 'g*')

    # unrestricted numerical range
    plt.subplot(224)
    plt.title('NR(L)')
    plt.plot(np.real(nr_unrestricted), np.imag(nr_unrestricted))
    plt.plot(np.real(e_unrestricted), np.imag(e_unrestricted), 'g*')
    return nr_restricted, e_restricted


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


def determine_polygon(adj, count_vect, singleton_adj, line_adj, poly_adj, pt):
    """
    Determines whether an an adjacency matrix has a singleton, line, or polygon numerical range, and saves
    four-part figures which include a visualization of the graph, the Lapalcian matrix, and the restricted
    and unrestricted numerical ranges.
    @param adj: adjacency matrix
    @param count_vect: counts the number of singletons, lines and polygons
    @param singleton_adj: list of adjacency matrices with singleton numerical ranges
    @param line_adj: list of adjacency matrices with line numerical ranges
    @param poly_adj: list of adjacency matrices with polygonal numerical range
    @param pt: used to allow for parallel computing with six or more vertices
    """
    global graph_num
    global title
    title = ""
    g = nx.DiGraph(adj)
    n = len(adj)

    # include graph in figure
    nx.draw_shell(g, with_labels=True, ax=plt.subplot(221))

    # form Laplacian
    diag = np.array([np.sum(adj[i, :]) for i in range(n)])
    lap = np.diag(diag) - adj

    # test for normality
    if np.linalg.norm(lap.transpose() @ lap - lap @ lap.transpose(), 'fro') == 0:  # normality test
        title = title + ' N_L '
    plt.subplot(221)

    # include graph number in figure
    if pt != -1:
        plt.title(title + "Graph Number %d.%d" % (pt, graph_num))
    else:
        plt.title(title + "Graph Number %d" % graph_num)
    plt.subplot(222)

    # include Laplacian in figure
    plt.matshow(lap, fignum=False, cmap='winter')
    for k in range(n):
        for j in range(n):
            c = lap[j, k]
            plt.text(k, j, str(c), va='center', ha='center', color='white')

    # find unrestricted numerical ranges and include both restricted and unrestricted in figure
    nr_restricted, eig_restricted = restricted_nr(lap)
    fig = plt.gcf()

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
        outfile = open('matrices/polygon_adjs/{}/{}_adjs.txt'.format(n, name), 'w+')
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

    line_list = infile.readlines()
    a = np.zeros((n, n))
    mat_line_count = 0
    count_vect = np.zeros(3)  # counts the number of singletons, lines, and polygons in the 0th, 1st, and 2nd positions
    singleton_adj = []
    line_adj = []
    poly_adj = []
    for row in line_list:
        if graph_num % 1000 == 0:
            print('Examined graph {}.{}'.format(pt, graph_num))
            print('We have {} singletons, {} lines, and {} polygons'.format(count_vect[0], count_vect[1],
                                                                            count_vect[2]))
            print('Total time elapsed ---- %.1f hours' % ((time.time() - start_time) / 3600))
        row = row.split(" ")
        if len(row) < n:  # we've read a full adjacency matrix
            determine_polygon(a, count_vect, singleton_adj, line_adj, poly_adj, pt)
            a = np.zeros((n, n))  # reset adjacency matrix
            mat_line_count = 0  # reset line count
        else:  # not the end of a matrix in infile, read the next line
            a[mat_line_count, :] = [eval(row[j]) for j in range(n)]
            mat_line_count = mat_line_count + 1

    write_out((singleton_adj, line_adj, poly_adj), n)

    report_file = open('reports/%d_report.txt' % n, 'w+')
    report_file.write('%d singletons \n%d lines \n%d polygons' % (count_vect[0], count_vect[1], count_vect[2]))


def main():
    find_poly_graphs(5)


if __name__ == '__main__':
    main()