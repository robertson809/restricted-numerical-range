from polygonal import qnr, normality
from find_polygonal import nr
import networkx as nx
from matplotlib import pyplot as plt
import numpy as np

def panel_display(l, name):
    # restricted numerical range
    f, e = qnr(l)
    # check if graph is polygonal
    a = (l - np.diag(np.diagonal(l))) * -1

    n = a.shape[0]
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
    plt.show()
    fig.savefig(f"demo_figures/{name}_quadriptych", dpi=400)

    plt.clf()


def restricted_nr(lap, title=None, name = ""):
    """
    Plots both the restricted and unrestricted numerical ranges
    @param lap: lapalcian matrix
    @return: the numerical range eigenvalue pair for the restricted
    and unrestricted numerical ranges
    """
    nr_restricted, e_restricted = qnr(lap)
    nr_unrestricted, e_unrestricted = nr(lap)
    
    plt.title(title)
    plt.plot(np.real(nr_restricted), np.imag(nr_restricted),'#03a5e0', linewidth=2.5)
    plt.plot(np.real(e_restricted), np.imag(e_restricted), '*',
                linestyle='None', marker='*', color='#0337e0', markersize = 8)
    plt.fill(np.real(nr_restricted), np.imag(nr_restricted), '#03a5e0')  # fill
    fig = plt.gcf()
    plt.show()
    fig.savefig(f"demo_figures/{name}_rnr", dpi=400)