# NR_AlgCon: Numerical Range and Algebraic Connectivity
#
# Author: Thomas R. Cameron
# Date: 6/1/2019
import numpy as np
from matplotlib import patches as mpatches
from matplotlib import pyplot as plt
from math import pi as pi
from math import sqrt as sqrt

###############################################
###             Numerical Range             ###
###############################################
#   Returns plot of the numerical range of a
#   matrix with its eigenvalues.
###############################################
def nr(a):
    """Returns plot of the numerical range of a matrix with its eigenvalues."""
    nv = 360
    m, n = a.shape
    if(m!=n):
        print('Warning: matrix is non-square')
        return
    else:
        e = np.linalg.eigvals(a)
        f = []
        for k in range(1,nv+1):
            z = np.exp(2*pi*1j*(k-1)/nv)
            a1 = z*a
            a2 = (a1 + np.transpose(np.conjugate(a1)))/2
            w, v = np.linalg.eig(a2)
            ind = np.argsort(w)
            w = w[ind]
            v = v[:,ind]
            v = v[:,n-1]
            f.append(np.dot(np.conjugate(v),np.dot(a,v))/np.dot(np.conjugate(v),v))
        f.append(f[0])
        fig = plt.figure()
        plt.plot(np.real(f),np.imag(f),'b',figure=fig)
        plt.plot(np.real(e),np.imag(e),'r*',figure=fig)
        return fig
        
###############################################
###             algCon                      ###
###############################################
#   Compute the algebraic connectivity.
###############################################
def algCon(l,q):
    """Compute the algebraic connectivity."""
    e = np.linalg.eigvalsh(0.5*np.dot(np.transpose(q), np.dot(l+np.transpose(l),q)))
    return np.amin(e)
    
###############################################
###             domNR                       ###
###############################################
#   Creates dominance graph of size n and plots
#   its numerical range along with information
#   regarding the foci and major and minor axis.
###############################################
def dom_nr(n):
    """Creates dominance graph of size n and plots its numerical range along with informationregarding the foci and major and minor axis."""
    # perfect dominance graph on n vertices
    l=np.zeros((n,n))
    for i in range(n):
        l[i,i] = n-(i+1)
        for j in range(i+1,n):
            l[i,j] = -1
    # orthonormal matrix q
    q = np.zeros((n,n-1))
    for j in range(n-1):
        q[0:j+1,j] = 1
        q[j+1,j] = -(j+1)
        q[:,j] = q[:,j]/np.linalg.norm(q[:,j])
    # projection transformation
    a = np.dot(np.transpose(q), np.dot(l,q))
    print(np.dot(np.transpose(a),a))
    # plt of numerical range and eigenvalues
    fig = nr(a)
    # algebraic connectivity
    alpha = algCon(l,q)
    print(alpha)
    # semi-major axis
    c = n/2.0
    smajor = c - alpha
    # semi-minor axis
    e = np.linalg.eigvals(0.5*np.dot(np.transpose(q), np.dot(l-np.transpose(l),q)))
    sminor = max(np.imag(e))
    # foci
    x = sqrt(smajor**2-sminor**2)
    f1 = c - x
    f2 = c + x
    # plot minor axes
    sline = np.linspace(0,sminor)
    zline = [c for k in range(len(sline))]
    plt.plot(zline,sline,'g',label="sminor = %.2f"%sminor,alpha=0.5,figure=fig)
    # plot major axes
    xline = np.linspace(f1,n/2.0)
    m = sminor/x
    dx = xline[1]-xline[0]
    yline = [m*dx*k for k in range(len(xline))]
    plt.plot(xline,yline,'r',label="smajor = %.2f"%smajor,alpha=0.5,figure=fig)
    # plot foci
    plt.plot(f1,0,'b*',label="foci1 = %.5f"%f1,alpha=0.5,figure=fig)
    plt.plot(f2,0,'b*',label="foci2 = %.5f"%f2,alpha=0.5,figure=fig)
    # plt legends
    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
               ncol=2, mode="expand", borderaxespad=0.)
    plt.show()
    
###############################################
###             specDomNR                   ###
###############################################
def specDomNR():
    # graph Laplacian
    l = np.array([[1.,0,-1],[-1,1.,0],[-1,0,1.]])
    # orthonormal matrix q
    n = 3
    q = np.zeros((n,n-1))
    for j in range(n-1):
        q[0:j+1,j] = 1
        q[j+1,j] = -(j+1)
        q[:,j] = q[:,j]/np.linalg.norm(q[:,j])
    # projection transformation
    a = np.dot(np.transpose(q), np.dot(l,q))
    print(np.dot(np.transpose(a),a))
    # plt of numerical range and eigenvalues
    fig = nr(a)
    # algebraic connectivity
    alpha = algCon(l,q)
    print(alpha)
    # semi-major axis
    c = n/2.0
    smajor = c - alpha
    # semi-minor axis
    e = np.linalg.eigvals(0.5*np.dot(np.transpose(q), np.dot(l-np.transpose(l),q)))
    sminor = max(np.imag(e))
    # foci
    x = sqrt(smajor**2-sminor**2)
    f1 = c - x
    f2 = c + x
    # plot minor axes
    sline = np.linspace(0,sminor)
    zline = [c for k in range(len(sline))]
    plt.plot(zline,sline,'g',label="sminor = %.2f"%sminor,alpha=0.5,figure=fig)
    # plot major axes
    xline = np.linspace(f1,n/2.0)
    m = sminor/x
    dx = xline[1]-xline[0]
    yline = [m*dx*k for k in range(len(xline))]
    plt.plot(xline,yline,'r',label="smajor = %.2f"%smajor,alpha=0.5,figure=fig)
    # plot foci
    plt.plot(f1,0,'b*',label="foci1 = %.5f"%f1,alpha=0.5,figure=fig)
    plt.plot(f2,0,'b*',label="foci2 = %.5f"%f2,alpha=0.5,figure=fig)
    # plt legends
    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
               ncol=2, mode="expand", borderaxespad=0.)
    plt.show()
    
###############################################
###             Imploding Start NR (ISNR)   ###
###############################################
def ISNR():
    # graph Laplacian
    l = np.array([[0.,0,0,0],[-1,1.,0,0],[-1,0,1.,0],[-1,0,0,1.]])
    # orthonormal matrix q
    n = 4
    q = np.zeros((n,n-1))
    for j in range(n-1):
        q[0:j+1,j] = 1
        q[j+1,j] = -(j+1)
        q[:,j] = q[:,j]/np.linalg.norm(q[:,j])
    # projection transformation
    a = np.dot(np.transpose(q), np.dot(l,q))
    print(a)
    # plt of numerical range and eigenvalues
    fig = nr(a)
    plt.show()

# Q NR of Complete Dominance
#dom_nr(3)
# Q NR of non-complete dominance with same eigenvalues
#specDomNR()
ISNR()