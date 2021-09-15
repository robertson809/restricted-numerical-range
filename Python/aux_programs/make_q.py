"""
Creates an 'q' matrix, that transports a vector from the basis R^n to the basis R^n / e, 
the all ones vector. 
"""

q = np.zeros((n,n-1))
for j in range(n-1):
    q[0:j+1,j] = 1
    q[j+1,j] = -(j+1)
    q[:,j] = q[:,j]/np.linalg.norm(q[:,j])