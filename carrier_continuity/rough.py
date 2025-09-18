import sympy as sm
import scipy.constants as con

V0,V1,V2=sm.symbols('V0 V1 V2')
Vprev0, Vprev1, Vprev2=sm.symbols('Vprev0:3')
n0,n1,n2=sm.symbols('n0 n1 n2')
nprev0, nprev1, nprev2=sm.symbols('nprev0:3')
p0,p1,p2=sm.symbols('p0:%d' % 3)
pprev0,pprev1,pprev2=sm.symbols('pprev0:%d'%3)
ni=sm.symbols('ni')
dt,dx=sm.symbols('dt dx')
G=sm.symbols('G')
q=sm.symbols('q')
k,T,un,up=sm.symbols('k T un,up')
def B(psi):
    return psi / (sm.exp(psi) - 1)

# Scharfetter-Gummel flux between node i and j (using node i on left, j on right)
def J(n_i, n_j, phi_i, phi_j,u):
    psi = (phi_j - phi_i)
    return (k*T*u/dx) * ( B(psi)*n_j - B(-psi)*n_i )
nF1=(J(n1,n2,V1,V2,un)-J(n0,n1,V0,V1,un))/q + G - n1*p1+ni**2
nF2=(J(nprev1,nprev2,Vprev1,Vprev2,un)-J(nprev0,nprev1,Vprev0,Vprev1,un))/q + G - nprev1*pprev1+ni**2
nF=(n1-nprev1)/dt - (nF1+nF2)/2
print("Residual nF =", nF)

print("\n--- Jacobian entries for nF ---")
print("∂nF/∂n0 =", sm.diff(nF, n0))
print("∂nF/∂n1 =", sm.diff(nF, n1))
print("∂nF/∂n2 =", sm.diff(nF, n2))
print("∂nF/∂p0 =", sm.diff(nF, p0))
print("∂nF/∂p1 =", sm.diff(nF, p1))
print("∂nF/∂p2 =", sm.diff(nF, p2))

pF1=-(J(p1,p2,V1,V2,up)-J(p0,p1,V0,V1,up))/q + G - n1*p1+ni**2
pF2=-(J(pprev1,pprev2,Vprev1,Vprev2,up)-J(pprev0,pprev1,Vprev0,Vprev1,up))/q + G - nprev1*pprev1+ni**2
pF=(p1-pprev1)/dt - (pF1+pF2)/2
print("\n\nResidual pF =", pF)

print("\n--- Jacobian entries for nF ---")
print("∂pF/∂n0 =", sm.diff(pF, n0))
print("∂pF/∂n1 =", sm.diff(pF, n1))
print("∂pF/∂n2 =", sm.diff(pF, n2))
print("∂pF/∂p0 =", sm.diff(pF, p0))
print("∂pF/∂p1 =", sm.diff(pF, p1))
print("∂pF/∂p2 =", sm.diff(pF, p2))
