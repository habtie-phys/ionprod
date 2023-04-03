import numpy as np
from scipy.io import loadmat
#
class ionprod:
    def __init__(self, rho, H, Ec, dE, I):
        self.rho = rho
        self.H = H
        self.Ec = Ec
        self.dE = dE
        self.I = I
    #
    def gety(self):
        E0 = self.Ec[:, np.newaxis]
        y = (2/E0)*(self.rho*self.H/6e-6)**(0.7)
        return y
    #
    def  getf(self):
        d = loadmat('pdata.mat')
        P = d['P']
        E0 = np.tile(self.Ec, (4, 1))
        a = np.arange(4).reshape(4, 1)
        ln = np.power(np.log(E0), a)
        #
        C = np.exp(P@ln)
        c1 = C[0,:]; c2 = C[1,:]
        c3 = C[2,:]; c4 = C[3,:]
        c5 = C[4,:]; c6 = C[5,:]
        c7 = C[6,:]; c8 = C[7,:]
        y = self.gety()
        y1 = y.T
        f1 = (c1*np.power(y1,c2))*np.exp(-c3*np.power(y1,c4))
        f2 = (c5*np.power(y1,c6))*np.exp(-c7*np.power(y1,c8))
        f = f1+f2
        return f
    #
    def getq(self):
        F = (self.I*self.Ec*self.dE)
        Q = F[:, np.newaxis]
        f = self.getf()
        q1 = f@Q
        q = q1[:,0]/(0.035*self.H)
        return q
#
def maxw(Q0, E0, E):
    nE = E.size-1
    dE = np.diff(E)
    Ec = E[0:nE]+0.5*dE
    I = Q0*Ec*np.exp(-Ec/E0)/(2.0*E0**3)
    return (Ec, dE, I)