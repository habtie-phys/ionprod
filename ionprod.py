import numpy as np
from scipy.io import loadmat
#
class ionprod:
    def __init__(self, rho, H, E):
        """"
        inputs:
            rho: is the atmospheric mass density (g/cm^3), (nZ,)
            H: scale height (cm), (nZ,)
            E: electron energies (keV), (nE+1,)
        """
        nE = E.size-1
        dE = np.diff(E)
        Ec = E[0:nE]+0.5*dE
        # 
        self.rho = rho
        self.H = H
        self.Ec = Ec
        self.dE = dE
    #
    def gety(self):
        y1 = 2/self.Ec
        y2 = (self.rho*self.H/6e-6)**0.7
        y = y2[:,np.newaxis]@y1[np.newaxis,:]
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
        f1 = (c1*np.power(y,c2))*np.exp(-c3*np.power(y,c4))
        f2 = (c5*np.power(y,c6))*np.exp(-c7*np.power(y,c8))
        f = f1+f2
        return f
    #
    def getA(self):
        """"
        gives matrix of energy deposition rate (keV/cm) (nZ, nE) 
        """
        f = self.getf()
        A = f*self.Ec*self.dE/self.H[:,np.newaxis]/0.035
        return A
    def getq(self, I):
        A = self.getA()
        q = A@I[:,np.newaxis]
        return q[:,0]
    
    def getVER(self, natmos, I):
        nN2 =  natmos['N2']
        nO2 =  natmos['O2']
        nO = natmos['O']
        # 
        V4728 = 0.628*nN2/(nN2+0.7*nO2+0.4*nO) # photon/keV
        B = 0.035*self.getA()*V4728[:,np.newaxis]
        ver = B@I[:, np.newaxis]
        return ver[:,0]

#
def maxw(Q0, E0, E):
    nE = E.size-1
    dE = np.diff(E)
    Ec = E[0:nE]+0.5*dE
    I = Q0*Ec*np.exp(-Ec/E0)/(2.0*E0**3)
    return I
def monoFlux(Q0, E0, E):
    nE = E.size-1
    dE = np.diff(E)
    Ec = E[0:nE]+0.5*dE
    # 
    indE = np.argmin(np.abs(Ec-E0))
    I = np.zeros(Ec.shape)
    I[indE] = Q0/(Ec[indE]*dE[indE])
    return I
