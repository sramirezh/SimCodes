#Make as c++

import numpy as np


def Initialise (Npart):
    print("\nCreating a new configuration!\n")
    Pos=np.random.rand(Npart,3)*L
    Vel=np.random.rand(Npart,3)
    return Pos,Vel

def FreeStream(Vel,Pos,Npart):
    for i in xrange(3):
        Pos[:,i]=Pos[:,i]*Vel[:,i]
    return Pos

def StochasticRotation(Vel,Pos):
    Disp=RandDisplacement()
    for i in xrange(3):
        Pos[:,i]=GridShift(Pos[:,i],Disp[:,i],L)
    
    
    
    return Vel

def Vcm(Vel):
    """
    Computes the center of mass velocity
    """
    n,m=np.shape(Vel)
    Vcm=[np.sum(Vel[:,0])/n,np.sum(Vel[:,1])/n, np.sum(Vel[:,2])/n]
    
    return Vcm
    

def Random():
    """
    Generates a random pair of numbers, phi [0,2pi] and tetha [-1,1]
    """
    phi=np.random.rand()*2.0*np.pi
    tetha=(np.random.rand()-0.5)*2.0
    return phi, tetha

def RandDisplacement():
    """
    Generates a random displacement vector
    """
    V=[np.random.rand()-0.5,np.random.rand()-0.5,np.random.rand()-0.5]
    
    return V

def RotationMatrix(phi,tetha,alpha):
    alpha=alpha*np.pi/180
    Rx=(1-tetha**2)**(0.5)*np.cos(phi)
    Ry=(1-tetha**2)**(0.5)*np.sin(phi)
    Rz=tetha
    c=np.cos(alpha)
    s=np.sin(alpha)
    D=np.matrix([[Rx**2+(1-Rx**2)*c,Rx*Ry*(1-c)-Rz*s, Rx*Rz*(1-c)+Ry*s],\
                 [Rx*Ry*(1-c)+Rz*s, Ry**2+(1-Ry**2)*c , Ry*Rz*(1-c)-Rx*s],\
                 [Rx*Rz*(1-c)-Ry*s, Ry*Rz*(1-c)+Rx*s, Rz**2+(1-Rz**2)*c ]])
    return D



def GridShift(X,xshift,Lx):
    """
    Performs a shift of the particles positions in the direction given taking into account PBC
    X: is the vector with the positions in that direction
    xshift: is the shift 
    Lx: Simulation box parameters in that direction. 
    """
    X=X+xshift
    IndexUp=np.where(X>Lx) #Replace with a for in c++
    IndexLo=np.where(X<0)
    
    #Applying BC
    for index in IndexUp:
        X[index]=X[index]-Lx
    for index in IndexLo:
        X[index]=Lx+X[index]
    
    return X    
    
def CellParticles(Pos,Npart,L):
    """
    Performs Linked lists
    Output 
    Head has the the number of the first (Last) particle found in the Cell
    List= For every particle position in the array, the entry is the particle linked by the algorithm
    WARNING, in head and List, the particles have i+1 index, in order to avoid two zeros
    """
    Cell=np.zeros((Npart,3),dtype=int)
    List=np.zeros(Npart,dtype=int)
    Head=np.zeros((L,L,L),dtype=int)
    for i in xrange(Npart):
        Cell[i,:]=np.floor(Pos[i,:])
        List[i]=Head[Cell[i,0],Cell[i,1],Cell[i,2]]
        Head[Cell[i,0],Cell[i,1],Cell[i,2]]=i+1
    return Head, List
    


alpha=90
deltat=0.5
Npart=80
rho=np.float(10)
L=int(np.ceil((Npart/rho)**(1./3.)))
Nx=L #Number of partitions in x,y,z direction.
a=L/Nx #Cell Size

Nrun=100

Pos,Vel=Initialise(Npart)
rho=Npart/L**3.0

Head,List=CellParticles(Pos,Npart,L)
phi,tetha=Random()
D=RotationMatrix(phi,tetha, alpha)


#For testing Purposes