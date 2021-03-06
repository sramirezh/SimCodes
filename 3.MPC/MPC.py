#Make as c++
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt


def EKinetic (Vel,Npart):
    Ek=0
    for i in xrange(Npart):
        Ek+=Vel[i,0]**2+Vel[i,1]**2+Vel[i,2]**2
    return 0.5*Ek

def Initialise (Npart, KbT):
    print("\nCreating a new configuration!\n")
    Pos=np.random.rand(Npart,3)*L
    Vel=np.random.uniform(low=0, high=1.0, size=(Npart,3))
    
    #Imposing zero momentum
    Vcm=CenterMass(Vel)
    Vel=Vel-Vcm
    
    #Scaling Velocities
    scale=np.sqrt(3/2*Npart*KbT/EKinetic(Vel,Npart))
    Vel=scale*Vel
    
    print "Checking if the velocity definition and scaling is correct \n"
    print "The Temperature set temperature is: %f \n" %(EKinetic(Vel,Npart)*2/(3*Npart))
    print "Momentum in (x,y,z):(%f,%f,%f) \n" %(np.sum(Vel[:,0]),np.sum(Vel[:,1]),np.sum(Vel[:,2]))
    
    #Correction of the center of mass movement
    
#    for i in xrange(Npart):
#        Vel[i,:]=Vel[i,:]/(np.linalg.norm(Vel[i,:]))
        
    return Pos,Vel



    
def FreeStream(Vel,Pos,Deltat):
    for i in xrange(3):
        Pos[:,i]=Pos[:,i]+Vel[:,i]*Deltat
    return Pos

def StochasticRotation(Vel,Pos):
    
    #Grid Displacement
    Disp=RandDisplacement()
    for i in xrange(3):
        Pos[:,i]=GridShift(Pos[:,i],Disp[i],L)
        
    #Cell Division
    Head,List=CellDivision(Pos,Npart,L)
    
    #Stochastic Rotation Method
    for i in xrange(L):
        for j in xrange(L):
            for k in xrange(L):
                particles=CellParticles(i,j,k,Head,List) #Particles in the cell
                n=np.size(particles)
                if n<2: break
                Vcm=CenterMass(Vel[particles]) #Particle velocity in the cm reference
                #Parameters and creation of the rotation matrix
                phi,tetha=Random()
                R=RotationMatrix(phi,tetha,alpha)
                #Vel[particles]=Vel[particles]+np.transpose((R-I)*np.transpose(Velcm))
                for p in particles:
                    Vel[p]=Vcm+np.transpose(R*np.reshape((Vel[p]-Vcm),(3,1)))
                
                
    #Grid Displacement back        
    for i in xrange(3):
        Pos[:,i]=GridShift(Pos[:,i],-Disp[i],L)
    return Vel

#def CenterMass(Vel):
#    """
#    Computes the center of mass velocity and expresses the velocities in the cm Frame
#    """
#    n,m=np.shape(Vel)
#    Velcm=np.zeros((n,m))
#    Vcm=[np.sum(Vel[:,0])/n,np.sum(Vel[:,1])/n, np.sum(Vel[:,2])/n] 
#    for i in xrange(3):
#        Velcm[:,i]=Vel[:,i]-Vcm[i]
#    
#    return Velcm,Vcm

def CenterMass(Vel):
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
    
def CellDivision(Pos,Npart,L):
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
    
def CellParticles(Indx,Indy,Indz,Head,List):
    """rho=Npart/L**3.0
    returns the particles index for the cell described by the Input indexes
    Cell[Indx,Indy,Indz]
    """
    Indexes=[]
    part=Head[Indx,Indy,Indz]
    while part!=0:
        Indexes.append(part-1)
        part=List[part-1]
    Indexes=np.array(Indexes)
    
    return Indexes




# =============================================================================
# Input Parameters 
# =============================================================================
alpha=90
Deltat=0.5
Npart=100000
KbT=1/3
rho=np.float(10) 
L=int(np.floor((Npart/rho)**(1./3.))) #Be careful to check what is the average density.
Nx=L #Number of partitions in x,y,z direction.
a=L/Nx #Cell Size
print "The Average density of the system is: %f" %(Npart/L**3.0)
Nrun=10
Tsample=0.1
Tsample=Nrun*Tsample
m=1 #Particle mass



#==============================================================================
# Starting the MPC algorithm
#==============================================================================


#Initialization of the system
Pos,Vel=Initialise(Npart,KbT)
VelInitial=np.copy(Vel)







for t in xrange(Nrun):
    print t
    Pos=FreeStream(Vel,Pos,Deltat)
    Vel=StochasticRotation(Vel,Pos)
    
    if t%Tsample==0:
        print EKinetic(Vel,Npart)
#    if t%plotF==0:
#        fig1.canvas.draw()
#        plt.hist(Vel[:,0],bins='auto', normed=1)
#        plt.pause(0.05)
#        fig1.clf()
    


plt.close('all')
fig = plt.figure()
ax = fig.gca()
n, bins, rectangles = ax.hist(Vel[:,0], 50, normed=True)

v=np.linspace(np.min(Vel[:,0]),np.max(Vel[:,0]))
pv=(m/(2*np.pi*KbT))**(1/2)*np.exp(-m*v**2/(2*KbT))
plt.plot(v, pv)

plt.xlabel("$v_x$", fontsize=16)  
plt.ylabel("$P(v_x)$", fontsize=16)  

# Ensure that the axis ticks only show up on the bottom and left of the plot.  
# Ticks on the right and top of the plot are generally unnecessary chartjunk.  
ax.get_xaxis().tick_bottom()  
ax.get_yaxis().tick_left()  

# Remove the plot frame lines. They are unnecessary chartjunk. 
ax = plt.subplot(111)  
ax.spines["top"].set_visible(False)  
ax.spines["right"].set_visible(False)  

plt.savefig("Vdist.eps", bbox_inches="tight")  

#plt.hist(Vel[:,0], bins='auto', normed=1)
#plt.figure(2)
#plt.hist(VelInitial[:,0],bins='auto', normed=1)

#==============================================================================
# For testing Purposes
#==============================================================================

#To see if the Cell process is correct

#part=Head[0,0,0]
#lista=[]
#while part!=0:
#    lista.append(part)
#    part=List[part-1]
#    
#lista=np.array(lista)-1
#import matplotlib.pyplot as plt
#plt.close("all")
#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#ax.scatter(Pos[:,0],Pos[:,1],Pos[:,2])
#ax.scatter(Pos[lista,0],Pos[lista,1],Pos[lista,2],c='r', marker='o',s=40)


#To see if the random numbers are ok 

#Nsamples=10000
#T=np.zeros(Nsamples)
#P=np.zeros(Nsamples)
#V=np.zeros(Nsamples)
#RandDisp=np.zeros((Nsamples,3))
#for i in xrange(1000):
#   T[i],P[i]=Random()
#   V[i]=np.random.uniform(low=0, high=1.0)
#   RandDisp[i,:]=RandDisplacement()
#    
