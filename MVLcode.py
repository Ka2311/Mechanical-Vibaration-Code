import numpy as np
import matplotlib.pyplot as plt

global x0,M,K,v0,p

i=int(input('press 1 for free_undamped\npress 2 for free_damped\npress 3 for forced_damped or undamped\n'))

M  = float(input('enter the mass->'))               # mass of the body
K  = float(input('enter spring constant->'))        # stiffness 
x0 = float(input('enter intial displacement->'))    # initial displacement
v0 = float(input('enter intial velocity->'))        # initial velocity
p=np.sqrt(K/M)   

t = np.arange(0, 2, 0.01)

if(i==2 or i==3):
 C  = float(input('enter value of damping constant(0 for undamped and non zero for damped)'))
 global zeeta                                       #damping ratio
 zeeta=C / (2*np.sqrt(K*M))
 
 if(zeeta>1):  
   global c1,c2,c3,c4,a,b,a1,b1                 #constants required during various calculations
   a=1/2*(x0 + (v0+x0*p*zeeta) / (p* np.sqrt(zeeta*2-1)) )
   b=x0-a
   c1 = a*((-zeeta + np.sqrt(zeeta*2-1) )* p)
   c2 = b*((-zeeta - np.sqrt(zeeta*2-1))* p)
   c3=a* ((-zeeta + np.sqrt(zeeta*2-1) )* p)**2
   c4=b* ((-zeeta - np.sqrt(zeeta*2-1))* p)**2
   a1=(x0-v0)/(1+p)
   b1=x0-a1
   print(c1)

 if(zeeta<1):
   global wd,a2,b2,d1,d2                          #constants required during various calculations'''
   wd=np.sqrt(1-zeeta*2) * p
   a2=(v0+x0*zeeta*p)/wd
   b2=x0
   d1=a2* np.sin(wd*t) + b2* np.cos(wd*t)
   d2=a2*np.cos(wd*t)-b2*np.sin(wd*t)



def free_undamped_vib():
  # Calculate displacement as a function of time
  y1 = x0 * np.cos(p*t) + v0/p * np.sin(p*t)
  # calculate velocity as function of time
  y2= v0/p * np.cos(p*t) -x0 * np.sin(p*t)
   # calculate acceleration as function of time
  y3=-v0/p *np.sin(p*t)- x0 * np.cos(p*t)

  # Plotting both the curves simultaneously
  plt.plot(t,y1,color='r',label='displacement')
  plt.plot(t,y2,color='g',label='velocity')
  plt.plot(t,y3,color='b',label='acceleration')
    
# Naming x-axis, y-axis and labelling the graph
  plt.xlabel('Time (s)')
  plt.ylabel('Displacement (m)')
  plt.plot(t,0*t,color='black')
  plt.title('Undamped Vibration')
  plt.legend()
  plt.show()



def free_damped_vib():
  
  if(zeeta>1):
    
    # Calculate displacement as a function of time
    y1=a* np.exp((-zeeta + np.sqrt(zeeta*2-1))* p* t) +b* np.exp((-zeeta - np.sqrt(zeeta*2-1))* p* t)
    # calculate velocity as function of time
    
    y2=c1*  np.exp((-zeeta + np.sqrt(zeeta*2-1))* p* t) +c2* np.exp((-zeeta - np.sqrt(zeeta*2-1))* p* t)
    # calculate accelaration as function of time
    
    y3=c3* np.exp((-zeeta + np.sqrt(zeeta*2-1))* p* t) +c4*np.exp((-zeeta - np.sqrt(zeeta*2-1))* p* t)
    # Plotting both the curves simultaneously
    plt.plot(t,y1,color='r',label='displacement')
    plt.plot(t,y2,color='g',label='velocity')
    plt.plot(t,y3,color='b',label='acceleration')
    
# Naming the x-axis, y-axis and labelling the graph
    plt.xlabel('Time (s)')
    plt.ylabel('Displacement (m)')
    plt.title('Overdamped')
    plt.plot(t,0*t,color='black')
    plt.legend()
    plt.show()

  elif(zeeta==1):
    a1=(x0-v0)/(1+p)
    b1=x0-a1
    # Calculate the displacement as a function of time
    y1=( a1 + b1* t)* np.exp(- (p* t)) 
    # calculate the velocity as function of time
    y2=(b1-b1*t*p-a1*p)*np.exp(- (p* t)) 
     # calculate the accelaration as function of time
    y3=(a1*(p**2)+b1*t*(p**2)-2*b1*p)*np.exp(- (p* t))
    # Plotting both the curves simultaneously
    plt.plot(t,y1,color='r',label='displacement')
    plt.plot(t,y2,color='g',label='velocity')
    plt.plot(t,y3,color='b',label='acceleration')
    

    plt.xlabel('Time (s)')
    plt.ylabel('Displacement (m)')
    plt.title('Critically damped')
    plt.plot(t,0*t,color='black')
    plt.legend()
    plt.show()

  else:
    wd=np.sqrt(1-zeeta*2) * p
    a2=(v0+x0*zeeta*p)/wd
    b2=x0
    # Calculate the displacement as a function of time
   
    y1=np.exp(-zeeta*p*t)* d1
    # calculate the velocity as function of time
    
    y2=np.exp(-zeeta*p*t) * (wd*d2 -d1*-zeeta*p)
    # calculate the accelaration as function of time
    y3=d1 *np.exp(-zeeta*p*t)*((zeeta*p)**2 -wd**2)-2*(wd*p*zeeta)*d2*np.exp(-zeeta*p*t)
    # Plotting both the curves simultaneously
    plt.plot(t,y1,color='r',label='displacement')
    plt.plot(t,0*t,color='black')
    plt.plot(t,y2,color='g',label='velocity')
    plt.plot(t,y3,color='b',label='acceleration')
    

    plt.xlabel('Time (s)')
    plt.ylabel('Displacement (m)')
    plt.title('Underdamped')
    plt.legend()
    plt.show()

def forced_damped_undamped_vib():
  
  F0 = float(input('Enter the magnitude of force'))
  W=float(input('enter harmonic frequency  applied->'))
  r = W/p
  global angle_fhi,A

  angle_fhi=np.arctan((2*zeeta*r)/(1-r**2))
  A=F0/(K* np.sqrt((1-r**2)**2 +(2*zeeta*r)**2))
  print(A)
  
  if(p==W and C==0):
    print('resonance condition')
    

  elif(zeeta>1):
    '''Adding complementary function with particular integral to find displacement function'''
  
    # Calculate the displacement as a function of time
   
    p_i= A * np.sin(W*t-angle_fhi)
    y1=p_i + a* np.exp((-zeeta + np.sqrt(zeeta*2-1))* p* t) +b* np.exp((-zeeta - np.sqrt(zeeta*2-1))* p* t)
     # Calculate the velocity as a function of time
    y2=A* np.cos(W*t-angle_fhi)*W + c1* np.exp((-zeeta + np.sqrt(zeeta*2-1))* p* t) +c2* np.exp((-zeeta - np.sqrt(zeeta*2-1))* p* t)
     # Calculate the acceleration as a function of time
    y3=A* np.cos(W*t-angle_fhi)*W*W + c3* np.exp((-zeeta + np.sqrt(zeeta*2-1))* p* t) +c4*np.exp((-zeeta - np.sqrt(zeeta*2-1))* p* t) 
    # Plotting both the curves simultaneously
    
    plt.plot(t,y1,color='r',label='displacement')
    plt.plot(t,y2,color='g',label='velocity')
    plt.plot(t,0*t,color='black')
    plt.plot(t,y3,color='b',label='acceleration')
    plt.legend()
    plt.show()

  elif(zeeta==1):
    
    a1=(x0-v0)/(1+p)
    b1=x0-a1
    # Calculate the displacement as a function of time
    y1= A * np.sin(W*t-angle_fhi) + ( a1 + b1* t)* np.exp(- (p* t)) 
    # calculate the velocity as function of time
    y2=A* np.cos(W*t-angle_fhi)*W + c1* np.exp((-zeeta + np.sqrt(zeeta*2-1))* p* t) + (b1-b1*t*p-a1*p)*np.exp(- (p* t)) 
     # calculate the accelaration as function of time
    y3=A* np.cos(W*t-angle_fhi)*W*W + c3* np.exp((-zeeta + np.sqrt(zeeta*2-1))* p* t) + (a1*(p**2)+b1*t*(p**2)-2*b1*p)*np.exp(- (p* t))
    
    plt.plot(t,y1,color='r',label='displacement')
    plt.plot(t,0*t,color='black')
    plt.plot(t,y2,color='g',label='velocity')
    plt.plot(t,y3,color='b',label='acceleration')
    plt.legend()
    plt.show()

  else:
    wd=np.sqrt(1-zeeta*2) * p
    
    a2=(v0+x0*zeeta*p)/wd
    b2=x0
    # Calculate the displacement as a function of time
    y1 = A * np.sin(W*t-angle_fhi) + np.exp(-zeeta*p*t)* d1
    y2=A* np.cos(W*t-angle_fhi)*W + np.exp(-zeeta*p*t) * (wd*d2 -d1*-zeeta*p)
    # calculate the accelaration as function of time
    y3=A* np.cos(W*t-angle_fhi)*W*W + -v0/p *np.sin(p*t)- x0 * np.cos(p*t)
    
    plt.plot(t,y1,color='r',label='displacement')
    plt.plot(t,0*t,color='black')
    plt.plot(t,y2,color='g',label='velocity')
    plt.plot(t,y3,color='b',label='acceleration')
    plt.legend()
    plt.show()



if(i==1):
  free_undamped_vib()
elif(i==2):
  free_damped_vib()
elif(i==3):
  forced_damped_undamped_vib()
  
