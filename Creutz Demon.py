import numpy as np
from random import random
import pylab as plt
import time

tic=time.time()
N=1000
e=np.ones(N)
w=np.ones(N)
s=np.zeros(N)
m=np.zeros(N)
pi=3.14159
m=20
p1=np.linspace(0,0.95,m)
p2=np.linspace(0,1,m)
ustore=np.zeros([m,m])
lstore=np.zeros([m,m])
estore=np.zeros([m,m])
l=1
while 100*np.abs(l-0.8)/0.8>3:
#while 100*np.abs(l)>5:
    for ii in range(1):
        p1[0]=0.5
        p2[0]=0.95
        for jj in range(1):
            Lstore=np.zeros(N)
            Ustore=np.zeros(N)
            e=np.ones(N)
            for i in range(N):
                e[i]=random()*0.01+p1[ii]
                w[i]=random()*2*pi    
                if random()<p2[jj]:
                    s[i]=1
                else:
                    s[i]=-1
            for i in range(N):
                Lstore[i]=s[i]*np.sqrt(1-e[i]**2)
                for j in range(N):
                    if j!=i:
                        Ustore[i]+=np.log(e[i]**2+e[j]**2-2*e[i]*e[j]*np.cos(w[j]-w[i]))/(2*pi)-4*np.log(2)/pi
            U=sum(Ustore)/2
            L=sum(Lstore)
            u=U/N**2
            l=L/N
            print(l)
            print(u)
            estore[ii,jj]=np.average(e)
            ustore[ii,jj]=u
            lstore[ii,jj]=l
winitial=np.zeros(N)
winitial[:]=w[:]
total=2*10**4
Dlim=2*10**2
tic=time.time()
Lit=np.zeros(total)
Uit=np.zeros(total)
e1s=[]
e2s=[]
D=0
Ukeep=np.zeros(N)
eaverage=np.zeros(total)
Dit=np.zeros(total)
count1=0
count2=0
evector=np.zeros(total)
for it in range(total):

    x=0
    y=0
    for i in range(N):
        x+=e[i]*np.cos(w[i])
        y+=e[i]*np.sin(w[i])
    evector[it]=np.sqrt((x/N)**2+(y/N)**2)

    Lstore=np.zeros(N)
    eaverage[it]=np.average(e)
    #Ukeep[:]=Ustore[:]
    Ui=sum(Ustore)/2

    test=0
    e1n=2
    

    while test==0:
        k1=0
        k2=0
        while k1==k2:
            k1=int(np.floor(random()*N))
            k2=int(np.floor(random()*N))

        s2=s[k2]
        s1=s[k1]
        e1=e[k1]
        e2=e[k2]
        ##### E1N IS THE NEW ECCENTRICITY FOR ORBIT 1 ##########
        
        ###########################################################
        ##########   POSITIVE KICKS ##############################
        ##########################################################
        
        e1n=np.mod(e1+random()*0.1,1)
        
        ############################################################
        ########### RANDOM KICKS ###################################
        ###########################################################
        
        #e1n=random()
        
        
        ##############################################################
        ################ POS NEG KICKS ##############################
        #############################################################
        
#        e1n=e1+random()*0.8-0.8/2
#        f=1
#        while e1n>1 or e1n<0:
#            f+=1
#            e1n=e1+random()*(0.8/f)-0.8/(2*f)
        
        ##############################################################
        ############### END OF KICK TYPES ############################
        
        lol=(-(s1/s2)*(np.sqrt(1-e1n**2)-np.sqrt(1-e1**2))+np.sqrt(1-e2**2))
        if (lol<=1 and lol>=0):
            e2n=np.sqrt(1-lol**2)
            s2n=s2
            s1n=s1
            test=1
#        lol1=(-(s1/s2)*(-np.sqrt(1-e1n**2)-np.sqrt(1-e1**2))+np.sqrt(1-e2**2))
#        if (np.abs(lol)<=1) and (e1n>=0 and e1n<=1):
#            e2n=np.sqrt(1-lol**2)
#            test=1
#            s2n=s2
#            s1n=s1
#            if lol<0:
#                s2n=-s2
#        elif (np.abs(lol1)<=1) and (e1n>=0 and e1n<=1):
#            e2n=np.sqrt(1-lol1**2)
#            s1n=-s1
#            s2n=s2
#            test=1
#            if lol1<0:
#                s2n=-s2
                     
        
    w1=w[k1]
    w2=w[k2]
    e[k2]=e2n
    e[k1]=e1n
    w[k2]=random()*2*pi
    w[k1]=random()*2*pi
    s[k2]=s2n
    s[k1]=s1n

    for i in range(N):
        Ukeep[i]=Ustore[i]
        Lstore[i]=s[i]*np.sqrt(1-e[i]**2)
        if i!=k1 and i!=k2:
            Ustore[i]=Ustore[i]-np.log(e[i]**2+e1**2-2*e[i]*e1*np.cos(w1-w[i]))/(2*pi)\
                      -np.log(e[i]**2+e2**2-2*e[i]*e2*np.cos(w2-w[i]))/(2*pi)\
                      +np.log(e[i]**2+e1n**2-2*e[i]*e1n*np.cos(w[k1]-w[i]))/(2*pi)\
                      +np.log(e[i]**2+e2n**2-2*e[i]*e2n*np.cos(w[k2]-w[i]))/(2*pi)
        else:
            Ustore[i]=0
            for j in range(N):
                if j!=i:
                    Ustore[i]+=np.log(e[i]**2+e[j]**2-2*e[i]*e[j]*np.cos(w[j]-w[i]))/(2*pi)-4*np.log(2)/pi

    Un=sum(Ustore)/2
    deltaE=Un-Ui

    if deltaE<=D and D-deltaE<=Dlim:
        D=D-deltaE
        Uit[it]=Un
        count1+=1
        l=sum(Lstore)/N
        e1s.append(e[k1]-e1)
        e2s.append(e[k2]-e2)
    else:
        e[k2]=e2
        e[k1]=e1
        w[k1]=w1
        w[k2]=w2
        s[k2]=s2
        s[k1]=s1
        for i in range(N):
            Ustore[i]=Ukeep[i]
        Uit[it]=Ui
        count2+=1

    Dit[it]=D
    Lit[it]=l
    
print(np.average(eaverage[100000:total]))
toc=(time.time()-tic)/60
print(toc)
