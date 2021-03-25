import numpy as np
from random import random
import pylab as plt
import time

tic=time.time()
N=10
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
#while 100*np.abs(l-0.8)/0.8>3:
while 100*np.abs(l)>5:
    for ii in range(1):
        p1[0]=0.5
        p2[0]=0.5
        for jj in range(1):
            Lstore=np.zeros(N)
            Ustore=np.zeros(N)
            e=np.ones(N)
            for i in range(N):
                e[i]=random()
                #e[i]=random()*0.01+p1[ii]
                w[i]=random()*2*pi    
                if random()<p2[jj]:
                    s[i]=1
                else:
                    s[i]=-1
            for i in range(N):
                Lstore[i]=s[i]*np.sqrt(1-e[i]**2)
                for j in range(N):
                    if j!=i:
                        Ustore[i]+=np.log(e[i]**2+e[j]**2-2*e[i]*e[j]*np.cos(w[j]-w[i])+0.01)/(2*pi)-4*np.log(2)/pi
            U=sum(Ustore)/2
            L=sum(Lstore)
            u=U/N**2
            l=L/N
            print(l)
            print(u)
            estore[ii,jj]=np.average(e)
            ustore[ii,jj]=u
            lstore[ii,jj]=l

umin=-0.8
E=np.arange(umin,-0.38+0.01,0.01)
E[:]=np.round(E[:]*10**2)/10**2
L=np.arange(-1,1+0.01,0.01)
L[:]=np.round(L[:]*10**2)/10**2
H=np.zeros([len(E),len(L)])
g=np.ones([len(E),len(L)])
Ukeep=np.zeros(N)
f=2.71828
tic=time.time()
for it in range(10**6):
    uold=sum(Ustore)/(2*N**2)
    lold=l
    k1=0
    k2=0
    while k1==k2:
        k1=int(np.floor(random()*N))
        k2=int(np.floor(random()*N))
    e1=e[k1]
    e2=e[k2]
    w1=w[k1]
    w2=w[k2]
    s1=s[k1]
    s2=s[k2]
    l1=s[k1]*np.sqrt(1-e[k1]**2)
    l2=s[k2]*np.sqrt(1-e[k2]**2)
    
    e1n=random()
    e2n=random()
    w1n=random()*2*pi
    w2n=random()*2*pi
    s1n=np.sign(random()-0.5)
    s2n=np.sign(random()-0.5)
    
    e[k1]=e1n
    e[k2]=e2n
    w[k1]=w1n
    w[k2]=w2n
    s[k1]=s1n
    s[k2]=s2n
    
    
    Lstore[k1]=s[k1]*np.sqrt(1-e[k1]**2)
    Lstore[k2]=s[k2]*np.sqrt(1-e[k2]**2)
    lnew=sum(Lstore)/N
    
    for i in range(N):
        Ukeep[i]=Ustore[i]
        if i!=k1 and i!=k2:
            Ustore[i]=Ustore[i]-np.log(e[i]**2+e1**2-2*e[i]*e1*np.cos(w1-w[i])+0.01)/(2*pi)\
                      -np.log(e[i]**2+e2**2-2*e[i]*e2*np.cos(w2-w[i])+0.01)/(2*pi)\
                      +np.log(e[i]**2+e1n**2-2*e[i]*e1n*np.cos(w[k1]-w[i])+0.01)/(2*pi)\
                      +np.log(e[i]**2+e2n**2-2*e[i]*e2n*np.cos(w[k2]-w[i])+0.01)/(2*pi)
        else:
            Ustore[i]=0
            for j in range(N):
                if j!=i:
                    Ustore[i]+=np.log(e[i]**2+e[j]**2-2*e[i]*e[j]*np.cos(w[j]-w[i])+0.01)/(2*pi)-4*np.log(2)/pi
    unew=sum(Ustore)/(2*N**2)
    
    oldi=int(np.floor((uold-umin)*10**2))
    newi=int(np.floor((unew-umin)*10**2))
    oldj=int(np.floor((lold+1)*10**2))
    newj=int(np.floor((lnew+1)*10**2))
    
    ratio=np.exp(g[oldi,oldj]-g[newi,newj])
    p=random()
    if p<=ratio:
        g[newi,newj]=g[newi,newj]+np.log(f)
        H[newi,newj]=H[newi,newj]+1
        l=lnew
    else:
        e[k2]=e2
        e[k1]=e1
        w[k1]=w1
        w[k2]=w2
        s[k2]=s2
        s[k1]=s1
        Lstore[k1]=l1
        Lstore[k2]=l2
        l=lold
        H[oldi,oldj]=H[oldi,oldj]+1
        g[oldi,oldj]=g[oldi,oldj]+np.log(f)
        for i in range(N):
            Ustore[i]=Ukeep[i]

toc=time.time()-tic
print(toc)
X,Y=np.meshgrid(L,E)
fig = plt.figure()
ax = plt.axes(projection='3d')

ax.plot_surface(X, Y, H)
ax.set_title('Surface plot')
plt.show()
