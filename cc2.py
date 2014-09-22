#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      seba
#
# Created:     13-09-2014
# Copyright:   (c) seba 2014
# Licence:     <your licence>
#-------------------------------------------------------------------------------
from  numpy import  *
import numpy as np
from scipy import *
from matplotlib import *
from numpy.linalg import norm, solve



def crear_matriz():
    matriz1 =[]
    file = open('associations.dat', 'r')
    i=0
    for linea in file:
        if linea[-1] == '\n':
		linea = linea[:-1]
        linea = linea.split(':')
        matriz1.append(float(linea[2]))
        i=i+1
    matriz=np.array(matriz1)
    matriz=np.reshape(matriz,(20,40))
    return matriz

def pregunta1(matriz):
    matriz_transpuesta=matriz.transpose()
    matriz_a=np.dot(matriz,matriz_transpuesta)
    return matriz_a

def initial():
    x=np.ones((20))
    return x

def power(x,a):
    val_prop_ant=0
    val_prop=0.1
    while val_prop!= val_prop_ant:
        val_prop_ant=val_prop
        u =x/norm(x)
        x=np.dot(a,u)
        val_prop=float(np.dot(np.dot(u.T,a),u))
    u=x/norm(x)
    print("el valor propios es:",val_prop)
    print("el vector propio es:",u)
    return val_prop,u

def inverpower(x,a,s):
    val_prop_ant=0
    val_prop=0.1
    while val_prop!= val_prop_ant:
        val_prop_ant=val_prop
        u =x/norm(x)
        if(np.linalg.det(a -s*np.eye(*a.shape))==0):
            val_prop=np.dot(u.T,x)
            break
        x=solve(a -s*np.eye(*a.shape), u)
        val_prop=np.dot(u.T,x)
    u=x/norm(x)
    val_prop=(val_prop**(-1))+s
    print("el valor propios es:",val_prop)
    print("el vector propio es:",u)

def rayleigh(x,a):
    val_prop_ant=0
    val_prop=1
    while val_prop!= val_prop_ant:
        val_prop_ant=val_prop
        u =x/norm(x)
        val_prop=float(np.dot(np.dot(u.T,a),u))
        x = solve(a -val_prop*np.eye(*a.shape), u)
    u=x/norm(x)
    print("el valor propios es:",val_prop)
    print("el vector propio es:",u)

def unshiftedqr(A,k):
    m = A.shape[0]
    Q = eye(m)
    Qbar = Q.copy(); R = A.copy()
    for j in range(k):
        Q,R = linalg.qr( dot(R,Q) )	# QR factorization
        Qbar = dot(Qbar,Q)
    lam = diag( dot(R,Q) )	# Rayleigh quotient
    print("el valor propios es:",lam[0])
    print("el valor propios es:",lam[1])
    print("el vector propio es:",Qbar[0])
    print("el vector propio es:",Qbar[1])
    return lam,Qbar

def matriz5():
    matriz=np.random.rand(5,5)
    return matriz
def matriz10():
    matriz=np.random.rand(10,10)
    return matriz

if __name__ == '__main__':

    matriz=crear_matriz()
    matriz_a=pregunta1(matriz)
    t= initial()
    a = np.array([[1, 0],[0, 2]])
    x = np.array([[1.],[1.]])
    s,v=power(x,a)
    print("esto es v",v)
    print("algoritmo ray")
    s=np.dot(s,np.dot(v,v.T))
    print("esto es s",s)
    inverpower(x,a,s)
    print("primer ray")
    rayleigh(x,a)
    a = a - np.dot(s,np.dot(v,v.T))
    print("segundo ray")
    rayleigh(x,a)
    print("qr")
    unshiftedqr(a,6)

