#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 22 05:41:28 2025
VITE-clustering demo for dist= 1.1, LiH molecula
@author: Mengzhen
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy as sp
from qiskit.circuit import QuantumCircuit
from molecular.molecular import LiH
from molecular.vqe import makewave0, dtheta, Lv, L
from molecular.utils import normalization, sf
from molecular.clustering import hopkins
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score

#%% demo for dist= 1.1, LiH molecula
n = 12  # number of qubit
roop = 100
for d in range(1):
    # dist = np.round(0.5 + 0.3 * d, 2)
    dist = 1.1
    Hp, Pauli_Op = LiH(dist)
    u, w = sp.sparse.linalg.eigs(Hp, k=15, which='SR')
    uu = sorted(set(np.around(u.real, 5)))[0:5]
    H2 = Hp.T.conjugate() * Hp
    wavefunction = QuantumCircuit(n)
    for i in range(n):
        wavefunction.h(i)
    wavefunction = makewave0(wavefunction, [], 2)
    N = wavefunction.num_parameters
    Estop=[]
    Eroop=[]
    E=[]
    dt=2e-1
    e = (np.around(min(uu)-0.2*(max(uu)-min(uu)),3)).real
    le = 500
    Ustep = []
    I = []
    de = (np.around((max(uu)-min(uu))/(20*5),5)).real
    xo=np.random.rand(N)*0.1
    for i in range(le):
        E.append(e)
        # H=(Hp-e*np.eye(2**n)).dot(Hp-e*np.eye(2**n))
        H = (Hp-e*np.eye(2**n)).T.conjugate().dot(Hp-e*np.eye(2**n))
        x0=xo.copy()
        # if use prior
        # ind = np.argmin(abs(x0_prior[:,0]-e)) # for prior
        # x0 = x0_prior[ind,-N:]           # for prior
        U=[L(x0,wavefunction)]
        ###########
        for j in range(roop):
            dx0=dtheta(x0, wavefunction, H)
            x0=(dx0)*dt+x0
            # print(x0)
            U.append(L(x0,wavefunction))
            F = Lv(x0,wavefunction)
            print(e,'\n Energy:',U[-1],
                  '\n fid:',np.round(max(sf(F,u,w))/sum(sf(F,u,w)),3),',',np.round(max(sf(F,u,w)),2),
                  '\n no:',np.argsort(sf(F,u,w))[::-1])
            if abs(U[-1]-U[-2])<1e-4:
                I.append(e)
                Ustep.append(j)
                Estop.append(x0)
                Eroop.append(x0)
                break
            # if j==4:
            #     E5.append(x0)
            # if j==9:
            #     E10.append(x0)
            if j==roop-1:
                Ustep.append(roop)
                Eroop.append(x0)
            # v=Lv(x0, wavefunction)
        e+=de
        if E[-1]>(max(uu)+0.2*(max(uu)-min(uu))).real:
            break 
    # I = pd.DataFrame(I)
    # Ustep = pd.DataFrame(Ustep)
    E = np.around(E,5)
    Eroop = pd.DataFrame(Eroop)
    # save
    # I.to_excel('00no_'+'I'+'_hd'+str(dist)+'.xlsx')
    # Ustep.to_excel('00no_'+'Ustep'+'_hd'+str(dist)+'.xlsx')
    # Eroop.to_excel('00no_'+'Eroop'+'_hd'+str(dist)+'.xlsx')



#%% Hopkins test
# read
# X_H = pd.read_excel('00no_Eroop_hd1.1.xlsx').iloc[:,-N:]
X_H = Eroop.iloc[:,-N:]
hopkins(X_H) # closer to 1, better

#%% KMeans
# X_H = pd.read_excel('00no_Eroop_hd1.1.xlsx').iloc[:,-N:]
data = np.array(X_H)
data = normalization(data)
X_H = pd.DataFrame(data)
clusters = 5
model_km = KMeans(n_clusters=clusters, n_init=50)
model_km.fit(X_H)
centers_km = model_km.cluster_centers_
classes_km = model_km.labels_

s_km = silhouette_score(X_H, classes_km)
print(s_km)
#%% show result
df_km = X_H
df_km.insert(0, 'class', classes_km)
# df_km.insert(0, 'energy', UN)
df_km.insert(0, 's', E)
ind = []
for i in range(clusters):
    if np.size(classes_km[classes_km==i],0)<5:
        ind+=(list(np.where(classes_km==i)[0]))
df_km = df_km.drop(ind)
df_km = df_km.reset_index(drop=True)
# X_H = df_km
# E_km = [ele for idx, ele in enumerate(E_km) if idx not in ind]

nc={}
j=0
nc[df_km['class'][0]]=j
for i in range(len(df_km['class'])):
    if df_km['class'][i] in nc.keys():
        df_km['class'][i]=nc[df_km['class'][i]]
    else:
        j+=1
        nc[df_km['class'][i]]=j
        df_km['class'][i]=j
#%% plot       
# sns.boxplot(x='class',y='s',data=df_km)
# for i in range(np.size(u)):
#     plt.axhline(y=u[i], color='r', linestyle=':')
# first 3 eigenvalue
df_km_top3 = df_km[df_km['class'].isin([0, 1, 2])]
sns.boxplot(x='class', y='s', data=df_km_top3)
for i in range(3):
    plt.axhline(y=u[i], color='r', linestyle=':')

