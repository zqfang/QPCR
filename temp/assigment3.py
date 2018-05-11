# -*- coding: utf-8 -*-
"""
Created on Wed Sep  7 21:15:06 2016

@author: bioninja
"""

import matplotlib as mpl
import matplotlib.pyplot as plt

import pandas as pd
import numpy as np 
from sklearn.cross_validation import train_test_split
from sklearn.svm import SVC
from sklearn.preprocessing import  Normalizer, MaxAbsScaler, MinMaxScaler, StandardScaler

# laod data in 
data = pd.read_csv("./Datasets/parkinsons.data")



y = data['status']
X = data.drop(['name','status'], axis=1)

# first question
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=7)

svc = SVC()
svc.fit(X_train, y_train)
score = svc.score(X_test, y_test)
print score

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=7)
best_score = 0
#second question
for C in np.arange(0.05,2, 0.05):
    for gamma in np.arange(0.001, 0.1, 0.001):
        
        svc = SVC(C=C,gamma=gamma)
        svc.fit(X_train, y_train)
        score = svc.score(X_test, y_test)
        if score > best_score:
            best_score = score 
            print "C, gamma, score", C, gamma, score




#normalizer
norm = Normalizer()
norm.fit(X)
T = norm.transform(X)

X_train, X_test, y_train, y_test = train_test_split(T, y, test_size=0.3, random_state=7)

for C in np.arange(0.05,2, 0.05):
    for gamma in np.arange(0.001, 0.1, 0.001):
        
        svc = SVC(C=C,gamma=gamma)
        svc.fit(X_train, y_train)
        score = svc.score(X_test, y_test)
        if score > best_score:
            best_score = score 
            print "C, gamma, score", C, gamma, score

#maxabs
norm = MaxAbsScaler()
norm.fit(X)
T = norm.transform(X)

X_train, X_test, y_train, y_test = train_test_split(T, y, test_size=0.3, random_state=7)

        
for C in np.arange(0.05,2, 0.05):
    for gamma in np.arange(0.001, 0.1, 0.001):
        
        svc = SVC(C=C,gamma=gamma)
        svc.fit(X_train, y_train)
        score = svc.score(X_test, y_test)
        if score > best_score:
            best_score = score 
            print "C, gamma, score", C, gamma, score

# MinMaxScaler
norm = MinMaxScaler()
norm.fit(X)
T = norm.transform(X)

X_train, X_test, y_train, y_test = train_test_split(T, y, test_size=0.3, random_state=7)

for C in np.arange(0.05,2, 0.05):
    for gamma in np.arange(0.001, 0.1, 0.001):
        
        svc = SVC(C=C,gamma=gamma)
        svc.fit(X_train, y_train)
        score = svc.score(X_test, y_test)
        if score > best_score:
            best_score = score 
            print "C, gamma, score", C, gamma, score


#StandardScaler

norm = StandardScaler()
norm.fit(X)
T = norm.transform(X)

X_train, X_test, y_train, y_test = train_test_split(T, y, test_size=0.3, random_state=7)

for C in np.arange(0.05,2, 0.05):
    for gamma in np.arange(0.001, 0.1, 0.001):
        
        svc = SVC(C=C,gamma=gamma)
        svc.fit(X_train, y_train)
        score = svc.score(X_test, y_test)
        if score > best_score:
            best_score = score 
            print "C, gamma, score", C, gamma, score