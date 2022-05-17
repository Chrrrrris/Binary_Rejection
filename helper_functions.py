#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 17 03:47:47 2022

@author: chris.w
"""
from astropy.io import fits
from astropy.table import Table
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline

def get_sigma(x, sigma):
    for n in range(len(x)):
        if x[n] <13:
            precision = 8*(10**(-3))
        elif x[n] < 18:
            precision = 16*(10**(-3))
        else:
            precision = 200*(10**(-3))
        sigma.append(np.sqrt(2*precision**2))
    sigma = np.array(sigma)
    return sigma

def rejection(x_s1, y_s1, outliers_5sigma, y_5sigma):
    # get the difference between the values predicted by the model and the actual bp-rp values
    for j in range(20):
        diff = []
        model = UnivariateSpline(x_s1, y_s1, k=4)
        predicted = model(x_s1)
        # get the difference between the values predicted by the model and the actual bp-rp values
        diff = y_s1 - predicted

        #get the standard deviation of the differences
        sigma = []
        sigma = get_sigma(x_s1, sigma)
        x_filtered = []
        y_filtered = []
        if j == 0:
            for k in range(len(diff)):
                if np.abs(diff[k]) >= 10*sigma[k]: #and 13<= x_s1[k]< 18:
                    outliers_5sigma.append(x_s1[k])
                    y_5sigma.append(y_s1[k])

                else:
                    x_filtered.append(x_s1[k])
                    y_filtered.append(y_s1[k])
            x_s1 = x_filtered
            y_s1 = y_filtered
        else:
            for k in range(len(diff)):
                if np.abs(diff[k]) >= 5*sigma[k]:
                    outliers_5sigma.append(x_s1[k])
                    y_5sigma.append(y_s1[k])
                else:
                    x_filtered.append(x_s1[k])
                    y_filtered.append(y_s1[k])
            x_s1 = x_filtered
            y_s1 = y_filtered
        
    return x_s1, y_s1,outliers_5sigma, y_5sigma, model

def segments(df_cluster):
    x_s1=[]
    y_s1=[]
    x_s2=[]
    y_s2=[]
    x_non_ms=[]
    y_non_ms=[]
    x = np.array(df_cluster['gmag'])
    y = np.array(df_cluster['bpmag']- df_cluster['rpmag'])
    x_turn = turnoff(x,y)
    #function needed to exclude non main sequence
    for i in range(len(x)):
        if (x[i]<=13 and x[i] >= x_turn):
            x_s1.append(x[i])
            y_s1.append(y[i])
        elif(x[i]>13):
            x_s2.append(x[i])
            y_s2.append(y[i])
        elif(x[i]<x_turn):
            x_non_ms.append(x[i])
            y_non_ms.append(y[i])
    return x_s1, y_s1, x_s2, y_s2, x_non_ms, y_non_ms

def turnoff(x,y):
    x_copy = x
    y_copy = y
    x_rej=[]
    y_rej=[]
    x_filtered = []
    y_filtered = []
    x_turn = x_copy[0]
    for j in range(2):
        diff = []
        model = UnivariateSpline(x_copy, y_copy, k=4)
        predicted = model(x_copy)

        # get the difference between the values predicted by the model and the actual bp-rp values
        diff = y_copy - predicted

        #get the standard deviation of the differences
        sigma = []
        sigma = get_sigma(x_copy,sigma)
        x_critical = model.derivative().roots()
        for i in range(len(x_critical)):
            if model.derivatives(x_critical[i])[2]>0:
                x_turn = x_critical[i]
        if(j == 0):
            for k in range(len(diff)):
                if np.abs(diff[k]) < 20*sigma[k] or x_copy[k] <= x_turn: #and 13<= x_s1[k]< 18:
                    x_filtered.append(x_copy[k])
                    y_filtered.append(y_copy[k])
                elif diff[k] >= 20*sigma[k] and x_copy[k] > x_turn:
                    x_rej.append(x_copy[k])
                    y_rej.append(y_copy[k])
        x_copy = x_filtered
        y_copy = y_filtered
    return x_turn