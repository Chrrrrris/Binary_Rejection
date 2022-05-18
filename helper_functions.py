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


def get_sigma(x):
    """
    Parameters
    ----------
    x : list
        x is a list of g_band apparent magnitude

    Returns
    -------
    sigma : np.array
        sigma is a numpy array of adjusted gaia photometric uncertainty of corresponding g_band apparent magnitude

    """
    sigma =[]
    for n in range(len(x)):
        if x[n] <13:
            precision = 10*(10**(-3))
        elif x[n] < 18:
            precision = 16*(10**(-3))
        else:
            precision = 100*(10**(-3))
        #error propagation for sigma
        sigma.append(np.sqrt(2*precision**2))
    sigma = np.array(sigma)
    return sigma

def rejection(x_s1, y_s1, outliers_5sigma, y_5sigma):
    """

    Parameters
    ----------
    x_s1 : list
        segment of g_band apparent magnitude
    y_s1 : list
        segment of bp-rp band apparent magnitude corresponding to the g_band magnitude
    outliers_5sigma : list
        a list of g_band apparent magnitude that is rejected with a bp-rp deviation above 5 sigma level
    y_5sigma : list
        a list of bp-rp apparent magnitude that is rejected with a deviation above 5 sigma level

    Returns
    -------
    x_s1 : list
        segment of "good points'" g_band apparent magnitude
    y_s1 : list
        segment of "good points'" bp-rp apparent magnitude corresponding to gmag
    outliers_5sigma : list
        a list of g_band magnitude corresponding to rejected stars
    y_5sigma : list
        a list of bp-rp magnitude corresponding to rejected stars
    model : 1-D spline object
        a spline model object fitted from x_s1 and y_s1

    """
    # get the difference between the values predicted by the model and the actual bp-rp values
    for j in range(20):
        #diff is a list of difference
        diff = []
        #initialize a spline model from x_s1 and y_s1
        model = UnivariateSpline(x_s1, y_s1, k=4)
        #generate a list of predicted bp-rp band corresponding to x_s1
        predicted = model(x_s1)
        # get the difference between the values predicted by the model and the actual bp-rp values
        diff = y_s1 - predicted

        #get the photometric uncertainties of stars with the g_mag
        sigma = get_sigma(x_s1)
        #generate lists of "good points"
        x_filtered = []
        y_filtered = []
        #reject significant outliers before fittings that perform 5 sigma rejection
        if j == 0:
            for k in range(len(diff)):
                if np.abs(diff[k]) >= 10*sigma[k]:
                    outliers_5sigma.append(x_s1[k])
                    y_5sigma.append(y_s1[k])

                else:
                    x_filtered.append(x_s1[k])
                    y_filtered.append(y_s1[k])
            x_s1 = x_filtered
            y_s1 = y_filtered
        #perform 5 sigma rejection
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
    """
    
    Parameters
    ----------
    df_cluster : panda dataframe
        panda dataframe of stars in a cluster with sorted g_mag

    Returns
    -------
    x_s1 : list
        a list of g_mag segments with g_mag<13
    y_s1 : list
        a list of bp-rp segments corresponding to g_mag<13
    x_s2 : list
        a list of g_mag segments with g_mag>13
    y_s2 : TYPE
        a list of bp-rp segments corresponding to g_mag>13
    x_non_ms : list
        a list of g_mag segments of non-main-sequence stars
    y_non_ms : list
        a list of bp-rp_mag segments of non-main-sequence stars

    """
    #generate lists of segments
    x_s1=[]
    y_s1=[]
    x_s2=[]
    y_s2=[]
    x_non_ms=[]
    y_non_ms=[]
    #x is the g_mag of stars in a cluster
    x = np.array(df_cluster['gmag'])
    #y is the bp-rp mag of stars in a cluster
    y = np.array(df_cluster['bpmag']- df_cluster['rpmag'])
    #get the g_mag of the turnoff point
    x_turn = turnoff(x,y)
    #function needed to exclude non main sequence
    for i in range(len(x)):
        # append g_mag and bp-rp mag to their corresponding segments.
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
    """
    

    Parameters
    ----------
    x : list
        g_mag apparent magnitude of stars
    y : list
        bp-rp apparent mag of stars.

    Returns
    -------
    x_turn : float
        the g_band mag of the turnoff point.

    """
    #get copies of g_mag and bp-rp mag
    x_copy = x
    y_copy = y
    #lists of significant outliers
    x_rej=[]
    y_rej=[]
    #lists of filtered stars
    x_filtered = []
    y_filtered = []
    #assume the turnoff is the g_mag of the brightest star
    x_turn = x_copy[0]
    #iterate twice: the first iteration rejects significant outliers; the second iteration fits the model and finds the turnoff,
    for j in range(2):
        #get the difference between fitted model and the observed bp-rp
        diff = []
        model = UnivariateSpline(x_copy, y_copy, k=4)
        predicted = model(x_copy)

        # get the difference between the values predicted by the model and the actual bp-rp values
        diff = y_copy - predicted

        #get the photometric uncertainties of stars with the g_mag
        sigma = get_sigma(x_copy)
        #critical points are those whose derivatives vanish
        x_critical = model.derivative().roots()
        for i in range(len(x_critical)):
            #find the point whose second derivative is positive
            if model.derivatives(x_critical[i])[2]>0:
                x_turn = x_critical[i]
        #reject significant outliers in the first iteration
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