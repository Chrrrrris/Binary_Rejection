#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 17 03:49:10 2022

@author: chris.w
"""

from astropy.io import fits
from astropy.table import Table
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline
import helper_functions

#edit here to change the directory of the fit file
with fits.open('./open_cluster/NearbyClustersGaia_GALEX_SDSS_Skymapper_2MASS_WISE_ASASSN.fits') as hdu:
    asn_table = Table(hdu[1].data)
# convert fit to pd.dataframe
df = asn_table.to_pandas()
#edit here to change the name of the cluster
df_cluster = df[(df['Cluster']=='alphaPer')]
df_cluster.dropna(subset = ['gmag','rpmag','bpmag'], inplace = True)
#sort the dataframe
df_cluster.sort_values(by=['gmag'], inplace=True)
'''get the segments; x_s1 is a list of g_mag segments with g_mag<13, y_s1 is a list of bp-rp segments corresponding to g_mag<13;
x_s2 is a list of g_mag segments with g_mag>13; y_s2 is a list of bp-rp segments corresponding to g_mag>13;
x_non_ms is a list of g_mag segments of non-main-sequence stars; y_non_ms is a list of bp-rp_mag segments of non-main-sequence stars'''
x_s1, y_s1, x_s2, y_s2, x_non_ms, y_non_ms = helper_functions.segments(df_cluster)
outliers_5sigma=[]
y_5sigma = []
#perform iterative spline rejections
x_s1, y_s1, outliers_5sigma, y_5sigma, model = helper_functions.rejection(x_s1, y_s1, outliers_5sigma, y_5sigma)
x_s2, y_s2, outliers_5sigma, y_5sigma,model1 = helper_functions.rejection(x_s2, y_s2, outliers_5sigma, y_5sigma)


#get the total filtered g_mag and bp-rp mag
x_filtered = x_s1+x_s2
y_filtered = y_s1+y_s2

#generate Dataframe of photometric data of filtered stars
DF = df_cluster[df_cluster['gmag']==x_filtered[0]]
for i in range(1, len(x_filtered)):
    df2 = df_cluster[df_cluster['gmag']==x_filtered[i]]
    DF = DF.append(df2, ignore_index=True)
#output csv
#Edit here for csv name
DF.to_csv('alphaPer.csv', index=False)