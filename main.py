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

with fits.open('./open_cluster/NearbyClustersGaia_GALEX_SDSS_Skymapper_2MASS_WISE_ASASSN.fits') as hdu:
    asn_table = Table(hdu[1].data)
df = asn_table.to_pandas()
df_cluster = df[(df['Cluster']=='alphaPer')]
df_cluster.dropna(subset = ['gmag','rpmag','bpmag'], inplace = True)
df_cluster.sort_values(by=['gmag'], inplace=True)
x_s1, y_s1, x_s2, y_s2, x_non_ms, y_non_ms = helper_functions.segments(df_cluster)
outliers_5sigma=[]
y_5sigma = []
x_s1, y_s1, outliers_5sigma, y_5sigma, model = helper_functions.rejection(x_s1, y_s1, outliers_5sigma, y_5sigma)
x_s2, y_s2, outliers_5sigma, y_5sigma,model1 = helper_functions.rejection(x_s2, y_s2, outliers_5sigma, y_5sigma)


x_filtered = x_s1+x_s2
y_filtered = y_s1+y_s2

for i in range(len(outliers_5sigma)):
    if(outliers_5sigma[i]>=13 and y_5sigma[i]< model(outliers_5sigma[i])):
        x_filtered.append(outliers_5sigma[i])
        y_filtered.append(y_5sigma[i])
    elif(outliers_5sigma[i]<13 and y_5sigma[i]< model1(outliers_5sigma[i])):
        x_filtered.append(outliers_5sigma[i])
        y_filtered.append(y_5sigma[i])
DF = df_cluster[df_cluster['gmag']==x_filtered[0]]
for i in range(1, len(x_filtered)):
    df2 = df_cluster[df_cluster['gmag']==x_filtered[i]]
    DF = DF.append(df2, ignore_index=True)
DF.to_csv('alphaPer.csv', index=False)