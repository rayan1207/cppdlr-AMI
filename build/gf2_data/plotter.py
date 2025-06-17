# -*- coding: utf-8 -*-
"""
Created on Mon Jun  9 14:26:09 2025

@author: Rayan
"""

import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt


L11_5 = np.loadtxt('data/SE_U4_L11_i1_b5.txt')
L15_5 = np.loadtxt('data/SE_U4_L15_i1_b5.txt')
L21_5 = np.loadtxt('data/SE_U4_L21_i1_b5.txt')
L25_5 = np.loadtxt('data/SE_U4_L25_i1_b5.txt')
LT5 = np.loadtxt('data/LT_SE_B5_U4.txt')

L11 = np.loadtxt('data/SE_U4_L11_i1_b50.txt')
L15 = np.loadtxt('data/SE_U4_L15_i1_b50.txt')
L21 = np.loadtxt('data/SE_U4_L21_i1_b50.txt')
L25 = np.loadtxt('data/SE_U4_L25_i1_b50.txt')
L31 = np.loadtxt('data/SE_U4_L31_i1_b50.txt')
LT = np.loadtxt('data/LT_SE_B50_U4.txt')


plt.figure(figsize=[12,5])
plt.subplot(1,2,1)
plt.title(r'$\beta=5$')
plt.plot(L11_5[:,0],L11_5[:,1],'.-',label='L=11')
plt.plot(L15_5[:,0],L15_5[:,1],'.-',label='L=15')
plt.plot(L11_5[:,0],L21_5[:,1],'.-',label='L=21')
plt.plot(L11_5[:,0],L25_5[:,1],'.-',label='L=25')
plt.plot(LT5[:,0],LT5[:,1],'.-',label='$L=\infty$')
plt.legend(ncol=2,fontsize=12)
plt.xlabel(r'i$\omega_n$',fontsize=15)
plt.ylabel(r'Im $\Sigma(k_{an})$',fontsize=15)
plt.xlim(-0.5,20)
plt.subplot(1,2,2)
plt.title(r'$\beta=50$')
plt.plot(L11[:,0],L11[:,1],'.-',label='L=11')
plt.plot(L11[:,0],L15[:,1],'.-',label='L=15')
plt.plot(L11[:,0],L21[:,1],'.-',label='L=21')
plt.plot(L11[:,0],L25[:,1],'.-',label='L=25')
plt.plot(L11[:,0],L31[:,1],'.-',label='L=31')
plt.plot(LT[:,0],LT[:,1],'.-',label='$L=\infty$')
plt.legend(ncol=2,fontsize=12)
plt.xlabel(r'i$\omega_n$',fontsize=15)
plt.ylabel(r'Im $\Sigma(k_{an})$',fontsize=15)
plt.xlim(-0.5,9)
plt.show()


U4_i1 = np.loadtxt('data/SE_U4_L15_i1_b10_DCA.txt')
U4_i20 = np.loadtxt('data/SE_U4_L15_i20_b10_DCA.txt')
U4_i20i = np.loadtxt('data/SE_U4_L15_i20_b10_IPT.txt')


U8_i1 = np.loadtxt('data/SE_U8_L15_i1_b10_DCA.txt')
U8_i20 = np.loadtxt('data/SE_U8_L15_i20_b10_DCA.txt')
U8_i20i = np.loadtxt('data/SE_U8_L15_i20_b10_IPT.txt')



U12_i1 = np.loadtxt('data/SE_U12_L15_i1_b10_DCA.txt')
U12_i20 = np.loadtxt('data/SE_U12_L15_i20_b10_DCA.txt')
U12_i20i = np.loadtxt('data/SE_U12_L15_i20_b10_IPT.txt')



U16_i1 = np.loadtxt('data/SE_U16_L15_i1_b10_DCA.txt')
U16_i20 = np.loadtxt('data/SE_U16_L15_i20_b10_DCA.txt')
U16_i20i = np.loadtxt('data/SE_U16_L15_i20_b10_IPT.txt')


plt.figure(figsize=[12,8])
plt.suptitle(r'$\beta=10, L=15*15$ (patch size),$ N_{k}=81$ (inside patch)',fontsize=19)
plt.subplot(2,2,1)
plt.title('U=4',fontsize=12)
plt.plot(U4_i1[:,0],U4_i1[:,1],'^-',label='single shot')
plt.plot(U4_i1[:,0],U4_i20[:,1],'x--',label='DCA+IPT')
plt.plot(U4_i20i[:,0],U4_i20i[:,1],'.-',label='IPT')
plt.ylabel(r'Im $\Sigma(k_{an})$',fontsize=15)
plt.legend(fontsize=15)
plt.subplot(2,2,2)
plt.title('U=8',fontsize=12)
plt.plot(U8_i1[:,0],U8_i1[:,1],'^-',label='single shot')
plt.plot(U8_i1[:,0],U8_i20[:,1],'x--',label='DCA')
plt.plot(U8_i20i[:,0],U8_i20i[:,1],'.-',label='IPT')

plt.subplot(2,2,3)
plt.title('U=12',fontsize=12)
plt.plot(U12_i1[:,0],U12_i1[:,1],'^-',label='single shot')
plt.plot(U12_i1[:,0],U12_i20[:,1],'x--',label='DCA')
plt.plot(U12_i20i[:,0],U12_i20i[:,1],'.-',label='IPT')
plt.ylabel(r'Im $\Sigma(k_{an})$',fontsize=15)
plt.xlabel(r'i$\omega_n$',fontsize=15)
plt.subplot(2,2,4)
plt.title('U=16',fontsize=12)
plt.plot(U16_i1[:,0],U16_i1[:,1],'^-',label='single shot')
plt.plot(U16_i1[:,0],U16_i20[:,1],'x--',label='DCA+IPT')
plt.plot(U16_i20i[:,0],U16_i20i[:,1],'.-',label='IPT')
plt.xlabel(r'i$\omega_n$',fontsize=15)
plt.savefig('../DCA_IPT_b10_iter20.pdf')
plt.show()



