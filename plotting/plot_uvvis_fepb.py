#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: sashagardner
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as stats
from scipy.special import wofz

plt.rcParams['figure.dpi'] = 600


#data for FePB experimental plot
FePB_exp = pd.read_csv("/Users/sashagardner/Documents/Cooksy/Carrano/webplot_data_experimental_spectra/FePB_solid.csv")

exp_nm= FePB_exp["wavelength"].values
exp_abs = FePB_exp["absorbance"].values

conv = 1239.84193 #for converting between energy and wavelength (hv)

exp_ex= np.divide(conv, exp_nm)

#data for calculated FeAB spectra

#conformer 13
read13=pd.read_csv("/Users/sashagardner/Documents/Cooksy/Carrano/Biometals_paper/Biometal_output_files/temp/ex_states_13.csv")
ex13 = read13["E_ex"].values
nm13 = read13["nm"].values #transition state energies eV
osc13 = read13["osc"].values #oscilator strength for transition states

#conformer 14
read14=pd.read_csv("/Users/sashagardner/Documents/Cooksy/Carrano/Biometals_paper/Biometal_output_files/temp/ex_states_14.csv")
ex14 = read14["E_ex"].values
nm14 = read14["nm"].values #transition state energies eV
osc14 = read14["osc"].values #oscilator strength for transition states

#conformer 15
read15=pd.read_csv("/Users/sashagardner/Documents/Cooksy/Carrano/Biometals_paper/Biometal_output_files/temp/ex_states_15.csv")
ex15 = read15["E_ex"].values
nm15 = read15["nm"].values #transition state energies eV
osc15 = read15["osc"].values #oscilator strength for transition states

#conformer 16
read16=pd.read_csv("/Users/sashagardner/Documents/Cooksy/Carrano/Biometals_paper/Biometal_output_files/temp/ex_states_16.csv")
ex16 = read16["E_ex"].values
nm16 = read16["nm"].values #transition state energies eV
osc16 = read16["osc"].values #oscilator strength for transition states

#change gamma and sigma to fit curve
gamma = 0.25
#sigma = 0.25


x = np.linspace(1.75, 5.5, 1000)

n13 = len(ex13)
n14 = len(ex14)
n15 = len(ex15)
n16 = len(ex16)
m = len(x)

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(10, 8))

#----13--------------------------------------

# plot experimental spectrum as area
ax1.fill_between(exp_nm, exp_abs, label='FePB experimental', color='gainsboro')

#plot pdf
pdf13 = np.zeros(shape=x.shape)
for i in range(n13):
    uvvis13 = 0.63*osc13[i] * stats.cauchy.pdf(x, ex13[i], gamma)
    #v = x - ex13[i]
    #uvvis13 = 0.64 * osc13[i] * np.real(wofz((v + 1j*gamma) / (sigma * np.sqrt(2)))) / (sigma * np.sqrt(2 * np.pi))
    pdf13 += uvvis13
    w = np.divide(conv, x) #convert to wavelength
ax1.plot(w, pdf13, label='calculated')

#subplot details
ax1.set_xlabel('Wavelength (nm)')
ax1.set_ylabel('Absorbance')
ax1.set_title('FePB enol cis-cis para-fac')
ax1.set_xticks(np.arange(250, 701, step=50))
ax1.set_ylim(0, 0.8)
ax1.legend(loc= 'upper center')
ax1.margins(x=0, y=0)


#letters on plots
ax1.text(0.9, 0.98, '(a)', fontsize= 'large', horizontalalignment='left', verticalalignment='top', transform=ax1.transAxes)


#---14---------------------------------------------

# plot experimental spectrum as area
ax2.fill_between(exp_nm, exp_abs, label='FePB experimental', color='gainsboro')

#plot pdf
pdf14 = np.zeros(shape=x.shape)
for i in range(n14):
    uvvis14 = 0.63*osc14[i] * stats.cauchy.pdf(x, ex14[i], gamma)
    #v = x - ex14[i]
    #uvvis14 = 0.64 * osc14[i] * np.real(wofz((v + 1j*gamma) / (sigma * np.sqrt(2)))) / (sigma * np.sqrt(2 * np.pi))
    pdf14 += uvvis14
    w = np.divide(conv, x) #convert to wavelength
ax2.plot(w, pdf14, label='calculated')

#subplot details
ax2.set_xlabel('Wavelength (nm)')
ax2.set_ylabel('Absorbance')
ax2.set_title('FePB cis-cis meta-fac')
ax2.set_xticks(np.arange(250, 701, step=50))
ax2.set_ylim(0, 0.8)
ax2.legend(loc= 'upper center')
ax2.margins(x=0, y=0)


#letters on plots
ax2.text(0.9, 0.98, '(b)', fontsize= 'large', horizontalalignment='left', verticalalignment='top', transform=ax2.transAxes)

#----15--------------------------------------------

# plot experimental spectrum as area
ax3.fill_between(exp_nm, exp_abs, label='FePB experimental', color='gainsboro')

#plot pdf
pdf15 = np.zeros(shape=x.shape)
for i in range(n15):
    uvvis15 = 0.63*osc15[i] * stats.cauchy.pdf(x, ex15[i], gamma)
    #v = x - ex15[i]
    #uvvis15 = 0.64 * osc15[i] * np.real(wofz((v + 1j*gamma) / (sigma * np.sqrt(2)))) / (sigma * np.sqrt(2 * np.pi))
    pdf15 += uvvis15
    w = np.divide(conv, x) #convert to wavelength
ax3.plot(w, pdf15, label='calculated')

#subplot details
ax3.set_xlabel('Wavelength (nm)')
ax3.set_ylabel('Absorbance')
ax3.set_title('FePB cis-trans para-mer')
ax3.set_xticks(np.arange(250, 701, step=50))
ax3.set_ylim(0, 0.8)
ax3.legend(loc= 'upper center')
ax3.margins(x=0, y=0)


#letters on plots
ax3.text(0.9, 0.98, '(c)', fontsize= 'large', horizontalalignment='left', verticalalignment='top', transform=ax3.transAxes)

#---16----------------------------------------------

# plot experimental spectrum as area
ax4.fill_between(exp_nm, exp_abs, label='FePB experimental', color='gainsboro')

#plot pdf
pdf16= np.zeros(shape=x.shape)
for i in range(n16):
    uvvis16 = 0.63*osc16[i] * stats.cauchy.pdf(x, ex16[i], gamma)
    #v = x - ex16[i]
    #uvvis16 = 0.64 * osc16[i] * np.real(wofz((v + 1j*gamma) / (sigma * np.sqrt(2)))) / (sigma * np.sqrt(2 * np.pi))
    pdf16 += uvvis16
    w = np.divide(conv, x) #convert to wavelength
ax4.plot(w, pdf16, label='calculated')

#subplot details
ax4.set_xlabel('Wavelength (nm)')
ax4.set_ylabel('Absorbance')
ax4.set_title('FePB cis-trans meta-mer')
ax4.set_xticks(np.arange(250, 701, step=50))
ax4.set_ylim(0, 0.8)
ax4.legend(loc= 'upper center')
ax4.margins(x=0, y=0)

#letters on plots
ax4.text(0.9, 0.98, '(d)', fontsize= 'large', horizontalalignment='left', verticalalignment='top', transform=ax4.transAxes)
        


plt.tight_layout()
plt.show()
