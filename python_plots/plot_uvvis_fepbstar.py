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
FePBstar_exp = pd.read_csv("/Users/sashagardner/Documents/Cooksy/Carrano/webplot_data_experimental_spectra/FePBstar_exp_2.csv")

exp_nm= FePBstar_exp["wavelength"].values
exp_abs = FePBstar_exp["absorbance"].values

conv = 1239.84193 #for converting between energy and wavelength (hv)

#data for calculated FeAB spectra

#conformer 21
read21=pd.read_csv("/Users/sashagardner/Documents/Cooksy/Carrano/Biometals_paper/Biometal_output_files/temp/ex_states_21.csv")
ex21 = read21["E_ex"].values
nm21 = read21["nm"].values #transition state energies eV
osc21 = read21["osc"].values #oscilator strength for transition states

#conformer 22
read22=pd.read_csv("/Users/sashagardner/Documents/Cooksy/Carrano/Biometals_paper/Biometal_output_files/temp/ex_states_22.csv")
ex22 = read22["E_ex"].values
nm22 = read22["nm"].values #transition state energies eV
osc22 = read22["osc"].values #oscilator strength for transition states

#conformer 23
read23=pd.read_csv("/Users/sashagardner/Documents/Cooksy/Carrano/Biometals_paper/Biometal_output_files/temp/ex_states_23.csv")
ex23 = read23["E_ex"].values
nm23 = read23["nm"].values #transition state energies eV
osc23 = read23["osc"].values #oscilator strength for transition states

#conformer 24
read24=pd.read_csv("/Users/sashagardner/Documents/Cooksy/Carrano/Biometals_paper/Biometal_output_files/temp/ex_states_24.csv")
ex24 = read24["E_ex"].values
nm24 = read24["nm"].values #transition state energies eV
osc24 = read24["osc"].values #oscilator strength for transition states


#change gamma and sigma to fit curve
gamma = 0.25
#sigma = 0.25


x = np.linspace(1.75, 5.5, 1000)
n21 = len(ex21)
n22 = len(ex22)
n23 = len(ex23)
n24 = len(ex24)

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(10, 8))

#----21---------------------------------------------------

# plot experimental spectrum as area
ax1.fill_between(exp_nm, exp_abs, label='FePB* experimental', color='gainsboro')

#plot pdf
pdf21 = np.zeros(shape=x.shape)
for i in range(n21):
    uvvis21 = 0.63*osc21[i] * stats.cauchy.pdf(x, ex21[i], gamma)
    #v = x - ex21[i]
    #uvvis21 = 0.64 * osc21[i] * np.real(wofz((v + 1j*gamma) / (sigma * np.sqrt(2)))) / (sigma * np.sqrt(2 * np.pi))
    pdf21 += uvvis21
    w = np.divide(conv, x) #convert to wavelength
ax1.plot(w, pdf21, label='calculated / 2')

#subplot details
ax1.set_xlabel('Wavelength (nm)')
ax1.set_ylabel('Absorbance')
ax1.set_title('FePB* enol cis-cis para-fac')
ax1.set_xticks(np.arange(250, 701, step=50))
ax1.set_ylim(0, 0.8)
ax1.legend(loc= 'upper center')
ax1.margins(x=0, y=0)

#letters on plots
ax1.text(0.9, 0.98, '(a)', fontsize= 'large', horizontalalignment='left', verticalalignment='top', transform=ax1.transAxes)


#---22-----------------------------------------------

# plot experimental spectrum as area
ax2.fill_between(exp_nm, exp_abs, label='FePB* experimental', color='gainsboro')

#plot pdf
pdf22 = np.zeros(shape=x.shape)
for i in range(n22):
    uvvis22 = 0.63*osc22[i] * stats.cauchy.pdf(x, ex22[i], gamma)
    #v = x - ex22[i]
    #uvvis22 = 0.64 * osc22[i] * np.real(wofz((v + 1j*gamma) / (sigma * np.sqrt(2)))) / (sigma * np.sqrt(2 * np.pi))
    pdf22 += uvvis22
    w = np.divide(conv, x) #convert to wavelength
ax2.plot(w, pdf22, label='calculated / 2')

#subplot details
ax2.set_xlabel('Energy (eV)')
ax2.set_ylabel('Absorbance')
ax2.set_title('FePB* cis-cis meta-fac')
ax2.set_xticks(np.arange(250, 701, step=50))
ax2.set_ylim(0, 0.8)
ax2.legend(loc= 'upper center')
ax2.margins(x=0, y=0)

#letters on plots
ax2.text(0.9, 0.98, '(b)', fontsize= 'large', horizontalalignment='left', verticalalignment='top', transform=ax2.transAxes)

#----23------------------------------------------------

# plot experimental spectrum as area
ax3.fill_between(exp_nm, exp_abs, label='FePB* experimental', color='gainsboro')

#plot pdf
pdf23 = np.zeros(shape=x.shape)
for i in range(n23):
    uvvis23= 0.63*osc23[i] * stats.cauchy.pdf(x, ex23[i], gamma)
    #v = x - ex23[i]
    #uvvis23 = 0.64 * osc23[i] * np.real(wofz((v + 1j*gamma) / (sigma * np.sqrt(2)))) / (sigma * np.sqrt(2 * np.pi))
    pdf23 += uvvis23
    w = np.divide(conv, x) #convert to wavelength
ax3.plot(w, pdf23, label='calculated / 2')

#subplot details
ax3.set_xlabel('Wavelength (nm)')
ax3.set_ylabel('Absorbance')
ax3.set_title('FePB* cis-trans para-mer')
ax3.set_xticks(np.arange(250, 701, step=50))
ax3.set_ylim(0, 0.8)
ax3.legend(loc= 'upper center')
ax3.margins(x=0, y=0)

#letters on plots
ax3.text(0.9, 0.98, '(c)', fontsize= 'large', horizontalalignment='left', verticalalignment='top', transform=ax3.transAxes)

#---24---------------------------------------------------

# plot experimental spectrum as area
ax4.fill_between(exp_nm, exp_abs, label='FePB* experimental', color='gainsboro')

#plot pdf
pdf24= np.zeros(shape=x.shape)
for i in range(n24):
    uvvis24 = 0.63*osc24[i] * stats.cauchy.pdf(x, ex24[i], gamma)
    #v = x - ex24[i]
    #uvvis24 = 0.64 * osc24[i] * np.real(wofz((v + 1j*gamma) / (sigma * np.sqrt(2)))) / (sigma * np.sqrt(2 * np.pi))
    pdf24 += uvvis24
    w = np.divide(conv, x) #convert to wavelength
ax4.plot(w, pdf24, label='calculated / 2')

#subplot details
ax4.set_xlabel('Wavelength (nm)')
ax4.set_ylabel('Absorbance')
ax4.set_title('FePB* cis-trans meta-mer')
ax4.set_xticks(np.arange(250, 701, step=50))
ax4.set_ylim(0, 0.8)
ax4.legend(loc= 'upper center')
ax4.margins(x=0, y=0)

#letters on plots
ax4.text(0.9, 0.98, '(d)', fontsize= 'large', horizontalalignment='left', verticalalignment='top', transform=ax4.transAxes)
        

plt.tight_layout()
plt.show()
