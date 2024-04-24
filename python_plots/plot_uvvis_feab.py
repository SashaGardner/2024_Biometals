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


#data for FeAb experimental plot
FeAB_exp = pd.read_csv("/Users/sashagardner/Documents/Cooksy/Carrano/webplot_data_experimental_spectra/FeAB_solid.csv")

exp_nm= FeAB_exp["wavelength"].values
exp_abs = FeAB_exp["absorbance"].values

conv = 1239.84193 #for converting between energy and wavelength (hv)

#data for calculated FeAB spectra

#conformer 1
read1=pd.read_csv("/Users/sashagardner/Documents/Cooksy/Carrano/Biometals_paper/Biometal_output_files/temp/ex_states_1.csv")
ex1 = read1["E_ex"].values
nm1 = read1["nm"].values #transition state energies eV
osc1 = read1["osc"].values #oscilator strength for transition states

#conformer 2
read2=pd.read_csv("/Users/sashagardner/Documents/Cooksy/Carrano/Biometals_paper/Biometal_output_files/temp/ex_states_2.csv")
ex2 = read2["E_ex"].values
nm2 = read2["nm"].values #transition state energies eV
osc2 = read2["osc"].values #oscilator strength for transition states

#conformer 3
read3=pd.read_csv("/Users/sashagardner/Documents/Cooksy/Carrano/Biometals_paper/Biometal_output_files/temp/ex_states_3.csv")
ex3 = read3["E_ex"].values
nm3 = read3["nm"].values #transition state energies eV
osc3 = read3["osc"].values #oscilator strength for transition states

#conformer 4
read4=pd.read_csv("/Users/sashagardner/Documents/Cooksy/Carrano/Biometals_paper/Biometal_output_files/temp/ex_states_4.csv")
ex4 = read4["E_ex"].values
nm4 = read4["nm"].values #transition state energies eV
osc4 = read4["osc"].values #oscilator strength for transition states


#change gamma and sigma to fit curve
gamma = 0.25
#sigma = 0.25




x = np.linspace(2.00, 5.50, 1000)

n1 = len(ex1)
n2 = len(ex2)
n3 = len(ex3)
n4 = len(ex4)

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(10, 8))

#----1-----------------------------------------

# plot experimental spectrum as area
ax1.fill_between(exp_nm, exp_abs, label='FeAB experimental', color='gainsboro')



#plot pdf
pdf1 = np.zeros(shape=x.shape)

for i in range(n1):
    uvvis1 = 1.55 * osc1[i] * stats.cauchy.pdf(x, ex1[i], gamma)
    #v = x - ex1[i]
    #uvvis1 = 1.56 * osc1[i] * np.real(wofz((v + 1j*gamma) / (sigma * np.sqrt(2)))) / (sigma * np.sqrt(2 * np.pi))
    pdf1 += uvvis1
    w = np.divide(conv, x) #convert to wavelength
ax1.plot(w, pdf1, label='calculated')


#subplot details
ax1.set_xlabel('Wavelength (nm)')
ax1.set_ylabel('Absorbance')
ax1.set_title('FeAB cis-cis C-fac')
ax1.set_xticks(np.arange(250, 651, step=50))
ax1.set_ylim(0, 0.6)
ax1.legend(loc= 'upper center')
ax1.margins(x=0, y=0)

#letters on plots
ax1.text(0.9, 0.98, '(a)', fontsize= 'large', horizontalalignment='left', verticalalignment='top', transform=ax1.transAxes)


#---2-----------------------------------------------------

# plot experimental spectrum as area
ax2.fill_between(exp_nm, exp_abs, label='FeAB experimental', color='gainsboro')

#plot pdf
pdf2 = np.zeros(shape=x.shape)
for i in range(n2):
    uvvis2 = 1.55*osc2[i] * stats.cauchy.pdf(x, ex2[i], gamma)
    #v = x - ex2[i]
    #uvvis2 = 1.56 * osc2[i] * np.real(wofz((v + 1j*gamma) / (sigma * np.sqrt(2)))) / (sigma * np.sqrt(2 * np.pi))
    pdf2 += uvvis2
    w = np.divide(conv, x) #convert to wavelength
ax2.plot(w, pdf2, label='calculated')

#subplot details
ax2.set_xlabel('Wavelength (nm)')
ax2.set_ylabel('Absorbance')
ax2.set_title('FeAB cis-cis N-fac')
ax2.set_xticks(np.arange(250, 651, step=50))
ax2.set_ylim(0, 0.6)
ax2.legend(loc= 'upper center')
ax2.margins(x=0, y=0)

#letters on plots
ax2.text(0.9, 0.98, '(b)', fontsize= 'large', horizontalalignment='left', verticalalignment='top', transform=ax2.transAxes)

#----3---------------------------------------------------

# plot experimental spectrum as area
ax3.fill_between(exp_nm, exp_abs, label='FeAB experimental', color='gainsboro')

#plot pdf
pdf3 = np.zeros(shape=x.shape)
for i in range(n3):
    uvvis3 = 1.55*osc3[i] * stats.cauchy.pdf(x, ex3[i], gamma)
    #v = x - ex3[i]
    #uvvis3 = 1.56 * osc3[i] * np.real(wofz((v + 1j*gamma) / (sigma * np.sqrt(2)))) / (sigma * np.sqrt(2 * np.pi))
    pdf3 += uvvis3
    w = np.divide(conv, x)#convert to wavelength
ax3.plot(w, pdf3, label='calculated')

#subplot details
ax3.set_xlabel('Wavelength (nm)')
ax3.set_ylabel('Absorbance')
ax3.set_title('FeAB cis-trans C-mer')
ax3.set_xticks(np.arange(250, 651, step=50))
ax3.set_ylim(0, 0.6)
ax3.legend(loc= 'upper center')
ax3.margins(x=0, y=0)

#letters on plots
ax3.text(0.9, 0.98, '(c)', fontsize= 'large', horizontalalignment='left', verticalalignment='top', transform=ax3.transAxes)

#---4-------------------------------------------------------

# plot experimental spectrum as area
ax4.fill_between(exp_nm, exp_abs, label='FeAB experimental', color='gainsboro')

#plot pdf
pdf4 = np.zeros(shape=x.shape)
for i in range(n4):
    uvvis4 = 1.55*osc4[i] * stats.cauchy.pdf(x, ex4[i], gamma)
    #v = x - ex4[i]
    #uvvis4 = 1.56 * osc4[i] * np.real(wofz((v + 1j*gamma) / (sigma * np.sqrt(2)))) / (sigma * np.sqrt(2 * np.pi))
    pdf4 += uvvis4
    w = np.divide(conv, x) #convert to wavelength
ax4.plot(w, pdf4, label='calculated')

#subplot details
ax4.set_xlabel('Wavelength (nm)')
ax4.set_ylabel('Absorbance')
ax4.set_title('FeAB cis-trans N-mer')
ax4.set_xticks(np.arange(250, 651, step=50))
ax4.set_ylim(0, 0.6)
ax4.legend(loc= 'upper center')
ax4.margins(x=0, y=0)

#letters on plot
ax4.text(0.9, 0.98, '(d)', fontsize= 'large', horizontalalignment='left', verticalalignment='top', transform=ax4.transAxes)


plt.tight_layout()
plt.show()

