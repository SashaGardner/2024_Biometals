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
FeABphoto_exp = pd.read_csv("/Users/sashagardner/Documents/Cooksy/Carrano/webplot_data_experimental_spectra/FeABphoto_dotted.csv")

exp_nm= FeABphoto_exp["wavelength"].values
exp_abs = FeABphoto_exp["absorbance"].values

conv = 1239.84193 #for converting between energy and wavelength (hv)


#data for calculated FeAB spectra

#conformer 9
read9=pd.read_csv("/Users/sashagardner/Documents/Cooksy/Carrano/Biometals_paper/Biometal_output_files/temp/ex_states_9.csv")
ex9 = read9["E_ex"].values
nm9 = read9["nm"].values #transition state energies eV
osc9 = read9["osc"].values #oscilator strength for transition states

#conformer 10
read10=pd.read_csv("/Users/sashagardner/Documents/Cooksy/Carrano/Biometals_paper/Biometal_output_files/temp/ex_states_10.csv")
ex10 = read10["E_ex"].values
nm10 = read10["nm"].values #transition state energies eV
osc10 = read10["osc"].values #oscilator strength for transition states

#conformer 11
read11=pd.read_csv("/Users/sashagardner/Documents/Cooksy/Carrano/Biometals_paper/Biometal_output_files/temp/ex_states_11.csv")
ex11 = read11["E_ex"].values
nm11 = read11["nm"].values #transition state energies eV
osc11 = read11["osc"].values #oscilator strength for transition states

#conformer 12
read12=pd.read_csv("/Users/sashagardner/Documents/Cooksy/Carrano/Biometals_paper/Biometal_output_files/temp/ex_states_12.csv")
ex12 = read12["E_ex"].values
nm12 = read12["nm"].values #transition state energies eV
osc12 = read12["osc"].values #oscilator strength for transition states


#change gamma and sigma to fit curve
gamma = 0.25
#sigma = 0.25



x = np.linspace(2.00, 5.50, 1000)

n9 = len(ex9)
n10 = len(ex10)
n11 = len(ex11)
n12 = len(ex12)


fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(10, 8))

#----9----------------------------------------------

# plot experimental spectrum as area
ax1.fill_between(exp_nm, exp_abs, label='FeAB* experimental', color='gainsboro')

#plot pdf
pdf9 = np.zeros(shape=x.shape)
for i in range(n9):
    uvvis9 = 1.55*osc9[i] * stats.cauchy.pdf(x, ex9[i], gamma)
    #v = x - ex9[i]
    #uvvis9 = 1.56 * osc9[i] * np.real(wofz((v + 1j*gamma) / (sigma * np.sqrt(2)))) / (sigma * np.sqrt(2 * np.pi))
    pdf9 += uvvis9
    w = np.divide(conv, x) #convert to wavelength
ax1.plot(w, pdf9, label='calculated')

#subplot details
ax1.set_xlabel('Wavelength (nm)')
ax1.set_ylabel('Absorbance')
ax1.set_title('FeAB* enol cis-cis C-fac')
ax1.set_xticks(np.arange(250, 651, step=50))
ax1.set_ylim(0, 0.6)
ax1.legend(loc= 'upper center')
ax1.margins(x=0, y=0)


#letters on plots
ax1.text(0.9, 0.98, '(a)', fontsize= 'large', horizontalalignment='left', verticalalignment='top', transform=ax1.transAxes)


#---10--------------------------------------------

# plot experimental spectrum as area
ax2.fill_between(exp_nm, exp_abs, label='FeAB* experimental', color='gainsboro')

#plot pdf
pdf10 = np.zeros(shape=x.shape)
for i in range(n10):
    uvvis10 = 1.55*osc10[i] * stats.cauchy.pdf(x, ex10[i], gamma)
    #v = x - ex10[i]
    #uvvis10 = 1.56 * osc10[i] * np.real(wofz((v + 1j*gamma) / (sigma * np.sqrt(2)))) / (sigma * np.sqrt(2 * np.pi))
    pdf10 += uvvis10
    w = np.divide(conv, x) #convert to wavelength
ax2.plot(w, pdf10, label='calculated')

#subplot details
ax2.set_xlabel('Wavelength (nm)')
ax2.set_ylabel('Absorbance')
ax2.set_title('FeAB* enol cis-cis N-fac')
ax2.set_xticks(np.arange(250, 651, step=50))
ax2.set_ylim(0, 0.6)
ax2.legend(loc= 'upper center')
ax2.margins(x=0, y=0)

#letters on plots
ax2.text(0.9, 0.98, '(b)', fontsize= 'large', horizontalalignment='left', verticalalignment='top', transform=ax2.transAxes)

#----11---------------------------------------------

# plot experimental spectrum as area
ax3.fill_between(exp_nm, exp_abs, label='FeAB* experimental', color='gainsboro')

#plot pdf
pdf11 = np.zeros(shape=x.shape)
for i in range(n11):
    uvvis11 = 1.55*osc11[i] * stats.cauchy.pdf(x, ex11[i], gamma)
    #v = x - ex11[i]
    #uvvis11 = 1.56 * osc11[i] * np.real(wofz((v + 1j*gamma) / (sigma * np.sqrt(2)))) / (sigma * np.sqrt(2 * np.pi))
    pdf11 += uvvis11
    w = np.divide(conv, x) #convert to wavelength
ax3.plot(w, pdf11, label='calculated')

#subplot details
ax3.set_xlabel('Wavelength (nm)')
ax3.set_ylabel('Absorbance')
ax3.set_title('FeAB* enol cis-trans C-mer')
ax3.set_xticks(np.arange(250, 651, step=50))
ax3.set_ylim(0, 0.6)
ax3.legend(loc= 'upper center')
ax3.margins(x=0, y=0)

#letters on plots
ax3.text(0.9, 0.98, '(c)', fontsize= 'large', horizontalalignment='left', verticalalignment='top', transform=ax3.transAxes)

#---12-------

# plot experimental spectrum as area
ax4.fill_between(exp_nm, exp_abs, label='FeAB* experimental', color='gainsboro')

#plot pdf
pdf12= np.zeros(shape=x.shape)
for i in range(n12):
    uvvis12 = 1.55*osc12[i] * stats.cauchy.pdf(x, ex12[i], gamma)
    #v = x - ex12[i]
    #uvvis12 = 1.56 * osc12[i] * np.real(wofz((v + 1j*gamma) / (sigma * np.sqrt(2)))) / (sigma * np.sqrt(2 * np.pi))
    pdf12 += uvvis12
    w = np.divide(conv, x) #convert to wavelength
ax4.plot(w, pdf12, label='calculated')

#subplot details
ax4.set_xlabel('Wavelength (nm)')
ax4.set_ylabel('Absorbance')
ax4.set_title('FeAB* enol cis-trans N-mer')
ax4.set_xticks(np.arange(250, 651, step=50))
ax4.set_ylim(0, 0.6)
ax4.legend(loc= 'upper center')
ax4.margins(x=0, y=0)

#letters on plots
ax4.text(0.9, 0.98, '(d)', fontsize= 'large', horizontalalignment='left', verticalalignment='top', transform=ax4.transAxes)
        

plt.tight_layout()
plt.show()
