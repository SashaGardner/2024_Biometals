import matplotlib.pyplot as plt
import numpy as np
plt.rcParams['figure.dpi'] = 600


# Example data (you should replace this with your actual data)
molecules1 = ['cc C-fac', 'ct C-mer', 'cc N-fac', 'ct N-mer']
series1_energies = [0, 12.9, 3.8, 8.7]
series2_energies = [0.7, 0, 4.8, 4.4]

molecules2 = ['cc para-fac', 'ct para-mer', 'cc meta-fac', 'ct meta-mer']
series3_energies = np.array([2.8, 3.7, 1.1, 0])
series4_energies = np.array([0, 1.5, 1.4, 1.8])
#series5_energies = [ , , , ]

# Create subplots with 1 row and 2 columns
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))


# Plot the first scatter plot in the first subplot
ax1.scatter(molecules1, series1_energies, label='FeAB', color='black', marker='x', s=100)
ax1.scatter(molecules1, series2_energies, label='FeAB* enol', facecolors='none', marker='o', edgecolors='black', s=125)
#ax1.scatter(np.arange(len(molecules1)), series3_energies, label='FeAB* keto', color='orange', marker='_', s=600)
ax1.set_title('Relative Energies of Aerobactin structures in kcal/mol')
#ax1.set_xlabel('Conformation')
ax1.set_ylabel('Relative Energy (kcal/mol)')
ax1.set_xticks([])
ax1.set_xticklabels([])
ax1.set_ylim(-1.0, 14.0)
ax1.legend()


cell_text1 = np.array([series1_energies,
             series2_energies])

rows1 = ('FeAB', 'FeAB*')

the_table = ax1.table(cellText=cell_text1,
                      rowLabels=rows1,
                      colLabels=molecules1,
                      loc='bottom',
                      cellLoc = 'center', rowLoc = 'center')

the_table.auto_set_font_size(False)
the_table.set_fontsize(9)

# Plot the second scatter plot in the second subplot
ax2.scatter(molecules2, series3_energies, label='FePB', color='black', marker='x', s=100)
ax2.scatter(molecules2, series4_energies, label='FePB* enol', facecolors='none', marker='o', edgecolors='black', s=125)
#axs[1].scatter(np.arange(len(molecules2)), series5_energies, label='FePB* keto', color='orange', marker='_', s=600)
ax2.set_title('Relative Energies of Petroactin structures in kcal/mol')
#ax2.set_xlabel('Conformations')
ax2.set_ylabel('Relative Energy (kcal/mol)')
ax2.set_xticks([])
ax2.set_xticklabels([])
ax2.set_ylim(-1.0, 7.0)
ax2.legend()

cell_text2 = np.array([series3_energies,
             series4_energies])

rows2 = ('FePB', 'FePB*')

the_table = ax2.table(cellText=cell_text2,
                      rowLabels=rows2,
                      colLabels=molecules2,
                      loc='bottom',
                      cellLoc = 'center', rowLoc = 'center')

the_table.auto_set_font_size(False)
the_table.set_fontsize(9)

plt.subplots_adjust(left=0.2, bottom=0.2)

# Adjust layout for better spacing
plt.tight_layout()

# Show the plot
plt.show()
