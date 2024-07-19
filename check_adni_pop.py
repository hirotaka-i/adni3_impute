import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.lines as mlines
df=pd.read_csv('work/pca.eigenvec', delim_whitespace=True, names=['FID', 'IID'] + [f'PC{i}' for i in range(1,5)])
t = pd.read_csv('work/ADNI3_1st_set.fam', delim_whitespace=True, header=None)
adni3_set1 = t[1].to_list()
t = pd.read_csv('work/ADNI3_2nd_set.fam', delim_whitespace=True, header=None)
adni3_set2 = t[1].to_list()
df.head()
df['Population']=[t.split('_')[2] if t.startswith('NA') else 
                  'ADNI3_set1' if t in adni3_set1 else
                  'ADNI3_set2' if t in adni3_set2 else
                  'Not Assinged' for t in df.IID]
df['Population'].value_counts()
colors = ['red', 'green', 'blue', 'cyan', 'magenta', 
        'black', 'purple', 'orange', 'pink', 'brown','yellow',
        'gray', 'navy', 'maroon', 'violet', 'turquoise', 'lime',
        'teal', 'indigo', 'coral', 'gold', 'darkred', 'darkgreen',
        'darkblue', 'lightgray', 'darkgray', 'beige', 'lightgreen', 'lightblue']

# figure
fig, ax = plt.subplots(1, 1, figsize=(10, 8))
legend_lines = []

for i, (j, group) in enumerate(df.groupby('Population')):
    if j == 'ADNI3_set1':
        ax.scatter(x='PC1', y='PC2', color=colors[i], label=j, s=10, alpha=1, data=group)
    elif j == 'ADNI3_set2':
        ax.scatter(x='PC1', y='PC2', color=colors[i], label=j, s=10, alpha=1, data=group)
    else:
        sns.kdeplot(x='PC1', y='PC2', n_levels=4, ax=ax, color=colors[i], alpha=0.8, data=group)
        legend_lines.append(mlines.Line2D([], [], color=colors[i], label=j))

# Add custom legend
ax.legend(handles=legend_lines, title='Population', loc='upper right')
plt.title('Whole cohort (with HAPMAP)')
fig.savefig('population_plot_HAPMAP.png', dpi=300, format='png')


dfpcaVal = pd.read_csv(f'work/pca.eigenval', header=None)
dfpcaVal.plot.line()
plt.title('Scree plot (with HAPMAP)')
plt.savefig(f'pca_scree_plot_HAPMAP.png')


# figure 2
fig, ax = plt.subplots(1, 1, figsize=(10, 8))
legend_lines = []

for i, (j, group) in enumerate(df.groupby('Population')):
    if j == 'ADNI3_set1':
        ax.scatter(x='PC1', y='PC3', color=colors[i], label=j, s=10, alpha=1, data=group)
    elif j == 'ADNI3_set2':
        ax.scatter(x='PC1', y='PC3', color=colors[i], label=j, s=10, alpha=1, data=group)
    else:
        sns.kdeplot(x='PC1', y='PC3', n_levels=4, ax=ax, color=colors[i], alpha=0.8, data=group)
        legend_lines.append(mlines.Line2D([], [], color=colors[i], label=j))

# Add custom legend
ax.legend(handles=legend_lines, title='Population', loc='upper right')
plt.title('Whole cohort (with HAPMAP)')
fig.savefig('population_plot_HAPMAP2.png', dpi=300, format='png')
