import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.lines as mlines
df=pd.read_csv('pca.eigenvec', delim_whitespace=True, names=['FID', 'IID'] + [f'PC{i}' for i in range(1,11)])
t = pd.read_csv('ADNI3_1st_set.fam', delim_whitespace=True, header=None)
adni3_set1 = t[1].to_list()
t = pd.read_csv('ADNI3_2nd_set.fam', delim_whitespace=True, header=None)
adni3_set2 = t[1].to_list()
df['Population']=[t.split('_')[2] if t.startswith('NA') else 
                  'ADNI3_set1' if t in adni3_set1 else
                  'ADNI3_set2' if t in adni3_set2 else
                  'Not Assinged' for t in df.IID]
print(df['Population'].value_counts())
colors = ['red', 'green', 'blue', 'cyan', 'magenta', 
        'black', 'purple', 'orange', 'pink', 'brown','yellow',
        'gray', 'navy', 'maroon', 'violet', 'turquoise', 'lime',
        'teal', 'indigo', 'coral', 'gold', 'darkred', 'darkgreen',
        'darkblue', 'lightgray', 'darkgray', 'beige', 'lightgreen', 'lightblue']


def plot_pcs(ax, pcx, pcy, data, colors):
    legend_handles = []
    for i, (pop, group) in enumerate(data.groupby('Population')):
        if pop in ['ADNI3_set1', 'ADNI3_set2']:
            scatter = ax.scatter(x=pcx, y=pcy, color=colors[i], label=pop, s=10, alpha=1, data=group)
            legend_handles.append(scatter)
        else:
            sns.kdeplot(x=pcx, y=pcy, n_levels=4, ax=ax, color=colors[i], alpha=0.8, data=group)
            legend_handles.append(mlines.Line2D([], [], color=colors[i], label=pop))
    ax.set_xlabel(pcx)
    ax.set_ylabel(pcy)
    ax.legend(handles=legend_handles, title='Population', loc='upper right')
    ax.set_title(f'{pcx} vs {pcy}')

fig, axes = plt.subplots(2, 2, figsize=(15, 15))
# Call the function for each subplot
# Call the function for each subplot
plot_pcs(axes[0, 0], 'PC1', 'PC2', df, colors)
plot_pcs(axes[0, 1], 'PC3', 'PC2', df, colors)
plot_pcs(axes[1, 0], 'PC1', 'PC4', df, colors)
# axes[1,1] is screeplot
dfpcaVal = pd.read_csv(f'pca.eigenval', header=None)
dfpcaVal.plot.line(ax=axes[1, 1])
plt.tight_layout()
plt.suptitle('Whole cohort (with HAPMAP)', fontsize=16, y=1.02)  # Adjust y for title positioning above subplots
fig.savefig('population_plot_HAPMAP_combined.png', dpi=300, format='png')


# ancestry separation
dfpop = df.pivot_table(index='Population', values=['PC1', 'PC2', 'PC3', 'PC4', 'PC5'], aggfunc=['mean', 'std'])

# Get the threshold table of mean +/- 6SD
def funcThres(x):
    lwl = x['mean'] - 6 * x['std']
    hgl = x['mean'] + 6 * x['std']
    return pd.Series({'lwl':lwl, 'hgl':hgl})

thres = dfpop.apply(funcThres, axis=1)

# function to infer ancestry for "OTH"
def funcInfPop(x):
    if not "ADNI" in x.Population:
        InfPop = 'REF'
    else:
        InfPop = 'ADMIX'
        for Population in ['EUROPE', 'AFRICA', 'ASIA' ]:
            if (thres.loc[Population, 'lwl']['PC1'] < x.PC1) & \
              (x.PC1 < thres.loc[Population, 'hgl']['PC1']) & \
              (thres.loc[Population, 'lwl']['PC2'] < x.PC2) & \
              (x.PC2 < thres.loc[Population, 'hgl']['PC2']) & \
              (thres.loc[Population, 'lwl']['PC3'] < x.PC3) & \
              (x.PC3 < thres.loc[Population, 'hgl']['PC3']) & \
              (thres.loc[Population, 'lwl']['PC4'] < x.PC4) & \
              (x.PC4 < thres.loc[Population, 'hgl']['PC4'])& \
              (thres.loc[Population, 'lwl']['PC5'] < x.PC5) & \
              (x.PC5 < thres.loc[Population, 'hgl']['PC5']):
                    InfPop = Population
    return InfPop

df['InfPop'] = df.apply(funcInfPop, axis=1)
print(df.InfPop.value_counts())

# save the population
df.loc[df.InfPop!='REF', 
       ['FID', 'IID', 'InfPop'] + [f'PC{i+1}' for i in range(10)]].to_csv('genetic_ancestry_all_pca.csv', index=False)
for continent in df.InfPop.unique():
    if continent != 'REF':
        t = df.loc[df.InfPop==continent, ['FID', 'IID']]
        print(t.shape)
        t.to_csv(f'genetic_ancestry_{continent}.txt', index=False, sep='\t')