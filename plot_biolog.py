import os
import numpy as np
import pandas as pd
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

def dissect_name(fname):
    fname = fname.split('_')
    strain = fname[0]
    plate = fname[1]
    time = int(fname[2][1:-5])
    return strain, plate, time


def plot_readings_substrate(df, strain, colname='od', version='v1', ind_df=None):
    import matplotlib.pyplot as plt
    import seaborn as sns
    sns.set_style('whitegrid')
    sns.set_context('paper', font_scale = 0.4)

    ind_df_strain = pd.DataFrame(columns = ['PLATE', 'WELL', strain])
    if ind_df is None:
        ind_df_strain['PLATE'] = ['PM1'] * 96 + ['PM2A'] * 96
        ind_df_strain['WELL'] = [str(chr(ord('A') + i//12)) + str(i%12 + 1) for i in range(96)] * 2
        ind_df_strain[strain] = [0] * 192
        fnameappend = '.png'
    else:
        ind_df_strain['PLATE'] = ind_df['PLATE']
        ind_df_strain['WELL'] = ind_df['WELL']
        ind_df_strain[strain] = ind_df[strain]
        fnameappend = '_labeled.png'

    colors = ['blue', 'red']
    
    for p in ['PM1', 'PM2A']: # by plate
        print('    .. processing ' + p + '..')
        ymax = max(df[(df['wavelength'] == 600) & (df['plate']==p)]['od'])
        fig, axes = plt.subplots(8, 12, figsize=(18, 12))

        for i in range(96): # by well
            r = i // 12
            c = i % 12
            w = str(chr(ord('A') + r)) + str(c+1)
            color_ind = ind_df_strain[(ind_df_strain['PLATE'] == p) & 
                                      (ind_df_strain['WELL'] == w)][strain].values[0]
            # print('    ..... plotting well ' + w + '..')
            df_temp = df[(df['wavelength'] == 600)
                         & (df['plate'] == p) 
                         & (df['well'] == w)]
            sns.lineplot(x = 'time', y = colname, data = df_temp, 
                         ax = axes[r,c], legend=False, color=colors[color_ind])
            axes[r,c].set_title(p + '_' + w)
            axes[r,c].set_ylim([0, ymax])

        fname = strain + '_' + p + fnameappend
        
        plt.savefig('data/biolog/summary/curves/' + version + '/' + fname, dpi=450, bbox_inches='tight')
        # plt.show()


def get_read_df_v2(strain):
    df = pd.DataFrame(columns = ['ID', 'strain', 'plate', 'time', 
                                 'well', 'wavelength', 'od'])
    for p in ['PM1', 'PM2A']:
        print('    .. processing ' + p + '..')
        read_df = pd.read_csv('data/biolog/reads/v2/' + strain + '_' + p + '.csv')

        for i in range(96):
            df_temp = pd.DataFrame(columns = ['ID', 'strain', 'plate', 'time', 'well', 'wavelength', 'od'])
            w = str(chr(ord('A') + i//12)) + str(i%12 + 1)
            read_df_temp = read_df[['Time_h', w]]
            df_temp['ID'] = [p+'_'+w] * len(read_df_temp)
            df_temp['strain'] = [strain] * len(read_df_temp)
            df_temp['plate'] = [p] * len(read_df_temp)
            df_temp['time'] = read_df_temp['Time_h'].values.tolist()
            df_temp['well'] = [w] * len(read_df_temp)
            df_temp['wavelength'] = [600] * len(read_df_temp)
            df_temp['od'] = read_df_temp[w].values.tolist()
            df = pd.concat([df, df_temp], ignore_index = True)

    return df


## obtain list of strains
## list of strains 
## v1 = ['ARW1R1', '13A', '13C1', '3-2', 'ARW1T']
## v2 = ['EA2', '13M1', 'EAB7WZ']
list_strain_v1 = ['ARW1R1']
list_strain_v2 = ['N2S']
ind_df = pd.read_csv('data/biolog/' + 'summary_biolog_od.csv')

for strain in list_strain_v2:
    # print('V1: reading ' + strain + ' ..')
    # df = pd.read_csv('data/biolog/reads_merged/' + strain + '.csv')

    print('V2: reading ' + strain + ' ..')
    # df = get_read_df_v2(strain)
    # df.to_csv('data/biolog/reads_merged/' + strain + '_v2.csv', index = False)
    df = pd.read_csv('data/biolog/reads_merged/' + strain + '_v2.csv')

    print('    plotting ' + strain + ' ..')
    # plot_readings_substrate(df, strain, colname='od', version='v2', ind_df=None)
    plot_readings_substrate(df, strain, colname='od', version='v2', ind_df=ind_df)

print('done')
