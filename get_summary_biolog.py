import os
import numpy as np
import pandas as pd

def dissect_name(fname):
    fname = fname.split('_')
    strain = fname[0]
    plate = fname[1]
    time = int(fname[2][1:-5])
    return strain, plate, time


def merge_readings_strain(strain):
    list_fname = os.listdir('data/biolog/reads/' + strain)
    list_fname = [fname for fname in list_fname if fname[0] != '~']

    df = pd.DataFrame(columns = ['ID', 'strain', 'plate', 'time', 
        'well', 'wavelength', 'od'])

    
    for fname in list_fname:
        print('  processing ' + fname + ' ..')
        strain, plate, time = dissect_name(fname)
        read = pd.read_excel('data/biolog/reads/' + strain + '/' + fname, 
                             skiprows = 7, 
                             nrows = 96, 
                             usecols = 'C:E')
        for wavelength in [450, 600, 750]:
            for i in range(96):
                val = read.loc[i, wavelength]
                well = chr(ord('A') + i // 12) + str(i % 12 + 1)
                df.loc[len(df)] = [plate+'_'+well, strain, plate, time, 
                    well, wavelength, val]

    return df


def plot_readings(df, strain, colname='od'):
    import matplotlib.pyplot as plt
    import seaborn as sns
    sns.set_style('whitegrid')
    wavelengths = [450, 600, 750]
    # plt.figure(figsize = (12, 24))
    fig, axes = plt.subplots(1, 3, figsize=(24, 6))
    # sns.set_context('paper', font_scale = 1.5)
    for row in range(3):
        df_temp = df[df['wavelength'] == wavelengths[row]]
        sns.lineplot(x = 'time', y = colname, hue = 'ID', data = df_temp, 
                     ax = axes[row], legend=False)
        axes[row].set_title(str(wavelengths[row])+' nm, ' + strain)
    plt.savefig('data/biolog/' + strain + '_' + colname + '.png')
    plt.show()


def get_fold_rate(df, strain):
    import copy
    df_out = copy.deepcopy(df)
    df_out['fold'] = np.nan
    df_out['rate'] = np.nan
    list_time = sorted(df['time'].unique())

    for wavelength in [450, 600, 750]:
        for id in df_out['ID'].unique():
            index_id_wavelength = (df_out['ID'] == id) \
                & (df_out['wavelength'] == wavelength)
            # update fold column
            index = index_id_wavelength
            index0 = index_id_wavelength & (df_out['time'] == 0)
            df_out.loc[index, 'fold'] = df_out.loc[index, 'od'] \
                / df_out.loc[index0,'od'].values[0]
            
            # update rate column
            for i in range(len(list_time) - 1):
                index_curr = index_id_wavelength \
                    & (df_out['time'] == list_time[i])
                index_next = index_id_wavelength \
                    & (df_out['time'] == list_time[i+1])
                df_out.loc[index_curr, 'rate'] = np.log(
                    (df_out.loc[index_next, 'od'].values[0] \
                    / df_out.loc[index_curr, 'od'].values[0])) \
                    / (list_time[i+1] - list_time[i])

    return df_out


def get_stats(df, strain):
    df_out = pd.DataFrame(columns = ['ID', 'strain', 'row', 'column', 'plate', 'well', 
        'wavelength', 'time_max', 'rate_max', 'rate_mean', 'fold_max'])
    list_time = sorted(df['time'].unique())
    print('  summarizing ' + strain + ' ..')

    for wavelength in [600]:
        for id in df['ID'].unique():
            id_split = id.split('_')
            row = id_split[1][0]
            column = int(id_split[1][1:])
            
            index_id_wavelength = (df['ID'] == id) \
                & (df['wavelength'] == wavelength)

            # maximal rate
            rate_max = df.loc[index_id_wavelength, 'rate'].max()
            index_max = np.argmax(df.loc[index_id_wavelength, 'rate'])
            time_max = df.loc[index_id_wavelength, 'time'].values[index_max]

            # average rate
            index0 = index_id_wavelength \
                & (df['time'] == min(list_time))
            index1 = index_id_wavelength \
                & (df['time'] == max(list_time))
            rate_avg = np.log(
                (df.loc[index1, 'od'].values[0] \
                / df.loc[index0, 'od'].values[0])) \
                / (max(list_time) - min(list_time))
            
            # maximal fold
            fold_max = df.loc[index_id_wavelength, 'fold'].max()

            df_out.loc[len(df_out)] = [id, strain, row, column,
                df.loc[index_id_wavelength, 'plate'].values[0], 
                df.loc[index_id_wavelength, 'well'].values[0], 
                wavelength, time_max, rate_max, rate_avg, fold_max]

            # sort by row and column
            df_out.sort_values(by=['plate', 'column', 'row'], inplace=True)

    return df_out


def filter_substrate(df_stat, strain, fold_threshold=1.05, wavelength=600):
    print('  filtering ' + strain + ' ..')
    # default filter
    df_out = df_stat[(df_stat['wavelength'] == wavelength)]
    df_out.loc[:,'filter_bin'] = np.zeros(len(df_out), dtype = bool)
    index = (df_stat['rate_mean'] > 0) & (df_stat['fold_max'] > 1)
    index_out = np.zeros(len(df_out), dtype = bool)

    # filter by comparing to control (A1)
    index_A1 = (df_out['well'] == 'A1')
    for plate in ['PM1', 'PM2A']:
        index_plate = (df_out['plate'] == plate)
        fold_max_A1 = df_out.loc[index_plate & index_A1, 'fold_max'].values[0]
        index_out = (df_out['fold_max'] > fold_max_A1*fold_threshold) \
            & index_plate | index_out

    # filter by mean growth rate
    index_out = (df_out['rate_mean'] > 0) & index_out

    df_out.loc[index_out, 'filter_bin'] = True
    df_out.loc[:,'filter_bin'] = df_out.loc[:,'filter_bin'].astype(int)
    
    return df_out


# obtain list of strains
list_strain = os.listdir('data/biolog/reads')
list_strain.remove('.DS_Store')
list_strain.remove('B3')
list_strain.remove('P2B')
list_strain.remove('19DW')
# list_strain.remove('summary_biolog_dye.xlsx')
# list_strain.remove('~$summary_biolog_dye.xlsx')

df = pd.DataFrame(columns = ['strain', 'plate', 'time', 'well', 
    'wavelength', 'od'])
fold_threshold_dict = {
    '3-2': 2,
    '4BL': 1.2,
    '4D': 1.2,
    '11-3': 1.2,
    '13A': 1.2,
    '13C1': 1.2,
    'ARW1R1': 1.2,
    'ARW1T': 1.2,
    'ARW7G5W': 1.2,
    'ARW7G5Y1': 1.2,
    'EAB7W2': 1.2
    }

# for strain in list_strain:
for strain in ['3-2']:
    print('reading ' + strain + ' ..')
    df = merge_readings_strain(strain)
    df = get_fold_rate(df, strain)
    df.to_csv('data/biolog/reads_merged/' + strain + '.csv', index = False)
    # df = pd.read_csv('data/biolog/reads_merged/' + strain + '.csv')

    df_stat = get_stats(df, strain)
    df_stat.to_csv('data/biolog/summary/stats/' + strain + '.csv', index = False)
    # df_stat = pd.read_csv('data/biolog/summary/stats/' + strain + '.csv')
    
    df_stat_filtered = filter_substrate(df_stat, strain, fold_threshold=fold_threshold_dict[strain])
    print(sum(df_stat_filtered['filter_bin']))
    df_stat_filtered.to_csv('data/biolog/summary/' + strain + '.csv', index = False)

    # plot_readings(df, strain)
    # plot_readings(df, strain, colname='fold')
    # plot_readings(df, strain, colname='rate')
print('done')
