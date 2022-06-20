# sort signature matrix base on specific order
import pandas as pd

df = pd.read_csv('anno_temp.tsv', sep='\t', header=None)
anno_temp = list(df.iloc[:,0])

patient, timepoint, relapse = [],[], []
for item in anno_temp:
    info = item.split('-')
    patient.append(info[1][1:])
    if info[2] == 'BL':
        timepoint.append('Baseline')
    elif info[2] == '2W':
        timepoint.append('2 weeks')
    elif info[2] == '3M':
        timepoint.append('3 months')
    elif info[2] == '6M':
        timepoint.append('6 months')
    elif info[2] == '12M':
        timepoint.append('12 months')
    elif info[2] == 'MeDIP':
        timepoint.append('MeDIP')
    else:
        print('timepoint error')

    if info[1] in ['003','004','005','006','009']:
        relapse.append('Yes')
    elif info[1] in ['007','008','010']:
        relapse.append('LOF')
    else:
        relapse.append('No')

df = {}
df['Patient'] = patient
df['Timepoint'] = timepoint
df['Relapse'] = relapse

df = pd.DataFrame(df)
df['Timepoint'] = pd.Categorical(df['Timepoint'], ["Baseline","2 weeks","3 months",
                                 "6 months", "12 months", "MeDIP"])
df['Relapse'] = pd.Categorical(df['Relapse'], ["Yes","No","LOF"])

# Relapse: Patient 3,4,5,6,9

df1_sort = df.iloc[:-7,:].sort_values(['Relapse','Patient','Timepoint'], axis=0)

df_sort = pd.concat((df1_sort, df.iloc[-7:,:]))  #for matrix with MeDIP data

df_sort.to_csv('metadata.tsv', index=False, sep='\t')

data = pd.read_csv('uveal.tsv', sep='\t')
data = data.drop(['bin_chr','bin_start','bin_end','i.bin_start','i.bin_end'],axis=1)
data = data.iloc[:,:46]

ind = df1_sort.index ##for matrix without MeDIP data

data = data.iloc[:,ind]
data.to_csv('uveal_sort.tsv', index=False, sep='\t')
