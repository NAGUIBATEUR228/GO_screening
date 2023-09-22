import os
import pandas as pd 
import numpy as np
import re
from datetime import datetime
from sys import argv
if len(argv)>1:
    need_to_BLAST=bool(int(argv[1]))
else: 
    need_to_BLAST=True

u1 = "GATGTCCACGAGGTCTCT"
u2 = "CGTACGCTGCAGGTCGAC"
#reverse-complement string
def rc(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    reverse_complement = ""
    for nt in seq:
        if nt in complement:
            reverse_complement = complement[nt] + reverse_complement
    return reverse_complement

ud={'u1':u1,'u2':u2,'u1_rc':rc(u1),'u2_rc':rc(u2)}

print('directory with script: C:\\Users\\zokmi\\Desktop\\study\\coursework\\')
directory=input('enter directory with data: ')
path=f'C:\\Users\\zokmi\\Desktop\\study\\coursework\\{directory}\\'
print(str(datetime.now())+' '+path)
ref=pd.read_csv(f'{path}..\\reference.txt')

dmp=path+'..\\dark_matter'
if need_to_BLAST:
    if not os.path.exists(dmp):os.mkdir(dmp)
    if not os.path.exists(dmp+'\\out'):os.mkdir(dmp+'\\out')
    f=open(f'{dmp}\\query_u.fasta',"w")
    for i,j in ud.items():
        f.write(f'>{i}\n{j}\n')
    f.close()
    os.system(f'makeblastdb -in dark_matter/query_u.fasta -dbtype nucl -out dark_matter/uref')
    
sumdm = pd.DataFrame({
    'exp': [], 'total_count': [], 'u_count': [], 'barcoded': [], 'u_barcoded': [], 'not_barcoded': [], 'u_nb': [], 'matched': [], 'u_matched': [], 'not_matched': [], 'u_nm': []
})
dirs=list()
for i in os.listdir(path):#list of directory and file names
    if os.path.isdir(os.path.join(path, i)):
        dirs.append(i)
print(str(datetime.now())+' '+str(dirs))#directories with fastq files

if need_to_BLAST:
    seqs=pd.Series(dtype=object)
    for i in dirs:
        p=path+i+'\\artem'
        content = os.listdir(p)
        files = []
        for file in content:
            if os.path.isfile(os.path.join(p, file)) and file.endswith('not_barcoded_raw_qual_count.csv'):
                files.append(file)
        #parsing not_barcoded files in folders
        for j in files:
            nb=pd.read_csv(f'{p}\\{j}')#not_barcoded table
            #print(len(nb[~nb['barcode'].isin(seqs)]))
            seqs=pd.Series(pd.concat([seqs,nb['seq']],ignore_index=True).unique())

    print(len(seqs))
if need_to_BLAST:
    f=open(f'{dmp}\\{directory}_dm.fasta',"w")#making fasta with BLAST queries
    for k in seqs:
        f.write(f'>{k}-unknownseqname\n{k}\n')
    f.close()
    os.system(f'cd {path}')
    print(str(datetime.now())+' '+f'{dmp}\\{directory}_dm.fasta')

    #BLAST command, searching not_matched sequences in database uref

        changed=False
if changed:
    if ch(check):
        sacc=check['sacc'].values
        conditions = [
            sacc == 'u1',
            sacc == 'u2',
            sacc == 'u2_rc',
            sacc == 'u1_rc'
        ]

        values = [
            check['qend'].values,
            -check['qstart'].values,
            check['qend'].values,
            -check['qstart'].values
        ]
        check['tof'] = np.select(conditions, values, default=None)
        check=check[(check.groupby(['qacc','sacc'])['tof'].transform('min'))==check['tof']]
    else:
        changed=False
print(ch(check))
check=pd.concat([check,check1],ignore_index=True)
sacc=check['sacc'].values
conditions = [
    sacc == 'u1',
    sacc == 'u2',
    sacc == 'u2_rc',
    sacc == 'u1_rc'
]

values = [
    check['qend'].values,
    check['qstart'].values - 1,
    check['qend'].values,
    check['qstart'].values - 1
]
check['ind'] = np.select(conditions, values, default=None)
print(str(datetime.now())+' pivot')
w=pd.pivot(data=check, values='ind', index='qacc',columns='sacc').reset_index()
w = w.where((pd.notnull(w)), None)
u=w.columns.values.tolist()
if sum([x in u for x in ud.keys()])==4:
    mask=(~pd.isna(w['u1']) & ~pd.isna(w['u2']) & (w['u1']<=w['u2'])) | (~pd.isna(w['u1_rc']) & ~pd.isna(w['u2_rc']) & (w['u2_rc']<=w['u1_rc']))
elif 'u1' in u and 'u2' in u:
    mask=(~pd.isna(w['u1']) & ~pd.isna(w['u2']) & (w['u1']<=w['u2']))
elif 'u1_rc' in u and 'u2_rc' in u:
    mask=(~pd.isna(w['u1_rc']) & ~pd.isna(w['u2_rc']) & (w['u2_rc']<=w['u1_rc']))
    if 'u1' in u: w.drop(columns=['u1'],inplace=True)
    if 'u2' in u: w.drop(columns=['u2'],inplace=True)
else: 
    import sys
    print('data have no barcode')
    sys.exit()
w=w[mask]
print(str(datetime.now())+' extraction')
v=np.vectorize(lambda x,y,z:x[y:z])
vecrc=np.vectorize(rc)
w['barcode']=w['qacc']
if 'u1' in u: 
    w['barcode']=v(w['qacc'].values,w['u1'].values,w['u2'].values)
if 'u1_rc' in u: 
    w['barcode']=np.where(w['barcode'].values==w['qacc'].values,
        vecrc(v(w['qacc'].values,w['u2_rc'].values,w['u1_rc'].values)),
        w['barcode'].values)
w.set_index('qacc', inplace=True)
print(str(datetime.now())+' filtered and barcoded')
for i in dirs:
    p=path+i+'\\artem'
    content = os.listdir(p)
    files = []
    for file in content:
        if os.path.isfile(os.path.join(p, file)) and file.endswith('not_barcoded_raw_qual_count.csv'):
            files.append(file)
    #parsing not_barcoded files in folders
    for j in files:
        name=re.split(r'_not_barcoded_raw_qual_count\.csv',j)[0]
        nb=pd.read_csv(f'{p}\\{j}')#not_barcoded table
        inds=nb['seq'][nb['seq'].isin(set(w.index))]
        to_merge=w.loc[inds, :].reset_index()
        print(str(datetime.now())+' '+j)
        full=pd.merge(to_merge,nb,left_on='qacc',right_on='seq',how='outer')
        full=pd.merge(full,ref,left_on='barcode',right_on='UPTAG_seqs',how='left')
        #stats
        exess=round(full[pd.isna(full['qacc'])]['count'].sum())
        exessu=round(len(full[pd.isna(full['qacc'])]['seq'].unique()))
        
        full=full[~pd.isna(full['qacc'])]
        print(str(datetime.now())+' aggregation')
        m=full[~pd.isna(full['Confirmed_deletion'])].copy()
        m['count']=m.groupby('barcode')['count'].transform('sum')
        m['n']=m.groupby('barcode')['count'].transform('count')
        m['notes']=m['UPTAG_notes']
        m=m[['Confirmed_deletion','barcode','n','count','notes']].drop_duplicates()
        m = m.groupby('Confirmed_deletion').agg(
                barcode= ('barcode',lambda x: '|'.join(x[~pd.isna(x)].unique())),
                n=('n', 'sum'),
                count=('count', 'sum'),
                notes=('notes', lambda x: '|'.join(x[~pd.isna(x)].unique()))
            ).reset_index()
        print(str(datetime.now())+' '+f'{p}\\{name}_output_dm_count.csv')
        m=m.sort_values (by = ['count'], ascending = [ False ])
        m.to_csv (f'{p}\\{name}_output_dm_count.csv', index= False )
        
        print(str(datetime.now())+' nm')        
        nm=full[pd.isna(full['Confirmed_deletion'])].copy()
        nm['qual']=nm['qual']*nm['count']
        nm['count']=nm.groupby('barcode')['count'].transform('sum')
        nm['n']=nm.groupby('barcode')['count'].transform('count')
        nm['qual']=nm.groupby('barcode')['qual'].transform('sum')
        nm['qual']=nm['qual']/nm['count']
        nm=nm[['barcode','n','count','qual']].drop_duplicates()
        print(nm['qual'].mean())
        print(nb['qual'].mean())

        nm.to_csv (f'{p}\\{name}_not_matched_dm.csv', index= False )

        sumdm.loc[ len(sumdm.index )] = [
        name,
        round(nb['count'].sum()),
        round(len(nb['seq'].unique())),
        round(full['count'].sum()), 
        round(len(full['barcode'].unique())),
        exess,
        exessu,
        round(m['count'].sum()), 
        round(len(m['barcode'].unique())),
        round(nm['count'].sum()), 
        round(len(nm['barcode'].unique()))
        ]


#sum information in percent
sumdm['perc_barcoded']=sumdm['barcoded']/sumdm['total_count']
sumdm['perc_matched']=sumdm['matched']/sumdm['total_count']
sumdm.to_csv (f'{path}sumdm.csv', index= False )
print(str(datetime.now())+' '+'DONE')
