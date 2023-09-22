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


# def rc(sequence):   
#     complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
#     vec=np.vectorize(complement.get) 
#     seq_array = np.array(list(sequence))
#     rev_seq_array = seq_array[::-1]
#     comp_seq_array = vec(rev_seq_array)
#     rev_comp_seq = ''.join(comp_seq_array)
#     return rev_comp_seq
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
# refcd = ref.groupby('Confirmed_deletion').agg(
#     UPTAG_notes= ('UPTAG_notes',lambda x: '|'.join(x[~pd.isna(x)].drop_duplicates())),
#     UPTAG_seqs=('UPTAG_seqs',lambda x:'|'.join(x[~pd.isna(x)].drop_duplicates())
# )
#             ).reset_index()

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
        #parsing not_matched files in folders
        for j in files:
            nb=pd.read_csv(f'{p}\\{j}')#not_matched table
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

    #BLAST command, searching not_matched sequences in database ref. output file end with blastout.txt. e-value 10, alignment initiating word size 6, search on the same strand, finds only the best hit. outfile contains different information about alignment in tsv format.

    os.system(f'blastn -query dark_matter/{directory}_dm.fasta -db dark_matter/uref -out dark_matter/out/{directory}_blastout.txt -evalue 0.001 -word_size 6 -strand plus -max_target_seqs 2 -outfmt \"6 qacc qlen sacc slen length nident evalue qstart qend sstart send\"')#-num_threads 8

    print(str(datetime.now())+' BLASTed')
    #subtracting -unknownseqname in outfile.
    with open(f'{dmp}\\out\\{directory}_blastout.txt','r') as f:
        old_data=f.read()
    new_data = 'qacc\tqlen\tsacc\tslen\tlength\tnident\tevalue\tqstart\tqend\tsstart\tsend\n'+old_data.replace('-unknownseqname\t','\t')
    with open(f'{dmp}\\out\\{directory}_blastout.txt','w') as f:
        f.write(new_data)

#tables:output_dm_count (matched seqs),not_matched_dm (nm seqs)
#join three tables: blast_seqs, reference, not_barcoded. make dict with orfs; filter where nan in orf

check=pd.read_csv(f'{dmp}\\out\\{directory}_blastout.txt',sep='\t')
check=check[check.groupby(['qacc'])['sacc'].transform('count')>=2]
mask=check.groupby(['qacc','sacc'])['sacc'].transform('count')!=1
check1=check[~mask]
check=check[mask]
#chooses the best alignment of a subject if two by max length
print(check, len(check))
print(check1)
print(check.groupby(['qacc'])['sacc'].value_counts())

def ch(check):
    a=set(check.groupby(['qacc'])['sacc'].value_counts().values)
    print(a,max(a))
    return max(a)!=1

changed=True
if changed:
    print(str(datetime.now())+' 1')
    if ch(check):
        print(str(datetime.now()))
        #check=check.loc[(check.groupby(['qacc','sacc']))['length'].idxmax()]
        check=check[(check.groupby(['qacc','sacc'])['length'].transform('max'))==check['length']]
    else:
        changed=False
if changed:
    print(str(datetime.now())+' 2')
    if ch(check):
        print(str(datetime.now()) )
        sacc=check['sacc'].values
        conditions = [
            sacc == 'u1',
            sacc == 'u2',
            sacc == 'u2_rc',
            sacc == 'u1_rc'
        ]

        values = [
            check['send'].values,
            -check['sstart'].values,
            check['send'].values,
            -check['sstart'].values
        ]
        print(str(datetime.now())+' select')
        check['tof'] = np.select(conditions, values, default=None)
        print(str(datetime.now())+' tof')
        #check=check.loc[(check.groupby(['qacc','sacc']))['tof'].idxmax()]
        check=check[(check.groupby(['qacc','sacc'])['tof'].transform('max'))==check['tof']]
    else:
        changed=False

if changed:
    print(str(datetime.now())+' 3')
    if ch(check):
        print(str(datetime.now()) )
        check=check[(check.groupby(['qacc','sacc'])['nident'].transform('max'))==check['nident']] 
    else:
        changed=False
if changed:
    print(str(datetime.now())+' 4')
    if ch(check):
        print(str(datetime.now()) )
        check=check[(check.groupby(['qacc','sacc'])['evalue'].transform('min'))==check['evalue']]
    else:
        changed=False
if changed:
    print(str(datetime.now())+' 5')
    if ch(check):
        print(str(datetime.now()) )
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
print(str(datetime.now())+' 6')
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
print(str(datetime.now())+' 7')
print(str(datetime.now())+' pivot')
w=pd.pivot(data=check, values='ind', index='qacc',columns='sacc').reset_index()
w = w.where((pd.notnull(w)), None)
u=w.columns.values.tolist()
print(str(datetime.now())+' normality')
if sum([x in u for x in ud.keys()])==4:
    mask=(~pd.isna(w['u1']) & ~pd.isna(w['u2']) & (w['u1']<=w['u2'])) | (~pd.isna(w['u1_rc']) & ~pd.isna(w['u2_rc']) & (w['u2_rc']<=w['u1_rc']))
    print('yez')
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
# refmerge=pd.merge(ref[['UPTAG_seqs','Confirmed_deletion']],refcd,left_on='Confirmed_deletion',right_on='Confirmed_deletion',how='left')
# print(len(ref)==len(refmerge))#true
for i in dirs:
    p=path+i+'\\artem'
    content = os.listdir(p)
    files = []
    for file in content:
        if os.path.isfile(os.path.join(p, file)) and file.endswith('not_barcoded_raw_qual_count.csv'):
            files.append(file)
    #parsing not_matched files in folders
    for j in files:
        name=re.split(r'_not_barcoded_raw_qual_count\.csv',j)[0]
        nb=pd.read_csv(f'{p}\\{j}')#not_barcoded table
        print(str(datetime.now())+' indicies')
        inds=nb['seq'][nb['seq'].isin(set(w.index))]
        print(str(datetime.now())+' to merge')
        to_merge=w.loc[inds, :].reset_index()
        print(str(datetime.now())+' '+j)
        full=pd.merge(to_merge,nb,left_on='qacc',right_on='seq',how='outer')
        full=pd.merge(full,ref,left_on='barcode',right_on='UPTAG_seqs',how='left')
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


        # print(str(datetime.now())+' aggregation')
        # full['qual']=full['qual']*full['count']
        # full['qual']=full.groupby('barcode')['qual'].transform('sum')
        # f = full.groupby('barcode').agg(
        #         count=('count', 'sum'),
        #         n=('count', 'count'),
        #         Confirmed_deletion=('Confirmed_deletion', lambda x: '|'.join(x[~pd.isna(x)].unique())),
        #         notes=('UPTAG_notes', lambda x: '|'.join(x[~pd.isna(x)].unique())),
        #         qual=('qual', 'sum')
        #     ).reset_index()
        # print(f['qual'].mean())
        # f['qual']=f['qual']/f['count']
        # f = f.sort_values (by = ['count'], ascending = [ False ])
        # print(str(datetime.now())+' split')
        # nm=f[f['Confirmed_deletion']==''][['barcode','n','count','qual']]
        # print(nm['qual'].mean())

        # nm.to_csv (f'{p}\\{name}_not_matched_dm.csv', index= False )

        # m=f[f['Confirmed_deletion']!=''][['Confirmed_deletion','barcode','n','count','notes']]
        # m = m.groupby('Confirmed_deletion').agg(
        #         barcode= ('barcode',lambda x: '|'.join(x[~pd.isna(x)].unique())),
        #         n=('n', lambda x:int(x.sum())),
        #         count=('count', lambda x: int(x.sum())),
        #         notes=('notes', lambda x: '|'.join(x[~pd.isna(x)].unique()))
        #     ).reset_index()
        # print(str(datetime.now())+' '+f'{p}\\{name}_output_dm_count.csv')
        # m=m.sort_values (by = ['count'], ascending = [ False ])
        # m.to_csv (f'{p}\\{name}_output_dm_count.csv', index= False )
        # exess=round(full[pd.isna(full['qacc'])]['count'].sum())
        # exessu=round(len(full[pd.isna(full['qacc'])]['seq'].unique()))

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