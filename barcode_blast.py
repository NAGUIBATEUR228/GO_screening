import os
import numpy as np
import pandas as pd
import re
from datetime import datetime
from sys import argv
if len(argv)>1:
    need_to_BLAST=bool(int(argv[1]))
else:
    need_to_BLAST=True

print('directory with script: C:\\Users\\zokmi\\Desktop\\study\\coursework\\')
directory=input('enter directory with data: ')
path=f'C:\\Users\\zokmi\\Desktop\\study\\coursework\\{directory}\\'
print(str(datetime.now())+' '+path)
#import reference barcode table, made by R scirpt
ref=pd.read_csv(f'{path}..\\reference.txt')
print(str(datetime.now()))
print(ref.head(5))
# refcd = ref.groupby('Confirmed_deletion').agg(
#     UPTAG_notes= ('UPTAG_notes',lambda x: '|'.join(x[~pd.isna(x)].drop_duplicates())),
#     UPTAG_seqs=('UPTAG_seqs',lambda x:'|'.join(x[~pd.isna(x)].drop_duplicates())
# )
#             ).reset_index()
# print(refcd.head(5))

dirs=list()
for i in os.listdir(path):#list of directory and file names
    if os.path.isdir(os.path.join(path, i)):
        dirs.append(i)
print(str(datetime.now())+' '+str(dirs))#directories with fastq files
if need_to_BLAST:
    #making fasta file for BLAST database
    f=open('./ref_db.fasta',"w")
    for i in range(len(ref['UPTAG_seqs'])):
        i1=ref['UPTAG_seqs'][i]
        i2=i1+'-refseqname'
        f.write(f'>{i2}\n{i1}\n')
    f.close()
    #BLAST command, making nucleotide database 'ref'
    os.system('makeblastdb -in ./ref_db.fasta -dbtype nucl -out ref')

sumblast = pd.DataFrame({
    'exp': [], 'all_count': [], 'm_count': [], 'percentage': [], 'qual': []
})
sumblastdm = pd.DataFrame({
    'exp': [], 'all_count': [], 'm_count': [], 'percentage': [], 'qual': []
})

barcodes=pd.Series(dtype=object)
for i in dirs:
    p=path+i+'\\artem'
    content = os.listdir(p)
    files = []
    for file in content:
        if os.path.isfile(os.path.join(p, file)) and (file.endswith('not_matched.csv') or file.endswith('not_matched_dm.csv')):
            files.append(file)
    #parsing not_matched files in folders
    for j in files:
        nb=pd.read_csv(f'{p}\\{j}')#not_matched table
        #print(len(nb[~nb['barcode'].isin(barcodes)]))
        barcodes=pd.Series(pd.concat([barcodes,nb['barcode']],ignore_index=True).unique())

print(len(barcodes))
if need_to_BLAST:
    f=open(f'{path}\\all.fasta',"w")#making fasta with BLAST queries
    for k in barcodes:
        f.write(f'>{k}-unknownseqname\n{k}\n')
    f.close()
    os.system(f'cd {path}')
    print(str(datetime.now())+' '+f'{p}\\all.fasta')
    #BLAST command, searching not_matched sequences in database ref. output file end with blastout3.txt. e-value 10, alignment initiating word size 6, search on the same strand, finds only the best hit. outfile contains different information about alignment in tsv format.
    os.system(f'blastn -query {path}\\all.fasta -db ref -out {path}\\all_blastout.txt -evalue 10 -word_size 6 -strand plus -max_target_seqs 1 -max_hsps 1 -outfmt \"6 qacc qlen sacc slen length nident evalue qstart qend sstart send\" -num_threads 10')#-num_threads 8
    #subtracting -unknownseqname in outfile.
    with open(f'{path}\\all_blastout.txt','r') as f:
        old_data=f.read()
    new_data = 'qacc\tqlen\tsacc\tslen\tlength\tnident\tevalue\tqstart\tqend\tsstart\tsend\n'+old_data.replace('-unknownseqname\t','\t')
    with open(f'{path}\\all_blastout.txt','w') as f:
        f.write(new_data)
    with open(f'{path}\\all_blastout.txt','r') as f:
        old_data=f.read()
    new_data = old_data.replace('-refseqname\t','\t')
    with open(f'{path}\\all_blastout.txt','w') as f:
        f.write(new_data)

check=pd.read_csv(f'{path}\\all_blastout.txt',sep='\t')
print(check)
#filtering distance between sequences
check['dist']=(check['slen']+check['qlen']-check['length']-check['nident'])
check=check[check['dist']<=2]
print(len(check))
print(str(datetime.now())+' BLASTed')
print(str(datetime.now())+' adding reference and not_matched data')

barcode=pd.DataFrame({'barcode':barcodes})
full=pd.merge(barcode,check,left_on='barcode',right_on='qacc',how='left')
full=pd.merge(full,ref,left_on='sacc',right_on='UPTAG_seqs',how='left')
full=full[~pd.isna(full['qacc'])].reset_index()[['barcode','Confirmed_deletion','UPTAG_seqs','UPTAG_notes']]
print(full)


for i in dirs:
    p=path+i+'\\artem'
    content = os.listdir(p)
    files = []
    for file in content:
        if os.path.isfile(os.path.join(p, file)) and file.endswith('not_matched.csv'):
            files.append(file)
    #parsing not_matched files in folders
    for j in files:
        nb=pd.read_csv(f'{p}\\{j}')#not_matched table
        name=re.split(r'_not_matched\.csv',j)[0]
        f=pd.merge(full,nb,left_on='barcode',right_on='barcode',how='right')
        m=f[~pd.isna(f['Confirmed_deletion'])].reset_index()[['barcode',    'Confirmed_deletion',   'UPTAG_seqs', 'UPTAG_notes',    'n',    'count',    'qual']]
        #print(str(datetime.now()))
        m = m.groupby('Confirmed_deletion').agg(
                barcode= ('barcode',lambda x: '|'.join(x[~pd.isna(x)].unique())),
                n=('n', 'sum'),
                count=('count', 'sum'),
                notes=('UPTAG_notes', lambda x: '|'.join(x[~pd.isna(x)].unique())),
                original_barcode=('UPTAG_seqs',lambda x: '|'.join(x[~pd.isna(x)].unique()))).reset_index()
        # def summary1(x):
        #     result = {
        #         'barcode': '|'.join(x[~pd.isna(x['barcode'])]['barcode'].drop_duplicates()),
        #         'n': x['n'].sum(),
        #         'count': x['count'].sum(),
        #         'notes': '|'.join(x[~pd.isna(x['UPTAG_notes'])]['UPTAG_notes'].drop_duplicates()) ,
        #         'original_barcode': '|'.join(x[~pd.isna(x['UPTAG_seqs'])]['UPTAG_seqs'].drop_duplicates())
        #     }
        #     return pd.Series(result)
        # m=m.groupby('Confirmed_deletion').apply(summary1)
        m=m.sort_values (by = ['count'], ascending = [ False ])
        print(str(datetime.now())+' '+f'{p}\\{name}_blasted.csv')
        m.to_csv (f'{p}\\{name}_blasted.csv', index= False )

        sumblast.loc[ len(sumblast.index )] = [
        name,
        round(nb['count'].sum()),
        round(m['count'].sum()),
        round(m['count'].sum())/round(nb['count'].sum()),
        nb['qual'].mean().item()
        ]

print(str(datetime.now()))
print(str(sumblast['percentage'].tolist())+' '+'DONE')
sumblast.to_csv (f'{path}sumblast.csv', index= False )

for i in dirs:
    p=path+i+'\\artem'
    content = os.listdir(p)
    files = []
    for file in content:
        if os.path.isfile(os.path.join(p, file)) and file.endswith('not_matched_dm.csv'):
            files.append(file)
    #parsing not_matched files in folders
    for j in files:
        nb=pd.read_csv(f'{p}\\{j}')#not_matched table
        name=re.split(r'_not_matched_dm\.csv',j)[0]
        f=pd.merge(full,nb,left_on='barcode',right_on='barcode',how='right')
        m=f[~pd.isna(f['Confirmed_deletion'])].reset_index()[['barcode',    'Confirmed_deletion',   'UPTAG_seqs', 'UPTAG_notes',    'n',    'count',    'qual']]
        m = m.groupby('Confirmed_deletion').agg(
                barcode= ('barcode',lambda x: '|'.join(x[~pd.isna(x)].unique())),
                n=('n', 'sum'),
                count=('count', 'sum'),
                notes=('UPTAG_notes', lambda x: '|'.join(x[~pd.isna(x)].unique())),
                original_barcode=('UPTAG_seqs',lambda x: '|'.join(x[~pd.isna(x)].unique()))).reset_index()
        # def summary1(x):
        #     result = {
        #         'barcode': '|'.join(x[~pd.isna(x['barcode'])]['barcode'].drop_duplicates()),
        #         'n': x['n'].sum(),
        #         'count': x['count'].sum(),
        #         'notes': '|'.join(x[~pd.isna(x['UPTAG_notes'])]['UPTAG_notes'].drop_duplicates()) ,
        #         'original_barcode': '|'.join(x[~pd.isna(x['UPTAG_seqs'])]['UPTAG_seqs'].drop_duplicates())
        #     }
        #     return pd.Series(result)
        # m=m.groupby('Confirmed_deletion').apply(summary1)
        m=m.sort_values (by = ['count'], ascending = [ False ])
        print(str(datetime.now())+' '+f'{p}\\{name}_blasted_dm.csv')
        m.to_csv (f'{p}\\{name}_blasted_dm.csv', index= False )

        sumblastdm.loc[ len(sumblastdm.index )] = [
        name,
        round(nb['count'].sum()),
        round(m['count'].sum()),
        round(m['count'].sum())/round(nb['count'].sum()),
        nb['qual'].mean().item()
        ]

print(str(datetime.now()))
print(str(sumblastdm['percentage'].tolist())+' '+'DONE')
sumblastdm.to_csv (f'{path}sumblastdm.csv', index= False )
