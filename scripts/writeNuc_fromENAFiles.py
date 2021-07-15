import gzip
from Bio import SeqIO
import os,sys


DIR_ENA_rel142=sys.argv[1]

filesList=[x for x in os.listdir(DIR_ENA_rel142) if x.split('.')[-2:]==['dat', 'gz']]

#f_GB=sys.argv[1]

def writeNucFromGB(f_GB, DIR_ENA_rel142):
    f_handle=gzip.open(DIR_ENA_rel142+f_GB, 'rt')
    f=SeqIO.parse(f_handle, 'embl')
    f_out=open(f_GB.split('.')[0]+'.Nuc.short.fasta', 'w')
    #f_out_Taxonomy=open(f_GB+'.Taxonomy.txt', 'w')
    for rec in f:
        if len(str(rec.seq))<=50000:
            f_out.write('>'+rec.id+'\n'+str(rec.seq)+'\n')
        else:
            pass
    f_out.close()
    
for f_GB in filesList:
    print f_GB
    writeNucFromGB(f_GB, DIR_ENA_rel142)
