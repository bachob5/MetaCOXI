import sys

'''
This script formats the MetaCOXI_Taxonomy_Metadata.tsv file into a qiime compatible taxonomy file format

#Execution:
python2.7 formatTaxonomy4Qiime.py MetaCOXI_Taxonomy_Metadata.tsv
'''

f_metadata=sys.argv[1]

def isNA(S):
    if S=='NA':
        return ''
    else:
        return S

def formatT4Qiime(fMeta):
    f=open(fMeta,'r')
    f_out=open(fMeta.split('.')[0]+'_qiimeFormat.tsv', 'w')
    f.readline()
    for line in f:
        l=line.strip().split('\t')
        Accession=l[0]
        kingdom='k__'+isNA(l[2])
        phylum='p__'+isNA(l[4])
        Class='c__'+isNA(l[6])
        Order='o__'+isNA(l[8])
        Family='f__'+isNA(l[10])
        g_temp=l[12]
        Genus='g__'+isNA(g_temp)
        Sp_temp=isNA(l[14])
        Species=''
        if Sp_temp!='':
            if Sp_temp.split(' ')[0]==g_temp:
                Species+='s__'+' '.join(Sp_temp.split(' ')[1:])
            else:
                Species+='s__'+Sp_temp
        else:
            Species+=''
        TaxonomyLine=[kingdom,phylum,Class,Order,Family,Genus,Species]
        f_out.write(Accession+'\t'+'; '.join(TaxonomyLine)+'\n')
    f.close()
    f_out.close()

GetFormattedFile=formatT4Qiime(f_metadata)
