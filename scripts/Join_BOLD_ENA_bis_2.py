import os,sys
from Bio import SeqIO


FinalBOLD_fastafile=sys.argv[1] 
FinalENA_fastafile=sys.argv[2] 


FinalBOLD_Taxfile=sys.argv[3] 
FinalENA_Taxfile=sys.argv[4]


def format_NCBI_TaxonomyRecs(FinalENA_Taxfile):
    D_Acc_Taxonomy={}
    f=open(FinalENA_Taxfile,'r')
    for line in f:
        l=line.strip().split('\t')
        Acc=l[0]
        A_Taxonomy=[]
        for t in l[1:]:
            A_Taxonomy.append(t.split('_')[1])
        D_Acc_Taxonomy[Acc]='\t'.join(A_Taxonomy)
    return D_Acc_Taxonomy


ENA_Acc_Taxonomy=format_NCBI_TaxonomyRecs(FinalENA_Taxfile)

def format_BOLD_TaxonomyRecs(FinalBOLD_Taxfile):
    D_Acc_Taxonomy={}
    f=open(FinalBOLD_Taxfile,'r')
    for line in f:
        l=line.strip().split('\t')
        Acc=l[0]
        A_Taxonomy=[]
        for t in l[1:]:
            t1=t.split('_')
            try:
                A_Taxonomy.append(t1[1])
            except IndexError:
                A_Taxonomy.append(t1[0])
        D_Acc_Taxonomy[Acc]='\t'.join(A_Taxonomy)
    return D_Acc_Taxonomy

BOLD_Acc_Taxonomy=format_BOLD_TaxonomyRecs(FinalBOLD_Taxfile)



dB=SeqIO.to_dict(SeqIO.parse(FinalBOLD_fastafile, 'fasta'))
dN=SeqIO.to_dict(SeqIO.parse(FinalENA_fastafile, 'fasta'))
f_out=open('BOLD_ENA_Dereplicated_COXI_Seqs.fasta','w')
f_out_Tax=open('BOLD_ENA_Dereplicated_COXI_Taxonomy.tsv','w')

L=[]
d={}
for rec in ENA_Acc_Taxonomy:
    try:
        recB=BOLD_Acc_Taxonomy[rec.split('.')[0]]
        #f_out.write(recB)
        L.append(rec)
        d[rec]=True
    except KeyError:
        try:
            f_out.write('>'+dN[rec].description+'\n'+str(dN[rec].seq)+'\n')
            f_out_Tax.write(rec+'\t'+ENA_Acc_Taxonomy[rec]+'\n')
        except KeyError:
            continue

for recB in BOLD_Acc_Taxonomy:
    try:
        f_out.write('>'+recB+'\n'+str(dB[recB].seq)+'\n')
        f_out_Tax.write(recB+'\t'+BOLD_Acc_Taxonomy[recB]+'\n')
    except KeyError:
        continue

f_out.close()
f_out_Tax.close()

print len(L)


