import os, sys
from Bio import SeqIO

CurrFolder=os.getcwd().split('/')[-1].split('_')[0]
print CurrFolder

PATH=sys.argv[1] 

filesList=[x for x in os.listdir(PATH) if x.split('_')[-1]=='SelectedGoodCOX1.fasta']

D_division_FPsynonyms=['product=BLTX405', 'product=BLTX534', 'product=BLTX384', 'Gene=COII', 'product=BLTX505', 'Gene=nad1', 'product=CO-like protein', 'product=similar to cytochrome oxidase I', 'Gene=PowCR01_000243500', 'product=Ophyra', 'Gene=cox1/cox2', 'product=placenta immunoregulatory factor PLIF', 'Gene=ND4', 'product=cellular protein AbCp-2', 'Gene=COB', 'product=unknown', 'product=BLTX358', 'product=putative cytochrome C oxidase polypeptide I', 'Gene=cox1/3', 'Gene=cox1/2']

def writeFasta_NoFPsyn(f_in, PATH):
    f_fas=SeqIO.parse(PATH+f_in, 'fasta')
    f_out=open(PATH+f_in.split('.fasta')[0]+'_rmFPsynonyms.fasta', 'w')
    f_out_finalIDs=open(PATH+f_in.split('.fasta')[0]+'_rmFPsynFinalIDs.txt', 'w')#these are the good ones
    RmIDs=[]
    for rec in f_fas:
        Des=rec.description
        Annotation=Des.split('=')[0].split()[-1]+'='+Des.split('=')[-1]
        if Annotation not in D_division_FPsynonyms:
            f_out.write('>'+rec.description+'\n'+str(rec.seq)+'\n')
            f_out_finalIDs.write(rec.id+'\n')
        else:
            RmIDs.append(rec.description)
    if len(RmIDs)>0:
        f_out_rmSyno=open(PATH+f_in.split('.fasta')[0]+'_rmFPsynonymsIDs.txt', 'w')#these are the false positive ones
        f_out_rmSyno.write('\n'.join(RmIDs)+'\n')
        f_out_rmSyno.close()
    f_out.close()
    f_out_finalIDs.close()
    
    
for fName in filesList:
    writeFasta_NoFPsyn(fName, PATH)

