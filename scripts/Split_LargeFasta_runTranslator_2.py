import os
from Bio import SeqIO

filesList=[x for x in os.listdir(os.curdir) if x.split('_')[-1]=='withMitocode.fasta']

CallTranslator='python newTranslator_CExtract.py'
Frames='1,2,3,-1,-2,-3'


def partitionIndexes(l, n):
    L=[]
    try:
        for i,j in enumerate(l):
            L.append([l[i], l[i+1]])
    except IndexError:
        if l[-1]+1==n:
            pass
        else:
            L.append([l[-1],n])
    return L


Commands=[]

SmallFiles_TranslatorCode=[]
for f_in in filesList:
    GenCode=f_in.split('_')[-2]     
    L_TranslatorCode=[]
    D=SeqIO.to_dict(SeqIO.parse(f_in, 'fasta')).items()
    LinesNUM=len(D)
    if LinesNUM>50000:
        os.system('mkdir temp_'+f_in.split('.')[0])
        DIR='temp_'+f_in.split('.')[0]+'/'
        l=range(0, LinesNUM+1, 50000)
        Partitions=partitionIndexes(l, LinesNUM)
        for p in Partitions:
            F_name=str(p[0])+'_'+f_in
            L_TranslatorCode.append(CallTranslator + ' ' + DIR +F_name + ' ' + GenCode + ' ' + Frames)
            f_out=open(DIR+F_name, 'w')
            sub_Items=D[p[0]:p[1]]
            for rec in sub_Items:
                Id=rec[0]
                SEQ=str(rec[1].seq)
                f_out.write('>'+Id+'\n'+SEQ+'\n')
            f_out.close()
        f_out_runTranslator=open(DIR+'TranslateDNA.sh','w')
        f_out_runTranslator.write('\n'.join(L_TranslatorCode)+'\n')
        f_out_runTranslator.close()
        os.system('chmod +x ' +DIR+'TranslateDNA.sh')
        Commands.append(DIR+'TranslateDNA.sh')
        #os.system(DIR+'TranslateDNA.sh')
    else:
        f_out_SmallFilesTranslator=open('TranslateDNA.sh','w')
        SmallFiles_TranslatorCode.append(CallTranslator + ' ' + f_in + ' ' + GenCode + ' ' + Frames)
        f_out_SmallFilesTranslator.write('\n'.join(SmallFiles_TranslatorCode))
        f_out_SmallFilesTranslator.close()
        os.system('chmod +x TranslateDNA.sh')
        Commands.append('./TranslateDNA.sh')
        #os.system('./TranslateDNA.sh')

f_parallel=open('run_parallelTranslator.sh','w')
C=list(set(Commands))
for c in C:
    f_parallel.write(c+' &\n') #this allows to run all 4 translators together
f_parallel.close()
os.system('chmod +x run_parallelTranslator.sh')
os.system('./run_parallelTranslator.sh')


