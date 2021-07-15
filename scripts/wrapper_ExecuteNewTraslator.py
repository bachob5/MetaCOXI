import os,sys


# Dir_FasFiles='ENA_r142/'

Dir_FasFiles=sys.argv[1]

for f in os.listdir(Dir_FasFiles):
    f_temp=f.split('_')
    if f_temp[-1]=='withMitocode.fasta':
        MitoCode=f_temp[0]
        if MitoCode!='0':
            os.system('python newTranslator_CExtract.py ' + Dir_FasFiles+f + ' ' + MitoCode +' ' + '1,2,3,-1,-2,-3')

