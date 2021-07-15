from Bio import SeqIO
from Bio.SeqIO import SeqRecord
#from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import Alphabet
from Bio.Alphabet import IUPAC
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA 
import os,sys
import gzip

DIR_ENA_rel142=sys.argv[1] #~/DBs/ENA/r142/



def get_good_Excluded_IDs(f_txt):
    L_Ids=[]
    try:
        f=open(f_txt,'r')
        for line in f:
            l=line.strip().split()
            Id=l[0].split(';')[0]
            L_Ids.append(Id)
        f.close()
    except IOError:
        print f_txt
    return L_Ids

def Join_2lists(L1, L2):
    L=L1+L2
    #print L
    return list(set(L))


def get_gb_entries(f_gb_withSelectedEntries, L_Ids):
    #Seq_Records=[]
    f_gb_original='_'.join(f_gb_withSelectedEntries.split('_')[1:6])
    f_handle=gzip.open(DIR_ENA_rel142+f_gb_original, 'rt')
    f_gb=SeqIO.parse(f_handle, 'embl')
    f_gb_prefix='_'.join(f_gb_withSelectedEntries.split('_')[:5])
    out_gb=open(f_gb_prefix+'_Good_Excluded_Entries.gb', 'w')
    for rec in f_gb:
        if rec.id in L_Ids:
            SeqIO.write(rec, out_gb ,'genbank')
    out_gb.close()
            #Seq_Records.append(SeqRecord(rec))



def getItWork_all():
    f_out=open('Info_Readme.txt','w')
    f_out.write('FileName\tNumber of Good and Excluded Entries\n')
    #L_singleMatchFiles=[x for x in os.listdir(os.curdir) if '_'.join(x.split('_')[1:])=='SingleMatchIDs_PfamProfile.txt']
    #L_ExcludedIDs=[y for y in os.listdir(os.curdir) if '_'.join(y.split('_')[1:])=='ExcludedIDs_PfamProfile.txt']
    L_MatchFiles=[x for x in os.listdir(os.curdir) if '_'.join(x.split('_')[-2:])=='MatchIDs_PfamProfile.txt']
    L_gbs=[x for x in os.listdir(DIR_ENA_rel142) if x.split('.')[-2:]==['dat', 'gz']]
    for MatchFile in L_MatchFiles:
        f_gb_temp=MatchFile.split('_')
        f_gb_prefix='_'.join(f_gb_temp[1:5]) #change the split code to '.' in case you use the genbank files without '_SelectedEntries_'
        f_gb_withSelectedEntries='_'.join(f_gb_temp[:6])+'.dat.gz'
        #good_IDs=get_good_Excluded_IDs(f_gb_prefix+'_ExcludedIDs_PfamProfile.txt')
        #Excluded_IDs=get_good_Excluded_IDs(f_gb_prefix+'_SingleMatchIDs_PfamProfile.txt')
        #MatchedIDs_temp=[y for y in L_MatchFiles if y.split('_')[-1]=='_MatchIDs_PfamProfile.txt' and f_gb_prefix in y]
        #assert len(MatchedIDs_temp)==1
        #MatchedIDs=get_good_Excluded_IDs(MatchedIDs_temp[0])
        MatchedIDs=get_good_Excluded_IDs(MatchFile)
        #print 'MatchedIDs', MatchedIDs
        #good_excluded_Ids=Join_2lists(good_IDs, Excluded_IDs)
        f_out.write(f_gb_withSelectedEntries + '\t'+str(len(MatchedIDs))+'\n')
        WriteEntries=get_gb_entries(f_gb_withSelectedEntries, MatchedIDs)
    f_out.close()

getItWork_all()
















