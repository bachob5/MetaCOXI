import sys, os

AllAccs=sys.argv[1]
PseudoEntries=sys.argv[2]


def get_IdsList(f_in):
    Ids={}
    f=open(f_in,'r')
    for line in f:
        Ids[line.strip()]=True
    f.close()
    return Ids


Ids_Dict=get_IdsList(AllAccs)
print len(Ids_Dict)
Pseudo_Dict=get_IdsList(PseudoEntries)
print len(Pseudo_Dict)


def get_finalIdsList(Ids_Dict, Pseudo_Dict):
    FinalIDs=[]
    for i in Ids_Dict:
        try:
            Pseudo_Dict[i]
        except KeyError:
            FinalIDs.append(i)
    f_out=open(AllAccs.split('.')[0]+'.txt','w')
    f_out.write('\n'.join(FinalIDs)+'\n')
    f_out.close()
    
W=get_finalIdsList(Ids_Dict, Pseudo_Dict)
