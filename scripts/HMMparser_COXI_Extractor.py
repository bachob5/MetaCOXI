import re
import os, sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def getGeneticCode(f_in):
    goodcode={}
    searchFile=open(f_in, 'r')
    GC=[]
    for line in searchFile:
        if line[0]=='#': continue
        L=line.split()
        if not L: continue
        NAME=L[0]
        name=NAME.split(';code')[0]
        code=NAME.split(';code')[1].split(';frame')[0]
        GC.append(code)
        score=L[7]
        try:
            goodcode[name][code]+=float(score)
        except KeyError:
            try:
                goodcode[name][code]=float(score)
            except KeyError:
                goodcode[name]={code:float(score)}
    temp=[[x,max(goodcode[x].items(),key=lambda x:x[1])[0]] for x in goodcode]
    #print 'goodcode', goodcode
    Dcont={}
    for i in goodcode:
        temp=[[y,x] for x,y in goodcode[i].items()]
        temp.sort(reverse=True)
        if temp[0][1] in GC:
            Dcont.setdefault(temp[0][1],[]).append(i)
    #print Dcont
    ClassCode=[]
    for i in Dcont:
        ClassCode.append([len(Dcont[i]), i])
    #print ClassCode
    ClassCode.sort(reverse=True)
    if len(ClassCode)>0:
        BestCode=int(ClassCode[0][1])
        GoodSeqs=Dcont[ClassCode[0][1]]
        ExcludSeqs=[]
        if len(Dcont)>1:
            Codes=Dcont.keys()
            Codes.remove(str(BestCode))
            for code in Codes:
                ExcludSeqs.append(Dcont[code])
        return BestCode, GoodSeqs, ExcludSeqs
    searchFile.close()

def HmmDomainsOut(f_in, codeN, GoodSeqs, ExcludSeqs): #Parser of domtblout option of hmmsearch
    QueryDomains={} #the query as key and a list of domains as values
    QueriesDomsPos={} #the domains as keys and a list of queries as values
    AllDoms={}
    AllDoms_temp={}
    PosQueries={}
    Poss={}
    Bach={}
    HmmPoss_temp={}
    QueryAllPos={}
    HmmPoss={}
    f=open(f_in, 'r')
    for i in xrange(3):
        f.next()
    for line in f:
        L=re.split(' +', line)
        #HmmPoss.setdefault(L[0],[]).append([L[3], int(L[-8]),int(L[-7])])
        if ';code'+str(codeN)+ ';' in L[0] and L[0].split(';')[0] in GoodSeqs:
            Start=int(L[-6])-1
            End=L[-5]
            KEY1 = L[0].split(';')
            PosS=int(KEY1[-2][1:])
            PosE=int(KEY1[-1][1:])
            ID, PosOrigin, PosInternal, Profile = KEY1[0], [PosS, PosE], [int(L[-6]), int(L[-5])], L[3]
            QueriesDomsPos.setdefault(L[0], []).append(Profile+'_' + str(Start) + '_' +End)
            if ';frame-' not in L[0]:
                MappedStartNuc=PosS+(Start*3)
                MappedEndNuc=PosS+((int(End)*3)+1)
                MappedStartProt=Start
                MappedEndProt=End
                NucStEnd=';'+';'.join([str(MappedStartNuc), str(MappedEndNuc)])
                ProtStEnd=';'+';'.join([str(MappedStartProt), str(MappedEndProt)])
                Poss.setdefault(KEY1[0], []).append(Profile)
                QueryDomains.setdefault(L[0]+ProtStEnd+NucStEnd,[]).append(Profile)
                HmmPoss_temp.setdefault(L[0]+ProtStEnd+NucStEnd,[]).append(Profile + ';' +';'.join([L[-8],L[-7]]))
            else:
                MappedStartNuc=PosE-(int(End)*3)
                MappedEndNuc=PosE-(Start*3)
                MappedStartProt=Start
                MappedEndProt=End
                NucStEnd=';'+';'.join([str(MappedStartNuc), str(MappedEndNuc)])
                ProtStEnd=';'+';'.join([str(MappedStartProt), str(MappedEndProt)])
                Poss.setdefault(KEY1[0], []).append(Profile)
                QueryDomains.setdefault(L[0]+ProtStEnd+NucStEnd,[]).append(Profile)
                HmmPoss_temp.setdefault(L[0]+ProtStEnd+NucStEnd,[]).append(Profile + ';' +';'.join([L[-8],L[-7]]))
    f.close()
    for j in HmmPoss_temp:
        HmmPoss[j]='_'.join(HmmPoss_temp[j])
    #print 'HmmPoss\n', HmmPoss
    #print 'QueryDomains', QueryDomains
    #print 'QueriesDomsPos', QueriesDomsPos
    Bach_temp={}
    OnlyQueries=set(QueryDomains.keys())
    #print '\n\nOnlyQueries', len(OnlyQueries), len(set(OnlyQueries))
    for i in OnlyQueries:
        temp=[]
        DomSucc=[]
        IDnoPos=i.split(';')
        PosNuc=[int(IDnoPos[-2]), IDnoPos[-1]]
        temp.append(PosNuc)
        DomNumber=len(QueryDomains[i])
        DomSucc.append([temp[0]]*DomNumber)
        #print DomSucc
        #for m in QueryDomains[i]:
        Bach_temp.setdefault(IDnoPos[0],[]).append(sum(DomSucc,[]))
    #print '\n\nBach_temp', Bach_temp
    for t in Bach_temp:
        Bach[t]=sum(Bach_temp[t],[])
        Bach[t].sort()
    temp={}
    for rec in QueryDomains:
        rec1=rec.split(';')
        ID = rec1[0] + ';'+';'.join([rec1[-2],rec1[-1]])
        temp[ID] = QueryDomains[rec]
    #print temp,'\n\nBach1', Bach
    for Id in Bach:
        tempL=[]
        for n in Bach[Id]:
            Poss=';'+str(n[0]) + ';' +n[-1]
            #if temp[Id+Poss] not in tempL:
            tempL.append(temp[Id+Poss])
        AllDoms_temp[Id] = tempL
        #AllDoms[Id] = temp[Id+Poss]
    #print 'AllDoms_temp', AllDoms_temp
    for i in AllDoms_temp:
        for d in AllDoms_temp[i]:
            AllDoms.setdefault(i,[]).append(d[0])
    #print 'Bach\n', Bach
    #print 'AllDoms\n', AllDoms
    PossDomsNoInclusions=DomPos_Inclusions(Bach,AllDoms)
    #print 'PossDomsNoInclusions\n', PossDomsNoInclusions
    PossNoInclusions=PossDomsNoInclusions[0]
    #print 'PossNoInclusions', PossNoInclusions, '\n'
    DomsNoInclusions=PossDomsNoInclusions[1]
    #print 'DomsNoInclusions', DomsNoInclusions
    FinalDoms=DomPos_OverlapCanc(PossNoInclusions)
    #print ControlDomsSuccession(DomsNoInclusions)
    HomoloSeq=list(set(sum(ControlDomsSuccession(DomsNoInclusions)[-2],[])))
    #print 'HomoloSeq', HomoloSeq
    DomsSucc=DomsNoInclusions[HomoloSeq[0]]
    #print DomsSucc
    DomSuccession=[]
    for index, DomX in enumerate(DomsSucc): #DomsSucc is a list of sorted domains
        DomSuccession.append(str(index) + '_'+ DomX)
    #print 'DomSuccession', DomSuccession
    FinalDict={}
    for dom in DomSuccession:
        tempL1=[]
        Ind=int(dom.split('_')[0])
        for ID in HomoloSeq:
            #print dom, Ind, ID
            try:
                Pos_intermediate= ';'+str(PossNoInclusions[ID][Ind][0]) + ';' +PossNoInclusions[ID][Ind][1]
            except IndexError:
                continue
            for OrgID in OnlyQueries:
                temp=OrgID.split(';')
                if ID== temp[0] and Pos_intermediate == ';' + temp[-2] + ';' + temp[-1]:
                    Id_Align=';'.join(temp[:5]) + '/' + str(int(temp[5])+1) + '-' + temp[6]
                    tempL1.append(Id_Align)
        FinalDict[dom] = list(set(tempL1)) 
    len(QueriesDomsPos.keys()), len(set(QueriesDomsPos.keys())), '\n', 
    len(AllDoms.keys()), len(set(AllDoms.keys())), '\n', 
    len(FinalDict.keys()), len(set(FinalDict.keys())), '\n', 
    len(Bach.keys()), len(set(Bach.keys())), '\n'
    #print 'Bach', Bach 
    return QueryDomains, QueriesDomsPos, DomsNoInclusions, FinalDict, PossNoInclusions, HmmPoss, PossDomsNoInclusions

def ControlDomsSuccession(D_QuDoms):
    """this function Controls if all queries have the same domains and in the same order. if this is true it returns 
    only Doms, otherwise it returns Doms and the queries to exclude listed in ExcludeSeq"""
    Doms={}
    ExcludeSeq=[]
    GoodSeq=[]
    N=[]
    for query in D_QuDoms:
        Doms.setdefault(';'.join(D_QuDoms[query]), []).append(query)
    for x in Doms:
        N.append(len(Doms[x]))
    N.sort(reverse=True)
    for x in Doms:
        if len(Doms[x]) == N[0]:
            GoodSeq.append(Doms[x])
        else:
            ExcludeSeq.append(Doms[x])
    return Doms, ExcludeSeq, GoodSeq, N

def checkExonIntr(Dict_OrgIDsDoms, Dict_IDsFragHmmPoss):
    Dict_OrgID_HmmPoss={}
    OrgIDs=Dict_OrgIDsDoms.keys()
    for Id in OrgIDs:
        for IDfrag in Dict_IDsFragHmmPoss.keys():
            if Id == IDfrag.split(';')[0]:
                Dict_OrgID_HmmPoss.setdefault(Id,[]).append({IDfrag: Dict_IDsFragHmmPoss[IDfrag]})
    return Dict_OrgID_HmmPoss
            

def getDomsSucc(Dict):
    '''this function writes a file with the domains succession and 
    return a list of these domain to give to the merger'''
    f_out = open('DomainsSuccession.txt', 'w')
    L=[]
    Doms = Dict.keys()
    for i in Doms:
        L.append([int(i.split('_')[0]), '_'.join(i.split('_')[1:])])
    L.sort()
    SortedDoms=[]
    for Num in L:
        SortedDoms.append(Num[-1])
    f_out.write('\n' .join(SortedDoms))
    f_out.close()
    return Doms


def DomPos_Inclusions(ID_Pos, ID_Doms):
    '''This function controls if two lists of domains positions are included one another. 
    It returns new ID_Pos (dictionary {ID:[[poss],[poss]]})and new ID_Doms (dictionary {ID:[dom1,dom2]}). 
    In this case the domain included in the other is eliminated. The user will not have any feedback on that 
    because he will obtain the same result as if he searches the sequence on pfam website.'''
    new_IDPos={}
    new_ID_Doms={}
#     equalPoss={}
#     equalPoss_diffDoms={}
#     l_equalDoms=[]
    for Id in ID_Pos:
        #print len(ID_Doms[Id]), len(ID_Pos[Id])
        temp=[]
        temp_Doms=[]
        while len(ID_Pos[Id]) >=2:
            Pos_dom1=set(range(ID_Pos[Id][0][0],int(ID_Pos[Id][0][1])))
            Pos_dom2=set(range(ID_Pos[Id][1][0],int(ID_Pos[Id][1][1])))
            if Pos_dom2.issubset(Pos_dom1):
                ID_Pos[Id].pop(1)
                ID_Doms[Id].pop(1)
                if ID_Pos[Id][0] not in temp:
                    temp.append(ID_Pos[Id][0])
                    temp_Doms.append(ID_Doms[Id][0])
            elif Pos_dom1.issubset(Pos_dom2):
                ID_Doms[Id].pop(0)
                ID_Pos[Id].pop(0)
                if ID_Pos[Id][0] not in temp:
                    #print ID_Pos[Id]
                    temp.append(ID_Pos[Id][0])
                    temp_Doms.append(ID_Doms[Id][0])
            else:
                if ID_Pos[Id][0] not in temp:
                    #print ID_Pos[Id]
                    #print ID_Doms[Id]
                    temp.append(ID_Pos[Id][0])
                    temp_Doms.append(ID_Doms[Id][0])
                    #ID_Pos[Id].pop(0)
                    #ID_Doms[Id].pop(0)
                else:
                    ID_Pos[Id].pop(0)
                    ID_Doms[Id].pop(0)
            #print ID_Doms[Id]
            #print temp_Doms
        if len(ID_Pos[Id]) ==1:
            #print ID_Pos[Id], ID_Doms[Id]
            if ID_Pos[Id][0] not in temp:
                temp.append(ID_Pos[Id][0])
                temp_Doms.append(ID_Doms[Id][0])
                ID_Pos[Id].pop(0)
                ID_Doms[Id].pop(0)
#             else:
#                 ID_Pos[Id].pop(0)
#                 ID_Doms[Id].pop(0)                
        elif len(ID_Pos[Id]) ==0:
            continue
        new_IDPos[Id]=temp
        new_ID_Doms[Id]=temp_Doms
    return new_IDPos, new_ID_Doms#, equalPoss, equalPoss_diffDoms, l_equalDoms


# IP={'gi|507925503|ref|XM_004673961.1|:1-2381':[[69, '289'], [306, '544'], [310,'500'],[552, '787'], [822, '1057'], [1128, '1363'], [1440, '1681'], [1740, '2380'], [2400,'2500'], [2441,'2449']]}
# DP={'gi|507925503|ref|XM_004673961.1|:1-2381':['PAN_1', 'Kringle', 'try', 'Kringle', 'Kringle', 'Kringle', 'Trypsin', 'DUF1986','bach','getthefuck']}
# DomPos_Inclusions(IP,DP)


def getExons(DomsNoInclusions, PossNoInclusions, HmmPoss):
    PossNoInclusions_Exons={}
    DomsNoInclusions_Exons={}
    IDs = PossNoInclusions.keys()
    for Id in IDs:
        temp_doms=[]
        temp_poss=[]
        ID_DomsPoss=PossNoInclusions[Id]
        ID_DomsNames=DomsNoInclusions[Id]
        while len(ID_DomsNames) >= 2:
            if ID_DomsNames[0] != ID_DomsNames[1]:
                #print ID_DomsNames[0], ID_DomsNames[1]
                temp_doms.append(ID_DomsNames[0])
                temp_poss.append(ID_DomsPoss[0])
                ID_DomsPoss.pop(0)
                ID_DomsNames.pop(0)
            elif ID_DomsNames[0] == ID_DomsNames[1]:
                #print ID_DomsNames[0], ID_DomsNames[1]
                Id0_HmmPoss = [x for x in HmmPoss.keys() if Id == x.split(';')[0] and ';'.join([str(ID_DomsPoss[0][0]),ID_DomsPoss[0][1]]) in x]
                Id1_HmmPoss = [y for y in HmmPoss.keys() if Id == y.split(';')[0] and ';'.join([str(ID_DomsPoss[1][0]),ID_DomsPoss[1][1]]) in y]
                #print 'Id1_HmmPoss', Id1_HmmPoss
                #print 'Id0_HmmPoss', Id0_HmmPoss
                Dom0_HmmPoss=HmmPoss[Id0_HmmPoss[0]].split(';')
                Dom1_HmmPoss=HmmPoss[Id1_HmmPoss[0]].split(';')
                print Dom0_HmmPoss, Dom1_HmmPoss
                assert Dom0_HmmPoss[0] == Dom1_HmmPoss[0]
                if int(Dom0_HmmPoss[-1]) <= int(Dom1_HmmPoss[1]): #Intron detected
                    #print [temp_poss[-1][-2:]], '\t', ID_DomsPoss[0], '\t',Dom0_HmmPoss, Dom1_HmmPoss, Dom0_HmmPoss[0], Id0_HmmPoss, Id1_HmmPoss
                    if len(temp_poss) ==0:
                        temp_doms.append(ID_DomsNames[0])
                        temp_poss.append(ID_DomsPoss[0] + ID_DomsPoss[1])
                    else:
                        if temp_poss[-1][-2:] == ID_DomsPoss[0]:
                            pass
                        else:
                            temp_poss[-1]=temp_poss[-1]+ID_DomsPoss[0]
                else:
                    if len(temp_poss) !=0:
                        if temp_poss[-1][-2:] == ID_DomsPoss[0]:
                            pass
                        else:
                            temp_doms.append(ID_DomsNames[0])
                            temp_poss.append(ID_DomsPoss[0])
                    else:
                        temp_doms.append(ID_DomsNames[0])
                        temp_poss.append(ID_DomsPoss[0])
                ID_DomsPoss.pop(0)
                ID_DomsNames.pop(0)
        if len(ID_DomsNames) ==1: #account for the last 1 or 2 domains
            if len(temp_doms)!=0:
                if temp_doms[-1] == ID_DomsNames[0]:
                    if temp_poss[-1][-2:] == ID_DomsPoss[0]:
                        pass
                    else:
                        if int(Dom0_HmmPoss[-1]) < int(Dom1_HmmPoss[1]):
                            if temp_poss[-1][-2:] == ID_DomsPoss[0]:
                                pass
                            else:
                                temp_poss[-1]=temp_poss[-1] + ID_DomsPoss[0]
                        else:
                            temp_doms.append(ID_DomsNames[0])
                            temp_poss.append(ID_DomsPoss[0])
                else:
                    temp_doms.append(ID_DomsNames[0])
                    temp_poss.append(ID_DomsPoss[0])
            else:
                temp_doms.append(ID_DomsNames[0])
                temp_poss.append(ID_DomsPoss[0])
        elif len(ID_DomsNames) ==0:
            break
        PossNoInclusions_Exons[Id] = temp_poss
        DomsNoInclusions_Exons[Id] = temp_doms
    return DomsNoInclusions_Exons, PossNoInclusions_Exons


def DomPos_OverlapCanc(ID_Pos):
    D={}
    ExcludedSites={}
    for i in ID_Pos:
        ExcluSitesLst=[]
        if len(ID_Pos[i]) == 1:
            D[i] = ID_Pos[i]
        else:
            LL = ID_Pos[i]
            Lfinale = []
            j = 0
            while j <= range(len(LL)) [-2]:
                if int(LL[j][-1]) >= LL[j+1][0]:
                    Diff = int(LL[j][-1]) - int(LL[j+1][-1])
                    #NewInterval_0 = LL[j]#, int(LL[j][-1])
#                    Lfinale.append(NewInterval_0)
                    NewInterval_1 = [int(LL[j][-1])+1, int(LL[j+1][-1])]
                    Lfinale.append(NewInterval_1)
                    ExcludedSites.setdefault(i, []).append(int(LL[j][-1]) + Diff)
                else:
                    Lfinale.append(LL[j])
                    Diff = int(LL[j][-1]) - LL[j+1][0]
                    ExcludedSites.setdefault(i, []).append([int(LL[j][-1])+1, LL[j+1][0]])
                j +=1
            Lfinale.append(LL[-1])
            NewPoss = Lfinale
            D[i] = NewPoss
    return D, ExcludedSites





'''AY363898.1;code2;frame3;S158;E425;1;79;161;396 COX1;258;335
AY363898.1;code2;frame3;S2;E155;0;51;2;156 COX1;206;255
'''


def write_IDsWithPoss(MappedSites):
    f_out_Ids=open(f_in.split('.')[0]+'_MatchIDs_PfamProfile.txt','w') #remember to change the split character
    D_Acc_AccWithPoss={}
    for AccInfo in MappedSites:
        Acc=AccInfo.split(';')[0]
        D_Acc_AccWithPoss.setdefault(Acc,[]).append(AccInfo)
    for acc in D_Acc_AccWithPoss:
        f_out_Ids.write(' '.join(D_Acc_AccWithPoss[acc])+ '\n')
    f_out_Ids.close()



filesList=[x for x in os.listdir(os.curdir) if x.split('_')[-1]=='DomTableOut']
print filesList
for f_in in filesList:
    print f_in
    try:
        C=getGeneticCode(f_in)
        Code =C[0]
        GoodSeqs=C[1]
        ExcludSeqs=C[2]
        #print Code
        A=HmmDomainsOut(f_in, Code, GoodSeqs, ExcludSeqs)
        MappedSites=A[-2]
        write_IDsWithPoss(MappedSites)
    except TypeError:
        continue

