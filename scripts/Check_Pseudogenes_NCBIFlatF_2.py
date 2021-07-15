from Bio import SeqIO
import sys,os

PATH=sys.argv[1] 
GB_filesList=[x for x in os.listdir(PATH) if '_'.join(x.split('_')[5:]) =='Good_Excluded_Entries.gb']

def get_newPoss(oldPos):
    if oldPos.startswith('<'):
        newPos=oldPos[1:]
        return newPos
    elif oldPos.startswith('>'):
        newPos=oldPos[1:]
        return newPos
    else:
        newPos=oldPos
        return newPos


def get_GbAnnotations(f_gb, PATH):
    D_Acc_Annotations={}
    handle_GB=SeqIO.parse(PATH+f_gb, 'genbank')
    L_pseudogEntries=[]
    L_EntriesWithFrameshifts=[]
    for rec in handle_GB:
        tempL=[]
        for seq_feature in rec.features:
            for p in seq_feature.qualifiers:
                if p.lower().startswith('pseudo'):
                    #print f_gb, p, rec.id
                    L_pseudogEntries.append(rec.id)
            if seq_feature.type=='CDS': # Add the control of Pseudogenes
                #print seq_feature.qualifiers
                if len(seq_feature.location.parts) > 1:
                    L_EntriesWithFrameshifts.append(rec.id)
                    tempJoin=[]
                    for j in seq_feature.location.parts:
                        if 'gene' in seq_feature.qualifiers:
                            GeneName=seq_feature.qualifiers['gene'][0]
                            startPos=get_newPoss(str(j.start))
                            endPos=get_newPoss(str(j.end))
                            tempJoin.append('Gene='+GeneName+';'+startPos+';'+endPos)
                        else:
                            if 'product' in seq_feature.qualifiers:
                                ProductName=seq_feature.qualifiers['product'][0]
                                startPos=get_newPoss(str(j.start))
                                endPos=get_newPoss(str(j.end))
                                tempJoin.append('product='+ProductName+';'+startPos+';'+endPos) #this means that there is no annotation named 'gene' in the entry
                    tempL.append(tempJoin)
                    #print f_gb, tempJoin
                elif len(seq_feature.location.parts) == 1: # Add the control of Pseudogenes
                    PosS = str(seq_feature.location.start)
                    newPosS=get_newPoss(PosS)
                    PosE=str(seq_feature.location.end)
                    newPosE=get_newPoss(PosE)
                    if 'gene' in seq_feature.qualifiers:
                        GeneName=seq_feature.qualifiers['gene'][0]
                        tempL.append(['Gene='+GeneName+';'+newPosS+';'+newPosE])
                    else:
                        if 'product' in seq_feature.qualifiers:
                            ProductName=seq_feature.qualifiers['product'][0]
                            tempL.append(['product='+ProductName+';'+newPosS+';'+newPosE])
        D_Acc_Annotations.setdefault(rec.id, []).append(tempL)
    if len(list(set(L_pseudogEntries)))>0:
        f_out=open(PATH+f_gb.split('.')[0]+'_Pseudogenes_Entries.txt','w')
        f_out.write('\n'.join(list(set(L_pseudogEntries)))+'\n')
        f_out.close()
    if len(list(set(L_EntriesWithFrameshifts)))>0:
        f_out_2=open(PATH+f_gb.split('.')[0]+'_FrameShifts_Entries.txt','w')
        f_out_2.write('\n'.join(list(set(L_EntriesWithFrameshifts)))+'\n')
        f_out_2.close()        
    return D_Acc_Annotations

for f in GB_filesList:
    D=get_GbAnnotations(f, PATH)
    
