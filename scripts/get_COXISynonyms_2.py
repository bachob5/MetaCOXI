import os

filesList=[x for x in os.listdir(os.curdir) if 'Good_Excluded_Entries_COXISynonyms.txt' in x]
print filesList

FolderName = os.getcwd().split('/')[-1]
print FolderName
L_Syns=[]
for FileName in filesList:
    f=open(FileName,'r')
    for line in f:
        L_Syns.append(line.strip())
    f.close()


f_out_allCOXIsynonyms=open(FolderName+'_AllCOXISynonyms.txt', 'w')
f_out_allCOXIsynonyms.write('\n'.join(list(set(L_Syns))))    
f_out_allCOXIsynonyms.close()
