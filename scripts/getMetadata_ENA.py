from Bio import SeqIO
import sys
import dateutil
from dateutil import parser
import datetime
from datetime import datetime 

Final_Complete=sys.argv[1]
AllFlatfiles_ENA142=sys.argv[2]


def get_GoodIds(GoodIds_f):
    D_Ids={}
    f=open(GoodIds_f,'r')
    for line in f:
        Acc=line.strip().split('\t')[0]
        D_Ids[Acc]=True
    f.close()
    return D_Ids

Accessions=get_GoodIds(Final_Complete)

def transform_geoCoor(L):
    lat_lon=[]
    try:
        if L[1]=='N':
            lat_lon.append(L[0])
        elif L[1]=='S':
            lat_lon.append('-'+L[0])
        if L[3]=='E':
            lat_lon.append(L[2])
        elif L[3]=='W':
            lat_lon.append('-'+L[2])
        return ' '.join(lat_lon)
    except IndexError:
        print 'NAAA '+' '.join(L)
        return ' '.join(L)


def get_datetime(Str_datetime):
    d_temp=datetime.utcnow()
    YearNow=d_temp.year
    Monthnow=d_temp.strftime('%b')
    Today=d_temp.day
    TodaysDate=str(Today)+'-'+Monthnow+'-'+str(YearNow)
    try:
        x1=dateutil.parser.parse(Str_datetime, yearfirst=True)
        Year=x1.year
        month_abre=x1.strftime('%b')
        Day=x1.day
        DATE=str(Day)+'-'+month_abre+'-'+str(Year)
        if DATE==TodaysDate:
            return 'NA'
        else:
            return DATE
    except ValueError:
        return Str_datetime


def get_Metadata_Flatfile(Flatfile):
    D_Acc_Metadata={}
    #D_Acc_host={}
    f=SeqIO.parse(Flatfile, 'genbank')
    for rec in f:
        acc=rec.id
        M=[] #['country', 'country_address' ,'host', 'collection_year' ,'collection_complete_date', 'geographical_coordinates']
        #Host=[] 
        for seq_feature in rec.features:
            if seq_feature.type=='source':
                try:
                    Location=seq_feature.qualifiers['country'][0].split(':')
                    country=Location[0]
                    M.append(country)
                    try:
                        country_address=Location[1]
                        M.append(country_address)
                    except IndexError:
                        M.append('NA')
                except KeyError:
                    M.append('NA\tNA')
                try:
                    host=seq_feature.qualifiers['host'][0]
                    M.append(host)
                    #Host.append(host)
                except KeyError:
                    M.append('NA')
                    #Host.append('NA')
                try:
                    c_date_temp=seq_feature.qualifiers['collection_date'][0].split('-')
                    year=c_date_temp[-1]
                    exactDate_temp='-'.join(c_date_temp)
                    exactDate=get_datetime(exactDate_temp)
                    M.append(year)
                    M.append(exactDate)
                except KeyError:
                    M.append('NA\tNA')
                try: 
                    geo_coor=seq_feature.qualifiers['lat_lon'][0]
                    lat_lon=transform_geoCoor(geo_coor.split(' '))
                    #print acc, lat_lon
                    M.append(lat_lon)
                except KeyError:
                    M.append('NA')
            D_Acc_Metadata[acc]=M
    return D_Acc_Metadata

Metadata_ENA142=get_Metadata_Flatfile(AllFlatfiles_ENA142)


def update_recs(Metadata_rel228, Metadata_rel229):
    Metadata_Updated={}
    AllAccs_temp=Metadata_rel228.keys()+Metadata_rel229.keys()
    AllAccs=list(set(AllAccs_temp))
    print len(AllAccs)
    for rec in AllAccs:
        try:
            Metadata_Updated[rec]=Metadata_rel229[rec]
        except KeyError:
            Metadata_Updated[rec]=Metadata_rel228[rec]
    return Metadata_Updated

#MetaData_Updated=update_recs(Metadata_rel228, Metadata_rel229)


def write_MetadataFile(Accessions, MetaData_Updated):
    D_Acc_Metadata2Write={}
    L_BOLD_Ids=[]
    for acc in MetaData_Updated:
        try:
            if Accessions[acc]==True:
                D_Acc_Metadata2Write[acc]=MetaData_Updated[acc]
        except KeyError:
            acc1=acc.split('.')[0]
            try:
                if Accessions[acc1]==True:
                    D_Acc_Metadata2Write[acc1]=MetaData_Updated[acc]
            except KeyError:
                continue
    for rec in Accessions:
        try:
            D_Acc_Metadata2Write[rec]
        except KeyError:
            L_BOLD_Ids.append(rec)
    return D_Acc_Metadata2Write, L_BOLD_Ids



Meta_D=write_MetadataFile(Accessions, Metadata_ENA142)
Metadata_NCBI=Meta_D[0]
Ids_BOLD=Meta_D[1]


#LinkToNCBISource='https://www.ncbi.nlm.nih.gov/nuccore/'
LinkToENASource='https://www.ebi.ac.uk/ena/browser/view/'
f_out=open('Metadata_ENA.tsv', 'w')
f_out.write('\t'.join(['Accession', 'country', 'country_address' ,'host', 'collection_year' ,'collection_complete_date', 'geographical_coordinates', 'RetrievalLink'])+'\n')
for acc in Metadata_NCBI:
    f_out.write(acc+'\t'+'\t'.join(Metadata_NCBI[acc])+'\t'+LinkToENASource+acc +'\n')
f_out.close()
 
f_out1=open('IDs2SearchBold.tsv','w')
f_out1.write('\n'.join(Ids_BOLD)+'\n')
f_out1.close()

