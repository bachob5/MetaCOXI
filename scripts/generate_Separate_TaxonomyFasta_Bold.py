from Bio import SeqIO
import sys

f_fasta=sys.argv[1]

def writeTaxonomyFasta(f_fas):
	f=SeqIO.parse(f_fas, 'fasta')
	f_out_fas=open(f_fas.split('.')[0]+'_SHortID.fasta', 'w')
	f_out_Taxo=open(f_fas.split('.')[0]+'_Taxonomy.tsv', 'w')
	for rec in f:
		tempID=rec.description.split(';')
		shortID=tempID[0]
		Taxonomy_B1='\t'.join(tempID[1:])
		if len(Taxonomy_B1)<6:
			print shortID
		Taxonomy_B2='\t'.join(tempID[2:])
		if len(Taxonomy_B2)<6:
			print shortID
		f_out_fas.write('>'+shortID+'\n'+str(rec.seq)+'\n')
		if shortID.split('.')[-1].startswith('COI'):
			f_out_Taxo.write(shortID+'\t'+Taxonomy_B1+'\n')
		else:
			f_out_Taxo.write(shortID+'\t'+Taxonomy_B2+'\n')
	f_out_fas.close()
	f_out_Taxo.close()


writeTaxonomyFasta(f_fasta)
