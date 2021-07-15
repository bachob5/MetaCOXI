for f in *_codeframes.fas ; do hmmsearch --cpu 16 --cut_tc --domtblout $f'_DomTableOut' -o $f_'Out' /home3/bbalech/LifeWatch/COXI_Extractor/Updates_030420/COX1_Profile/COXI_Profile.hmm $f; done
