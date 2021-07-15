import requests

TAXA=['Acanthocephala', 'Acoelomorpha', 'Brachiopoda', 'Bryozoa', 'Chaetognatha', 'Cnidaria', 'Ctenophora', 'Cycliophora', 'Echinodermata', 'Entoprocta', 'Gastrotricha', 'Gnathostomulida', 'Hemichordata', 'Kinorhyncha', 'Nematoda', 'Nematomorpha', 'Nemertea', 'Onychophora', 'Phoronida', 'Placozoa', 'Platyhelminthes', 'Porifera', 'Priapulida', 'Rhombozoa', 'Rotifera', 'Sipuncula', 'Tardigrada', 'Xenacoelomorpha','Chordata','Mollusca','Annelida','Arthropoda']

TAXA=['Arthropoda']

Arthro=['Arachnida', 'Branchiopoda', 'Cephalocarida', 'Chilopoda', 'Collembola', 'Copepoda', 'Diplopoda', 'Diplura', 'Hexanauplia', 'Ichthyostraca', 'Malacostraca', 'Merostomata', 'Oligostraca_class_incertae_sedis', 'Ostracoda', 'Pauropoda', 'Pentastomida', 'Protura', 'Pycnogonida', 'Remipedia', 'Symphyla', 'Metatiron', 'Insecta']

Insecta=['Archaeognatha', 'Blattodea', 'Coleoptera', 'Dermaptera', 'Embioptera', 'Ephemeroptera', 'Hemiptera', 'Mantodea', 'Mecoptera', 'Megaloptera', 'Neuroptera', 'Notoptera', 'Odonata', 'Orthoptera', 'Phasmatodea', 'Plecoptera', 'Psocodea', 'Raphidioptera', 'Siphonaptera', 'Strepsiptera', 'Thysanoptera', 'Trichoptera', 'Zoraptera', 'Zygentoma', 'Lepidoptera', 'Hymenoptera', 'Diptera']

Diptera=['Acartophthalmidae', 'Acroceridae', 'Agromyzidae', 'Anisopodidae', 'Anthomyiidae', 'Anthomyzidae', 'Apioceridae', 'Apsilocephalidae', 'Apystomyiidae', 'Asilidae', 'Asteiidae', 'Atelestidae', 'Athericidae', 'Aulacigastridae', 'Australimyzidae', 'Austroleptidae', 'Axymyiidae', 'Bibionidae', 'Blephariceridae', 'Bolbomyiidae', 'Bolitophilidae', 'Bombyliidae', 'Brachystomatidae', 'Braulidae', 'Calliphoridae', 'Camillidae', 'Canacidae', 'Canthyloscelidae', 'Carnidae', 'Cecidomyiidae', 'Celyphidae', 'Ceratopogonidae', 'Chamaemyiidae', 'Chaoboridae', 'Chironomidae', 'Chloropidae', 'Chyromyidae', 'Clusiidae', 'Coelopidae', 'Conopidae', 'Corethrellidae', 'Cryptochetidae', 'Ctenostylidae', 'Culicidae', 'Curtonotidae', 'Cylindrotomidae', 'Cypselosomatidae', 'Deuterophlebiidae', 'Diadocidiidae', 'Diastatidae', 'Diopsidae', 'Diptera_family_incertae_sedis', 'Ditomyiidae', 'Dixidae', 'Dolichopodidae', 'Drosophilidae', 'Dryomyzidae', 'Empididae', 'Ephydridae', 'Evocoidae', 'Fanniidae', 'Fergusoninidae', 'Glossinidae', 'Gobryidae', 'Helcomyzidae', 'Heleomyzidae', 'Helosciomyzidae', 'Hesperinidae', 'Hilarimorphidae', 'Hippoboscidae', 'Homalocnemiidae', 'Huttoninidae', 'Hybotidae', 'Inbiomyiidae', 'Ironomyiidae', 'Iteaphila-group', 'Keroplatidae', 'Lauxaniidae', 'Limoniidae', 'Lonchaeidae', 'Lonchopteridae', 'Lygistorrhinidae', 'Marginidae', 'Megamerinidae', 'Mesembrinellidae', 'Micropezidae', 'Milichiidae', 'Mormotomyiidae', 'Muscidae', 'Mycetophilidae', 'Mydidae', 'Mystacinobiidae', 'Mythicomyiidae', 'Nannodastiidae', 'Natalimyzidae', 'Nemestrinidae', 'Neminidae', 'Neriidae', 'Neurochaetidae', 'Nothybidae', 'Nymphomyiidae', 'Odiniidae', 'Oestridae', 'Opetiidae', 'Opomyzidae', 'Oreogetonidae', 'Oreoleptidae', 'Pachyneuridae', 'Pallopteridae', 'Pantophthalmidae', 'Paraleucopidae', 'Pediciidae', 'Pelecorhynchidae', 'Periscelididae', 'Perissommatidae', 'Phoridae', 'Piophilidae', 'Pipunculidae', 'Platypezidae', 'Platystomatidae', 'Polleniidae', 'Pseudopomyzidae', 'Psilidae', 'Psychodidae', 'Ptychopteridae', 'Pyrgotidae', 'Rhagionidae', 'Rhiniidae', 'Rhinophoridae', 'Richardiidae', 'Ropalomeridae', 'Sarcophagidae', 'Scathophagidae', 'Scatopsidae', 'Scenopinidae', 'Sciaridae', 'Sciaroidea_incertae_sedis', 'Sciomyzidae', 'Sepsidae', 'Simuliidae', 'Somatiidae', 'Sphaeroceridae', 'Stratiomyidae', 'Strongylophthalmyiidae', 'Syringogastridae', 'Syrphidae', 'Tabanidae', 'Tachinidae', 'Tanyderidae', 'Tanypezidae', 'Tephritidae', 'Teratomyzidae', 'Thaumaleidae', 'Therevidae', 'Tipulidae', 'Trichoceridae', 'Ulidiidae', 'Vermileonidae', 'Xenasteiidae', 'Xylomyidae', 'Xylophagidae']


BaseURL='http://www.boldsystems.org/index.php/API_Public/combined?taxon=' 
format='&format=tsv'


def Extract_Taxon_Data(BaseURL, taxon, format):
    #BLOCK_SIZE=1024
    response=requests.get(BaseURL+taxon+format, stream=True)
    handle=open(taxon+'_FullRecords_BOLD.tsv.tsv', "wb")
    for chunk in response.iter_content(chunk_size=512):
        if chunk:  # filter out keep-alive new chunks
            handle.write(chunk)
    handle.close()


for t in TAXA:
    T=Extract_Taxon_Data(BaseURL, t, format)

for t in Arthro:
    T=Extract_Taxon_Data(BaseURL, t, format)


for t in Insecta:
    print 'Processing ', t
    T=Extract_Taxon_Data(BaseURL, t, format)
    print t, 'DONE'

for t in Diptera:
    print 'Processing ', t
    T=Extract_Taxon_Data(BaseURL, t, format)
    print t, 'DONE'

