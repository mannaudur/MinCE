## Fyrstu skref í MinCE þegar ég er kominn út til Caltech (mögulega, ég veit ekki)

- Þýða mince.py og triangulate.py yfir í C++
    - Passa að nota Bifröst fyrir hakkafallið, basically bara tengja það inn út sketch
    - Gæti verið gott að styðjast við b-tré fyrir dist_all()-fallið
        - Þarf að gera það hraðara einhvern veginn
- Þýða solve_features.py og deBruijn_solver.py yfir í C++
    - Gæti þurft að hugsa formattið á atom-skránum upp á nýtt, json er stórt og ekki vinur C++
    - Þarf að finna fljótlegri leið til að greina það fljótt
        - Google er vinur minn en ég gæti líka skoðað Burrows-Wheeler transform eða eitthvað álíka
        - Gæti líka sniðið þetta eins og tré og reynt að finna niðurstöður þannig
            - Tréð búið til úr þeim nóðum sem á að finna
            - Sérhver k-mer úr leitarsafninu þýddur yfir í hakkagildi og sendur inn í tréð
        - Er tvítekning að þýða k-merana yfir í hakkagildi?
            - Þess þarf ekki, því þetta er algjörlega aðskilið rissunum
            - Betra að vinna bara með k-mera mengið?
                - Pottþétt til einhver góð tól til þess að greina slíkt
                    - Augljóslega myndi ég nota kallisto í það
- Þýða extract_features.py og deBruijn_extractor.py yfir í C++
    - Nota mmh3
        - *Nema ég ætli að vinna bara með k-merana og tengja þetta þannig við kallisto*
    - Gæti þurft að hugsa formattið á atom-skránum upp á nýtt, json er stórt og ekki vinur C++
    - Ef ég ætla að gera þetta almennilega, þarf að skoða allt reikniritið upp á nýtt
        - Tryggja áreiðanleika niðurstaðna
            - Líklega ómögulegt nema með ítrekuðum tilraunum
        - Reyna að spara minni og tíma hvar sem ég get
            - Skoða reikniritið í heild sinni, rissa það upp og greina frá A til Ö
- Byrja að skilgreina Leiden-undiratóm og útfæra á nýja safninu

# Jákvæð atriði við útfærsluna á bottom up cliques:
- Einkvæm vörpun, hver rissa er bara í einni klíku og klíkurnar eru þ.a.l. ekki 'overlapping' - einfaldar margt
- Minnstu klíkurnar fá mest væri, svo þegar maður skoðar radíus í kringum besta match er líklegt að maður fái líkar rissur saman í klíku (þarf að orða þetta betur, þyrfti að hugsa smá um þetta)
- Ef rissa er í fjarlægð 1 frá fullt af 0-klíkum, þá munu þær 0-klíkur fylgja með í deBruijn lausnarforritinu, svo maður missir ekki af þeim
- Þéttari klíkur, þ.e. þær með lágan radíus eins og 0-klíkur eða 1-klíkur, verða minni og gefa þ.a.l. nákvæmari niðurstöður úr deBruijn forritinu. Stærri klíkur innihalda fjarskyldari meðlimi, svo þar ætti einnig að vera auðveldara að velja þessar raðir. Verst er að blanda saman skyldum og óskyldum tegundum í stórum klíkum. Þessi útfærsla sneiðir framhjá því.
- Ef ske kynni að stakur einstaklingur A missir af sinni nálægustu klíku, t.d. ef hann er í 1 fjarlægð frá þéttri 0-klíku, þá mun hann parast í aðra klíku í næstu skrefum, kannski 1-klíku eða 2-klíku, en ef engin klíka getur tekið hann að sér mun honum verða skeytt inn í 0-klíkuna í lokin. Síðasta tilfellið er alls ekkert vandamál, á meðan slíkar klíkur stækka ekki of mikið, en hið fyrra gæti leitt til þess að raðirnar fyrir A í 2-klíkunni sem hann raðast í gætu verið of ónákvæmar til þess að greina afbrigði hans sérstaklega. Þetta mun þó bara leiða til false-positive, þar sem breytileikinn verður meiri innan 2-klíkunnar og auðveldara verður að greina í sundur innan hennar. Þess vegna er séns að við fengjum FP í annað hvort nálægu 0-klíkunni (finnum óvart A en áttum ekki að finna) eða 2-klíkunni sem A raðast í (áttum ekki að finna A en finnum A). Það er mjög ólíklegt að við fáum FP í báðum, ef við eigum ekki að vera að leita að A - nema við bætum þriðju klíkunni í dæmið sem er með sambærilegt samband við 0-klíkuna og 2-klíkan sem A er í - og því er þetta ásættanlegur annmarki. 

# To do list:

- Pipe to sketch all .fna files for database
    - *DONE using Sketch.cpp*
- Create rough cliques with union-find
    - *DONE using LoadFirstCliques.cpp and CliquesFromHL.cpp*
- Create finer cliques using f.x. annoy or other clustering algorithms
    - *NOT DONE*
    - Should use MDS on union-find cliques
    - Modify algorithm to place overlapping members in all cliques connected to them
    - Fixed size for cliques in respect to number of members?
- Pipe to churn .clique files into .tsv and .fasta files for bitmatrix work
    - *DONE using bitmats.cpp*
    - Needs small adjustments for multiple reference files at once...
- DCSS algorithm to choose columns from bitmatrices
    - *NOT DONE*
    - Needs to be modified to pick multiple columns with some redundancy
        - This could perhaps be specified with k = 10*members
- Algorithm to sort out which member of clique is close enough to input fastq file
    - *NOT DONE*
    - Needs to be modular with respect to cliques, multiple inputs should not scale linearly
- Modify results stage of Mince.cpp to best fit use-case of metagenomic input files
    - *NOT DONE*