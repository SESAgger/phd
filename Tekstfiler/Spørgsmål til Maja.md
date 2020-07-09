# Spørgsmål

lav i autogen hvor man kan vælge hvor lang tid ny.sh kører.

Brug Ensembles variant effect annotater: Viser alle annotationer og hvilken effect det kan have. Læs om. 

CFA=Biobank nr. kan matches i summary til Agria med fænotype.

Kan vi filtrere vores DNA baseret på fragment-længde? Fleste cfDNA fragmenter er 150-180bp. Ville mindske mængden af gDNA i prøver WP3

## WP2

* WP2 efter 40, så sekvenser

  * [ ] Analyser

    

WP1/HERR2 ligger i fryser.

LSA+OSA+MMT: Ligger noget data. OSA fra USA 									

### Spørgsmål

* Laver vi den CT eller hvad?


### Mangler svar

- Ændrer vi active region threshhold?
- Skal vi bruge active region til specifikt gen? (F.eks. HER2)
- Bruger vi neural network?
- Undervisningsvirksomhed
- Hvor mange klinikdage.
  - 15 pr 6 mdr.
- WP 2, collection of 40 samples form 1 of the tumor types/breeds
  - Which one(s)?
- Mammatumor longitudinel
  - ​	Aggressive tumorer
    - mets, str, vækst
  - Hvordan rent praktisk?
    - Udtag bp+cyto/histo ved op
      - gå ud fra deltager
    - Screen for ctDNA eller afvent histologi og screen først derefter
- Hvad giver vi gratis?
  - Opfølgende besøg
    - Evt. RTG/UL
    - Ikke andet gratis
- Truth sets? Nope
- FCR: Germline SNP+indels??
  - HaplotypeCaller -> GVCF-> Consolidation ->Joint-Call -> vcf
- Hvorfor gatk3?
  - Fordi bare
- Hvilken bwa? sampe/**mem**/andet

Kommandoer:

* I bwa mem -t hvordan vælges antallet af threads?
* AddOrReplaceGroup:
  * -Xmx16g???
    * scilifelab
  * Funktion af sort_order?
    * coordinate er krav fra BuildBam index, som bruges til at lave fastlookup af .bam fil
  * 
    * rglb=exome lib, hvordan vælges? hvilke muligheder? bruges det til germline wgs?
  * 
* CollectHsMetrics
  * kun til exome data? Ja

### Besvarede

- [Er Bayesian Analysis relevant?](http://phdcourses.dk/Course/64284#.XDR55S2ZM1I)
- Forslag til kurser?
  - Bioinformatik
  - Statistik
  - Kurser i Uppsala eller Broad? (Tjek SciLifeLab jævnligt)
- Hvad gør vi med patienter med både maligne og benigne? Især hvis cfDNA fra tumor ==findes==  i blodet?
- Parallelprøveanalyse
- Foretrukken Bib-program?
  - endnote
- Hvilket sekvenseringsmaskiner benyttes?
  - Illumina
  - Evt. pacbio
- Britisk eller amerikansk?
  - US
- Skal UU og SLU også på consent form?
  - Nope
- MuTect2 kun β, what's up with that? Ikke anbefalet til production, regner vi med at det bliver prod, eller siger vi bare fuck it??
  - De lyver
- Bruger vi firecloud? 
  - GUI, brug command-line
- Hvordan bruges pladsen? Er det udelukkende på hvilken mappe man står i?
  - kaldes i sbatch.sh fil
- Hvorfor SAM-filer? Virker ikke som om gatk teamet selv bruger det?
  - Hvis muligt uden, så prøv/gør
- dog_dbSNP.recode.vcf? 
  - Muligvis skal der gøres mere ved den



## Hvad mangler ansøgningen?

- Budget and funding
  - Alt
- Activities and confirmations
  - Expected work obligations
- Supervisors
  - Fødselsdatoer

bwa mem -t 8 /proj/snic2018-3-256/nobackup/Resources/canFam3.fa /proj/uppstore2017228/KLT.05.GRUS/QJ-1615/171205_ST-E00215_0251_AHF257CCXY/Sample_QJ-1615-CFA007614/QJ-1615-CFA007614_S1$







Ejer skal godkende at vi må tage prøver specifikt til vores projekt pga. tidsskrifter

Møde hver uge i starten, find dag

HVad vil jeg lave, hvad har jeg lavet, isssues?

Uppsala ved svoc

Capture til hund, ikke illumina, kun humant. Nimblegen, dyrt, dårlig hitrate

Kan vi tage prøver fra alle tumorer. 

RNAlater, sp'rg Annemarie om rørkøb. 

Indel måske ikke brugbart, dårlig mapping, drop?

dbSNP, sammenlign veronicas, med DATA_merge fil

Veronikas ligger på snic et sted.

UPP har hjælpefunkt, spørg om online





breed comparisons

do assay herr2

RNA from tumors already submitted. 

Ask Eva where stuff is

300kr pr sample

funds for danish samples

Review of all tumor normal data. 



**DO**
Webpage
	DUG
	ASK QUESTIONS
	Send cheekswaps out
	Darwins arc webpage, look.

Questionaire for dogs

FCR

GWAS MMT

GWAS cross-breeds OSA. 

Nomenclature: stat-significance, nominal-significance

What software for gwas?

