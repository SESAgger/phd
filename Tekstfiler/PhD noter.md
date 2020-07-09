# PhD Noter

## Møde med Kerstin 8-5-19

#### Spørgsmål til Maja

* Hvad har vi adgang til?

  * RNA
  * Exomes
* Gene-set enrichment analysis

  * Bruger vi det?
* Hvad er Lycka, Leon og ??
* Yderligere fænotypisk info? 
* Hvor er info for resten af de 813 dyr


### Progress

#### Flat-coated Retrievers

* I have submitted most cases, some of them fails due to something some ressource being used, will try again after service window
* I have some results from the first run. I have looked at them superficially, looks good.
* If possible use Golden Retrievers as controls if possible, as they are very closely related?
* Histiocytic disease and behavioural disease are big enough to use as outcomes. 
* 

| *Phenotype* |      | Genes | Sex  | Breed | AgE  |
| ----------- | ---- | ----- | ---- | ----- | ---- |
|             |      |       |      |       |      |

Use GR to sass out HS relevant genes. Simplified they will be $FCR_m$ but $GR_{wt}$. Just add to control??  

Focus on CDK5RAP2 genes

#### Plan going forward

| Thing                         | Start   | End     | Time     |
| ----------------------------- | ------- | ------- | -------- |
| Run pipeline                  | 17/4/19 | 11/5/19 | ~1 week  |
| Run pipeline on Dobson        |         |         | < 1 week |
| Genome-wide association study | 13/5/19 | ???     |          |
|                               |         |         |          |
|                               |         |         |          |

#### Samples

* Possibly VetGruppen and Majbritts maternity-leave covers (very competent).

#### Uppsala trip

* Preferably somewhere between July and October or between March and May
* At least 1 month total, I would prefer longer, max 4 months at a time. 

### Things I want to do

#### Oncology related

* Expression patterns of (m)RNA in different cancers. 
* KIT, TP53, rb etc.





Onsdag: Hav plan. Lav summary, slides til snak. Kan vi udvide prøveindsamling (Christian), hvornår uppsala, hvor længe. Noget jeg gerne vil
mødes 9-10:30K9 pipeline

Hver work mappe har en forskellig funktion, så f.eks. d8 er quality-recalibration. 3b er gate realign. Usikker på om det er en forskellig mappe pr. run, så d8 er for CFA0001 og den svarer til d1? i CFA0002. Den nyeste er altid den uden et tal bagved (så .html og ikke .html.1)



ALTID

```
1		0a/44bd52	4618		fastqc (QJ-1615-CFA008402_S5_L003)					COMPLETED	0	2019-04-17 13:37:58.187	1h 15m 9s		1h 15m 5s		139.4%		299.7 MB	3.2 GB	39.8 GB		1.1 MB
2		8f/fcc0ea	16979		bwa (QJ-1615-CFA008402_S5_L003)							COMPLETED	0	2019-04-29-08:38:30.498	2d 20h 38m 9s	2d 20h 38m 6s	4655130.0%	4.1 GB		4.4 GB	322.3 GB	278.7 GB
3		12/b2b58c	30207		gatk_realign (QJ-1615-CFA008402_S5_L003)		COMPLETED	0	2019-04-26 10:01:49.371	5h 55m 35s		5h 55m 31s		154135.9%	3.9 GB		7.5 GB	76.9 GB		33.9 GB
4		2f/0a5e20	8856		mark_duplicates (QJ-1615-CFA008402_S5_L003)	COMPLETED	0	2019-05-02-10:12:24.756	2h 55m 15s		2h 55m 10s		781.5%		3.8 GB		7.4 GB	81.5 GB		48.2 GB

```



### gVCFcombine

### fastqc

**Time: 1h 15m 5s**

**Functions** 
fasqc -q

### bwa

**Time: 2d 20h 38m 9s**

**Functions**
bwa mem -t
samtools sort —threads
samtools index

### gatk_realign

**Time: 4h 55m 44s**

**Functions**
RealignerTargetCreator
IndelRealigner

### mark_duplicates

**Time: 2h 55m 15s**
**Functions**
MarkDuplicates

### quality_recalibration

**Time ???**
**Functions**
BaseRecalibrator
PrintReads



## Rækkefølge:

fastqc: 23/4 16 til 23/4 17 -> 1 hour:  fastqc -q 
bwa: 29/4 8:45 til 2/5 kl 5.15 -> 2 days, 20 hours, 30 min: readGroup, bwa mem -t, samtools sort, samtools index
gatk realign: 2/5 kl 5.15 til 2/5 kl 10 - > 4 hours, 45 min: Realign targetCreator, Indel realigner
mark duplicates: 2/5 10:15 til 13 - 3 hours, 45 min: MarkDuplicates
quality recalibration: 2/5 13.15 til ???: BaseRecalibrator x 2, printreads



# <span style='text-transform: none'>cfDNA</span>

One issue with cfDNA is contamination with gDNA, this can be fixed by either immediately processing blood after phlebotomy or by stabilizing the WBC. In addition, gDNA are often longer than cfDNA, theoretically making it possible to filter by length [^Norton 2013].

It has been shown that more mutated DNA fragments in cfDNA compared to CTC [^Karachaliou 2015]



###Bibliography

[^Norton 2013]: Norton, S. E. E., Lechner, J. M. M., Williams, T., & Fernando, M. R. R. (2013). A stabilizing reagent prevents cell-free DNA contamination by cellular DNA in plasma during blood sample storage and shipping as determined by digital PCR. *Clinical Biochemistry*, *46*(15), 1561–1565. https://doi.org/10.1016/j.clinbiochem.2013.06.002
[^Karachaliou 2015]: Karachaliou, N., Mayo-de-Las-Casas, C., Molina-Vila, M. A., & Rosell, R. (2015). Real-time liquid biopsies become a reality in cancer treatment. *Annals of Translational Medicine*, *3*(3), 36. https://doi.org/10.3978/j.issn.2305-5839.2015.01.16

###Litteraturstudie histiocytisk sygdom

Histiocytisk sygdom er sjældent humant, derfor vil hunde være en oplagt model da der er høj forekomst i visse racer. FCR har typisk den solitære form, hvorimod BS ofte har dissemineret. I hunde er der set variationer i MAPK og ERK pathways, værd at fokusere på disse. Det vil også passe godt sammen med BRAF V600E mutation i mennesker som ofte er påvirket evt. andre mutationer i MAPK pathwayen (men findes aldrig sammen, hvilket antyder at det er pathways-funktion der påvirkes). Mange af de gener der er forandret ekspression af er fælles for begge typer. Mus relativt dårlige modelle (as usual), ikke naturligt forekommende og genetisk mere lig. 

Mutations often shared between FCR and BMD

20 timer 

#### Noter

* Genetiske variationer ifb. histiocytær sygdom
  - Mus
    - BRAF V600E
  - Rotte
  - Menneske
    - BRAF V600E
  - Hund
    - ~~BRAF V600E ortolog~~: CanFam 3.1 16:8296227-8296345. Ikke set i hematopoietiske cancer i et studie.
      - MAPK pathway: Tjek med DAVID? Forskellige resultater afhængig af artikel. I hvert fald værd at lede efter. 
    - ERK pathway
    - SINEs/LINEs
  - FCR specifikt: Crossmatch med BS?
    - PPBP
    - SpiC
    - VCAM1
    - ENPEP
    - ITGAD
    - GTSF1
    - Col3a1
    - CD90
    - LUM
* Germline genetiske variationer i FCR
* Interessante germline mutationer
  * MDR, farver, uspecifikke arvelige sygdomme, renal dysplasi

Targets: 

* CDK5RAP2

  BRAF V600E equiv. CanFam 3.1 16:8296227-8296345

  MAPK pathway

  ERK pathway

  PPBP

  SpiC

  VCAM1

  ENPEP

  ITGAD

  GTSF1

  Col3a1

  CD90

  LUM

  Short interspersed nuclear elements (SINEs)

  Long interspersed nuclear elements (LINEs)

### Forslag til kurser

- Statistical methods in bioinformatics 

- Epidemiological methods in medical research 

- Next-Generation-Sequencing Analyse

- (Bioinformatik for Sundhed og Informatik)

- Algoritmer i bioinformatik

- **Etikkursus ** 


### Tjek ud

- Læs UPPMAX


### Noter



Delprojekter skal være færdigt hurtigere, til publikationer.

Mere ondartede?

12 mdr?

Send sverige eller her?

- Sverige



Lab arbejde:

​	Digital PCR (biorat) kig kopinr oncogener?

​	HER2 ab behandling

​	Se på duplikationsanalyse, sensitiv og meget specifik, hvor mange kopier er i prøven releativt til ref-gen.

​	Design assay? (Freeware)

​	Data skal vise variation i genkopier i det relevante gen. 
