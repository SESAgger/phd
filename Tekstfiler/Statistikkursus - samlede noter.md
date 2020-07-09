**Table of Contents**

_________________

[TOC]





# [Introduction to statistical methods for high-dimensional data, linear models and regularization methods](/Users/jqc305/Library/Mobile Documents/3L68KQB4HG~com~readdle~CommonDocuments/Documents/Phd/Bioinformatics kursus/day1_copy.pdf)



## [Klassisk statistik](#Klassisk statistik)

$outcome_i=\alpha+\beta_1 \cdot gen_i + \beta_2 \cdot køn + \beta_3 \cdot alder_i+\varepsilon$

Dvs. i praksis er ens outcome påvirket af diverse miljø-påvirkninger oven i de ^6 gener

Lineær regression er ikke brugbart til bioinformatik, da det ikke kan håndtere antallet af variable 

### Mønstergenkendelse

Vores programmer kræver en god model/algoritme, hvis den skal kunne fungere og billedet skal normaliseres på en sådan måde at den opfylder de basale krav som er nødvendige for en model/algoritme kan virke (f.eks. skal et blot justeres så punkterne rent faktisk ligger i et grid.)

### Datasæt

Datasæt har klassisk været "lange", dvs. mange cases relativt til antal variable. Bioinformatik giver "brede" datasæt med få cases og mange variable. Det er nødvendigt at reducere antallet af variable såfremt man vil benytte de "gamle" metoder, hvilket sjældent er nødvendigt. Man laver en sammenligning pr. variabel. 

##Håndtering af falske resultater

Pga. de mange test ender man med mange p-værdier, hvilket øger risikoen for falsk positive. Hvis risikoen stiger med antallet af test (humant kan man lave omkring 10<sup>6</sup> tests).

P-værdien skal sænkes, da det er håbløst at bruge $α=0.05$
$$
P(at\ least\ 1\ false\ positive)=1-(1-α)^P
$$
Husk $α=signifikans$ og $1-\beta=power$

Der er forskellige metoder til at håndtere dette problem og det er vigtigere at beskrive hvordan der er taget højde for problemet end den præcise metode. Det gælder også ift. hvorvidt man har taget højde for andre studier lavet på samme data, som i princippet også øger antallet af test og derfor bør der tages højde for disse.



###Family-wise error rate

----------------------------------

hvor $R \ = \ afviste \ hypotese$r og den eneste observerede variabel​. Resten er ikke observerede.

$FWER$[^a]: Sandsynligheden for at lave mindst en type I fejl.

$ FWER\ =\ P(V > 0)\ =\ 1-P(V=0) ≤ m \cdot α $
$m \cdot α$ = er den øvre grænse.

Den 2. ulighed holder kun ved uafhængighed, hvilket sjældent er sandt. 

Man skal huske at FWER fokuserer på p-værdier og i mange tilfælde er de egentlig ikke så interessante i sig selv. Det vi er interesseret i er om *der er* falsk positive.

####Bonferroni korrektion

|             | $H_0$ er sand | … is false | Total |
| ----------- | ------------- | ---------- | ----- |
| Afvist      | V             | S          | R     |
| Ikke afvist | U             | T          | m-R   |
| Total       | $m_0$         | $m-m_0$    | m     |



Den mest konservative metode, men er uafhængig af både dependence og distribution.

$FWER$ er stadig sandt, derfor skal signifikansen sættes til $\frac{α}{m}$ dog skal man huske at $m=antallet\ af\ tests$ så hvis generne ikke er uafhængige er $m$. 
Det betyder af vi afviser $H_i$ hvis:
$m \cdot p_i ≤ α \Leftrightarrow p_i ≤ \frac{α }{m}$
Da ligningen forudsætter at de $m$ test er uafhængige, vil $\frac{α}{m}$ være underestimeret, hvilket fører til en urimelig lav signifikans-grænse.

Et alternativ er Šídák's test, hvor signifikansniveauet er givet som:
$1-(1- α)^m = α^*  \Leftrightarrow α = 1- \sqrt[m]{1- α ^*}$

####Holm's korrektion

I praksis mindre konservativ, man regner og sorterer p-værdierne og finder $\hat{k}$, som er den laveste forkert-afviste $H_0$.
For $p_{(1)}≤p_{(2)}≤…≤p_{(m)}$
find
$\hat{k} = min\{k: p_{(k)}>\frac{α }{m+1-k}\}$ = mindste p-værdi uden signifikans
Hvis denne eksisterer, afvis hypoteser svarende til $p_{(1)},…,p_{(\hat{k})}$
Derfor:
$k ≤ m -(m_0-1)$
så
$p_{(k)}≤\frac{\alpha}{m+1-k}≤\frac{\alpha}{m+1- (m -(m_0-1))}≤\frac{\alpha}{m_0}$

Siden der findes $m_0$ sande hypoteser, så er sandsynligheden for at en af dem er signifikant max $\alpha$ og $FWER$ er dermed kontrolleret, men problemet med afhængige variable er stadigvæk til stede og dermed et forhøjet antal falsk negative.

###False discovery rate

FWER er teorien en god løsning, den er let at udregne og kontrollerer godt for falsk positive. Til gengæld er poweren så tilstrækkelig lav at det i praksis er svært at få brugbare resultater ud. Derfor kan man istedet vælge at bruge FDR.

|             | $H_0$ er sand | … is false | **Total** |
| ----------- | ------------- | ---------- | --------- |
| Afvist      | V             | S          | **R**     |
| Ikke afvist | U             | T          | **m-R**   |
| Total       | $m_0$         | $m-m_0$    | **m**     |

Den eneste kolonne der er kendt er **Total**-kolonnen.
Da: $V \ = \ falsk \ positive$ og $R \ = \ afviste \ hypoteser$ vil
$Q = \frac{V}{R}$ = proportionen af falsk positive resultater. Derfor vil FDR være gennemsnittet af Q for alle test:
$FDR = \mathbb{E}(Q)= \mathbb{E}({\frac{V}{R}})$

Bemærk at vi ikke ved hvilken $m_i$ der er falsk, eller hvordan/hvor meget den er falsk, FDR fokuserer på et andet problem end FWER, men det kan dog være mere kompliceret at udregne $\widehat{\mathbb{E}(Q)}$. En mulighed er at lave x antal test f.eks. $10^5$ og ud fra disse estimere hvor sandsynligt det er at få et givent resultat:

| Ingen sammenhæng                                             | Sammenhæng                                                   |
| ------------------------------------------------------------ | ------------------------------------------------------------ |
| ![image-20190415115100881](/Users/jqc305/Library/Application Support/typora-user-images/image-20190415115100881.png) | ![image-20190415115131672](/Users/jqc305/Library/Application Support/typora-user-images/image-20190415115131672.png) |

#### Benjamini-Hochberg step-up procedure

For $FDR=\alpha$

Udregn og sorter p-værdier:

For $p_{(1)}≤p_{(2)}≤…≤p_{(m)}$
find
$\hat{k} = max\{k: \frac{m}{k} \cdot p_{(k)} ≤ \alpha\}$

hvis $\hat{k}$ eksister så afvis hypoteser baseret på:

$p_{(1)}≤…≤p_{(\hat{k})}$
$$
\begin{array}{c}
\tilde{p} _{(1)} &=&min\{\tilde{p}_{2},mp_{(1)}\}\\
\vdots && \vdots \\
\tilde{p} _{(m-1)} &=&min\{\tilde{p}_{2},\frac{m}{m-1}p_{(m-1)}\}\\
\tilde{p}_{(m)} & = & p_{(m)}

\end{array}
$$
Desuden vil $p_i≤p_i \ for \ Holm's$ så FWER er også kontrolleret. Hvis id(i) for $m_0$ uafhængige test og sorteret så $m_0$ sand test kommer først, så vil $FDR$ være kontrolleret ved niveau $q$ som er den justerede p-værdi/$FDR$'s pendant til p-værdi og er det første punkt hvor p er signifikant:
$$
\begin{array}{l}
E(\frac{V}{R}) & = & \sum\limits_{r=1}^m \mathbb{E}[\frac{V}{r}\mathbb{1}_{R=r}]& = &\sum\limits_{r=1}^m \frac{1}{r}\mathbb{E}[V\mathbb{1}_{R=r}]\\
&=& \sum\limits_{r=1}^m \frac{1}{r}\mathbb{E}\Big[\sum\limits_{i=1}^{m_0}\mathbb{1}_{p_i≤\frac{qr}{m}}\mathbb{1}_{R=r}\Big]& = & \dots\\
&=& \sum\limits_{r=1}^m \frac{m_0}{r}\Big[\mathbb{1}_{p_i≤\frac{qr} {m}}\mathbb{1}_{R=r}\Big] & = & \dots\\
&=& \sum\limits_{r=1}^m \frac{m_0}{r}\Big[\sum\limits_{i=1}^{m_0}\mathbb{1}_{p_i≤\frac{qr}{m}}\mathbb{1}_{R=r}\Big]& = & \dots\\
&=& q\frac{m_0}{m}≤ q
\end{array}
$$


q er defineret som:
$q \ value(p_i) = \min\limits_{t≥p_i} \widehat{FDR}(t)$
Derfor er q-værdien for en individuel hypotese den mindste FDR, hvor testen kan kaldes signifikant.

Hvis $m H_0$ er sande, så: $FDR_{control}=FWER_{control}$
FDR giver generelt mere power og færre falsk positive end ukorrigeret testning. Den er også mere robust end FWER og kan bedre håndtere afhængige test, som er tilfældet for genetiske test, som aldrig er helt uafhængige (tænk codons).

##Sæt ind skema med sammenligning af metoder

## Multipel regression

$y_i=\beta_0+\beta_1\cdot x_{1i}+\beta_2\cdot x_{2i}+\beta_3\cdot x_{3i}+\beta_p\cdot x_{pi}+\varepsilon_i$

I vores tilfælde skal vi udvide til en generaliseret lineær regressions model ved at transformere med en passende funktion.

$g = link\ funktionen$ som kan mappe $\mathbb{E}(population)$ til den lineære prediktor.
$$
\begin{array}{l}
g(\mathbb{E}[Y_i|X_{1i},\dots , X_{pi}]) & = & \beta_0+\beta_1\cdot x_{1i}+\dots +\beta_p\cdot x_{pi}\\
\mathbb{E}[Y_i|X_{1i},\dots , X_{pi}] &= & g^{-1}(\beta_0+\beta_1\cdot x_{1i}+\dots +\beta_p\cdot x_{pi})\\
\mathbb{E}[Y_i|X_{1i},\dots , X_{pi}] &= & g^{-1}(X \beta) \\
\mathbb{E}[Y_i|X_{1i},\dots , X_{pi}] &= & X \beta &(for \ gaussisk \ data)\\
Dermed\\
Y=\boldsymbol{X}\beta + \varepsilon 
\end{array}
$$
Vi ender altså med en simpel lineær regression.
Least squares estimator:

$\hat{\beta}= (\boldsymbol{X}^t \boldsymbol{X})^{-1}\boldsymbol{X}^ty$
For minimering af $\beta$ for den bedste linje:
$(y-\boldsymbol{X}\hat{\beta})^t(y-\boldsymbol{X}\hat{\beta})$

###Penalized regression

For  $\boldsymbol{Y}=\boldsymbol{X}\beta+\varepsilon$  vil der for hver $\beta$ være en 'straf' og derfor ønsker man at have så få $\beta≠0$ som muligt.

####[Lasso](#Lasso)

$Z_n(\beta)=\frac{1}{n}(\boldsymbol{Y}-\boldsymbol{X}\beta)^t(\boldsymbol{Y}-\boldsymbol{X}\beta)+\lambda_n||\beta||$
$\hat{\beta}=arg\ min_{\beta\in\mathbb{R}^M}Z_n(\beta)$

Formålet er her at finde den bedste balance mellem antallet af variable og en god model, man optimerer mod en "sparse model" som stadig kan give præcise estimater.

Lasso vælger ikke variable konsekvent, dvs. for flere afhængige variable vil der være variationer af hvilken af de variable vil blive udvalgt til en given udregning. Det betyder også at variable kan blive fravalgt hvis deres værdi bliver for høj, da det ikke kan betale sig at inkludere dem i modellen. Det betyder også at alle $\beta$ er minimeret og biased mod 0. Vi bruger derfor de som ≠ 0 og laver normal multipel regression på disse for at "debiase".

####[Ridge](#Ridge)

Da Lasso har tendens til at undervurdere størrelsen af $\beta$ kan man bruge en anden metode: Ridge. Denne metode minimerer den "straffede" lineære regression ved:

$Z_n(\beta)=\frac{1}{n}(\boldsymbol{Y}-\boldsymbol{X}\beta)^t(\boldsymbol{Y}-\boldsymbol{X}\beta)+\xi_n||\beta||^2$

Dvs. den tilføjer en (lille) værdi til diagonalerne i matricet, det giver bias, men mindsker variationen. Forskellen mellem Lasso og Ridge er at for Lasso $\beta=0$ og for Ridge $\min{\{\beta}\}$

#### [Elastic net](#Elastic net)

Elastic net er en kombination af Lasso og Ridge, hvor de to indgår som specielle tilfælde hvor $\xi=0|\lambda=0$

$Z_n(\beta)=\frac{1}{n}(\boldsymbol{Y}-\boldsymbol{X}\beta)^t(\boldsymbol{Y}-\boldsymbol{X}\beta)+\lambda||\beta||+\xi_n||\beta||^2$

Benytter en 2-trins procedure, først udregnes $\lambda$ og derefter bruges denne i glmnet oftest vil man bruge 1. sd[$\lambda$]. Dette vil give en sparse-model, men som stadig er til at arbejde med, hvis man vælger $\lambda$ direkte, vil man ende med en meget kontrolleret/konservativ model, så man benytter en lidt "dårligere" model, men med mere fleksibilitet og dermed mere brugbar.

###[Principal component analyse](#Principal Component Analysis)

Ideen bag PCA er at dimensions-reducere kovariaterne og derved "strømline" variansen efter akserne. Dvs. for 3D vil man forsøge at lægge det meste af variancen langs x-aksen, den næste langs y-aksen og den sidste langs z-aksen og derved maksimere informationen pr. vektor.

Det betyder dog at hver komponent består af data fra flere forskellige kovarianter, hvilket kan gøre resultatet svært at tolke. Det man dog altid kan sige er at hvis f.eks. $PCA_1=0.9\beta_1+0.1\beta_2$ så er $PCA_1≈\beta_1$ og man kan udtale sig om $\beta_1$.

Metoden vil enten ikke påvirke størrelsen af en $PCA$ eller mindske den til 0. Man vælger de k største komponenter og bruger dem i et elastic net. Det betyder i praksis at man benytter de k største komponenter og udsmider de p-k mindste komponenter. Igen vil mængden af information pr. variabel blive maksimeret.

Man udfører den ved først at udregne kovariant-matricen for de uafhængige variable og derefter udregne eigenværdier og tilhørende vektorer, som vil svare til akserne. Man kan derefter sortere efter eigenværdier hvilket vil:

1. Give os hvilke PCA der har størst effekt.
2. Give os hvilke kovarianter som har størst effekt på de PCA der har størst effekt.

##### Algoritme:

1. Reducer dimensionerne ved at vælge en enhedsvektor (e)
2. Udskift datapunkter [$x$)] med sin procition $u^tx$
3. Hvis $varians(x)=\sum$ så $varians(u^tx)=u^t\sum u$
4. Find $u$ hvor $u^t\sum u$ er maksimeret, hvilket svarer til eigenvektoren med højeste værdi.

<img src="/Users/jqc305/Library/Application Support/typora-user-images/image-20190416111405014.png" width="50%"/>

# [Network biology](/Users/jqc305/Library/Mobile Documents/3L68KQB4HG~com~readdle~CommonDocuments/Documents/Phd/Bioinformatics kursus/dag2.pdf)

### Ordbog

Nodes: Punkterne som netværket består af
Edges: Sammenkoblingerne mellem punkterne
Degree: Hvor mange edges en node har. angiver kun det samlede antal og siger intet om hvilke typer edges det er.
Centrality: Hvor vigtig er en node for et netværk, det er både i form af antal af edges, men også om den er vigtig for at holde netværket forbundet. Denne er dog også afhængigt af hvad der er undersøgt, f.eks. vil tp53 være centralt, fordi det er undersøgt meget.
Clustering coefficient: Et mål for clustering af degrees.
Robustness: Kan netværket holde sammen, hvis en del fjernes, desuden: Er det relevant for biologi. Dvs. er der et flow som ødelægges eller er det falsk clustering. Eksempelvist vil at fjerne en del af en pathway have biologisk betydning, hvorimod fjernelse af en tilfældig clustering ikke har en biologisk effekt.
**Funktionel interaktion:** a påvirker b
	Fysisk interaktion: En direkte interaktion: a binder til b.

I netværk omtales punkterne som nodes og sammenkoblingerne som edges.

Netværk gør det muligt at overskue meget kompliceret data på en visuel og intuitiv måde, man kan lave netværk af alting. Sammenkoblinger kan være undirected eller rettet (i enten den ene eller anden retning).
$a \leftrightarrow b$ (Uspecifik binding af to proteiner) eller $a \leftarrow b = b \rightarrow a$ (a reguleres af b).
Edges kan også være weighted, så hvor sandsynligt er det at den pågældende edge dannes (f.eks. bindingsaffinitet).

Generelt kan man sige at hvis man ser linking ved f.eks. fusion af gener eller hvis de altid findes sammen, så er det sandsynligt at de samarbejder.

## Databaser

Brug altid det system som alle siger er det næstbedste (for alle vil altid sige at deres eget er det bedste). Husk at være **kildekritisk** og at data i forskellige databaser kan variere i kvalitet, navngivning, formater og derfor ikke nødvendigvis er sammenlignelige.

#### [STRING](https://string-db.org)

En database over gener, proteiner etc. 

### Quality scores

Forskellige assays for vurdering af interaktioner. De vil generelt finde forskellige subset af virkeligheden. Kvalitetsscoring varierer også mellem metoder hvilket gør det endnu sværere. 

### Cytoscape

Visualiseringsprogram 

#[Genome-wide associations tests and extensions](/Users/jqc305/Library/Mobile Documents/3L68KQB4HG~com~readdle~CommonDocuments/Documents/Phd/Bioinformatics kursus/day3real.pdf)

Formålet er at finde regioner i genomet, som er associeret med en fænotype, variationer som findes oftere i individer med bestemt fænotype og håndtere mange SNPs på samme tid. 

GWAS består af flere trin:

* Matchede, stratificerede case-control grupper
* SNP genotypning
* Find haplotyper
* Find  associations-signaler 
  * $\chi^2$eller lign. test
  * Ukorrigeret p-værdi$< 10^{-7}$ eller FDR korrektion
* Fine-mapping af aossiciations-signal
  * Genotyp SNPs i regioner hvor association er fundet
  * Fine mapping af LD i regioner
  * Find haplotyper
  * Undersøg effekt af stratificering
* Replikér resultater
* Biologisk validering af association
  * Hvilke varianter øger risikoen
  * Hvad er den funktionelle konsekvens af varianten
  * Hvad er mekanismen bag den øgede risiko.

Alleler kan findes i 0-2 kopier i et individ og er oftest relateret til omkringliggende nukleotider, og alle andre nukleotider i det pågældende gen. Hvert gen svarer til genotypen for det SNP.

$Y_i=\boldsymbol{X}\beta+\varepsilon_i=\alpha\sum\limits_{j=1}^m\beta x_{im}+\varepsilon_i$ men problem med mange prediktiorer, derrfor skal lasso, elastic net, rridge eller lign. benyttes. Desuden skal der også tages højde for interaktioner, manglende information og udefrakommende information (køn, alder, etc.)

### Standard metode

* Kør en DNA-prøve på SNP chip og mål genotyper for millioner af SNPs

* Find SNPs med signifikant association mellem allel og fænotype.

* Find kromosom-regioner hvor en haplotype er signifikant associeret med fænotype. Dette kan enten være kvantitativt eller kvalitativt.

* Lav analyse for hver SNP

  

  ![image-20190417125208524](/Users/jqc305/Library/Application Support/typora-user-images/image-20190417125208524.png)

Pga. multiple test $(\~10^6)$ vil det nødvendige signifikans niveau være: $q=5 \cdot \frac{10^{{-2}}}{10^6}=5\cdot 10^{-8}$





##Eksempel

| id   | bmi  | g1   | g2   | g3   | g4   | g5   | g6   | g7   |
| ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
| 1    | 23   | 0    | 2    | 0    | 1    | 0    | 0    | 1    |
| 2    | 31   | 0    | 1    | 0    | 0    | 2    | 0    | 1    |
| 3    | 26   | 0    | 1    | 1    | 0    | 0    | 0    | 0    |
| 4    | 35   | 2    | 1    | 1    | 0    | 2    | 1    | 0    |

$bmi_i=\boldsymbol{X}\beta+\varepsilon_i=\alpha\sum\limits_{j=1}^m\beta x_{im}+\varepsilon_i$

#### Ordbog

SNP: Enkelt nukleotid polymorfisme (i en signifikant del af befolkningen > 1%). De fungerer som markører for nærliggende gener.
INDEL: Insertions and delitions.
MAF: Minor Allele Frequency. Den mindst forekommende allel i populationen (non-wild type).

#R-kode

## Klassisk statistik

```R
#For at tage højde for multiple test kan man benytte:
p.adjust(x, method=a) #hvor x er vektor af p-værdier. method kan bl.a. være "bonferroni", "fdr" eller "holm"

#Eksempel
library("MESS")
data(superroot2)

pval <- by(superroot2,
           superroot2$gene,
           FUN=function(dat) {anova(lm(log(signal) ~ array + color + plant, data=dat))[3,5]} )
## Vælg de 8 mindste p-værdier:
sort(pval)[1:8]
##Efter holms korrektion
sort(adjust(pval, method= 'holm'))[1:8]
## Efter fdr korrektion
sort(adjust(pval, method= 'fdr'))[1:8]
```



## Lasso

```R
library('glmnet')
library('broom')
res <- glmnet(x,y,family="gaussian/binomial/andet") #standardisering bliver gjort automatisk, kan fravælges såfremt predictors har samme enheder
plot(res,lwd=2)
```

<img src="/Users/jqc305/Library/Application%20Support/typora-user-images/image-20190415145040539.png" width="50%" align="left"/>

```R
coef(res,s=2)#Svarer til λ=2
pick <- which(coef(res,s=2)!= 0) #finder hvilke rækker som ≠ 0, den første er altid intercept, så derfor -1
select <- pick[-1]-1
select
broom::tidy(lm(y ~ x[,select]))
#Tjek selectiveInference pakken ud om et stykke tid.
```

## Ridge

```R
library('glmnet')
res <- glmnet(x,y,alpha=0)#alpha=0 er ridge og alpha=1 er lasso, standard er 1, den er ikke binær
plot(res,lwd=2)
```

<img src="/Users/jqc305/Library/Application Support/typora-user-images/image-20190416104212899.png" width="50%" align="left"/>



## Elastic net

```R
res <- cv.glmnet(x,y) #Kryds-validerer k=10 gange
plot(res,lwd=2)
```

<img src="/Users/jqc305/Library/Application Support/typora-user-images/image-20190416105020356.png" width="50%"/>

## Principal Component Analysis

```R
library(MASS)
data(biopsy)
names(biopsy)
predictors <- biopsy[complete.cases(biopsy),2:10]
fit <- prcomp(predictors, scale=TRUE)
summary(fit) #Viser os hvor stor en andel af Variansen som en given PCA står for.
plot(fit) #Viser os et histogram over varians pr. PCA. Giver en god visuel forståelse af predictornes effekt
biplot(fit) #Viser hvilke componenter som primært påvirker en given PCA. Kan dog lettere læses af daten selv
fit
```

## Genome Wide Association Study

```R
library("data.table")   
library("MESS")# Get packages
geno <- fread("hapmap1.ped")  # Read file
for(j in seq_along(geno)){    # Convert 0 to NA+          
  set(geno, i=which(geno[[j]]==0), j=j, value=NA)+ }
pheno <- fread("qt.phe")      # Read phenotypes
geno$V1 <- NULL ; geno$V2 <- NULL ; geno$V3 <- NULL ;> geno$V4 <- NULL ; geno$V5 <- NULL ; geno$V6 <- NULL> # Collapse pairs of alleles to genotype
geno2 <- geno[, lapply(1:(ncol(.SD)/2),+ 
                       function(x) sum(.SD[[2*x-1]], .SD[[2*x]])-2),+               												by = 1:nrow(geno),+
                       .SDcols = grep('^V', names(geno), value = TRUE)]
keep <- which(apply(geno2, 2, sd) > 0) # Remove those with no var
geno2 <- geno2[,keep, with=FALSE]
geno2$nrow <- NULL
                       
mres <- mfastLmCpp(pheno$V3, as.matrix(geno2))
pval <- 2*pt(-abs(mres[,3]), df=nrow(pheno)-2)
head(sort(pval))
	[1] 5.272596e-09 1.095443e-06 1.630357e-05 2.885730e-05 3.580539e-05
	[6] 3.715969e-05
head(names(geno2)[order(pval)])
	[1] "V10602" "V81525" "V37137" "V12225" "V18546" "V53636"
```





Standard analyse til FCR

[^a]:Family-wise error rate
