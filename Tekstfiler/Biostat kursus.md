# Mandag

The `p.adjust()` function takes a vector of p-values as input and returns a vector of p-values adjusted for multiple comparisons. The `method` argument sets the correction method, and the values `"holm"` for Holm’s correction (the default), `"bonferroni"` for Bonferroni correction, and `"fdr"` (for false discovery rate) are among the possibilities. The Holm and Bonferroni methods both seek to control the overall false positive rate with Holm’s method being less conservative than the Bonferroni method. The `fdr` method controls the false discovery rate which is a less stringent criterion than the overall false positive rate.

We illustrate the `p.adjust()` function with a typical microarray gene expression experiment found as the `superroot2` data frame in the `MESS` package. The expression levels of 21,500 different genes are examined and for each gene we wish to determine if there is any difference in gene expression between two plant types (mutant and wild type) after correcting for dye color and array. A standard analysis of variance model is used to test the hypothesis of no difference in gene expressions between the two types of plants for each of the 21,500 genes. The `by` function is combined with `anova` and `lm`to calculate the 21,500 p-values.



sort(pval)[1:8]

sort(adjust(pval, method= 'holm'))[1:8]

sort(adjust(pval, method= 'fdr'))[1:8]

$bonf_1=holm_1$ fordi de starter med det samme niveau (m)

1. How many genes would you expect to be deemed ``significant’’ by chance with a standard p-value cut-off of 0.05?

   1. 1075

2. In the code below what is going on in the line starting with 

3. Run the code to obtain the p values.

4. Would you conclude that there are any significantly differentially expressed genes?
   1. No

5. Does this conclusion change for any of the methods?

   1. No, all > 0.1

6. Do the Bonferroni correction yourself and see how that changes your results.
   1. Changes as 21499 p-values.
   2. Why are the smallest p-values for the FDR correction all identical?

   1. Always chooses the smallest one and brings it along, pga.  * 0?

7. When might these multiple testing methods be inadequate?

   1. 

8. Why might it not always be a good idea to focus on p values (think about prediction and discovery, and the experimental setup).

   1. Vi vil hellere fokusere på fdr, for vi tror ikke rigtig på p værdier. Ikke interessant



Man kunne lave en "stop regning" hvis ingen boef er sign. for de andre vil ikke give bedre.

## Sæt 2

1. Fit a lasso model to the genotype/phenotype data.
2. Is it necessary to standardize the input data before running the analysis? [Hint: look at the `standardize` argument]
   1. Nej, den gør det automatisk
      Logical flag for x variable standardization, prior to fitting the model sequence. The coefficients are always returned on the original scale. Default is `standardize=TRUE`. If variables are in the same units already, you might not wish to standardize. See details below for y standardization with `family="gaussian"`.
3. Why would it normally make sense to standardize the columns of the  predictors? Explain what might happen if we do not and how the penalty  will influence the different predictors.
   1. Hvis ikke normalized risikere man at f.eks. amylase og alder begge er predictors, så fucker man en af dem op.
4. Use cross-validation to obtain a reasonable estimate for the penalty parameter.
5. Use the estimated penalty parameter and extract the corresponding list of coefficients. How many predictors are selected?
   1. 9

Above we just considered the outcome continuous (even though it is a  series of 0s and 1s). A better model would be to use a binomial model  like logistic regression. To analyze a dichotomous outcome such as  case/control status we use `family=”binomial”`.

1. Try to do that and compare the results. What should/shouldn’t you be looking for here?
2. Run the same analysis using ridge regression and compare to the lasso results.
3. Although none of the parameters are set to zero for ridge regression  would you still think it would be possible to at least get information  about a sparse solution? How? [Hint: this is an *ad hoc* question/answer so just state a general idea]
4. Run the same analysis using elastic net and compare to the previous results.



# Principal component analysis

1. What is the interpretation of the results that you have obtained?
   1. At PC1->PC9 er faldende i 'variationsandel'
2. Try using `plot(fit)` and `biplot(fit)` and  explain what you see. Make sure you produce these two plots and  understand what they show Especially the biplot is important.
   1. Viser hvor stor andel som de forskellige predictorer står for
   2. plot(fit) giver en god visuel forstøelse for PC delene. Hvis meget skarp skift i plot, vælg kun dem før a aftager (a*x+b)f
3. Which variable(s) provide(s) any real contribution to the second principal component? 
   1. V9 kan ses på graf (0.905)
   2. Ellers kig på fit dataen i sig selv.
4. Look at the loading factor for the result. Do you see any patterns in the composition of the principal components?



# Tirsdag

## Exercise 1

### 1.1 Single protein query

We will first retrieve a STRING network for human insulin receptor (INSR). Go to <https://string-db.org/>, open the **Protein by name** search interface, and type **INSR** in the field **Protein Name**. You can either specify **Homo sapiens** in the **Organism** field leave it on **auto-detect**. Click **SEARCH**.  If you specified the organism, you will immediately receive a protein  network; otherwise you will first be presented with a disambiguation  page on which you can specify that you meant the human protein.

### 1.2 Visual representations

The STRING web interface provides several different visual representations of the network. The **Settings**  tab below the network view allows you to change between different  visual representations of the same network. Try changing between the **confidence** and **evidence** views; do not forget to press the **UPDATE** button.

*Which information is shown for the edges in each representation?  Why are there sometimes multiple lines connecting the same two proteins  in the evidence representation?*

* Evidence =  line color indicates the type of interaction evidence
  * Evidence for flere forskellige typer af interaktioner?
* Confidence = line thickness indicates the strength of data support
  Better for complex networks, makes the picture less messy. Removes evidence type though



## 1.3 Evidence viewers

A key feature of the STRING web interface is the evidence viewers.  One should not rely purely on the confidence scores; it is important to  inspect the actual evidence underlying an interaction before relying on  it, for example, for designing experiments.

*Which types of evidence support the interaction between insulin receptor (INSR) and insulin receptor substrate 1 (IRS1)?*

​	Experiments, databaser og textminin.

Further detail on the evidence of an interaction can be seen in a  popup by clicking on the corresponding edge in the network. Click on the  edge between INSR and IRS1 to view its popup; you may need to move the  nodes to make this easier.

*Which type of evidence gives the largest contribution to the confidence score?*

​	Association in Cur db, derefter Exp and then co-mentioning in PubMed

Click on the **Show** button to view the experimental evidence for the interaction.

*Which types of experiments support this interaction?*

​	Pull-down
​	2-hybrid assays
​	affinity chromatography technology
​	etc.	

### 1.4 Query parameters

The **Settings** tab also allows you to modify detailed parameters for the search, such as the types of evidence to use (**active interaction sources**), the **minimum required interaction score**, and the **max number of interactors to show**.

Change the **minimum required interaction score** to high confidence (0.700).

*Does this change the set of proteins shown? Does it change the interactions shown?*

​	Næh

Turn off all evidence types except experiments.

*Does this change the set of proteins shown in the network?*

​	Yes, fjerner en og færre interaktionstyper ses

​	Kan være relevant ift. fysiske interactione
​	Kan også være problematisk hvis man både tager eksperimenter og også papaers om de experimenter.

Increase the **max number of interactors to show** to 20.

*How many interaction partners of INSR do you get in the network?*

Change the **minimum required interaction score** back to 0.400.

*How many INSR interactors do you now get?*
	9

Turn on all evidence types back on.

*Does this change the set of proteins shown in the network?*
	A lot of proteins now

## Exercise 2

### 2.1 Chemical compound query 

Bruger 1-5 fordi dårlig gold-standard. If no. then people use it as propability, men det er det ikke.

To retrieve for a chemical compound, go to <http://stitch-db.org/>, type **aspirin** in the **Item Name** field, type **Homo sapiens** in **Organism**, and click the **SEARCH**  button. You should now see a network of proteins associated with  aspirin, including PTGS1 and PTGS2, which are the direct targets of  aspirin.

### 2.2 Binding assay data

Click the interaction between aspirin and PTGS1 and then the **Show** button next to **Experimental/Biochemical Data**  in the popup. The evidence viewer now shows all the available binding  assay data for this interaction, including the binding affinities (IC50  values).

## Exercise 3

### 3.1 Disease query

Go to <https://diseases.jensenlab.org/>, type **Parkinson** in the search field, and click the **search**  button. The web interface will now show the search results, which  include all diseases and protein names starting with the search term.  Click **Parkinson’s disease** to get to the results page showing proteins associated with the disease.

Like STRING and STITCH, the DISEASES database integrates several  types of evidence, in this case automatic text mining, manually curated  knowledge, and experimental evidence from genome-wide association  studies.

*Are there any proteins that are supported by all three types of evidence?*
	SNCA f.eks



Text-mining lettere at kvalitetsvurderer. Computere fejltolker ofte, dobbelttjek.

### 3.2 Validation of text mining

The DISEASES database too allows you to inspect the underlying  evidence for an association. Since the predominant source of evidence is  automatic text mining, it is always wise to read the underlying text to  manually validate the results. Click on **SNCA** in the **Text mining** table to view the text based on which it was associated with the disease. Click **View abstract** for a given entry to see the complete abstract rather than only the title.

*Do the abstracts all mention both the protein and the disease? Do they all use the same name for the protein?*
	Nej

## Exercise 4

### 4.1 Complete virus query

Go to <http://viruses.string-db.org>, and select **Complete Virus** from the menu on the left.  In the **Virus** dropdown, enter “Measles”, and the **Host** dropdown can be left as auto-detect to detect the host with the most interactions, in this case, Homo sapiens.  Click **Search** to retrieve the network.

### 4.2 Inspect virus evidence

Click on the edge connecting the measles virus P/V protein and the human STAT2 protein.

*What types of evidence support an interaction between these proteins?  List two publications that the evidence comes from.*
	Exp and CO-mentioning
	"Morbillivirus Experimental Animal Models: Measles Virus Pathogenesis Insights from Canine Distemper Virus. 
*da Fontoura Budaszewski R, von Messling V*
Viruses. 8(10) (2016)."

"Host-Pathogen Interactions in Measles Virus Replication and Anti-Viral Immunity. 
*Jiang Y, Qin Y, Chen M*
Viruses. 8(11) (2016)"

### 4.3 Single virus protein query

Click on the logo at the top of the page to go back to the main search screen.  Select **Virus by Single Protein** from the left, and then enter “P” as the **Virus Protein Name** and “bacteriophage lambda” as the **Virus**.  The **Host** can again be left to auto-detect E. coli.  Protein P is responsible for the bi-directional replication of phage DNA.

*Which host proteins does P interact with?  What types of evidence supports these interactions?*![image-20190409100329377](/Users/jqc305/Library/Application Support/typora-user-images/image-20190409100329377.png)





# Cytoscape

## Exercise 1

In this exercise, we will perform some simple queries to retrieve  molecular networks based on a protein, a small-molecule compound, a  disease, and a topic in PubMed.

### 1.1 Protein queries

Go to the menu **File → Import → Network from Public Databases**. In the import dialog, choose **STRING: protein query** as **Data Source** and type your favorite protein into the **Enter protein names or identifiers** field (e.g. SORCS2). You can select the appropriate organism by typing the name (e.g. Homo sapiens). The **Maximum number of interactors**  determines how many interaction partners of your protein(s) of interest  will be added to the network. By default, if you enter only one protein  name, the resulting network will contain 10 additional interactors. If  you enter more than one protein name, the network will contain only the  interactions among these proteins, unless you explicitly ask for  additional proteins.

Unless the name(s) you entered give unambiguous matches, a  disambiguation dialog will be shown next. It lists all the matches that  the stringApp finds for each query term and selects the first one for  each. Select the right one(s) you meant and continue by pressing the **Import** button.



*How many nodes are in the resulting network? How does this  compare to the maximum number of interactors you specified? What types  of information does the **Node Table** provide?*

### 1.2 Compound queries

Go to the menu **File → Import → Network from Public Databases**. In the import dialog, choose **STITCH: protein/compound query** as **Data Source** and type your favorite compound into the **Enter protein or compound names or identifiers**  field (e.g. imatinib). You can select the organism and number of  additional interactors just like for the protein query above, and the  disambiguation dialog also works the same way.



*How is this network different from the protein-only network with respect to node types and the information provided in the **Node Table**?*

### 1.3 Disease queries

Go to the menu **File → Import → Network from Public Databases**. In the import dialog, choose **STRING: disease query** as **Data Source** and type a disease of interest into the **Enter disease term**  field (e.g. Alzheimer’s disease). The stringApp will retrieve a STRING  network for the top-N proteins (by default 100) associated with the  disease.

The next dialog shows all the matches that the stringApp finds for  your disease query and selects the first one. Make sure to select the  intended disease before pressing the **Import** button to continue.



*Which additional attribute column do you get in the **Node Table** for a disease query compared to a protein query? (Hint: check the last column.)*

### 1.4 PubMed queries

Go to the menu **File → Import → Network from Public Databases**. In the import dialog, choose **STRING: PubMed query** as **Data Source** and type query representing a topic of interest into the **PubMed Query**  field (e.g. jet-lag). You can use any query that would work on the  PubMed website, but it should obviously a topic with related genes or  proteins. The stringApp will query PubMed for the abstracts, find the  top-N proteins (by default 100) associated with these abstracts, and  retrieve a STRING network for them.



*Which attribute column do you get in the **Node Table** for a PubMed query compared to a disease query? (Hint: check the last columns.)*

### 1.5 New search interface

The types of queries described above can alternatively be performed  through the new Cytoscape search interface. Click on the drop-down menu  with an icon on it, located on the left side below the **Network** tab in the **Control Panel**.  Select one of the four possible STRING queries and directly enter your  query in the text field. To change settings such as organism, click the ☰  button next to the text field. Finally, click the 🔍 button to retrieve  a STRING network for your query.

## Exercise 2

In this exercise, we will work with a list of 78 proteins that  interact with TrkA (tropomyosin-related kinase A) in neuroblastoma cells  10 min after stimulation with NGF (nerve growth factor) ([Emdal et al., 2015](http://stke.sciencemag.org/content/8/374/ra40)). An adapted table with the data from this study is available [here](https://jensenlab.org/assets/stringapp/Emdal2015SciSignal.tsv).

### 2.1 Protein network retrieval

Start Cytoscape or close the current session from the menu **File → Close**. Go to the menu **File → Import → Network from Public Databases**. In the import dialog, choose **STRING: protein query** as the **Data Source** and paste the list of UniProt accession numbers from the first column in the table into the **Enter protein names or identifiers** field.

Next, the disambiguation dialog shows all STRING proteins that match  the query terms, with the first protein for each query term  automatically selected. This default is fine for this exercise; click  the **Import** button to continue.



*How many nodes and edges are there in the resulting network? Do the proteins all form a connected network? Why?*

​	Mange
​	Nej, kun mange af dem.

Cytoscape provides several visualization options under the **Layout** menu. Experiment with these and find one that allows you to see the shape of the network easily. For example, you can try the **Degree Sorted Circle Layout**, the **Prefuse Force Directed Layout**, and the **Edge-weighted Spring Embedded Layout**.



*Can you find a layout that allows you to easily recognize  patterns in the network? What about the Edge-weighted Spring Embedded  Layout with the attribute ‘score’, which is the combined STRING  interaction score.*

### 2.2 Discrete color mapping

Cytoscape allows you to map attributes of the nodes and edges to  visual properties such as node color and edge width. Here, we will map  drug target family data from the [Pharos](https://pharos.nih.gov/idg/targets) database to the node color.

Select **Style** from the top menu in the left panel (it is between **Network** and **Select**). Click the **◀** button to the right of the property you want to change, in this case **Fill Color**, and set **Column** to the node column containing the data that you want to use (i.e. **target family**).  This action will remove the rainbow coloring of the nodes and present  you with a list of all the different values of the attribute that are  exist in the network.

*Which target families are present in the network?*

​	Kinases
​	GPCR

To color the corresponding proteins, first click the field to the right of an attribute value, i.e. **GPCR** or **Kinase**, then click the **⋯**  button and choose a color from color selection dialog. You can also set  a default color, e.g. for all nodes that do not have a target family  annotation from Pharos, by clicking on the white button in the first  column of the same row.

*How many of the proteins in the network are kinases?*

Note that the retrieved network contains a lot of additional  information associated with the nodes and edges, such as the protein  sequence, tissue expression data (Node Table) as well as the confidence  scores for the different interaction evidences (Edge Table). In the  following, we will explore these data using Cytoscape.

### 2.3 Data import

Network nodes and edges can have additional information associated  with them that we can load into Cytoscape and use for visualization. We  will import the data from the [text file](https://jensenlab.org/assets/stringapp/Emdal2015SciSignal.tsv).

To import the node attributes file into Cytoscape, go to **File → Import → Table from File**.  The preview in the import dialog will show how the file is interpreted  given the current settings and will update automatically when you change  them. To change the default selection chosen by Cytoscape, click on the  arrow in the column heading. For example, you can decide whether the  column is imported or not by changing the **Meaning** of  the column (hover over each symbol with the mouse to see what they  mean). This column-specific dialog will also allow you to change the  column name and type.

Now you need to map unique identifiers between the entries in the  data and the nodes in the network. The key point of this is to identify  which nodes in the network are equivalent to which entries in the table.  This enables mapping of data values into visual properties like Fill  Color and Shape. This kind of mapping is typically done by comparing the  unique Identifier attribute value for each node (Key Column for  Network) with the unique Identifier value for each data value (key  symbol). As a default, Cytoscape looks for an attribute value of ‘ID’ in  the network and a user-supplied Key in the dataset.

The **Key Column** for Network can be changed using a  combo box and allows you to set the node attribute column that is to be  used as key to map to. In this case it is **query term**  because this attribute contains the UniProt accession numbers you  entered when retrieving the network. You can also change the Key by  pressing the key button for the column that is to be used as key for  mapping values in the dataset.

If there is a match between the value of a Key in the dataset and the  value the Key Column for Network field in the network, all  attribute–value pairs associated with the element in the dataset are  assigned to the matching node in the network. You will find the imported  columns at the end of the Node Table.



### 2.4 Continuous color mapping

Now, we want to color the nodes according to the protein abundance  (log ratio) compared to the cells before NGF treatment. From the left  panel top menu, select **Style** (it is to the right of **Network**). Then click on the **◀** button to the right of the property you want to change, for example **Fill Color**. Next, set **Column** to the node column containing the data that you want to use (10 min log ratio). Since this is a numeric value, we will use the **Continuous Mapping** as the **Mapping Type**,  and set a color gradient for how abundant each protein is. The default  Cytoscape color gradient blue-white-red already gives a nice  visualization of the negative-to-positive abundance ratio.

*Are the up-regulated nodes grouped together?*
	Ish

To change the colors, double click on the color gradient in order to bring up the **Continuous Mapping Editor**  window and edit the colors for the continuous mapping. In the mapping  editor dialog, the color that will be used for the minimum value is on  the left, and the max is on the right. Double click on the triangles on  the top and sides of the gradient to change the colors. The triangles on  the top represent the values at which the data will be clipped;  anything above the right triangle will be set to the max value. This is  useful if you have a small number of values that are significantly  higher than the median. To have three colors, you need to add a new  triangle (for the white color) by pressing the Add button and set the  Handle position value to 0. As you move the triangles and change the  color, the display in the network pane will automatically update — so  this is all easier to do than to explain! If at any point it doesn’t  seem to work as expected, it is easiest to just delete the mapping and  start again.

*Can you improve the color mapping such that it is easier to see which nodes have a log ratio below -2 and above 2?*
	Set min og max til ±2

### 2.5 Functional enrichment

Next, we will retrieve functional enrichment for the proteins in our  network. After making sure that no nodes are selected in the network, go  to the menu **Apps → STRING Enrichment → Retrieve functional enrichment**  and keep the default p-value of 0.05. A new STRING Enrichment tab will  appear in the Table Panel on the bottom. It contains a table of enriched  terms and corresponding information for each enrichment category.



*Which are the three most statistically significant terms?*

To explore only specific types of terms, e.g. GO terms, and to remove  redundant terms from the table, click on the filter icon in the **Table panel** (leftmost icon). Select the three types of GO terms, enable the option to **Remove redundant terms** and set **Redundancy cutoff**  to 0.2. In this way, you will see only the statistically significant GO  terms that do not represent largely the same set of proteins within the  network. You can see which proteins are annotated with a given term by  selecting the term in the **STRING Enrichment** panel.



*Do the functional terms assigned to a protein correlate with whether it is up- or down-regulated?*

Next, we will visualize the top-3 enriched terms in the network by  using split charts. Click the settings icon (rightmost icon) and set the  **Number of terms** to chart to 3; you can optionally also **Change Color Palette** before clicking **OK** to confirm the new settings. Click the colorful chart icon to show the terms as the charts on the network.

 

*Do these terms give good coverage of the proteins in network?*
	Jae

Finally, save the list of enriched terms and associated p-values as a text file by going to **File → Export → Export STRING Enrichment**.

## Exercise 3

We are going to use the stringApp to query the [DISEASES](https://diseases.jensenlab.org)  database for proteins associated with Parkinson’s disease and with  Alzheimer’s disease, retrieve STRING networks for both, created a  combined network for the two neurodegenerative diseases, identify  clusters in the network, and color it to compare the diseases.

### 3.1 Disease network retrieval

Close the current session in Cytoscape from the menu **File → Close**. Use the menu **File → Import → Network from Public Databases** and the **STRING: disease query** option. Retrieve a network for **Parkinson’s disease** and another for **Alzheimer’s disease**.

 

### 3.2. Working with node attributes

Browse through the node attributes table and find the disease score  column. Sort it descending by values to see the highest disease scores.  You can highlight the corresponding nodes by selecting the rows in the  table, bringing up the context menu (right-click the selected rows) and  choosing the ‘Select nodes from selected rows’ option. You can also use  one of the icons in the menu to zoom into the selected node.

Look for an example for a node with a disease score of 5 and one with a disease score below 4.

Rename the **disease score** columns in the two networks to **PD score** and **AD score**, respectively, by right-clicking the name and choosing the **Rename column** option.

 

### 3.3 Merging networks

Cytoscape provides functionality to merge two or more networks,  building either their union, intersection or difference. To merge the  two disease networks go to the menu **Tools → Merge → Networks**. Select the two disease networks in the **Available Networks** list and move them to the **Networks to Merge** list by clicking the **>** button. Make sure the **Union** button is selected and click the **Merge** button.

*How many nodes and edges are there in the merged network compared to the two individual disease networks?*

​	HSA har 45 osa har 45

​	87 sammen

Because the merged network was not created by the stringApp, but  rather by Cytoscape’s merge tool based on two separately retrieved  STRING networks, we now have two problems. First, Cytoscape does not  know the merged network is a STRING network, and most menu points in the  stringApp menu are thus grayed out; fix this by going to the menu **Apps → STRING → Set as STRING network**.  Second, because the two disease networks were retrieved separately, the  merged network does not contain any interactions between proteins  involved only in Parkinson’s disease and proteins involved only in  Alzheimer’s disease, even if the proteins interact according to STRING.  To solve this, first go to the menu **Apps → STRING → Change confidence menu**, set the **New confidence cutoff** to 1, and press **OK**;  this will remove all STRING interactions, leaving only the proteins.  Then bring up the same dialog and lower the confidence cutoff back down  to 0.4; the stringApp will now query the server again to retrieve  interactions among all the proteins.

 

*How many nodes and edges does the network contain compared to before?*

87 

### 3.4 Network clustering

Next, we will use the MCL algorithm to find identify clusters of  tightly connected proteins within the combined network. Go to the menu **Apps → clusterMaker → MCL Cluster**. Set the **Granularity parameter (inflation value)** to 4 and choose the **score** attribute (i.e. the overall STRING confidence score) as **Array Sources**, select the option **Create new clustered network**, and click **OK** to start the clustering. The app will now run the algorithm and automatically create a network showing the clusters.



*How many clusters have at least 4 nodes?*

Kun 2

### 3.5 Selection filters and style bypass

Finally, we want to color the nodes based on which disease(s) the  proteins are involved in. To do so we will make use of selection filters  to select nodes based on their attributes and the visual style bypass  to explicitly specify the colors of individual nodes.

First, to select all proteins associated with Alzheimer’s disease, go to the **Select** tab, click the **+** button, and choose **Column Filter**. Select **Node: AD score**  in the drop-down menu and the 100 nodes associated with Alzheimer’s  disease will be highlighted yellow in the network. Next, go to the **Style** tab and click the square just left of **Fill color** to set a bypass color for the selected nodes (e.g. red). Repeat this process for the node attribute column **Node: PD score** and choose different color (e.g. blue).

To select the nodes that are shared between the diseases, go to the **Select** tab and create two column filters, one for **Node: AD score** and one for **Node: PD score**. Since the filters by default are combined using the **Match all (AND)** rule, only the nodes that have both scores will be selected. Go to the **Style** tab to set a third bypass color for this last group.

*How are the two diseases distributed across the clusters?*
	Meget lidt overlap 

## Exercise 4

In this exercise, we will retrieve virus-host networks for two  closely related viruses, merge them into a single network, and then will  retrieve the functional enrichment for the host proteins in this  network.

### 4.1 Virus queries

Go to the menu **File → Import → Network from Public Databases**. In the import dialog, choose **STRING: protein query** as the **Data Source**.  As of version 1.4 of the STRING app, 236 virus species are included in  the species dropdown menu. Since most viruses are small (they have a  median of 9 proteins in their genomes) it is reasonable to import **all proteins of this species**  for a given virus, so select this checkbox underneath the species  dropdown. For this example we will query all proteins of “Human  papillomavirus type 16 (HPV 16)”. Simply type HPV 16 and select the  species from the resulting shorter dropdown menu.

*How many virus proteins are encoded for by this virus? What node information is imported along with the names of the proteins?*
	8
	Sekvens + alt muligt andet

### 4.2 Expand with host interactors

To retrieve interactions with host proteins, go to **Apps → STRING → Expand network**. In the resulting dialog, enter the number of desired host proteins, and select the host species from **Type of interactors to expand network by**.  All host species for which we have interactions with the currently  imported virus genes, will be shown in the dropdown menu. The **selectivity of interactors**  can also be specified – we recommend a default value of 0.5, but you  can move the slider towards 0 to decrease the number of network-sepcific  interactors or towards 1 to increase it. In this example, we will  import 10 human proteins, and keep the default selectivity.

The resulting network will be automatically re-styled such that the  nodes representing virus proteins are red and host proteins are  green-blue. These attributes can be changed from the Cytoscape Style  menu.

*Which human protein has the highest interaction score to one of  the virus proteins? What cellular functions is this protein involved in?  (Hint: open the results panel under **Apps → STRING → Show results** panel.)*
	Cadherin -> cell-binding

Additional viruses or hosts can be added to the network by iterating  on this procedure, but this will only add proteins that interact with  the proteins that are already in the network. This will work fine when  adding new hosts, since all virus proteins are already in the network.  However to add new viruses, we recommend merging the expanded networks  for each virus.

### 4.3 Add specific host proteins

If a specific host protein is desired, it can also be included in the network from the **Apps → STRING → Query for additional nodes**  menu option. In this example, p53 is not one of the proteins that was  included in the network in the previous step, however it is known that  the HPV E6 protein mediates ubiquitination of p53. To include this  protein, choose “Homo sapiens” for the species (you may have to scroll  up in the list), and enter “tp53” into the text area box in the dialog,  then click **Import**.

*Which HPV proteins does p53 interact with?*
	E6+E7

Note that p53 will be added to the network in the previous step if  more proteins are imported or the selectivity is set to a lower value.  Choosing a lower selectivity will include more hub proteins (such as  p53) that are connected to many proteins, and that do not interact  specifically with proteins in your network. Conversely, choosing a  higher selectivity will include more proteins that are more specific to  your network, but these interactions will have lower confidence (since  any higher confidence hub proteins will be filtered out). Further, be  aware that changing the selectivity parameter will change the enrichment  results in step 4.5, since different proteins will be included in the  host network.

### 4.4 Merge two host-virus networks

Let us now compare the networks for HPV 16 and HPV 1a. Create a new  host-virus network for “Human papillomavirus type 1a (HPV 1a)” by  repeating steps 4.1 and 4.2. Merge the two networks using **Tools → Merge → Networks**. Move both the HPV 16 and HPV 1a networks into the **Networks to merge** box and otherwise use the defaults for the merge. In the resulting network, use the menu option **Apps → STRING → Set as STRING network**  to manipulate the network as a STRING network again. To show any  interactions between host nodes that were present in one source network  but not the other, first set the confidence to 1, then set the  confidence to the desired confidence (0.4) to retrieve any missing  interactions.

The resulting network can be styled to give the nodes of each species  a distinct color so that the proteins of the two viruses can be  distinguished from each other.

*How many host proteins interact with E6 from both HPV species?*

### 4.5 Functional enrichment

We will now examine the human proteins to see what pathways are enriched in this network.

Next, we will retrieve functional enrichment for the human proteins. Go to the menu **Apps → STRING Enrichment → Retrieve functional enrichment**  and keep the default p-value of 0.05. Homo sapiens will be selected by  default in the species dropdown. It is currently only possible to  retrieve enrichment for host proteins. A new STRING Enrichment tab will  appear in the Table Panel. It contains a table of enriched terms and  corresponding information for each enrichment category. Use the filter  button in the top left of the STRING Enrichment panel to show only **KEGG Pathways**. Click on the draw charts icon to the right of the filter icon to plot the enrichment values on the network.

*Which two KEGG pathways have the lowest p-values? Which host  proteins are associated with the KEGG pathways “cell cycle”? (Hint:  click on the associated row in the enrichment table to select the  proteins with this term.)*





# Exercises for day 3

-------------------------------------------------------------

These exercises will primarily use functions in R but we will also briefly touch upon the `plink`  software package. The exercises are build up over two types of  problems: introductory text that you should run through to make sure you  understand how to use and a few additional questions for you to  explore.

We will be needing the `MESS`, `lme4`, `coxme`, and `glmnet` packages, which can all be installed using the `install.packages()` function.

```
install.packages(c("MESS", "lme4", "coxme", "glmnet"))
```



# Standard GWAS analysis

In this exercise we will try to run a standard genome-wide  association study. As much as possible we will be using standard  functionality in R and try to use the techniques we looked at this  Monday. There are a few specialized packages that wrap some of  fucntionality into complete functions but we will try to do as much of  the analyses by hand to make sure we understand what is going on. Also,  we are considering only a part of a full GWAS dataset since we want to finish the analysis today.

Here we will look at a dataset concerned with the genetic  associations to body-mass index (BMI). The data can be downloaded  directly from the web page using

```R
load(url("http://www.biostatistics.dk/teaching/bioinformatics/data/gwasgt.rda"))
load(url("http://www.biostatistics.dk/teaching/bioinformatics/data/gwaspt.rda"))
```

and it contains information from two chromosomes: the first 11283  # brug denne til at gruppere på chromosomer (columns in the genotype file are from chromosome 1 and the rest for  chromosome 2.

1. How many genes are available in our dataset
   1. 32019 v
      How many individuals are we analyzing? What else is available?
   2. 1324 v 
   3. BMI, gender, age
2. What type of analysis would me do if we had just observed a single gene variant?
   1. lin reg
3. Run a simple analysis corresponding to the model you proposed in 2.
4. Make a full GWAS analysis as shown at the lectures. That means you  should run an analysis for each column in your genotype dataset. [Hint:  you can use the `mfastLmCpp()` function from the `MESS` package to run an analysis between an outcome and each predictor.]
5. Create a Manhattan plot with the resulting *p*

BMI= alpha+ beta1*age + beta2\*gender +beta3*genes + epsilon

~BMI=BMI-\^beta1\*age-\^beta2*gender

Assumes no confounder between age, gender og genotype (hvis genotype begrænser alder, så er der en confounder)



Other assumption: Genotype er linær sammenhæng. Så forskellen mellem 0-1 mutant og 1-2 mutant er ens. Svarer til kodominans. Hvis ikke sandt -> factor(beta3), kategorical instead of numbers. so recessive=no effect 0-1 but effect 1-2
dominant=effect 0-1 but no effect 1-2

-values. Note: on the 

*y*



-axis of a Manhattan plot we typically plot 

−log10(*p*)



 instead of the 

*p*

-value itself.

Try to analyze the data using the lasso instead. Compare the results  to the previous marginal findings. Is something consistant? Is  something disappearing? [Note: we cannot get *p*

1. -values directly out of `glmnet()` so we should instead compare the set of interesting findings.]
2. Try to analyze the data using the elastic net instead with at least  half of the penalty weight put on ridge regression. What should be the  advantage of that in this situation compared to the regular lasso?
3. Compute the genomic control factor and correct the analyses if necessary.

Then do principal components (see monday) and use these for regression. Use lm or glm. Use relevant predictors as binary. Rotate principal components to maximize fit$x



glm(binary(outcome)~fit[,1], familiy=binomial)

Instead of using original predictors which are many, use principal components, to get a maximum of information in the fewest predictors. så outcome=PC1+PC2 f.ex. decide how many based on cumulative proportion. But no biological value in and of itself. You can conclude that 1 PC is mostly driven by Vx, then you can conclude that effect<sub>PC</sub>≈effect<sub>Vx</sub>

I(outcome=='category') -> category=true=1



# Handling unobserved population admixture: 

## PCA in GWAS

A group of researchers are nervous that there is unmeasured  population admixture in the data. This essentially means that *the full  population in reality contains two or more “hidden” sub-populations.* 

1. Why is population admixture a potentially big problem in GWAS studies?
   1. Unknown and different confounders
   2. Obscures Genetic component
   3. DIvergent allele frequencies

We can use principal component analysis based on all SNPs to produce a  combined measure (of the SNPs). *If there is population substructure  then the allele frequencies are likely to be different among different  sub-populations and across many SNPs and this effect is likely to be  captured by the principal component since they are computed regardless  of the outcome.*

Use the data from the previous exercise for this exercise. You may  have to look at a single chromosome only, if your computer is slightly  slow

1. Make a full analysis where you also take potential population  structure into account. How does this change the results? [ Hint:  compute a principal component or two and check if they contribute  anything to the analyses. The **`prcomp()`**  function that we  have used to far is rather slow on large datasets like this. It will  take 5 minutes or so to compute the PCAs using **`prcomp()`**. I can really recommend the function `batchpca()` from the package `onlinePCA`]
2. Which markers should ideally be used for the principal components?  What risks for the study power are there for this approach on the GWAS  findings?



# Heritability (slow)

In this exercise we will try to estimate the proportion of variance explained. To do this we will use the `lmekin()` function from the `coxme` package. Use the large dataset from the exercises above for the following.

1. We should start by computing the allele frequencies for each of the SNPs. Use the `colMeans()` function.
2. The mean of each SNP is 2*p*



 where 

*p*



 is its allele frequency, and it has a standard deviation of 

2*p*(1−*p*)‾‾‾‾‾‾‾‾‾√

1. . Standardize the values for each gene by subtracting the correct mean and dividing by its standard deviation. Call the result `W` in R. See the code below for inspiration.

2. Run the following code

   ```R
   MAF <- 1-colMeans(genotypes)/2  # Minor allele frequency: 1-gns/2. SÅ 1-det halve gennemsnit
   W <- t((t(genotypes)/2 - MAF)/sqrt(2*MAF*(1-MAF)))#transpose(genotype) og / med 2, træk derefter allelfrekvensen fra og / med 2*allelfrekvensen*(1-allelfrekvensen). Transpose resultatet. Så lav rækker til kolonner
   y <- phenotypes$BMI #liste af BMI
   N <- length(y) #antal af BMI værdier
   A <- tcrossprod(W) / ncol(genotypes)#tcrossprod() takes the cross-product of the transpose of a matrix/divideret med antallet af kolloner i genotyper(1)
   library("coxme")
   id <- 1:N # lav en liste af 1-1329 
   res <- lmekin(y ~ 1 + (1|id), varlist=list(id=A))#The lmekin function fits a linear mixed effects model, with random effects specified in the same structure as in the coxme function.
   #Så BMI som funktion a 1+1 eller id. med varianslisten som hedder A
   #varlist= the variance family to be used for each random term. If there are multiple terms it will be a list of variance functions. The default is coxmeFull. Alternatively it can be a list of matrices, in which case the coxmeMlist function is used.
   ```

   The first couple of lines compute the correlation matrix we should use for the calculations. The last line fits the model.

3. Explain what goes on in the code above.

4. Estimate the proportion of variance explained. What is the heritability in this case?



# Find linkage disequilibrium (LD) blocks

The following text explains how to make a heat map.

Heat map plots are used to visualize data from a two-dimensional  array using colors to indicate the values in the array. In R, the  function `heatmap()` plots a heat map, and it requires a numeric matrix as its first argument.

There are several formulas for computing linkage disequilibrium  between two genes. A simple measure that is fast to compute and easy to  interpret is the correlation coefficient. `cor()` computes the correlation matrix for all pairwise correlations.

1. Use the same data as before. Pick 50 SNPs on each side around a  region of interest - for example one of the SNPs found to be associated  with the outcome. The columns in the genotype dataset are ordered  according to their position on the genome. If we use all the SNPs then  the computations will take too long.
2. Compute the correlation matrix for the 100 genes/SNPs.
3. Use `heatmap(x, Rowv=NA, Colv=NA)` (`x` is the matrix of correlation values computed above). What do you see in the graph?
4. What happens if you do not include the arguments `Rowv=NA, Colv=NA`?

The color scheme of the heat map plot can be set with the `col` option. By default it uses heat colors obtained from the `heat.colors` function, but that can be changed to give other color schemes. If the options `ColSideColors` and `RowSideColors`  are provided, they should contain a vector of color names for a  horizontal/vertical side bar that annotates the columns and rows of the  data matrix, respectively. This can be used to indicate known groupings  by color code as shown below.



# Rare variants

We want to see if there is any reason to consider a burden/SKAT test  in order to find association between rare variants and the outcome.

1. Compute the minor allele frequency for each of the genes. What is the distribution of allele frequencies? You can use the `colMeans()` function to compute the mean of the columns in the genotype data.
2. Identify the SNPs found to be in the same clusters using the heatmap  from above and make a simple burden test for each of the clusters  identified. An example of this can be analyzed in R as follows (needs to  be modified):

```
library("lme4")
# Include, say, V6 and V7 as random variable
y <- phenotypes$BMI
fullmodel <- lmer(y ~ 1 + (1|V6) + (1|V7), data=genotypes, REML=FALSE)
smallmodel <- lm(y ~ 1, data=genotypes)
logLik(fullmodel)
logLik(smallmodel)
```

1. How should we compare the results? Has the power improved?



# plink

This exercise is to try the software program `plink`. You need to download and install the programme before you can run the exercise. `plink`  is super versatile and has many more options than what we can consider  here. Instead we will try to make small extensions to the analyses we  made above.

1. Download the `plink` program from

   `http://zzz.bwh.harvard.edu/plink/download.shtml`

   and install it.

The `plink` executable file should be placed in either the  current working directory or somewhere in the command path. All  commands involve typing plink at the command prompt followed by a number  of options (all starting with `--option`) to specify the data files methods to be used.

`plink` needs two standard files: a ped-file and a map-file. The necssary files can be downloaded from the course website.

1. Check the files and make sure you understand their formats. The  map file contains chromosome number, the SNP name, the position in  morgans on the chromosome, and the position in base pairs on the  chromosome.
2. Run the following commands and discuss the output:
   - Type `plink --noweb --file hapmap1` and look at the  output. It takes two files (a ped file and a map file) as input. They  can be downloaded from the course website.
   - Type `plink --noweb --missing --file hapmap1` and look at the output. What information is contained in the missing files.
   - Type `plink --noweb --freq --file hapmap1` and look at  the output. What information is contained in the frequency files. How  does the results compare to the results you computed manually earlier in  R?
   - Type `plink --noweb --assoc --file hapmap1` and look at the output. What information is contained in the association files.
   - Type `plink --noweb --assoc --adjust  --file hapmap1` and look at the output. What information is contained in the adjusted association files.

# Torsdag

## Exercise 2: Building a suffix array

  Suffix arrays are data structures used in sequence alignment and read mapping, allowing for quick localization of exact sequence matches.  

1.  Given the reference sequence "BANANA$", construct the suffix array        for it. (Hint: 1. Create a list of all suffixes; 2. Sort them        alphabetically; 3. Write down the list of starting indices in order.)   
2.  Consider the read "NA", and find all instances of it in your        reference using the constructed index structure from (1).        (Hint: Use binary search.) How many steps do you need to find it?   
3.  Consider if you have to find it without using a suffix array.        How many steps do you need in this naive approach?   

0. BANANA$
1. ANANA$
2. NANA$
3. ANA$
4. NA$
5. A$
6. $

6
5
3
1
0
4
2



## Exercise 3: Transcriptome analysis with Galaxy

  Galaxy contains tools for mapping and analyzing RNA-seq data. This is an opportunity to get some experience in handling these data, including quality control and mapping to the reference genome. Go to the 

Galaxy workbench

 and have a short look into the tool list on the left side to see which tools are available. 

  This exercise is based on data from the 

Human Body Map project 2

, and is loosely based on this exercise 

from usegalaxy.org

. In the exercise you will 

(I)

 clean, 

(II)

 map and 

(III)

 find differentially expressed genes in transciptomic RNA-seq data.  In the exercise we will work with our 

local Galaxy server

.

  First, register and login, so that you can re-access your Galaxy history at any time. You will need this for the next exercise.

  Second, download and save the four raw sequencing data sets (subsets of the full data sets to save time and disk space):
 [adrenal_1.fastqsanger](https://rth.dk/~hull/teaching/bioinformatics_1/transcriptome_analysis_with_Galaxy/Galaxy2-adrenal_1.fastq.fastqsanger)
 [adrenal_2.fastqsanger](https://rth.dk/~hull/teaching/bioinformatics_1/transcriptome_analysis_with_Galaxy/Galaxy3-adrenal_2.fastq.fastqsanger)
 [brain_1.fastqsanger](https://rth.dk/~hull/teaching/bioinformatics_1/transcriptome_analysis_with_Galaxy/Galaxy4-brain_1.fastq.fastqsanger)
 [brain_2.fastqsanger](https://rth.dk/~hull/teaching/bioinformatics_1/transcriptome_analysis_with_Galaxy/Galaxy5-brain_2.fastq.fastqsanger)

 Remember that to download a file you just have to right click on it and select save destination as. Upload the data sets by clicking on `Get data` → `Upload File from your computer` in the Galaxy Tools list on the left side. Click `Choose local file` and select the four files. Then set the data type to `fastqsanger` (not `fastqcsanger`) in the dropdown menu below the text field, set the genome to `hg38`, click `Start` and `Close`. The dataset should now show up on the right side in your `History`. The background color of each dataset or task in the `History` describes the status of the job: *gray* means job is pending; *yellow* means job is running, and *green*  means job is done. Some tasks may run for a while but you do not have  to wait until its finished. Instead you can continue to prepare the next  steps and queue up several jobs. Galaxy will automatically start your  jobs in the history as soon as all its input data it available. Now you  are ready to start your analysis.

  

### I. Cleaning

1. Click on the `eye icon` to the right of the RNA-seq data in your history to see the content of the data files. 

   **(a)** Do we work with single-end or paired-end reads? 

   ​	

      **(b)** What is the length of the reads?   

2. 

3.  The first step in the bioinformatics analysis of RNA-seq data is its quality control. Run Fast Quality Check (`Local` → `FastQC`). Under `Short read data from your current history` select `multiple datasets`. Select all four RNA-seq datasets and press `Execute`.  In your history an entry for FASTQC raw data and FASTQC Webpage should  appear for each data set (total 8 entries). When FASTQC is done (green  color) check the FASTQC webpage which will show you different quality  criteria. 

      **(a)** How many sequences are there in each of the data sets and how long are they? 

   Brain 1: Total Sequences	37992, Sequence length	50, PQ: 0, %GC	54

   Brain 2: Total Sequences	37992, Sequence length	50, PQ: 0, %GC	55

   Adr 1: Total Sequences	50121, Sequence length	50, PQ: 0, %GC	58

   Adr 2: Total Sequences	50121, Sequence length	50, PQ: 0, %GC	59

   **(b)** How many sequences are flagged as poor quality and how is the sequence quality distributed over the sequence length? 
      **(c)** How large is the average GC content and the GC content distribution? 
      **(d)** Does the RNA-seq data need cleaning?   

4. 

5.  In RNA-seq experiments the sequence quality is expected to  decrease with the length. But there can be other problems with the read  qualities. An easy way of improving the quality of your RNA-seq data is  cutting the end of the reads. A base quality of 20 is usually considered  good enough. 

      **(a)** Is it enough to clean the reads by cutting the reads after certain length? 

      Use the tool `Local` → `Trimmomatic`. This tool has  several ways of cleaning up the data. You can build a cleaning pipeline  by adding several steps of Trimmomatic operations. To add a step click  on `+ Insert Trimmomatic Operation`. The operation is set by clicking on `Select Trimmomatic operation to perform`.  Make a pipeline which Step 1. Cuts bases at the beginning of a read if  the quality of the base is below 20. Step 2. Performs a sliding window  trimming requiring an average score of 20 over 4 bases. Step 3. Drop  reads with a length below 20 nucleotides. Set the `Single-end or paired-end reads` to `Paired-end (two separate input files)`.  Run the pipeline for two sets of two data sets. This will create four  new files for each pair of data sets. Two files contains the reads where  there is no partner after cleaning (named R1 unpaired and R2 unpaired),  and two files contains the reads which has a partner, named R1 paired  and R2 paired. In the rest of the exercise we will only use the paired  reads. Run FASTQC on the paired files. 

   ​    **(b)** Compare with the FASTQC results from the untrimmed data.  Answer the same questions as in 2.a to c. Which differences do you  observe? Did the quality of the data improved?   
   



### II. Mapping

1. Now that the sequences are cleaned you are ready to map them to  the reference genome. Choose Bowtie2 for the mapping because its pretty  fast (`Local` → `Bowtie2`). The different mapping programs  can perform very different on the same data set with different  parameter settings. Select the analysis mode `2: Full parameter list` and look into the different parameters. 

      **(a)** Enumerate some parameter settings that impact running time and mapping quality and explain why? 

      For now keep the default settings. Choose the correct library type and  reference genome and run Bowtie2 for both adrenal gland and brain by  clicking `Execute`!   

2. 

3. After the sequence mapping has been finished, check the output of the Bowtie2 job inside your `History`. lick on `Bowtie2 on data ..`  step. Thus will show you a box with information about the run. However  in some cases there may be a few error messages drowning the relevant  output. To see all the output click on the `i` icon and find the field `Tool Standard Error`. Click on the `stderr` link. 

   **(a)** How many reads have been mapped exactly one time and how many multiple times? How many reads could not be mapped? 

   "40692 reads; of these: 40692 (100.00%) were paired; of these: 13040 (32.05%) aligned concordantly 0 times 26860 (66.01%) aligned concordantly exactly 1 time 792 (1.95%) aligned concordantly >1 times"

      **(b)** Any idea why some sequences have not been aligned to the reference genome?    

4. 

5. The output of Bowtie2 is a BAM file. 

      **(a)**  Shortly describe what a BAM file is. 

      The content of the BAM file can be visualized by the UCSC genome  browser. Find the link to the UCSC genome browser in the bottom of the  expanded job description in the `History`. In the UCSC genome  browser the Bowtie2 run should appear as one of the Custom tracks (in  the top of all tracks). Navigate in the browser to see the read profiles  `chr19:2,950000-3,550000` (The sequences of the RNA-seq experiment are not covering the entire genome.). 

      **(b)** How could you find this region if you had to do it yourself? 
      **(c)** Where do the sequences map and how are they distributed?   

6. 

7.  Another mapping tool supported by Galaxy is TopHat. 

      **(a)** What is TopHat doing and how is it related to Bowtie2? (*Hint:* Check the tool description in Galaxy.) 

      TopHat has a longer run time than Bowtie2, hence, you might have to  wait for the results, even beyond the class. To avoid the waiting time  pre-aligned files (using Galaxy) are available (Tophat accepted_hits  files, aligned using default parameters): `accepted_hits`. These files are available here: [Adrenal](https://rth.dk/internal/index.php/s/rBD9apDiYd15ogr/download?path=%2F&files=Tophat_adrenal_paired_R1_R2_accepted_hits.bam) and [Brain](https://rth.dk/internal/index.php/s/rBD9apDiYd15ogr/download?path=%2F&files=Tophat_brain_paired_R1_R2_accepted_hits.bam). Download the two files and upload them to Galaxy using `Get data` → `Upload File`. Make sure to set the genome to `hg38` and the data type to `bam` before pressing `Start` to start the upload. These files contains the paired mapped reads. 

      **(b)** Why is it a good idea to use Tophat? Which mapper is most suitable for this kind of data and why? 
      **(c)** View the Tophat predictions in the genome browser and discuss what you see.   
      

### III. Differential expression

  In this part of the exercise we will search for transcripts that are  differentially expressed between the adrenal gland and the brain. To do  this you will the cufflinks pipeline.  

1. The first step after mapping the reads is to assemble them into  transcripts. This can be done with the Cufflinks. Run cufflinks on the  uploaded Tophat files. 

      **(a)** Why is it necessary to assemble the reads? 
      **(b)** Which output files do you get? What information does these files contain? 
      **(c)** What does FPKM mean? 
      **(d)** Try viewing the `assembled transcripts` in the  genome browser. Do the results look as you would expect? Also compare  the Tophat and Cufflinks predictions with the gencode annotation.   

2. 

3.  The second step is merge the assembled transcripts, to get a meta  assembly. This uses the information from the two conditions to build a  better transcript assembly. Run `cuffmerge` using the two `assembled transcripts` data sets. Once the job finish you can view the result in the UCSC genome browser.   

4. 

5.  The last step is to find out if there are any differentially expressed genes between the two tissues. This can be done using `Cuffdiff`. Set the `transcripts` parameter to the output from the `cuffmerge` step. Set the brain dataset as one condition, and the adrenal dataset as the second condition. 

      **(a)** Are there any differentially expressed genes? 
      **(b)** What does the log2 (fold_change) tell you? 

      Try viewing some of the differentially expressed transcripts in the genome browser.  





# Fredag

# Exercises for day 5

These exercises will primarily use functions in R but we will also briefly touch upon the `freebayes`  software package. The exercises are build up over two types of  problems: introductory text that you should run through to make sure you  understand how to use and a few additional questions for you to  explore.

We will be needing the `MESS`, `energy`, `PTAk`, `pscl` and `compositions` packages, which can all be installed using the `install.packages()` function.

```
install.packages(c("MESS", "energy", "compositions", "PTAk", "pscl"))
```



# Distance correlations

Here we will look at a dataset concerned with the genetic  associations to body-mass index (BMI) that we looked at earlier. The  data can be downloaded directly from the web page using

```
load(url("http://www.biostatistics.dk/teaching/bioinformatics/data/gwasgt.rda"))
load(url("http://www.biostatistics.dk/teaching/bioinformatics/data/gwaspt.rda"))
```

We are interested in looking for *other* relationships that we did not find previously.

1. Use the same data as before. Pick 10 SNPs from columns 791 to 810  which are located around one of the SNPs of interest. If we use all the  SNPs then the computations will take too long.

2. Compute the Pearson correlation matrix for the 20 genes/SNPs selected above. You can use the `cor()` function for this.

3. Load the `energy` package to get access to the `dcor()`  function. Compute the distance correlation for each of the SNPs against  SNP 802. Does the results differ from the findings from `cor()`? [ Hint: you can use something like

   ```
   sapply(791:810, function(i) { dcor(genotypes[,803], genotypes[,i])})
   ```

   to do the computations. Make sure you understand what goes on here. ]

   ​	For 791:810, brug funktionen dcor med række 803 som x og genotypes (1-i) som y.

4. Use `dcor()` to compare the SNPs to the BMI outcome  for the SNPs that were found on Wednesday to be related to the outcome.  Do you find any apparent relationships between the SNPs and the  phenotypes, that you did not find using the traditional regression  model?

5. Compute t tests for each of the distance correlations. What are your findings? [Hint: use the `dcor.ttest()` function. ]

6. Why does it make sense to compare the outcome from the regression model to results from `dcor()`?



# Zero-inflated and hurdle models

Zero-inflated and hurdle regression models with Poisson and negative-binomial models can be modeled in R using the `pscl` package. There are two primary functions we will be using: `hurdle()` and `zeroinfl()`. Both functions work similary to the `glm()` function that we have used previously.

1. Load the data from webpage by running the following command

   ```
   load(url("http://www.biostatistics.dk/microbiome.rda"))
   ```

The data consists of two conditions (so have to add up to 100), A and B, the relative abundance   from 15 taxa and 35 samples from each condition. There is also  information on age an sex of the individuals the provided the samples.

1. Fit a zero-inflated Poisson model for the raw relative abundance of taxa 1 (otu 1) to compare the two conditions using the `zeroinfl()` function. Hint: you can use the `subset` argument to restrict the analysis to only one taxa.

   Also, the `zeroinfl()` function uses the Poisson model (for the non 0) by  default so it is necessary to round the relative abundance to the  closest integer using, say, the `round()` function. value=normalised

   What is the conclusion? What are the directions of the effect of the conditions? Is this what you would expect?

   B er negativ, ses ikke på rå-tælning

2. Fit the same model with argument `dist="negbin"` (still for the raw relative abundance of taxa 1). Which of the two models fit best?

3. Now try the same two models but include additional information on  the gender and age. You may run into problems with the model fit /  convergence. You can specify individual models for the proper  distribution and the logistic regression model by specifying individual  models. This can be done in the model formula by giving two models  separate split with a vertical bar such as `y~x | x+z`. (så enten model eller model x på anden model)

4. Try fitting the same models with a hurdle model

2 models, summary -> start with the bottom part =binomial = 0 part of model, upper part - poisson = non-0 part. B has a lower risk of being 0, no sign. log ≈ ±20 -> somethings wrong.

two models kan have different models, defined with binom model | poisson model (check altid hvilken der har hvilken)

no b zero, 8 for a… messes up 0-part. Poisson is log'ed, so result needs to be exp'ed, så B 25% mere

log(inter)=niveau af a, log(cond)=prop forskel mellem a og b. log(inter)*log(cond)=niveau af b.

# Compositional data analysis

Use the `microbiome` data from exercise 2 in this  exercise. The data consists of two conditions, A and B, the relative  abundance from 15 taxa and 35 samples from each condition. There is also  information on age an sex of the individuals the provided the samples.  We start by extracting the intensities as a matrix with 15 columns and  70 rows.

We start be extracting the values as a matrix with 70 rows and 15 columns

**check link type**

```R
exprmatrix <- matrix(microbiome$value, nrow=70)
```

1. Load the `compositions` package and compute the  isometric-log-ratio for each row of the expression matrix. Store the  result in some object for later use. [ Hint: You can use the `ilt()` function.]

2. What is the dimension of the resulting matrix? What is returned?

   1. 70x15= samme som oprindelig data, den som er sat til 1 er et random midtpunkt

3. We can now analyse the data in different ways: 

   - We could use each of the columns in the transformed matrix  independently of each other one by one and use those as outcomes. The  predictors could be condition, age, and gender.
   - We could use the full set of transformed values to predict the  condition based on the microbiome composition. This could be a regular  logistic regression model or a penalized regression model.

    Do one or both of these analyses.



# PARAFAC

For this exercise we will be using the `PTAk` package and the array `x` that is grabbed from an R dataset on the web. We wish to decompose 5 small metabolite arrays that are found in `x` and possibly get information on the basis functions.

```
library("PTAk")
load(url("http://www.biostatistics.dk/parafac.rda")) # Loads x
```

The individual metabolite scans can be accessed through the last index of the array. For example will `x[,,1]` give the 120×120

 matrix of values for individual 1. To get an idea of the data try plotting a 3D plot of the output from person 1:

```
persp(x[,,1])
```

If you have the `rgl` package installed you can create a 3D plot that can be rotated:

```
library("rgl")
persp3d(x[,,1], col="lightblue")
```

There are two peaks on the array. The same is true for the other  individuals although the relative height of the peaks vary from person  to person.

Let us create the matrices used to create the tensor product:

```
res <- CANDPARA(x)
```

1. Look at the output from the analysis. What does `res` contain?
2. Try to plot the “basis” functions that make up the different  components. For example, the first component of the first matrix can be  extracted with `res[[1]]$v[1,]`. Plot the three components from the first matrix and interpret them.

We can create our matrix approximation with the `REBUILD()` function.

```
appx <- REBUILD(res)
```

1. Compare the approximated tensor with the original tensor. How well  would you say the matrix was approximated (you should just be  eye-balling here).
2. By default `CANDPARA` uses 3 components but this can be changed with the `dim` argument. Try to approximate the tensor with just 2 components. How much loss is there in the tensor approximation now?



# Gene calling

This exercise is outside the scope of what we can manage within the  current version of the course, so it is mostly included here as a super  brief introduction to gene variant calling using the `freebayes` software. `freebayes`  is not distributed as a binary package so you need to download and  compile it for it to work. If you are running a windows computer then  this might give you more of a headache than you’d like. Alternatively,  you can use a docker container, with the relevant image.

To install the program you should first get a copy of the source

```
git clone --recursive git://github.com/ekg/freebayes.git
```

then you need to compile it. FreeBayes requires g++ and the standard C  and C++ development libraries. Additionally, cmake is required for  building the BamTools API.

```
cd freebayes ; make
```

In its simplest operation, freebayes requires only two inputs: a  FASTA reference sequence, and a BAM-format alignment file sorted by  reference position. See the web-site for an example.