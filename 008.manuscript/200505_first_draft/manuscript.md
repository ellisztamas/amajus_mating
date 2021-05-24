
"How's your father? Male fitness due to differences in flower colour in a snapdragon hybrid zone"
=================================================================================================



Thomas James Ellis,^1,2^, David Luke Field,^2,3^ and Nicholas H. Barton,^2^

**1** Gregor Mendel Institute of Molecular Plant Sciences; **2** IST Austria; **3** University of Perth

# 1. Introduction {#MPSection:BD1E0B64-DC41-42C4-86C6-C51FD1A35516}

Mating patterns are important because they determine which alleles unite in the same individuals to be tested by selection. When individuals with particular combinations of genotypes are more or less likely to mate with each other than expected by chance, this can bring together combinations of beneficial alleles in the offspring. It can also bring recessive deleterious alleles together and cause inbreeding depression. There may also be variation between individuals in how many mates, and hence offspring, they have, causing heritable changes in the population. We would like to understand who mates with whom, and characterise variation in reproductive success.\[Efficient inference of paternity and sibship inference given known maternity via hierarchical clustering\](\#Ellis2018);;

For many plants mating is mediated by animal pollinators. Two aspects of pollinator behaviour are important. Firstly, pollinators make decisions about which plants to visit based on the phenotypes of the plants. They may preferentially visit flowers that are especially attractive or rewarding, causing differences in reproductive success between plants. There is an extensive literature exploring differences in the female component of reproductive success measured through seed set. However, although the majority of flowering plants are hermaphrodites \[ref?\], the male component of fitness, reflected in successful pollen export to other plants, has received much less attention, because this requires accurate assignment of paternity in seeds. Since every offspring receives exact half its genetic material from the mother and half from the father, by neglecting variation in male fitness we also potentially miss substantial variation in overall reproductive success.

Secondly, it is more energy efficient for pollinators to transition between plants that are close to one another than to move at random between plants. This means that plants are more likely to mate with close neighbours than more distant plants. The distribution of pollen dispersal distances is described by the pollen dispersal kernel. In plants, it is common to find that the dispersal kernels that show strong leptokurtosis, or overdispersion; most dispersal is between close neighbours, but there is a 'fat tail' of long-distance dispersal events that cause the dispersal kernel to decay more slowly than exponentially. Although these events are relatively rare, they contribute disproportionally to the spread of migrants, allowing populations to expand rapidly in space, and to homogenise distant populations. This is especially important in a hybrid zone, where divergent populations meet and exchange migrants. Here, the extent to which alleles at individual loci are able to introgress across the hybrid zone is determined by the balance between the homogenising effect of dispersal, and selection against migrant alleles keeping populations spatially distinct. Given estimates of dispersal distance and distribution of alleles in space we can infer the strength of selection acting on those alleles. Accurate estimates of the shape of the dispersal kernel are this invaluable in understanding the evolutionary processes acting in natural populations.

Pedigree reconstruction has emerged as a powerful tool to investigate both both differences in mating success and dispersal ;;;;\[Wild pedigrees: the way forward\](\#Pemberton2008);;\[Parentage and sibship exclusions: higher statistical power with more family members\](\#Wang2007)\[Towards unbiased parentage assignment: combining genetic, behavioural and spatial data in a Bayesian framework\](\#Hadfield2005). The idea is to reconstruct mating events by identifying the parents of a set of offspring asssuming Mendelian inheritance based on shared alleles at loci for which offspring and candidate parents have been genotyped. This is most successful when as much informative data as possible can be included in the analysis ;;\[A Bayesian Framework for Parentage Analysis: The Value of Genetic and Other Biological Data\](\#Neff2001). One way to achieve this is to jointly infer sibling relationships and the paternity of those sibships . Likewise, inference of parentage jointly with other biological parameters that influence mating, such as those related to dispersal or mating success increases the accuracy to infer both the pedigree and biological parameters ;. It follows from this that an even more powerful approach would be to unify these approaches and jointly infer parentage, sibships and biological parameters. However, we currently lack software which is able to do this efficiently for modern SNP data.

Here, we examine a hybrid zone between two subspecies of the snapdragon, ;*Antirrhinum majus* that differ for flower colour. Flower colour in *A. majus* is controlled primarily by the genes *Rosea*, controlling the distribution of magenta anthocyanin pigments, and *Sulfurea*, controlling the distribution of yellow aurone pigments. In the Vall de Ribès in the Spanish Pyrenees, the yellow-flowered *A. m. striatum* and the magenta-flowered *A. m. pseudomajus* meet and form a narrow hybrid zone where loci at *Rosea* and *Sulfurea* recombine to give rise to six flower colour morphs. Alleles at these loci nevertheless form sharp clines consistent with selection against recombinant genotypes over many generations. Furthermore, parentage analysis indicates genotype-by-environment interactions for fitness, with the *A. m. striatum*- and *A. m. pseudomajus*\[Patterns of floral colour neighbourhood and their effects on female reproductive success in an Antirrhinum hybrid zone\](\#TASTARD2011)-like genotypes having an advantage in their home territory. As there are no obvious difference in habitat or other phenotypes across the hybrid zone, and given that pollinators are known to discriminate between flowers based on colour, non-random foraging behaviour by pollinators is an appealing explanation for these differences in fitness. Statement about female fitness (in the core only?) from ;, but we still don't know much about male fitness. Parentage estimates reflect the whole life cycle, which includes variation in mating success, seedling establish and survival, and we are not able to distinguish pollen and seed dispersal. As such, this is not ideal to investigate mating patterns directly to properly characterise male fitness and pollen dispersal.

In this study we reconstruct mating events between wild snapdragons in the hybrid zone. We use a panel of open-pollinated seedlings of known maternity collected in a natural hybrid zone and dense sampling of possible sires. This differs from the pedigree because there is no seed dispersal or establishment, so parentage reflects mating patterns through pollen movement only. We build on previous methods to jointly infer paternity, sibship relationships and the pollen dispersal kernel. Use these results to examine the extent to which pollinators cause assortative mating and differences in male fitness beyond what would be expected due to local population structure. Finally we use simulations to verify that these results are unlikely to be due to false positive paternity assignment due to insufficient statistical power or missing fathers. These results imply that pollinators are an important force mediating selection on flower colour.

# 2. Materials and Methods {#MPSection:9C0EFB03-BC84-4C7F-9624-F43C6FF827FE}

## Study population {#MPSection:62D1C788-9DB5-4F35-9595-6BF0BC122890}

<div id="MPFigureElement:D7553A62-2F56-4DE4-9A24-69630214E618">

![](Figure_1.png)

Figure 1: Map of the hybrid zone. The map shows the distribution of parental- and hybrid-like plants along the lower (South) and upper (North) roads. Maternal plants were selected from the lower road between -628 and 890 metres.

</div>

The study population grows along two parallel roads running East-West close to Ribès in the Spanish Pyrenees (<span class="kind elementIndex">**Figure 1** </span>). The 'lower' road is at 1150-1200m, whilst the 'upper' road climbs 1250-1500m, and is 500-1000m north of the lower road. Hybrids are mostly confined to a \~1km 'core' hybrid zone, with *A. m. striatum*- and *A. m. pseudomajus*-like plants becoming dominant to the West and East respectively (<span class="kind elementIndex">**Figure 1** </span>). We surveyed as many flowering plants as we could find in June and July of 2012 (n=2128), and collected information on flower number and location using a Trimble GeoXT datalogger. Antirrhinum grows in disturbed habitat such as roadsides and railways and much of the habitat between the two roads is forest and pasture, so it is likely that we sampled most of the plants that flowered. *A. majus* has a sporophytic self-incompatibility system, and self-pollinated seeds are very rare.

<div id="MPFigureElement:6AA05495-92BA-4E8A-A43B-22B3FBF46705">

![](Figure_2.png)

Figure 2: Flower colour phenotypes. Combinations of alleles at the Rosea and Sulfurea loci controlling anthocyanin and aurone pigmentation give rise to six flower phenotypes. A. m. pseudomajus and A. m. striatum show the phenotypes displayed top left and bottom right respectively.

</div>

We collected two to three leaves for DNA extraction and dried these in silica gel (Fischer Scientific) and one mature healthy flower for phenotyping. Following the visual scoring system developed by \[Evolutionary Paths Underlying Flower Color Variation in Antirrhinum\](\#Whibley2006);; we assigned plants to one of six flower phenotypes based on their inferred genotypes at the loci *Rosea* and *Sulfurea* (<span class="kind elementIndex">**Figure 2** </span>). ROS/ROS plants have intense magenta pigmentation, whilst ros/ros plants have none, and heterozygotes are intermediate. The SULF allele is fully dominant over sulf, and SULF/SULF and SULF/sulf plants cannot be easily distinguished by human eyes, and will henceforth be referred to collectively as SULF/+. sulf/sulf plant are yellow throughout the face of the flower, whilst SULF/+ plants express aurones only in a small patch at the centre of the flower.

In August 2012 we collected a single, mature, wild-pollinated fruit from each of 60 mothers (figure \[fig:map-2012\]). In order to minimise disturbance to the population we only sampled from plants which had set a minimum of five mature fruits. These mothers were chosen to represent an even sample of pigmentation genotypes, spread as evenly as possible across the core of the hybrid zone where hybrids are most dense.

## Genotyping {#MPSection:F4F6DA47-3093-47CE-BF77-1D9E53C39F41}

We grew seeds in 5cm plug trays filled with potting compost (Gramaflor) in a greenhouse under Sylvania GroLux lights on a 16-hour cycle. We sowed three seeds per plug for 50-70 plugs per maternal family and thinned seedlings to a single seedling per plug after cotelydons had appeared. We transferred approximately 1cm\^2\^ of fresh tissue from 1419 seedlingsto 96-well DNA-extraction plates (LGC Genomics, Berlin) and allowed tissue to dry using the sample bag and silica gel provided. For parental tissue from the hybrid zone we transferred approximately 1cm\^{2} tissue dried in the field to the same plates. DNA extractions of the plated tissue samples were carried out by LGC Genomics \[Is this in another publication yet I can cite?\].

We genotyped tissue samples at 71 SNPs by KASPR sequencing (LGC Genomics). These SNPs are a subsample of a panel used for a wider survey of the hybrid zone (David Field, unpublished data). The total SNP panel is a mixture of \[how many?\] diagnostic (showing a gradient in allele frequency) and (how many?) parentage (with as even a gradient in allele frequency as possible) SNPs. For parentage loci we chose only biallelic loci with a minor allele frequency greater than 0.3 in each of inner four pools closest to the centre of the cline, selected to maximise mapping distance between loci. Diagnostic SNPs were either linked to pigmentation loci, or else showed sharp clines across the hybrid zone. To select a subsample of SNPs for this study, we selected markers which were at least 2cM apart to maximise mapping distance between individuals. We removed 474 offspring and four adults that had missing data at more than 7.5% of the SNPs. We also pruned 7 SNPs that showed more than 10% missing data, or less than 15% heterozygosity. This left us with a set of 984 offspring from 60 maternal families, with between two and 29 offspring per family (mean=16.4).

## Pollen dispersal kernel {#MPSection:A4553086-BE4E-4AC6-983B-5B0136605857}

A useful function for describing plant dispersal kernels is the generalised normal distribution (GND). This is a generalisation of the exponential family of probability distributions and includes the exponential and standard normal distributions as special cases, but allows for fat and thin tails \[Why Trees Migrate So Fast: Confronting Theory with Dispersal Biology and the Paleorecord\](\#Clark1998);\[A generalized normal distribution\](\#Nadarajah2004);\[Using genetic markers to estimate the pollen dispersal curve\](\#Austerlitz2003);;;;;;. The GND describes the probability of observing dispersal distance *d* through scale parameter *a* and shape parameter *b*:

<div id="MPEquationElement:D5113F89-58C2-45D6-9845-2A076BDC6CBE">

GND(d\_{mj},a,b) = K \\exp \\left(\\frac{-d}{a} \\right)\^b

</div>

where K = b/({2a\\Gamma(\^1/\_b)}) and is a normalising constant and \\Gamma is Euler's gamma function. The function takes the form of the standard exponential and normal distributions when *b*=1 and *b*=2 respectively. Values of *b*&lt;1 indicate leptokurtic distributions with tails that decay more slowly than would be expected under the exponential distribution. The mean dispersal distance is given by the standard deviation of the GND:

<div id="MPEquationElement:15FAE91C-B1F3-49E4-A5A4-3168CA8B8AA6">

\\sqrt\\frac{a\^2 \\Gamma(\^3/\_b)}{\\Gamma(\^1/\_b)}

</div>

We model the pollen dispersal kernel based on the distances between known mothers of full-sibling families and the plants inferred to be the fathers of those families. Leptokurtosis in this distribution may be caused by a long tail of genuine mating events with distant fathers. However, when unrelated individuals are incorrectly identifed as being fathers of a full sibship, we expect these individuals to be drawn at random from the population since the marker panel is designed to show no spatial structure,, whereas true fathers will (on average) be close to maternal plants. This will inflate the apparent kurtosis in the data and underestimate *b*. To accommodate this we model the probability that candidate pollen donor individual *j* is the sire of a full sibship with mother *m* given distance d\_{mj} between them as a mixture between a GND and a uniform distribution:

<div id="MPEquationElement:8AEE4BB4-7F7D-4F88-AB68-6D9E9F3178A9">

\\Pr(d\_{mj}|a,b,\\lambda) = \\lambda GND(d\_{mj},a,b) + \\frac{1-\\lambda}{N}

</div>

where *N* is the number of candidate fathers and parameter 0 \\leq \\lambda \\leq 1 determines the proportion of the probability mass due to 'real' dispersal. Conveniently, 1-\\lambda also provides an estimate of the false positive rate if we had run paternity analysis using genetic data alone.

## Joint estimates of paternity, sibships and dispersal {#MPSection:29656B57-0A74-4750-BBCE-7D43BBFEEA7A}

### Allowing for covariates in paternity inference {#MPSection:53CF8519-341F-4860-90D7-981DA062C081}

We have previously described a Python package ;;*FAPS*\[Efficient inference of paternity and sibship inference given known maternity via hierarchical clustering\](\#Ellis2018) which performs joint analysis of paternity and sibship relationships for sibling arrays based on SNP data, and allows for integrating out uncertainty in relationships ;\[Towards unbiased parentage assignment: combining genetic, behavioural and spatial data in a Bayesian framework\](\#Hadfield2005). Here we extended the software to allow for additional non-genetic information to be included ;.

*FAPS* is based on a matrix **G**, with a row of each offspring and a column for each candidate father, where element g\_{ij} is the probability that candidate *j* is the father of offspring *i*. Rows in **G** sum to one, and describe a multinomial distribution of probabilities of paternity over all candidate fathers. The final column of **G** is the probability that the father of each offspring was missing from the sample of candidates, based on the probability of observing offspring alleles from population allele frequencies, and an estimate of the proportion of possible pollen donors in the population that had been sampled. *FAPS* then uses hierarchical clustering to partition offspring into a sets of plausible full sibships. The probability that *j* is the father of putative full sibship *k* is then \\prod\_i g\_{ij}.

Similarly, we can describe the probability z\_{j} that *j* is the father of sibship *k* based on non-genetic information through a suitable function of those data as z\_j\\prod\_i g\_{ij}. Here, we model z\_j as Pr(d\_{mj}|a,b,\\lambda), but this formulation is general to any suitable function relating non-genetic data to probabilities of paternity. The likelihood of the whole array is then the product of probabilities for each full sibship within a partition, summed over all possible partition structures. When there are multiple maternal families, the likelihood of the whole dataset is the product of per-family likelihoods over each maternal family.

### Inference via MCMC {#MPSection:C463125F-8581-4DB7-8976-D512E5286B23}

We can then infer paremeters of interest by finding combinations of values that maximise the likelihood of the data. Here, there are four parameters of interest, while using the algorithm in *FAPS* to integrate out uncertainty in paternity and sibship relationships: (1) the proportion of missing fathers, (2) the shape and (3) scale parameters of the generalised normal distribution describing dispersal, and (4) the mixture parameter \\lambda describing the strength of signal from the GND distribution of dispersal.

We used the Metropolis-Hastings Markov-chain Monte Carlo algorithm to infer the posterior distribution of these parameters. We ran four independent chains beginning from distant areas of the parameter space. At each iteration, each parameter was peturbed by a a factor drawn from a normal distritution with a fixed standard deviation for each parameter (\\sigma = 0.025 for proportion of missing fathers and \\lambda, 0.05 for shape and 2 for scale). We ran each chain for 60,000 iterations to ensure mixing, and subsequently removed the first 1000 iterations of each as burn-in. After checking chain convergence (supplementary figures 2-4) we thinned subsequent iterations to retain 250 posterior draws from each chain for further analyses.

We used beta prior distributions for proportion of missing fathers and \\lambda (Beta(\\alpha=3, \\beta=15) and Beta(\\alpha=1.1, \\beta=1.1) respectively), and Gamma priors for dispersal shape and mean dispersal distance (\\Gamma( shape=200, scale=200) and \\Gamma(shape=2, scale=200) respectively). We used a prior for mean dispersal distance instead of the scale parameter of the GND because this has a more intuitive interpretation.

## Flower colour and male fitness {#MPSection:F32655F7-E452-44D2-B0AF-D9BE5A5CCC43}

We next investigated whether certain flower-colour genotypes sire more offspring than would be expected given local genotype frequencies, and whether this varies across the hybrid zone. We did this by comparing 'observed' mating events, based on genetic and dispersal data, to 'expected' mating events simulated based on dispersal only for mothers in each of five bins of approximately 300m (n=3, 6, 21, 27, 3 mothers per bin; supplementary figure 1). We calculated observed and expect mating probabilities for each of 1000 sets of dispersal parameters and proportions of missing fathers from the trimmed MCMC results.

To calculate observed mating probabilities we used the likelihoods from ;*FAPS*\[Efficient inference of paternity and sibship inference given known maternity via hierarchical clustering\](\#Ellis2018) that an individual candidate father had sired at least one offspring with a maternal plant based on genotype data, his distance from the mother, and dispersal parameters, integrating over uncertainty in sibling relationships ;. For expected mating probabilities under random mating we calculated the probability of mating between each pair of mothers and candidate fathers based on the distance between them and dispersal parameters. We then summed these values over all candidates of each *Rosea* and *Sulfurea* genotype, and over mothers within each spatial bin. We normalised likelihoods for each maternal genotype to sum to one to give relative probabilities than mothers in each bin received pollen from males of each flower colour genotype.

## Power analyses {#MPSection:13E33082-AFF8-41AF-888C-849D20292F00}

One explanation for apparent leptokurtosis in dispersal is incorrect assignment of paternity to candidates far from maternal plants, either because the true father was not sampled or due to insufficient statistical power of the marker set. This would inflate the tail of the dispersal kernel even if true dispersal were not leptokurtic (i.e. cause the shape parameter of the GND to be less than one, even when it is really greater or equal to one).

To test this, we first examined how often we should expect incorrect paternity assignment using this marker set, and the extent to which incorrect assignment and missing fathers would inflate apparent dispersal. We used FAPS to simulate offspring based on observed maternal and candidate male genotypes \[Efficient inference of paternity and sibship inference given known maternity via hierarchical clustering\](\#Ellis2018);;. For each of the 60 mothers we sampled mates in proportion to their distance, assuming that the pollen dispersal distances are exponentially distributed with a mean of 10, 50, 100 and 200m. Since clustering offspring into sibships is the most computationally demanding part of the analysis, we instead consider paternity of individuals rather than sibships; this is valid because he are primarily interested in mating events, not of sibling relationships themselves. We simulated one offspring for each of *N* fathers for each mother, where *N* is the most likely number of full sibship families in observed families rounded to the nearest integer (supplementary figure 2), to give a total of 433 offspring. We then calculated (1) the posterior probability of paternity for each true sire in the absence of information about dispersal, (2) mean distance between mothers and known true sires and (3) weighted mean distance between mothers and possible all sires, weighted by probability of paternity for all candidates with a probability &gt;0.001 of having sired at least offspring. To test dependency on sampling effort we then removed 10%, 20%, 30% or 40% of true sires at random from the dataset and recalculated weighted-mean distances between mothers and sires. We repeated this procedure on 100 replicate simulated datasets.

We then used a second set of simulations to investigate the extent to which a joint analysis of paternity and dispersal could mitigate the inflation of leptokurtosis. We generated offspring as described above assuming exponential dispersal of 50m only, and performed a simplified joint-analysis of paternity and dispersal by grid interpolation. We multiplied the matrix of paternity probabilities from genetic data by a matrix of dispersal probabilities from \\Pr(d\_{mj}|a,b,\\lambda)using combinations of scale (5 \\leq a \\leq 150) and shape (0.1 \\leq b \\leq 2) parameter values, with \\lambdafixed at the posterior mean of 0.82. We calculated the likelihood of each combination of parameters by summing joint probabilities of paternity over each candidate father, and multiplying those likelihoods over all 433 offspring. We then recorded the combination of dispersal parameters with the highest (log) likelihood. We also repeated this procedure for datasets excluding true sires as described above. We repeated these analyses on 100 replicate simulated datasets.

# 3. Results {#MPSection:D30BF7C0-46E3-417E-93EA-8E16C3A912B3}

## Dispersal is leptokurtic {#MPSection:06F797FC-6FED-4805-BAEE-E05112C39AB9}

FAPS**grouped the 937 offspring into between 434 and 437 full-sibling families, and identified 316 mating events between the 60 mothers and candidate fathers with posterior probabilities greater than 0.99 (supplementary figure 2). Mixture parameter \\lambdawas strongly weighted towards signal coming from generalised normal dispersal rather than random draws from the population (mean=0.82, 95% credible intervals = 0.77, 0.93), indicating that most of the mating events reflect real mating events (supplementary figure 4). The data were compatible with broad range of plausible values for the proportion of missing fathers (mean=0.21, 95% credible intervals = 0.06, 0.42; supplementary figure 4), but the posterior distribution was very similar to the prior distribution (supplementary figure 3), suggesting that the data are not informative about missing fathers.

Pollen dispersal distances are characterised by many dispersal events between nearby plants, and a long tail of more distant mating events. 50% of the high-probability fathers were within 40m of the mother; distances to remaining fathers decayed slowly up to a maximum of 2398m. Consistent with this, the full pollen dispersal kernel inferred jointly with sibships and paternity consistently showed shape parameters less than one (posterior mean=0.40, 95% credible intervals: 0.329-0.483). These results indicate strong leptokurtosis in the pollen dispersal kernel.

<div id="MPFigureElement:E0424EB3-1085-42CE-96B3-C56587B263CD">

![](Figure_3.png)

Figure 3: Pollen dispersal kernel for posterior mean parameters (shape=0.40, scale=8.24). Plots show histograms of distances between mothers and candidate fathers, weighted by the posterior probability of paternity. Red curves show GND curve for the same parameter values. Inset: the same data censored to dispersal within 40m.

</div>

## Parents sire more offspring on either side of the hybrid zone {#MPSection:E1404EBD-6262-4CC6-9AB7-683A568A075C}

*A. m. striatum-* and *A. m. pseudomajus*-like fathers were more likely to sire offspring with mothers towards the East and West of the hybrid zone respectively than would be expected under random mating and plausible dispersal distances (<span class="kind elementIndex">**Figure 4** </span>; suplementary figure 6). In the Western-most spatial bin, maternal plants were 5.2- and 3.0-fold more likely to receive pollen from ros/ros and sulf/sulf pollen donors respectively, but 2.1- and 1.6-fold less likely to receive pollen from ROS/ROS and S/+ donors, than would be expect chance given local genotype frequencies. This pattern reverse around the centre of the hybrid zone, and by the 300-600 bin, maternal plants were 1.3- and 1.1-fold more likely to receive pollen from ROS/ROS and SULF/+ plants respectively, but 1.7- and 1.3-fold less likely to receive pollen from ros/ros and sulf/sulf donors. In the Eastern-most spatial bin there were no significant differences between observed and expected siring probabilities between paternal genotypes, but this bin contained only three plants. There were no clear differences in observed and expected probabilities of receiving pollen from ROS/ros pollen donors across the hybrid zone. These results indicate that parental genotypes have increased pollen export on the sides of the hybrid zone where they are dominant.

<div id="MPFigureElement:14A31887-A9DE-4A5D-A296-0ED46F98B19D">

![](Figure_4.png)

Figure 4: Pollen export of flower colour genotypes across the hybrid zone. Plots show differemves between observed and expected probabilities that maternal plants receive pollen from each donors of each Rosea and Sulfurea genotype in 300m bins. Figures show the mean and 95% credible intervals based on dispersal parameters from 1000 posterior draws.

</div>

## Leptokurtosis unlikely to be an artifact {#MPSection:1B5A6CEB-F146-4192-953D-2F6C8B6114F3}

In simulations using with complete sampling of candidate fathers, true fathers had mean posterior probabilities of 0.993, with a minimum in any simulation of 0.984, indicating that markers have high statistical power to resolve relationships. However, when no information about dispersal is included in paternity analysis, the remaining uncertainty about paternity causes apparent distance between mates to be inflated by approximately 40m, regardless of the true dispersal distance (<span class="kind elementIndex">**Figure 5** </span>A, left-most panel). This bias increases to several hundred metres when sampling of fathers is incomplete as paternity is incorrectly assigned to unrelated individuals further from the mother, inflating the tail of the dispersal kernel (<span class="kind elementIndex">**Figure 5** </span>A).

In contrast, when dispersal is inferred jointly with paternity, the inferred shape parameter values of the dispersal kernel are close to the true value of one (<span class="kind elementIndex">**Figure 5** </span>B). In fact, when sampling of father is complete, the shape parameter is somewhat inflated (i.e. the pollen dispersal kernel is inferred to decay more rapidly than the exponential distribution). This inflation decreases as the proportion of missing fathers increases, and begins to underestimate the shape parameter of the dispersal kernel when many fathers are unsampled.

<div id="MPFigureElement:2FEB7E7D-5C82-4DBA-B320-0F34321434FD">

![](Figure_5.png)

Figure 5: Simulation results under exponentially-distributed dispersal distances. A: Mean distances between mothers and true sires, and between mothers and all candidate sires weighted by probability of paternity in the absence of dispersal information. Subpanels show the proportion of missing fathers in each dataset. B: Inferred GND shape parameter from a joint analysis of paternity and dispersal.

</div>

# 4. Discussion {#MPSection:ACE6A524-9318-4E80-A0CF-84A447DAE2FA}

In this study we performed joint analysis of paternity, sibships and pollen dispersal in an *A. majus* hybrid zone. We find that dispersal is leptokurtic, with half of seeds being sired by pollen donors within 40m of the maternal plant, and remaining donors being spread up to more than 2km away. Based on this dispersal we find

## Can we trust the dispersal kernel? {#MPSection:F55B305B-95D1-47E9-B2E1-0DEBCF207069}

## Selection for parental phenotypes {#MPSection:E56A728E-D3DF-4972-8C89-91628CA25080}

Male fitness and maintenance of the hybrid zone

Matches pedigree. Implies GxE

Why might parents do better? ‘Good’ colours that contrast with background. Frequency-dependence. Odour.

# 5. Bibliography {#MPSection:9BFEA772-C48A-421C-88F7-4EED2D0F68BC}

<div id="MPBibliographyElement:AEC5F8AE-819A-45B1-9C2E-1C05990D227F" xmlns="http://www.w3.org/1999/xhtml">

Austerlitz, F., Dick, C. W., Dutech, C., Klein, E. K., Oddou-Muratorio, S., Smouse, P. E., & Sork, V. L. (2004). Using genetic markers to estimate the pollen dispersal curve. *Molecular Ecology*. Wiley Online Library.

Clark, J. S. (1998). [Why Trees Migrate So Fast: Confronting Theory with Dispersal Biology and the Paleorecord](https://doi.org/10.1086%2F286162). *The American Naturalist*, *152*(2), 204–224. [doi:10.1086/286162](https://doi.org/10.1086%2F286162)

Ellis, T. J., Field, D. L., & Barton, N. H. (2018). [Efficient inference of paternity and sibship inference given known maternity via hierarchical clustering](https://doi.org/10.1111%2F1755-0998.12782). *Molecular Ecology Resources*, *18*(5), 988–999. [doi:10.1111/1755-0998.12782](https://doi.org/10.1111%2F1755-0998.12782)

Hadfield, J., Richardson, D., & Burke, T. (2006). Towards unbiased parentage assignment: combining genetic, behavioural and spatial data in a Bayesian framework. *Molecular Ecology*. Wiley Online Library.

Nadarajah, S. (2005). A generalized normal distribution. *Journal of Applied Statistics*. Taylor & Francis.

Neff, B. D., Repka, J., & Gross, M. R. (2001). [A Bayesian Framework for Parentage Analysis: The Value of Genetic and Other Biological Data](https://doi.org/10.1006%2Ftpbi.2001.1520). *Theoretical Population Biology*, *59*(4), 315–331. [doi:10.1006/tpbi.2001.1520](https://doi.org/10.1006%2Ftpbi.2001.1520)

Pemberton, J. M. (2008). [Wild pedigrees: the way forward](https://doi.org/10.1098%2Frspb.2007.1531). *Proceedings of the Royal Society B: Biological Sciences*, *275*(1635), 613–621. [doi:10.1098/rspb.2007.1531](https://doi.org/10.1098%2Frspb.2007.1531)

TASTARD, E., FERDY, J.-B., BURRUS, M., THÉBAUD, C., & ANDALO, C. (2011). [Patterns of floral colour neighbourhood and their effects on female reproductive success in an Antirrhinum hybrid zone](https://doi.org/10.1111%2Fj.1420-9101.2011.02433.x). *Journal of Evolutionary Biology*, *25*(2), 388–399. [doi:10.1111/j.1420-9101.2011.02433.x](https://doi.org/10.1111%2Fj.1420-9101.2011.02433.x)

Wang, J. (2007). [Parentage and sibship exclusions: higher statistical power with more family members](https://doi.org/10.1038%2Fsj.hdy.6800984). *Heredity*, *99*(2), 205–217. [doi:10.1038/sj.hdy.6800984](https://doi.org/10.1038%2Fsj.hdy.6800984)

Whibley, A. C. (2006). [Evolutionary Paths Underlying Flower Color Variation in Antirrhinum](https://doi.org/10.1126%2Fscience.1129161). *Science*, *313*(5789), 963–966. [doi:10.1126/science.1129161](https://doi.org/10.1126%2Fscience.1129161)

</div>
