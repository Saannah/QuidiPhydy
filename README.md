# QuidiPhydy

## Inferring introductions for phylodynamic studies of trait evolution, outbreaks, phylogeography, and zoonotic spillovers

*Package development*

The functions found in this package were previously developed for a phylogeography study (Murall et al. 2021) and a zoonotic spillover study (Naderi et al. 2023) of SARS-CoV-2. Due to the large size of the SARS-CoV-2 phylogenetic trees and due to the need to analyze many large trees at once, we developed these functions to run character optimization functions and tree crawling in a more automated fashion. We have modified the code developed for these projects so that they are more generally applicable and will make up a new R package called **QuidiPhyDy** (https://github.com/Saannah/QuidiPhydy).  
The package has been tested locally and on the automatic testing platform for CRAN. We plan to further extend the package to add more features and functions.  

So far, this package only contains one function: `transitionsfinder()`.  
This transitions finding algorithm infers and counts the transmission events between a change of states of a trait along a phylogenetic tree, i.e. an introduction or spillover event is infered to have happened between the most recent common ancestor (MRCA) from the outgroup to the ingroup and the MRCA within the ingroup. This algorithm applies acncestral state reconstruction, ASR, (impmented in the *ape* package, Paradis et al.), using methods such maximum likelihood (ML) or maximum parsimony, and state transition models, such as equal rates, to infer state transitions. Once applied to the tree, the algorithm crawls the tree to count all found transitions, and returns a summary table with all possible events (based on the liberal or conservative choices made by the user when choosing ASR methods). 

We also added the feature to distinguish between primary and secondary transmissions. First, all branches along with a state change has happened will be identified. Then, the mosst basal ones are marked as primary transmission events, such that nested events aren't reported as independent occurences.

This function takes the following as input:

``` tree ```
This is a tree object of type nexus. Preferably it should be a time tree so that the branch lengths correspond to time (this way the inference is of when an introduction or spillover event happened). The option for the function to accept other common tree types (e.g. Newick or divergence) is in the works.

``` trait```
This is an array that describes the traits or states of the tips of the tree. It should be an array whose length is the same as the number of the tips of the tree, and the order of the states should correspond to the order that tips of the tree appear in the tree object. This is pretty easy to arrange, since you can get the tip labels of most regular tree objects in R on its tree$tip.label attribute. Currently, the function only supports only two discrete states for the tips. Thus, here should be exactly two unique values in this vector. In future versions, incorporation of n discrete states and continuous traits will be added. 

```fromto = c("out", "in")```
This input determines the direction in which you wish to find transmission events between the outgroup and the ingroup. It should be an array of two, and should be compatible with the trait variable (should include only the two unique values in trait).  
For instance, in a spillover study `fromto` would be `c("animal", "human")` or for a phylogeographic study it would be `c("not-focal country", "focal country")`. 

```method```
Currently this argument has been set to default = "ML". This is the method with which the ancestral state reconstruction will be done. In the future, we hope to add the option to do this with a other inference methods such as maximum parsimony or empirical Bayes methods.

The function outputs a list. This is what this list looks like:

```
intros = list( numberOfEvents = length(unlist(MRCAout)),
                 MRCAout = MRCAout, MRCAin = MRCAin, reconstruction = all,
                 all_transitions = trans_nodes)
```
 ```numberOfEvents ```The number of primary events that were found in the desired directions determined by the user.

``` MRCAout``` The node numbers corresponding to parent node of the primary transitions branches found.

 ```MRCAin``` The node numbers corresponding to child node of the primary transitions branches found.

```reconstruction``` A table describing the ancestral state reconstruction probabilities for every internal node.

```all_transitions``` A list of all infered transition events found, identified by their parent node of the transition branches. This is useful for identifying/counting secondary events. The difference between this list and MRCAout would give the list of secondary events.

If you plan to use this function and need more help, please feel free to contact Sana at sana.naderi@mcgill.ca.  

**References**  
- Murall et al. (2021) *A small number of early introductions seeded widespread transmission of SARS-CoV-2 in Quebec, Canada.* Genome Medicine. 13:169 [https://link.springer.com/article/10.1186/s13073-021-00986-9]  
- Naderi et al. (2023) *Zooanthroponotic transmission of SARS-CoV-2 and host-specific viral mutations revealed by genome-wide phylogenetic analysis.* eLife. 12:e83685 [https://elifesciences.org/articles/83685]  
- Paradis et al. *ape* package https://cran.r-project.org/web/packages/ape/index.html  


