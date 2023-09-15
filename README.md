# QuidiPhydy

*Package development*

The core algorithm that we used for inference of transmission counts on the phylogenetic tree was previously developed (Naderi et al.). During this project, we modified the code so that it is a more generally usable function in R. Additionally, we formatted the code into an R package called QuidiPhyDy (https://github.com/Saannah/QuidiPhydy) which has been tested locally and also on the automatic testing platform for CRAN. We plan to further extend the package to add more features and functions to it. So far, it only contains one function: Spillover().

We also added the feature to distinguish between primary and secondary transmissions. First, all branches along with a state change has happened will be identified. Then, the mosst basal ones are marked as primary spillover events, such that nested events aren't reported as independent occurences.

This function takes the following as input:

``` tree ```
This is a tree object of type nexus. Preferably it should be a time tree so that the branch lengths correspond to time (this way the inference makes more sense in terms of time). I will add the option for the function to accept trees of other common types (Newick) in the future.

``` trait```
This is an array that describes the traits or states of the tips of the tree. It should be an array whose length is the same as the number of the tips of the tree, and the order of the states should correspond to the order that tips of the tree appear in the tree object. This is pretty easy to arrange, since you can get the tip labels of most regular tree objects in R on its tree$tip.label attribute. Currently, the function only supports only two states for the tips in this version. So there should be exactly two unique values in this vector.

```fromto = c("human", "animal")```
This input determines the direction in which you wish to find transmissions. It should be an array of two, and should be compatible with the trait variable (should include only the two unique values in trait).

```method```
Currently this argument has been set to default = "ML". This is the method that the ancestral state reconstruction will be done with. In the future, we might add the option to do this with a maximum parsimony method as well.

The function outputs a list. This is what this list looks like:

```
intros = list( numberOfEvents = length(unlist(MRCAout)),
                 MRCAout = MRCAout, MRCAin = MRCAin, reconstruction = all,
                 all_transitions = trans_nodes)
```
 ```numberOfEvents ```The number of primary spillovers that were found in the desired directions determined by the user.

``` MRCAout``` The node numbers corresponding to parent node of the primary transmission branches found.

 ```MRCAin``` The node numbers corresponding to child node of the primary transmission branches found.

```reconstruction``` A table describing the ancestral state reconstruction probabilities for every internal node.

```all_transitions``` A list of all transmissions found, identified by their parent node of the transmission branches. This is useful for identifying/counting secondary spillovers. The difference between this list and MRCAout would give the list of secondary spillovers.

If you plan to use this function and need more help, please feel free to contact me at sana.naderi@mcgill.ca.
