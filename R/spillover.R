
#' Transitions finder
#'
#' This function searches for the transitions between states of a trait along a given phylogenetic tree. 
#' It implements ancestral state reconstruction and finds the most recent common ancestor (MRCA) both in and out of the clade of interest.
#' This function can be used to count and locate the number of spillover events (e.g. from animal to human), introduction events (e.g. into a new country or outbreak institution), or for estimating when trait evolution occurred. 
#' 
#' @param tree Phylogenetic tree object of type Nexus
#' @param trait An array describing traits in the tips of the tree. The order should correspond to the order of the tips of the tree. Current version only supports 2 traits
#' @param fromto 2-element array specifying the direction in which transitions should be searched for, should match the values in "trait"
#' @param method Ancestral state reconstruction method, current version supports maximum likelihood which is the default
#'
#' @return transitions object
#' @export
#' @author Sana Naderi and Carmen Lia Murall
#' @references
#'   Murall et al. (2021) *A small number of early introductions seeded widespread transmission of SARS-CoV-2 in Quebec, Canada.* Genome Medicine. 13:169
#'   Naderi et al. (2023) *Zooanthroponotic transmission of SARS-CoV-2 and host-specific viral mutations revealed by genome-wide phylogenetic analysis.* eLife. 12:e83685
#'
#' @examples 
#' fromto = c("animal", "human")
#' fromto = c("country1", "country2")
#' fromto = c("general population", "institution outbreak")
#' fromto = c("without mutation", "with mutation")
#'
transitionsfinder = function(tree,
                       trait,
                       fromto,
                       #species,
                       #number,
                       method = "ml"){


  #tree = ape::read.nexus(paste0("/Users/sana/Documents/AARMS/project/data/",species,"_",number,"_timetree.nexus"))
  #trait = rep("human", length(tree$tip.label))
  #trait[which(stringr::str_detect(tree$tip.label, paste0("/", species, "/")))] = "animal"
  labels = data.frame(id = tree$tip.label, trait = trait)


  #method = "ml"
  ml_asr_model = "ER"
  root_state = NULL
  transition_fromto = fromto


  ################################# data loading and cleaning starts here
  strainlist = tree$tip.label #list of all strain id's we have
  labels= labels[which(labels[,1] %in% strainlist),]
  ################################# data loading and cleaning ends here


  ###### check if tree is rooted
  if (!ape::is.rooted(tree)){
    tree = ape::root(tree, outgroup = 1, resolve.root = TRUE)
  }


  ################################# converting the tree to a fully dichotomous tree, starts here
  dichotomous_tree = ape::multi2di(tree, random = FALSE)   #format change of tree
  dichotomous_tree$edge.length[dichotomous_tree$edge.length<=0] = 1e-8  #set very small values to zero and negative values (needed for ace to run properly)
  ################################# converting the tree to a fully dichotomous tree, starts here



  ################################# creating the labels for reconstruction, starts here
  # "labels" is the matrix that includes the tip labels for the tree
  # the dimension of labels is ntips x 2. the first column has sample IDs and the second column has the states corresponding to each sample.

  states = unique(labels[,2]) #the unique values of the states DO NOT CHANGE THE ORDER OF THIS VECTOR IN THE REMAINDER OF THE CODE
  n_states = length(states)


  mpr_tipstates = rep(0, length(dichotomous_tree$tip.label))

  for (i in 1:n_states){ #loop to assign the tip labels for the traits
    s = states[i]
    whos = labels[which(labels[,2] == s),1]
    mpr_tipstates[match(whos, tree$tip.label)] = i
  }

  dichotomous_tree$tip.label -> names(mpr_tipstates)
  ################################# creating the labels for reconstruction, ends here



  ################################# ML reconstruction ends here
  # labels vector format is same as the MPR analysis, the same vector was used

  if (method == "ml"){
    ml_recon = ape::ace(mpr_tipstates, dichotomous_tree, type = "discrete", model="ER")
  }
  #
  ################################# ML reconstruction ends here


  #### NOTE:
  # the output of the mpr function corresponds to the edge attribute of the tree
  # the output of the ace fucntion increments from the top of the tree (starting from the root)



  comparison = ape::comparePhylo(tree, dichotomous_tree) #compares the original polytomous tree to the dichotomous tree
  if(is.null(comparison$NODES)){
    e1 = length(tree$tip.label)+1
    e2 = (length(tree$tip.label) + (tree$Nnode))
    ns = seq(e1, e2)
    seq = paste0("(", as.character(ns), ")")
    comparison$NODES = data.frame(tree = seq, dichotomous_tree = seq)
  }

  comp = matrix(0, dim(comparison$NODES)[1], dim(comparison$NODES)[2])


  for (i in 1:dim(comparison$NODES)[1]){ #loop for extracing the node number out of the compare$node table returned by comparePhylo()
    for (j in 1:dim(comparison$NODES)[2]){
      string = comparison$NODES[i,j]
      str = stringr::str_extract(string = comparison$NODES[i,j], pattern = "(?<=\\().*(?=\\))")
      comp[i,j] = as.numeric(str) #replacing the string number to numeric value, optional, might have to change later
    }
  }



  ################################# cleaning up ace output starts here

  # the columns of ml_recon$lik.anc have the same name as tip state numbers. meaning that column names will be numbers from 1 to "n_state"


  if (method == "ml"){


    ml_asr = ml_recon$lik.anc #the matrix containing likelihoods for everynode
    ml_colnames = colnames(ml_asr)

    ace_dataframe = matrix(0, dim(comp)[1], 2)
    lik_animal = vector()
    for (i in 1:dim(comp)[1]){
      row = comp[i,2] - length(dichotomous_tree$tip.label)
      s = which.max(ml_asr[row,]) #the states are numbered from 1 to n, so the column number is the state itself
      ace_dataframe[i,1] = i + length(tree$tip.label) #node number in the actual tree
      ace_dataframe[i,2] = s #the state of that node
      lik_animal = c(lik_animal, ml_asr[row,2])
    }

    ace_dataframe = data.frame(node = ace_dataframe[,1], state = ace_dataframe[,2])
  }
  ################################# cleaning up ace output ends here



  # converting transition_fromto to numbers [1 to n_states]
  # THE ORDER OF THE "states" VARIABLE MATTERS!!!
  from = transition_fromto[1]
  to = transition_fromto[2]
  n_from = which(states == from) #numerical value for the given states (from)
  n_to = which(states == to) #numerical value for the given states (to)



  ########## transitions:


  ################################# transition nodes starts here

  #i put this conditions to later be able to manipulate data types that are different in type among different methods, might remove later
  if (method == "ml"){
    all = data.frame(node = c(1:length(tree$tip.label), ace_dataframe$node), state = c(mpr_tipstates, ace_dataframe$state))
  }



  ## add a column to all dataframe
  all$labels = states[1]
  all$labels[which(all$state == 2)] = states[2]

  #sort edge table by the second column (children)
  edge = tree$edge
  edge = data.frame(parent = edge[,1], child = edge[,2])
  edge = edge[order(edge$child),]


  # finding *all* transition events in the tree
  o = vector()
  for (i in 1:dim(edge)[1]){
    s = all[edge$child[i], 2]
    if(s == n_to){
      o = append(o,i)
    }
  }
  transitions_childeren = edge[o,]
  z = vector()
  for (i in 1:dim(transitions_childeren)[1]){
    s = all[transitions_childeren$parent[i], 2]
    if(s == n_from){
      z = append(z,i)
    }
  }
  transitions_parents = transitions_childeren[z,]
  #trans_nodes includes ALL events where an OUT node has an IN child. this includes embedded transitions and singletons
  trans_nodes = sort(unique(transitions_parents$child))
  trans_nodes_store = trans_nodes




  if(length(trans_nodes) == 0){
    intros = list( numberOfEvents = length(unlist(MRCAout)),
                   MRCAout = MRCAout, MRCAin = MRCAin, reconstruction = all,
                   all_transitions = trans_nodes)
    return(intros)
  }


  #function for finding the parents of a node up to the root
  # inputs: tree, node number outputs: a vector of all parents of the node up to the root
  parent_finder = function(tree, nodes){ #function starts here
    parents_list = list()
    for (node in nodes){
      edge = tree$edge
      parents = vector()
      flag = 0
      root = length(tree$tip.label) + 1

      parents = append(parents, node)
      while(flag == 0){
        if (node == root){
          flag = 1
          next
        }
        ind = which(edge[,2] == node)
        node = edge[ind, 1]
        parents = append(parents, node)
        if (node == root){
          flag = 1
        }
      }
      parents_list = append(parents_list, list(parents))
    }
    return(parents_list)
  } #function ends here

  parents_list = parent_finder(tree, trans_nodes)

  #function for finding transition parents of transition nodes
  # this function gets a list of vectors representing nodes ancestors, returns transition events among them


  checked_list = vector()
  transitions_obj = list()
  for (i in 1:length(parents_list)){ #loop starts here
    p = unlist(parents_list[i])
    if (p[1] %in% checked_list){
      next
    }

    parents = vector()
    for (j in 1:length(p)){
      if (p[j] %in% trans_nodes){
        parents = append(parents, p[j])
      }
    }
    checked_list = append(checked_list, parents)
    transitions_obj = append(transitions_obj, list(parents))

  }#loop ends here

  independent_transitions = vector()
  for (i in 1:length(transitions_obj)){ #loop starts here
    p = unlist(transitions_obj[i])
    independent_transitions = append(independent_transitions, p[length(p)])
  }#loop ends here
  independent_transitions = unique(independent_transitions)
  ################################# transition nodes ends here


  new_independent_transitions = independent_transitions
  new_trans_nodes = trans_nodes




  # finding the OUT parents of "new transitions"
  new_independent_transitions_parents = vector()
  for (i in 1:length(new_independent_transitions)){
    ind = which(edge[,2] == new_independent_transitions[i])
    new_independent_transitions_parents = append(new_independent_transitions_parents, edge[ind,1])
  }

  MRCAout = list(new_independent_transitions_parents)
  MRCAin = list(new_independent_transitions)

  intros = list( numberOfEvents = length(unlist(MRCAout)),
                 MRCAout = MRCAout, MRCAin = MRCAin, reconstruction = all,
                 all_transitions = trans_nodes)
  return(intros)

}
