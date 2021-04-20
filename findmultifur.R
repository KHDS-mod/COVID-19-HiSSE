##
## This is simply a script to see how many nodes are multifucating. Source this
## and it will produce some variables in the R image.
##
suppressPackageStartupMessages({
  library("phytools")
  library("diversitree")
})

tr = read.tree("RAW_PHYLOGENY/nextstrain_ncov_global_timetree.nwk")
ntips = length(tr$tip.label)
internalnodes = {
    x = unique(c(tr$edge))
    x[x > ntips]
}
nodeinfo = data.frame(edge_no = integer(length(internalnodes)),
		      nchild = integer(length(internalnodes)),
		      nchildtip = integer(length(internalnodes)))
for (j in seq_along(internalnodes)) {
    nodeinfo[['nchild']][j] = sum(tr$edge[,1] == internalnodes[j])
    nodeinfo[['nchildtip']][j] = sum(tr$edge[,1] == internalnodes[j] & tr$edge[,2] <= ntips)
    nodeinfo[['edge_no']][j] = internalnodes[j]
}

nnode_multifur = sum(nodeinfo[['nchild']] > 2)
nnode_multifur_contains_tip = sum(nodeinfo[['nchild']] > 2 & nodeinfo[['nchildtip']] >= 1)
nnode_multifur_allchildis_tip = sum(nodeinfo[['nchild']] > 2 & nodeinfo[['nchildtip']] == nodeinfo[['nchild']])
nnode_multifur_allchildis_tip_proportion = sum(nodeinfo[['nchildtip']][nodeinfo[['nchild']] > 2 & nodeinfo[['nchildtip']] == nodeinfo[['nchild']]]) / ntips

ntips_contained_in_multifur = sum(nodeinfo[['nchildtip']][nodeinfo[['nchild']] > 2 & nodeinfo[['nchildtip']] >= 1])
ntips_with_fully_tipsy_mother = sum(nodeinfo[['nchildtip']][nodeinfo[['nchild']] > 2 & nodeinfo[['nchildtip']] == nodeinfo[['nchild']]])

proportion_ntips_contained_in_multifur = ntips_contained_in_multifur/ntips *100
proportion_ntips_with_fully_tipsy_mother = ntips_with_fully_tipsy_mother/ntips *100

dist_tipchildren_of_multifur_mothers = as.data.frame(table(nodeinfo[['nchildtip']][nodeinfo[['nchild']] > 2 & nodeinfo[['nchildtip']] >= 1]))
colnames(dist_tipchildren_of_multifur_mothers) = c('number_of_tip_children', 'count')
dist_tipchildren_of_multifur_mothers[['%count']] = dist_tipchildren_of_multifur_mothers[['count']] / tr$Nnode * 100

dist_children_of_multifur_mothers = as.data.frame(table(nodeinfo[['nchild']][nodeinfo[['nchild']] > 2 & nodeinfo[['nchildtip']] >= 1]))
colnames(dist_children_of_multifur_mothers) = c('number_of_children', 'count')
dist_children_of_multifur_mothers[['%count']] = dist_children_of_multifur_mothers[['count']] / tr$Nnode * 100


dist_tipchildren_of_fully_tipsy_mothers = as.data.frame(table(nodeinfo[['nchildtip']][nodeinfo[['nchild']] == nodeinfo[['nchildtip']]]))
colnames(dist_tipchildren_of_fully_tipsy_mothers) = c('number_of_tip_children', 'count')
dist_tipchildren_of_fully_tipsy_mothers[['%count']] = dist_tipchildren_of_fully_tipsy_mothers[['count']] / tr$Nnode * 100


