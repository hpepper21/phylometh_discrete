getwd()
setwd("/Users/elipepper/Documents/GitHub")

r <- getOption("repos")
r["CRAN"] <- "http://cran.cnr.berkeley.edu/"
options(repos = r)
install.packages("yearn")
yearn::yearn(ape)
yearn::yearn(geiger)
yearn::yearn(phytools)
yearn::yearn(phangorn)
yearn::yearn(corHMM)
library(rotl)
library(datelife)
gecko.trees <- datelife_search(input="Gekkonidae", get_spp_from_taxon=TRUE)
#gek.id<-tnrs_match_names("Gekkonidae")$ott_id
#gek.tree <- tol_subtree(ott_id=gek.id)
gek.tree <- gecko.trees[[which.max(lapply(gecko.trees, Ntip))]]
library(devtools)
plot.phylo(gek.tree, cex=0.3)
print(paste("The Gekkonid tree has ", Ntip(gek.tree), " terminals and ", 
            Nnode(gek.tree), " internal nodes out of ",Ntip(gek.tree)-2,
            " possible, which means it is ", 
            round(100*(Nnode(gek.tree)-1)/(Ntip(gek.tree)-3), 2),
            "% resolved", sep=""))

my.dat<-read.delim("BayesTraits_binary.data.txt", header=FALSE, stringsAsFactors=FALSE)
colnames(my.dat)<-c("taxon", "nocturnal")
#w/ 0 being diurnal & 1 being nocturnal 
head(my.dat)

CleanData<-function(phy,data){
  data.vector <- data$nocturnal
  names(data.vector) <- gsub("_", " ", data$taxon)
  
  cleaned <- treedata(phy, data.vector)
  return(cleaned)
}
clean.discrete<-CleanData(phy = gek.tree,data = my.dat)
pruned.dat <- clean.discrete$data
clean.phy<-clean.discrete$phy
VizualizeData<-function(phy,data){
  plot(phy)
  print(data)
}
VizualizeData(phy = clean.phy,data = pruned.dat)



treedata(clean.phy, pruned.dat)

#using parsimony to look at ancestral states:
cleaned.discrete.phyDat <- phangorn::phyDat(pruned.dat, type = "USER", levels = c("0","1"), return.index = TRUE)
anc.p <- phangorn::ancestral.pars(clean.phy, cleaned.discrete.phyDat)
pdf(file="discrete.pars.tree.pdf", width=10, height=10)
pars.tree <- plotAnc(clean.phy, anc.p, 1, cex.pie = 0.1, pos = "bottomleft", show.tip.label = F, no.margin = F)
tiplabels(text = clean.phy$tip.label, cex = 0.2, frame = "none", adj = c(0.05,0.5), offset = 0.75)
dev.off()

#plotting the likelihood estimates:
anc.ml <- ancestral.pml(pml(clean.phy, cleaned.discrete.phyDat), type = "ml")
pdf(file="discrete.ml.tree", width = 10, height = 10)
ml.tree <- plotAnc(clean.phy, anc.ml, 1, cex.pie = 0.1, pos = "bottomleft", show.tip.label = F, no.margin = F)
tiplabels(text = clean.phy$tip.label, cex = 0.2, frame = "none", adj = c(0.05,0.05), offset = 0.75)
dev.off()

#using the ace() to obtain 95% CI on ancestral state values
anc.states <- ace(pruned.dat, clean.phy, type = "discrete", CI = T)
anc.states$rates #gives you the mle of the transition rates
#the big difference b/w anc.p & anc.ml is that only anc.ml methods can provide CIs on the ancestral states

#anc.p is attempting to find the distribution of ancestral states w/in
#a given tree which minimizes the total number of character state changes 
#that would be necessary to explain the states observed at the tips of the tree;
#under parsimony, the rate of change is assumed to be relatively low;
#another debilitating assumption is that ancestors could have had one state (diurnal)
#or another (nocturnal) but not both at the same time;

#anc.ml is attempting to find the parameter values that maximize the 
#probability of the data given the hypothesis; this method treats the 
#character states at internal nodes as parameters; the rate
#of change is estimated from the data and then secondarily 
#used to infer the probability of an ancestral area; the ml method
#might be able to fit the biological reality better than parsimony
#in most cases, due to the fact that ml can model multiple changes
#along a single branch & takes branch-lengths into account;

###Biological Questions:
library(corHMM)
##How can you estimate transition rates b/w states?
#using rayDISC() to estimate transition rates & 
#ancestral states for two binary characters (diurnal-0 vs. nocturnal-1)
pruned.dat=data.frame(species=rownames(pruned.dat), trait=pruned.dat[,1], stringsAsFactors = FALSE)
trans.rate <- rayDISC(phy = clean.phy, data = pruned.dat, ntraits = 1, rate.mat = NULL, model = "ARD", node.states = "marginal")
##How could you examine if transition rates are equal?
ER.trans.rate <- rayDISC(phy = clean.phy, data = pruned.dat, ntraits = 1, rate.mat = NULL, model = "ER", node.states = "marginal")
##running the tree and discrete dat through the Lewis (2001) MKV model
#via lewisMkv()
lewis.mkv <- lewisMkv(phy = clean.phy, data = pruned.dat, include.gamma = FALSE, include.beta = FALSE)
##how could you test order of state evolution?




