#######################################################
# Samples from a mixed Markov random field
# jonashaslbeck@gmail.com
# January,2015
#######################################################


#setwd("G:\\_THESIS\\mMRF\\src")
#sourceCpp("mMRFSampler_Ccore.cpp")



mMRFsampler <- function(
  n, #number of samples
  type, #type of data/from which distribution comes the data
  lev, # Number of levels
  graph, #graph structure
  nIter = 1000, #number of samples for each node
  thresh #thresholds, for every node (& category)
){
  
  lev <- as.numeric(lev)
  # some checks on the input
  stopifnot(isSymmetric(graph))
  stopifnot(length(type)==nrow(graph))
  stopifnot(length(type)==length(lev))
  stopifnot(sum(diag(graph)!=0)==0)
  
  # check: is the within-gaussian covariance matrix positive definite?
  g_cov <- graph[type=="g",type=="g"] #select within-gaussian cov matrix
  
  if(length(g_cov)>1)
  {
    g_covm <- matrix(g_cov, nrow(g_cov), ncol(g_cov))
    diag(g_covm) <- 1 # we are working with unit variance
    stopifnot(is.positive.definite(g_covm))
  }
  
  
  graphe <- f_set_specified_effects(graph, type, lev, thresh) #create model.parameter.matrix
  
  
  nNodes <- ncol(graph)  # number of nodes
  Data <- matrix(0, n, nNodes)  # create empty data matrix
  inde <- as.numeric(colnames(graphe))
  colnames(graphe) <- rownames(graphe) <- NULL
  
  
  #transform thresh into a matrix (C doesnt take lists)
  thresh_m <- matrix(0, nrow=length(thresh), ncol=max(lev))
  for(t in 1:nNodes) {
    thresh_m[t,1:length(thresh[[t]])] <- thresh[[t]]
  }
  
  #transform types into integers
  type_c <- numeric(nNodes)
  type_c[type=="c"] <- 1
  type_c[type=="g"] <- 2
  type_c[type=="p"] <- 3
  type_c[type=="e"] <- 4
  
  #CALL C CORE
  c_out <- mMRFCsampler(Data, n, nNodes, type_c, lev, nIter=100, thresh_m, graphe, inde)
  
  return(c_out)
  
} # end of function



graph <- matrix(0,3,3)
graph[1,2] <- graph[2,1] <- 1
thresh <- list(c(0,0,0), c(0,0,0), c(.1))

mMRFsampler(n=10, type=c("c", "c", "g"), lev=c(3,3,1), graph = graph, thresh=thresh)




