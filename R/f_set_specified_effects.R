##### Set_specified_effects.R

## Wy do we need this function?
# When categorical variables are involved, we have more parameters in the model
# than are provided by a (weighted) adjacency matrix with dim=#variables
# therefore we need no (arbitrarily) set some effects that reflect the conditional dependence
# that is reflected by the present edge (weight)

## What does the function do?
# Input: (weighted) adjacency matrix + effect specifications
# Output: Model-parameter-matrix


############### development


# whats difficult?
# input: categorical variables with different k; + a continuous variable
# output: coherent model-parameter-matrix

#specified rules:
# - interactions between categorical variables: equal categories = nonzero effect;
#   unequal group sizes: k=1,2; j=1,2,3; then; effects: 11; 22; 23;
#   which means the last category of the variable with k<j picks up effects until j
# - interactions categorical * continuous: splithalf (if odd, round down) categories;
#   lower "half" has an effect, upper "half" doesnt



f_set_specified_effects <- function(
  graph, #(weighted) adjacency matrix ;matrix
  type,  # vector indicating the type of the variables ;vector
  levels, #how many levels do variables have (>1 only for categorical) ;vector
  thresh #threshold for each variable; for each category for categoricals ;list
)
{

  ## create empty matrix with correct dimensions

  model.par.matrix <- matrix(0, sum(levels), sum(levels))
  nVar <- nrow(graph)

  r <- 1 #current row; we need to specify that extra, as we use two indices (variables and cat) in one

  for(v in 1:nVar) { #loop variables


    # TYPE OF VAR CHECK
    if(type[v]=="c") {

      # A) CATEGORICAL

      subm.mpm <- matrix(0,levels[v], sum(levels[-v])) #fill a dummy matrix; then include slip in thresholds and merge with full matrix


      for(k in 1:levels[v]) { #loop categories

        c <- 1 #current column index

        for(v2 in (1:nVar)[-v]) #loop other variables
        {

          if(graph[v,v2]!=0) { #check whether edge present between variables; if not, we just set all parameters to 0 (means: we dont to anything)

            #checking type of interaction (type of 2nd variable)
            if(type[v2]=="c") {

              #REWRITE, this only works when k<j, but not when j<k!!


              subm.mpm[k,c:(c+levels[v2]-1)] <- (k==(1:levels[v2])) * graph[v,v2] #check whether j=k

              if(k==levels[v]) { # fill up when k<j or j<k

                if(levels[v]>levels[v2]) #k<j
                {
                  subm.mpm[levels[v2]:k, levels[v2]] <- graph[v,v2]
                } else { #j<k
                  subm.mpm[k, (c+levels[v]-1):(c+levels[v2]-1)] <- graph[v,v2]
                }
              }


            } else { #same for all continuous, no further specification

              subm.mpm[k,c] <-  (k <= levels[v]/2) * graph[v,v2] # smaller "half" has effect, other not; MISSING: multiply actual edge weight

            }
          }#end if: edge present?

          c <- c + levels[v2]
        }#end: other variables

      }#end: categories

      #fill in thresholds


      p1 <- sum(levels[0:(v-1)]) #cols before
      if(v+1 > length(levels)) { p2 <- 0 } else {
        p2 <- sum(levels[(v+1):length(levels)]) } #cols after

      dummy <- 1 #necessary in the case where the submatrix is at the right edge of the subm.mpm.t
      if((ncol(subm.mpm)-p2+1)>ncol(subm.mpm)) { dummy <- 0 }

      subm.mpm.t <- cbind(subm.mpm[,0:p1],
                          matrix(rep(0,levels[v]^2),levels[v],levels[v]),
                          subm.mpm[,((ncol(subm.mpm)-p2+1)*dummy):(ncol(subm.mpm)*dummy)])


      #merge with
      model.par.matrix[r:(r+levels[v]-1),] <- subm.mpm.t #include weights from (weighted) adj matrix


    } else {

      # B) CONTINUOUS: we only need to get the continuous interactions; the rest we get by symmetry

      subv.mpm <- numeric(sum(levels))
      c <- 1 #current column index

      for(v4 in (1:nVar)){ #loop other variables

        # only categorical
        if(type[v4]!="c") { subv.mpm[c] <-graph[v,v4] }
        c <- c + levels[v4]
      }#end: other variables

      #merge
      model.par.matrix[r,] <- subv.mpm

    }

    r <- r + levels[v] #move down the target matrix

  } #end for variables



  ## "mirror" matrix
  pars <- ncol(model.par.matrix)
  for(rows in 1:pars)
  {
    for(cols in 1:pars) {
      if(model.par.matrix[rows,cols]!=0) {
        model.par.matrix[cols,rows] <- model.par.matrix[rows,cols]
      }
    }
  }



  # set diagonal to 0
  diag(model.par.matrix) <- unlist(thresh)


  # add labels indicating which row/col corresponds to which variable; this helps to feed the parameters in the sampler
  dummy_matrix <- cbind(levels, 1:length(levels))
  ind.e <- unlist(apply(dummy_matrix,1,function(x) { rep(x[2],x[1])}))
  colnames(model.par.matrix) <- rownames(model.par.matrix) <- ind.e


  return(model.par.matrix)

  ### end of FUNCTION
}





## testing function

type <- c("c","c", "g", "g", "c", "p", "e")
levels <- c(2, 4, 1, 1, 6, 1, 1)
graph <- matrix(0,7,7)
graph[2,5] <- 1
graph[3,4] <- 0

graph[6,7] <- 1
graph[6,1] <- 1

graph <- pmax(graph, t(graph))
thresh <- list(rep(0,2), rep(0,4), rep(0,1), rep(0,1), rep(0,6), rep(0,1), rep(0,1))

graph

f_set_specified_effects(graph, type, levels, thresh)





graph <- matrix(1,4,4)
diag(graph) <- 0
levels <- c(2,3,1,1)
type <- c("c", "c", "g", "g")
thresh <- list(rep(0,2), rep(0,3), rep(0,1), rep(0,1))
f_set_specified_effects(graph, type, levels, thresh)






