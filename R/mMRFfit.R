


#### FUNCTION ####

mMRFfit <- function(
  data, #data matrix, col=variables
  type, #data type for col 1:ncol; c=categorical, g=gaussian, p=poisson, e=exponential
  levels, #number of categories of categorical variables, continuous variables have level=1
  d, #maximal clique size
  rule.w="OR", #parameter-aggregation of categorical variables
  rule.b="OR"  #parameter-aggregation across symmetric estimation
)

{

  # step 1: sanity checks & info from data
  stopifnot(ncol(data)==length(type))
  stopifnot(ncol(data)==length(levels))
  n <- nrow(data)
  nNode <- ncol(data)
  
  # step 2: prepare data
  data <- as.data.frame(data)
  colnames(data) <- paste("V",1:nNode, sep="") #necessary for formula input
  emp_levels <- numeric(nNode)
  for(z in 1:nNode)
  {
    emp_levels[z] <- 1
    if(type[z]=="c") {
      emp_levels[z] <- length(unique(data[,z]))
      data[,z] <- as.factor(data[,z]) #turn categoricals into factors  
      }
  }

  #indexing-dummy that helps us to put parameters involving categorical variables in the right place
  dummy_par.sort <- logical()
  dummy_levels <- numeric()
  for(i in 1:nNode)
  {
    lev.e <- 1
    tar.e <- TRUE
    if(emp_levels[i]>1)
    {
      tar.e <- c(FALSE,rep(TRUE,emp_levels[i]-1))
      lev.e <- emp_levels[i]-1
    }
    dummy_par.sort <- c(dummy_par.sort, tar.e)
    dummy_levels <- c(dummy_levels, lev.e)
  }
  dummy_par.sort <- as.logical(dummy_par.sort)
  dummy_matrix.int <- cbind(dummy_levels, 1:length(dummy_levels))
  dummy.ind <- unlist(apply(dummy_matrix.int,1,function(x) { rep(x[2],x[1])}))

  dummy_matrix <- cbind(emp_levels, 1:length(emp_levels))
  ind <- as.numeric(unlist(apply(dummy_matrix,1,function(x) { rep(x[2],x[1])})))
  

  # step 3: create storage for parameters
  model.par.matrix <- matrix(0, sum(emp_levels), sum(emp_levels))
  m_lambdas <- matrix(0,nNode,2) #storing lambda threshold and actual lambda


  # step 4: estimation
  for(v in seq_len(nNode))
  {
    # step 4.1: compute design matrix (adding interactions as a function of d)
    if(d>(nNode-1)) {
      stop("Order of interactions can be maximal the number of predictors!")
    } else if (d==1){ form <- as.formula(paste(colnames(data)[v],"~ (.)"))
    } else { form <- as.formula(paste(colnames(data)[v],"~ (.)^",d)) }


    X <- model.matrix(form, data=data)[,-1]

    #define link function
    if(type[v]=="c") {
      fam <- "multinomial"
    } else if(type[v]=="g" | type[v]=="e") {
      fam <- "gaussian"
    } else if(type[v]=="p") {
      fam <- "poisson"
    }

    # step 4.2: call glmnet
    y <- data[,v]
    fit <- cv.glmnet(X, y, family=fam, alpha=1, nfold=10, type.measure = "deviance", type.multinomial="ungrouped")
    coefs <- coef(fit, s="lambda.min")

    
    #to do:  taking out intercept only for binary!
    
    
    #list to matrix
    coefsm <- matrix(do.call(rbind,lapply(coefs, as.numeric)),nrow=emp_levels[v])[,-1]


    # step 4.3: save lambda + save & apply tau threshold
    m_lambdas[v,2] <- bound <- sqrt(d) * sqrt(sum(coefsm^2)) * sqrt(log(nNode)/n)
    m_lambdas[v,1] <- fit$lambda.min
    coefsm[abs(coefsm)<bound]<-0 # apply tau threshold


    # step 4.4: write into model.par.matrix

    #select corresponding row in model par matrix & fill in
    #get correct dummy
    dummy_par.sort.v <- dummy_par.sort[ind!=v]

    #select corresponding row in model par matrix & fill in
    if(emp_levels[v]==1) {
      exp.n.c <- length(model.par.matrix[ind==v,ind!=v][dummy_par.sort.v])
      model.par.matrix[ind==v,ind!=v][dummy_par.sort.v] <- coefsm[1:(exp.n.c)]

    } else {

      for(L in 1:emp_levels[v])
      {
        exp.n.c <- length(model.par.matrix[ind==v,ind!=v][,dummy_par.sort.v][L,])
        model.par.matrix[ind==v,ind!=v][,dummy_par.sort.v][L,] <- coefsm[L,1:(exp.n.c)]
      }
    }

    } # end variable-loop


  # step 5: derivates on model parameter matrix

  # 5.1: aggregate on within categories

  #select only colums where paramater are actually estimated (glmnet estimates k-1 not k parameters)
  m.p.m <- model.par.matrix[,dummy_par.sort]

  # averaging over  columns
  m.p.m.1 <-  t(apply(m.p.m, 1, function(x) {

    out <- (numeric(0))
    for(i in 1:nNode)
    {
      out.n <- mean(abs(x[dummy.ind==i]))
      if(rule.w=="AND") {
        out.n <- out.n * (sum(x[dummy.ind==i]==0)<1) #the second term = 0 when not all coefficients are nonzero
      }
      out <- rbind(out, out.n)
    }
    out <- matrix(out, nrow=1)
  }))

  # averaging over rows
  m.p.m.2 <-  apply(m.p.m.1, 2, function(x) {
    out <- numeric()
    for(i in 1:nNode)
    {

      out.n <- mean(abs(x[ind==i]))
      if(rule.w=="AND") {
        out.n <- out.n * (sum(x[ind==i]==0)<1) #the second term = 0 when not all coefficients are nonzero
      }
      out <- rbind(out, out.n)
    }
    out <- matrix(out, ncol=1)
  })


  ### 5.3: aggregate across diagonal

  if(rule.b=="AND") {
    m.p.m.2_nonzero <- m.p.m.2!=0
    m.p.m.2_nonzero <- m.p.m.2_nonzero * t(m.p.m.2_nonzero)
    m.p.m.2 <- m.p.m.2 * m.p.m.2_nonzero
  }

  #dichotomize
  m.p.m.2_d <- (m.p.m.2!=0)*1

  #make symmetric
  m.p.m.2_d <- pmax(m.p.m.2_d, t(m.p.m.2_d))
  m.p.m.2 <- pmax(m.p.m.2, t(m.p.m.2))
  model.par.matrix_d <- (model.par.matrix!=0)*1

  # step 6: output

  output_list <- list("wadj"=m.p.m.2, "adj"=m.p.m.2_d, "wpar.matrix" = model.par.matrix, "par.matrix" = model.par.matrix_d, "lambda"=m_lambdas)
  
  return(output_list)
  
} # end function


