#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]


NumericVector mMRFCsampler(NumericMatrix Data, int n, int nNodes, NumericVector type_c,
                   NumericVector levels, int nIter, NumericMatrix thresh_m,
                   NumericMatrix graphe, IntegerVector inde) {


NumericVector indicatorstate;
NumericVector potcat;
NumericVector potcon;

for(int p=0; p<n; p++) { // loop cases


    // generate initial states
    for (int node=0; node<nNodes; node++) {
      if(type_c[node] == 1)
      Data(p,node) =  (double)round(R::runif(0.5, levels[node]+.5));
      if(type_c[node] == 2)
      Data(p,node) =  (double)R::rnorm(1,1);
      if(type_c[node] == 3)
      Data(p,node) =  (double)R::rpois(1);
      if(type_c[node] == 4)
      Data(p,node) =  (double)R::rexp(1);
    }


    for(int iter = 0; iter<nIter; iter++) { // loop iterations

      for(int node=0; node<nNodes; node++) { // loop nodes

        // CATEGORICAL
        if(type_c[node]==1)
        {

          NumericVector potential_stor(levels[node]); //storage for potential of each category

          for(int l = 0; l<levels[node]; l++) { // loop categories


            //create storage for interaction terms
              // gives an indicator how many cat vs. cont variables (to create approporate storage)
            NumericVector ind_type(nNodes);
            for(int i=0; i<nNodes; i++)
            {
              if(type_c[i]==1)
              ind_type[i] = 1;
              else
              ind_type[i] = 0;
            }
              //create storage objects
            potcat = rep(0,(sum(ind_type)-1)); //we know we deal type_c[node] = categorical
            potcon = rep(0,(nNodes - sum(ind_type)));
            int cat = 0; //"manual" indices, because depends on sequence of con-cat-....
            int con = 0;

            for(int node2 = 0; node2<nNodes; node2++) { // loop other nodes

              if(node2!=node) {

              if(type_c[node2]==1)
              {

                  // which categoryy is present
                  indicatorstate = rep(0,levels[node2]);
                  for(int k=0;k<levels[node2]; k++) {
                    if(Data(p,node2)==(k+1)) indicatorstate[k]=1;
                    }


                  // cut out correct part of graphe
                  IntegerVector ind_rows(inde.size());
                  for(int k=0;k<inde.size(); k++) {if(inde[k]==(node+1)) ind_rows[k]=1;}

                  IntegerVector ind_cols(inde.size());
                  for(int k=0;k<inde.size(); k++) {if(inde[k]==(node2+1)) ind_cols[k]=1;}

                  NumericMatrix n_m1(sum(ind_rows), ind_cols.size());

                  int rowc = 0;
                  for(int r = 0; r<ind_cols.size(); r++) {
                    if(ind_rows[r]==1) {
                      n_m1(rowc,_) = graphe(r,_);
                      rowc++;
                    }
                  }

                NumericVector n_v1 = n_m1(l,_); //cut out level/state
                NumericVector n_v2 = n_v1[ind_cols==1]; //cut out columns (node 2)

                potcat[cat] = sum(n_v2 * indicatorstate); // interaction term of node2 categorical variable
                cat++;


              } else { //if node2=(con)

                // cut out correct part of graphe (same as above)
                IntegerVector ind_rows(inde.size());
                for(int k=0;k<inde.size(); k++) {if(inde[k]==(node+1)) ind_rows[k]=1;}

                IntegerVector ind_cols(inde.size());
                for(int k=0;k<inde.size(); k++) {if(inde[k]==(node2+1)) ind_cols[k]=1;}

                NumericMatrix n_m1(sum(ind_rows), ind_cols.size());

                int rowc = 0;
                for(int r = 0; r<ind_cols.size(); r++) {
                  if(ind_rows[r]==1) {
                    n_m1(rowc,_) = graphe(r,_);
                    rowc++;
                  }
                }

                NumericVector n_v1 = n_m1(l,_); //cut out level/state

                NumericVector n_v2;
                n_v2 = n_v1[ind_cols==1]; //cut out columns (node 2)


                NumericVector dummy_postcon; //apparently necessary to get from vector->double
                dummy_postcon = n_v2 * Data(p,node2); // interaction term of node2 (continuous variable)
                potcon[con] = dummy_postcon[0]; // interaction term of node2 (continuous variable)
                con++;

              } // end: if node2: cond vs. con

              } // end: if node2!=node

            } // end: other nodes

            potential_stor[l] =  exp(thresh_m(node,l) + sum(potcat) + sum(potcon));

        } // end: loop categories

          //sample stae proportional to potentials
          double samplesum = sum(potential_stor);
          double randomNumber = R::runif(0, samplesum);
          double newsum = 0;
          for (int ind=0; ; ind++) {
            newsum += potential_stor[ind];
            if (randomNumber <= newsum) {
              Data(p,node) = ind+1;
              break;
            }
          }

              } else { // CONTINUOUS (most code just copied from above)


            //create storage for interaction terms
              // gives an indicator how many cat vs. cont variables (to create approporate storage)
            NumericVector ind_type(nNodes);
            for(int i=0; i<nNodes; i++)
            {
              if(type_c[i]==1)
              ind_type[i] = 1;
              else
              ind_type[i] = 0;
            }
              //create storage objects
            potcat = rep(0,(sum(ind_type)));
            potcon = rep(0,(nNodes - sum(ind_type)-1)); //we know we deal type_c[node] = continuous
            int cat = 0; //"manual" indices, because depends on sequence of con-cat-....
            int con = 0;


              for(int node2 = 0; node2<nNodes; node2++) { // loop other nodes

                if(node2!=node) {

                  //node2 = cat
                  if(type_c[node2]==1) {


                  // which categoryy is present
                  indicatorstate = rep(0,levels[node2]);
                  for(int k=0;k<levels[node2]; k++) {
                    if(Data(p,node2)==(k+1)) indicatorstate[k]=1;
                    }


                  // cut out correct part of graphe
                  IntegerVector ind_rows(inde.size());
                  for(int k=0;k<inde.size(); k++) {if(inde[k]==(node+1)) ind_rows[k]=1;}

                  IntegerVector ind_cols(inde.size());
                  for(int k=0;k<inde.size(); k++) {if(inde[k]==(node2+1)) ind_cols[k]=1;}

                  NumericMatrix n_m1(sum(ind_rows), ind_cols.size());

                  int rowc = 0;
                  for(int r = 0; r<ind_cols.size(); r++) {
                    if(ind_rows[r]==1) {
                      n_m1(rowc,_) = graphe(r,_);
                      rowc++;
                    }
                  }

                NumericVector n_v1 = n_m1(0,_); //to use the same structure is above; there is of course only 1 row
                NumericVector n_v2 = n_v1[ind_cols==1]; //cut out columns (node 2)

                potcat[cat] = sum(n_v2 * indicatorstate); // interaction term of node2 categorical variable
                cat++;

                    //node2 = con
                    } else {


                // cut out correct part of graphe (same as above)
                IntegerVector ind_rows(inde.size());
                for(int k=0;k<inde.size(); k++) {if(inde[k]==(node+1)) ind_rows[k]=1;}

                IntegerVector ind_cols(inde.size());
                for(int k=0;k<inde.size(); k++) {if(inde[k]==(node2+1)) ind_cols[k]=1;}

                NumericMatrix n_m1(sum(ind_rows), ind_cols.size());

                int rowc = 0;
                for(int r = 0; r<ind_cols.size(); r++) {
                  if(ind_rows[r]==1) {
                    n_m1(rowc,_) = graphe(r,_);
                    rowc++;
                  }
                }

                NumericVector n_v1 = n_m1(0,_); //cut out level/state
                NumericVector n_v2;
                n_v2 = n_v1[ind_cols==1]; //cut out columns (node 2)


                NumericVector dummy_postcon; //apparently necessary to get from vector->double
                dummy_postcon = n_v2 * Data(p,node2); // interaction term of node2 (continuous variable)
                potcon[con] = dummy_postcon[0]; // interaction term of node2 (continuous variable)
                con++;


                    } // end if: node2: cat vs. cont

                  } // end if: node2=node

                } // end: loop: other nodes (in node=continuous)

        double natpar;
        natpar = thresh_m(node,0) + sum(potcat) + sum(potcon);

        if(type_c[node]==2) { //gauss
          if (exp(natpar)>(10^70)) Rcpp::stop("mu > 10^70 for Gaussian node");
          Data(p,node) = R::rnorm(natpar,1);
        } else if(type_c[node]==3) { //pois
          if (exp(natpar)<=0) Rcpp::stop("Lambda <= 0 for poisson node.");
          if (exp(natpar)>(10^70)) Rcpp::stop("Lambda > 10^70 for Poisson node");
          Data(p,node) = R::rpois(exp(natpar));
        } else if(type_c[node]==4) { //exp
          if (1/(-natpar)<=0) Rcpp::stop("Lambda <= 0 for exponential node.");
          if (1/(-natpar)>(10^70)) Rcpp::stop("Lambda > 10^70 for Exponential node");
          Data(p,node) = R::rexp(1/(-natpar)); //"bug" in Rcpp!! it doesnt take the input as lambda, but as 1/lambda
        }


        } // end: if categorical vs. continuous

      } // end: loop nodes

  } // end: loop interactions

} // end: loop cases


return Data;

} // end of function


