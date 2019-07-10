## QPuB - Quantification of peptide using Bayesian approach

   QPuB employs Bayesian statistical inference based on Markov chain Monte Carlo sampling to learn the posterior distribution of the converison factors for the peptide products
   
* **Execution**

  Execution of QPuB proceeds through two simple steps
  
    * manually download and unzip the .zip file in your local computer
    
    * run the following command in your terminal
    
    ```R
       Rscript runQPuB.r -fol ifolder -pf paramfile -tim tp -titr titration
    ```
 * **Dependencies**
 
   1. **R.utils** -for processing the command line arguments.
   
   2. **tictoc** -provides timing function
   
   3. **mcmcse** to compute the effective sample size.


tmvtnorm to generate the truncated multivariate normal distri-
bution.
corpcor to make a matrix postive definite.
dqRNGkind to generate fast pseudo random numbers.
coda to produce an mcmc object for further anaylsis of the
Markov chain.
base to execute Base R functions
matrixcalc to check whether a matrix is positive definite or not
ggdmc to generate univariate truncated normal distributon.
grDevices for cairo pdf.
matrixStats to compute standard deviation of columns of a ma-
trix.
