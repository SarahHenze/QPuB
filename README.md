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
 
   1.  **R.utils** for processing the command line arguments.
   
   2. **tictoc** provides timing function
   
   3. **mcmcse** to compute the effective sample size.
   
   4. **tmvtnorm**  to generate the truncated multivariate normal distri-
bution.
   5. **corpcor** to make a matrix postive definite. 
   
   6. **dqRNGkind** to generate fast pseudo random numbers.
   
   7. **coda** to produce an mcmc object for further anaylsis of the
Markov chain.

   8. **base** to execute Base R functions

   9. **matrixcalc** to check whether a matrix is positive definite or not

   10. **ggdmc** to generate univariate truncated normal distributon.

   11. **grDevices** for cairo pdf.

   12. **matrixStats** to compute standard deviation of columns of a ma-
trix.
