## QPuB - Quantification of peptide using Bayesian approach

   QPuB employs Bayesian statistical inference based on Markov chain Monte Carlo sampling to learn the posterior distribution of the converison factors for the peptide products
   
* **Execution**

  Execution of QPuB proceeds through two simple steps
  
    * manually download and unzip the .zip file in your local computer
    
    * after going to the **QPuB** folder, run the following command from the terminal 
    
    ```R
       Rscript runQPuB.r -fol <INPUT FOLDER> -pf input.txt -tim <file contains timepoints> -titr <file contains titration>
    ```
    
 *  **Output**
 
    **boxplot\_chain.pdf** &nbsp; boxplot corresponding to the distributions of conversion factors  
	 
   	**chain\_XYZ.RData**   &nbsp; stores trace of the Markov chain or the time series of the parameter draws  
			
   	**chain\_backscaled.RData** &nbsp;   stores the Markov chain after dividing the products by the scaling factor pre-calculated according to the substrate concentration     
	 
 	**conc\_means\_Nrep.csv, conc\_sd\_Nrep.csv** &nbsp; means and standard deviations of absolute concentrations of the products at different time points for the **Nrep**th replicate, respectively  
	
	**massdeviation.png** &nbsp; plot of total mass deviation of the products over Monte Carlo iterations                                                                                                   
 
	**massdev.RData**  &nbsp;  stores mass deviation of the products over Monte Carlo iterations   
 
	**relation.png**   &nbsp; plot of relation between the estimated conversion factors and the peptide lengths                                                                                          
 
 	**runQPuB\_ROUT.txt** &nbsp; contains the initial parameters that were feed to the algorithm and acceptance rates as the chain progresses. The output of any print command directly goes into the file. 
 
 
 	**chain\_Niter.pdf**   &nbsp; trace plots of the Markov chain at Niter\(^{th}\) iteration   
 
 	**residuals\_M.txt**   &nbsp; residual plot at M\(^{th}\) iteration. The plot display mass deviation for individual amino acids of the substrate.
 
 	**statistics.csv**    &nbsp; summary statistic for the Markov chain  
      
    
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



