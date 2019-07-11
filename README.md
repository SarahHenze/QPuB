# QPuB

[![Bayesian inference][bayestat-image]][bayestat-url]
[![Peptide quantification][peptide-image]][peptide-url]
[Immunoproteomics](immunoprot-url)

QPuB (**Q**uantifcation of **p**eptides **u**sing **B**ayesian inference) employs [Bayesian statistical inference](https://en.wikipedia.org/wiki/Bayesian_inference) based on [Markov chain Monte Carlo (MCMC)](https://en.wikipedia.org/wiki/Markov_chain_Monte_Carlo) sampling to learn the posterior distributions of the converison factors for peptide products. 

## Getting started

The following instructions will help you to download QPuB and execute the same to estimate the conversion factors for your peptide products

### Prerequisites

In order to run QPuB, users are strongly advised to install [**R** ≥ 3.5.0](https://www.r-project.org/) and the following packages:

* [**R.utils**](https://cran.r-project.org/web/packages/R.utils/index.html)
   
* [**tictoc**](https://cran.r-project.org/web/packages/tictoc/index.html) 
   
* [**mcmcse**](https://cran.r-project.org/web/packages/mcmcse/index.html) 
   
* [**tmvtnorm**](https://cran.r-project.org/web/packages/tmvtnorm/index.html)  

* [**corpcor**](https://cran.r-project.org/web/packages/corpcor/index.html) 
   
* [**dqRNGkind**](https://www.rdocumentation.org/packages/dqrng/versions/0.2.1/topics/dqRNGkind) 
   
* [**coda**](https://cran.r-project.org/web/packages/coda/index.html) 

* [**base**](https://www.rdocumentation.org/packages/base/versions/3.6.1) 

* [**matrixcalc**](https://www.rdocumentation.org/packages/matrixcalc) 

* [**ggdmc**](https://www.rdocumentation.org/packages/ggdmc)
 
* [**grDevices**](https://www.rdocumentation.org/packages/grDevices) 

* [**matrixStats**](https://www.rdocumentation.org/packages/matrixStats) 

To install the above packages, start your **R** shell, and execute the following command:

```R
 > install.packages("package_name")
 ```
REMEMBER, for Linux the above command will work provided the user has the root access.

## Running QPuB
   
   Once all the prerequisites are met, QPuB is needed to be download either manually from its [Github repository](https://github.com/QuantSysBio/QPuB). The next step is to navigate to **QPuB** sub-folder inside the **QPuB** and run the following command from the terminal 
    
 ```R     
       
       Rscript runQPuB.r -fol <INPUT FOLDER> -pf input.txt -tim <file contains timepoints> -titr <file contains titration>
   
 ```
## Examples

   The **examples** folder contains two toy examples of endopeptidase digestion. 
   In order to run the examples, execute the following command for a date without noise 
   
   ### Example 1: Endopeptidase digestion without noise
   ```R     
       
       Rscript runQPuB.r -fol examples/toy_endo3_nonoise -pf input.txt 
   
 ``` 
   ### Example 2: Endopeptidase digestion with noise  
  ```R     
       
       Rscript runQPuB.r -fol examples/toy_endo3_nonoise -pf input.txt 
   
 ``` 
 ## Output of QPuB
 
  QPuB generates the following set of output files in the respective sub-folder **OUTPUT_example** inside the  [Output](https://github.com/QuantSysBio/QPuB/Output) folder
   
 | File | Description |
| ------ | ------ |
| **boxplot\_chain.pdf** | boxplot corresponding to the distributions of conversion factors |
| **chain\_Niter.pdf**  | trace plots of the Markov chain at **Niter**<sup>th</sup> iteration     |  
| **chain\_XYZ.RData** | tores trace of the Markov chain or the time series of the parameter draws   |
| **chain\_backscaled.RData** | tores the Markov chain after dividing the products by the scaling factor pre-calculated according to the substrate concentration    |
| **massdev.RData**  | stores mass deviation of the products over Monte Carlo iterations    |  
| **runQPuB\_ROUT.txt**  | contains the initial parameters that were feed to the algorithm and acceptance rates as the chain progresses. The output of any print command directly goes into the file.    |  
| **residuals\_M.txt** | t residual plot at **M**<sup>th</sup> iteration. The plot display mass deviation for individual amino acids of the substrate     |
| **conc\_means\_Nrep.csv, conc\_sd\_Nrep.csv** | stores the Markov chain after dividing the products by the scaling factor pre-calculated according to the substrate concentration   |
| **statistics.csv**  | summary statistic for the Markov chain     |
| **massdeviation.png**  | plot of total mass deviation of the products over Monte Carlo iterations  |
| **relation.png**  | plot of relation between the estimated conversion factors and the peptide lengths    |  
 
   
 ## Documentation
 
      A PDF version of the documentation will be soon added to the repository
 
 ## Developers
               
  * [**Sarah Henze**](https://www.mpibpc.mpg.de/person/59990/84522)
  
  * [**Debdas Paul**](https://www.mpibpc.mpg.de/person/97709/2169)  
  
  * [**Juliane Liepe**](https://www.mpibpc.mpg.de/person/52238/15851745)
 
 ## References 
    
   [1] Henze, S., Paul, D., Mansurkhodzhaev, A., Henklein, P., Textoris-Taube, K., Henning, U., Mishto,
M., and Leipe, J. (2019). Quantification of in vitro peptide hydrolysis and protease-catalyzed
peptide splicing using bayesian inference. Submitted.

   [2] Paul, D*., Henze, S*., Mishto, M., Liepe, J (2019). QPuB - Quantification of peptides using Bayesian inference. In preperation.  *Equal contributions
  
  ## LICENSE:
  
  This project is licensed under the GNU General Public License v3.0 - see the [LICENSE.md](LICENSE.md) file for details
  
<!-- Markdown link & img dfn's -->
[bayestat-image]: https://en.wikipedia.org/wiki/Markov_chain_Monte_Carlo#/media/File:Bayes_icon.svg
[peptide-image]: https://en.wikipedia.org/wiki/Peptide#/media/File:Tetrapeptide_structural_formulae_v.1.png
[bayestat-url]: https://en.wikipedia.org/wiki/Bayesian_inference
[peptide-url]: https://en.wikipedia.org/wiki/Peptide
[immunoprot-url]: https://en.wikipedia.org/wiki/Immunoproteomics
