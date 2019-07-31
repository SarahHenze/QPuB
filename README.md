
<img src="logoqpub.png" width="200">

:black_medium_small_square: [Bayesian inference][bayestat-url] 
:black_medium_small_square:[Peptide quantification][peptide-url]

QPuB (**Q**uantifcation of **p**eptides **u**sing **B**ayesian inference) employs [Bayesian statistical inference](https://en.wikipedia.org/wiki/Bayesian_inference) based on [Markov chain Monte Carlo (MCMC)](https://en.wikipedia.org/wiki/Markov_chain_Monte_Carlo) sampling to learn the posterior distributions of the converison factors for peptide products. 

## Getting started

The following instructions will help you to download QPuB and execute the same to estimate the conversion factors for your peptide products

### Prerequisites

In order to run QPuB, users MUST install [**R** ≥ 3.5.0](https://www.r-project.org/) and the following packages:

* [**R.utils**](https://cran.r-project.org/web/packages/R.utils/index.html)
   
* [**tictoc**](https://cran.r-project.org/web/packages/tictoc/index.html) 

* [**sys**](https://cran.r-project.org/web/packages/sys/index.html) 
   
* [**mcmcse**](https://cran.r-project.org/web/packages/mcmcse/index.html) 

* [**mvtnorm**](https://cran.r-project.org/web/packages/mvtnorm/index.html) 
   
* [**tmvtnorm**](https://cran.r-project.org/web/packages/tmvtnorm/index.html)  

* [**corpcor**](https://cran.r-project.org/web/packages/corpcor/index.html) 
   
* [**dqrng**](https://www.rdocumentation.org/packages/dqrng/index.html) 
   
* [**coda**](https://cran.r-project.org/web/packages/coda/index.html) 

* [**matrixcalc**](https://www.rdocumentation.org/packages/matrixcalc) 

* [**ggdmc**](https://www.rdocumentation.org/packages/ggdmc)

* [**matrixStats**](https://www.rdocumentation.org/packages/matrixStats) 

To install the above packages one by one, start your **R** shell, and execute the following command:

```R
 install.packages("package_name")
 ```
where the "package_name" should be replaced by the name of packages that are listed above. Alternatively, one can copy-paste the following command in the R terminal 

```R
install.packages(c("R.utils", "tictoc", "mcmcse", "mvtnorm", "tmvtnorm", "corpcor", "dqrng", "coda", "matrixcalc", "ggdmc", "sys", "matrixStats"))
 ```
and press enter.

REMEMBER, for Linux, the above command might need the root access.


## Running QPuB
   Once all the prerequisites are met, QPuB is needed to be downloaded from its Github repository. Next,
   to use QPuB on your data, follow these steps:

**1.** Make sure everything is properly installed on your computer (see
Chapter 3 of the documentation).

**2.** Make sure the data files and the inputfolder have the right
structure (see Sections 2.3.1 and 2.3.2  of the [Documentation](Documentation.pdf) for more details)

**3.** Open the terminal/command prompt (How to:
[Linux](https://www.wikihow.com/Open-a-Terminal-Window-in-Ubuntu),
[Mac](https://www.wikihow.com/Open-a-Terminal-Window-in-Mac),
[Windows](https://www.wikihow.com/Open-Terminal-in-Windows)).

**4.** Go to your working directory:

```sh 
   $ cd hostname:˜/workingdirectory
```

**5.** The QPuB main script is executed using flags:

```sh
  $ Rscript <name or path to runQPuB.r> -infol -outfol -titr
```

QPuB is equipped with facilities to accept and parse command line
arguments, where

|                          |                                                                               |
| :----------------------- | :---------------------------------------------------------------------------- |
| <span>**name or path to runQPuB.r**</span>   | name/path of the <span>**runQPuB.r**</span>. This input is mandatory.         |
| <span>**-infol**</span>  | name/path of the <span>**input folder**</span>. This input is mandatory.      |
| <span>**-outfol**</span> | userdefined name of the **output folder**. This input is optional.            |
| <span>**-titr**</span>   | name of the <span>**titration data**</span> csv-file. This input is optional. |

Flags can be specified as `-infol inputfolder` or `--infol=inputfolder`.

Depending on your working directory, you also have to provide the path
to the runQPuB.r file:

1.  working directory is QPuB directory: `runQPuB.r`

2.  working directory is somewhere else: `path\QPuB\runQPuB.r`

The first argument will always be taken as the name of the R script:
`Rscript runQPuB.r`.

Depending on the location of the input folder, you also have to provide
the path to the folder:

1.  inputfolder is in same parent directory as QPuB: `-infol
    inputfolder`

2.  inputfolder is somewhere else: `-infol path\inputfolder`
## Examples

   The **examples** folder contains two toy examples of endopeptidase digestion. 
   In order to run the examples, execute the following commands in the terminal assuming your current directory is the QPuB-master
   
   ### Example 1: Endopeptidase digestion without noise
   ```sh     
       
       Rscript ./QPuB/runQPuB.r -infol examples/toy_nonoise
   
 ``` 
   ### Example 2: Endopeptidase digestion with noise  
  ```sh     
       
        Rscript ./QPuB/runQPuB.r -infol examples/toy_noise
   
 ``` 
   ### Example 3: gp100<sub>40-52</sub> digestion by 26S proteasomes
  ```sh     
       
        Rscript ./QPuB/runQPuB.r -infol examples/Prot_K386
   
 ``` 
 ## Output of QPuB
 
  QPuB generates the following set of output files
   
 | File | Description |
| ------ | ------ |
| <span>**boxplot\_chain.pdf**</span>                                                  | boxplot corresponding to the distributions of conversion factors                                     |
| <span>**chain.RData**</span>                                                         | Markov chain: time series of all parameters                                                          |
| <span>**chain\_backscaled.RData**</span>                                             | Markov chain after backscaling (see Section 2.3.5 of the [Documentation](Documentation.pdf)                                  |
| <span>**conc\_five\_rp.csv, conc\_median\_rp.csv, conc\_ninetyfive\_rp.csv** </span> | 0.05, 0.5 and 0.95 quantiles of concentration kinetics of all products for the rp<sup>th</sup> replicate |
| <span>**runQPuB\_ROUT.txt**</span>                                                   | Progress report of the algorithm. The output of every print command directly goes into this file.    |
| <span>**chain\_j.pdf**</span>                                                        | Trace plots and distributions of the Markov chain at j<sup>th</sup> iteration                            |
| <span>**residuals\_j.txt**</span>                                                    | Residual plots at j<sup>th</sup> iteration                                                               |
| <span>**statistics.csv**</span>                                                      | Summary statistics of the Markov chain                                                               |


   
 ## Documentation
 
   A PDF version of the documentation is available at [Documentation](Documentation.pdf)
 
 ## Developers
               
  * [**Sarah Henze**](https://www.mpibpc.mpg.de/person/59990/84522)
  
  * [**Debdas Paul**](https://www.mpibpc.mpg.de/person/97709/2169)  
  
  * [**Juliane Liepe**](https://www.mpibpc.mpg.de/person/52238/15851745)
 
 ## References 
    
   [1] Henze, S., Paul, D., Mansurkhodzhaev, A., Henklein, P., Textoris-Taube, K., Henning, U., Mishto,
M., and Leipe, J. (2019). Quantification of in vitro peptide hydrolysis and protease-catalyzed
peptide splicing using bayesian inference. Submitted.

   [2] Henze, S*., Paul, D*., Mishto, M., Liepe, J (2019). QPuB - Quantification of peptides using Bayesian inference. In preperation.  *Equal contributions
  
  ## LICENSE
  
  This project is licensed under the GNU General Public License v3.0 - see the [LICENSE.md](LICENSE.md) file for details
  
<!-- Markdown link & img dfn's -->
[bayestat-image]: https://en.wikipedia.org/wiki/Markov_chain_Monte_Carlo#/media/File:Bayes_icon.svg
[peptide-image]: https://en.wikipedia.org/wiki/Peptide#/media/File:Tetrapeptide_structural_formulae_v.1.png
[bayestat-url]: https://en.wikipedia.org/wiki/Bayesian_inference
[peptide-url]: https://en.wikipedia.org/wiki/Peptide
[immunoprot-url]: https://en.wikipedia.org/wiki/Immunoproteomics
