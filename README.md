
<img src="logoqpub.png" width="200"> -- **adding Bayesian flavour in peptide quantification**

:black_medium_small_square: [Bayesian inference][bayestat-url] 
:black_medium_small_square:[Peptide quantification][peptide-url]

QPuB (**Q**uantifcation of **p**eptides **u**sing **B**ayesian inference) employs [Bayesian statistical inference](https://en.wikipedia.org/wiki/Bayesian_inference) based on [Markov chain Monte Carlo (MCMC)](https://en.wikipedia.org/wiki/Markov_chain_Monte_Carlo) sampling to learn the posterior distributions of the converison factors for peptide products. 

## Getting started

The following instructions will help you to download QPuB and execute the same to estimate the conversion factors for your peptide products. QPuB is designed to be invoked using the command line arguments as well from RStudio but with little modifcations in the code.

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
where the **package_name** should be replaced by the name of packages that are listed above. Alternatively, one can copy-paste the following command in the R terminal 

```R
install.packages(c("R.utils", "tictoc", "mcmcse", "mvtnorm", "tmvtnorm", "corpcor", "dqrng", "coda", "matrixcalc", "ggdmc", "sys", "matrixStats"))
 ```
and press enter.

REMEMBER, for Linux, the above command might need the root access.


## Running QPuB
   Once all the prerequisites are met, QPuB is needed to be downloaded from its Github repository. Click [here](https://github.com/QuantSysBio/QPuB/archive/master.zip) for to download the .zip file. Next, in order 
   to run QPuB on your data, follow the steps below:

**1.** Make sure everything is properly installed on your computer (see
Chapter 3 of the documentation).

**2.** Make sure the data files and the inputfolder have the right
structure (see Sections 2.3.1 and 2.3.2  of the [Documentation](Documentation.pdf) for more details)

**3.** Open the terminal/command prompt (How to:
[Linux](https://www.wikihow.com/Open-a-Terminal-Window-in-Ubuntu),
[Mac](https://www.wikihow.com/Open-a-Terminal-Window-in-Mac),
[Windows](https://www.wikihow.com/Open-Terminal-in-Windows)).

**4.** Navigate to the directory from where you want to run the QPuB, For example if the name of the directory is "ABC" then 
       the command takes the following generic form:
```sh 
   $ cd <path to the ABC folder>
```
For more details on how to navigate to a particular directory using the **cd** command, see [here](https://en.wikipedia.org/wiki/Cd_(command)). 

**5.** Once you are in your working directory (the directory from where you wish to run QPuB), say "ABC", the command to run the QPuB has the following generic form

```sh
  $ Rscript <path to runQPuB.r relative to the folder ABC> -infol <path to input_folder> -outfol <output_folder> -titr <titration_file>
```
**NOTE** The commands described above are providing a generic structure and must not be copy-pasted directly in the terminal.

For execution, the fields enclosed within **<>** including the symbol **<>** MUST be replaced by appropriate folder and file names along with their relative paths (if required) as explained in the **Examples** below. For output folder you should not provide the path, but only the name of the folder. The same applies for the titration data file. Note that For now, ./ and ../ notation does not work. Please use the full path. Also avoid leading and trailing
slashes.

The table below provides the meaning of the flags used in the above shell command.

|         Flags                 |                 Meaning                                                              |
| :----------------------- | :---------------------------------------------------------------------------- |
| <span>**-infol**</span>  | flag for the <span>**input_folder**</span>. This input is mandatory.      |
| <span>**-outfol**</span> | flag for the <span>**output_folder**</span>. This input is optional, and the name can be defined by the user            |
| <span>**-titr**</span>   | flag for the <span>**titration_file**</span>. It is a csv-file. This input is optional. |

Flags can also be specified using `--` as prefix. For example, `--infol=inputfolder`. If you do not wish to provide the output folder and/or titration data, please do not use the respective flags. Please follow the examples below for illustrations.


## Examples

   The **examples** folder contains two toy examples of endopeptidase digestion, one exopeptidase digestion, and a gp100<sub>40-52</sub> digestion by 26S proteasomes. In order to run these examples, execute the following commands in the terminal **assuming your current directory is the QPuB-master**
   
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
       
        Rscript ./QPuB/runQPuB.r -infol examples/toy exo
   
 ``` 
   ### Example 4: gp100<sub>40-52</sub> digestion by 26S proteasomes
  ```sh     
       
        Rscript ./QPuB/runQPuB.r -infol examples/Prot_K386
   
 ``` 
 ## Running QPuB from [RStudio](https://www.rstudio.com/)
 
One can execute QPuB from RStudio, but in that case one first needs to open the main Rscript **runQPuB.r** and then disable/commenting out the following code snippet in that scripts 

```R
args = commandArgs(trailingOnly=FALSE, asValue=TRUE)
keys <- attachLocally(args)

if (!exists('infol')){
stop("This is the QPuB package. \n
At least one argument must be supplied: -infol.\n 
-infol inputfolder (COMPULSORY): folder containing all of the following input files \n AND folder 'data' with at least one csv file with input data. \n 
-outfol outputfolder (OPTIONAL): if not provided, default name will be 'OUTPUT_inputfolder(i)'. \n
-titr titrationfile (OPTIONAL): csv file with titration data \n
For documentation see: https://github.com/QuantSysBio/QPuB.\n", call.=FALSE)
}
```
and enable/uncommenting the following snippet and editing the variables according to the desired example. For **Prot_K386** with an output folder **my_fancy_name**, the snippets looks like the following:

```R
# file <- 'runQPuB.r'
# infol <- 'examples/Prot_K386'
# outfol <- 'my_fancy_name'
# titr <- '190423_K386_titration_substrate_charge_3_K386.csv'
```
However, in order to avoid complications espcially crashing of RStudio for large datasets, users are strongly recommended to use the command line as described in the following sections.

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

Output generated while running the examples are illustrated in Chapter 5 of the [Documentation](Documentation.pdf).
    
## Troubleshooting

In this chapter, we listed a non-exhaustive list of issues that the user
might get into while running the QPuB. However, users are <span>strongly
recommended</span> to keep track of the **runQPuB\_ROUT.txt** from
<span>the very beginning of the run</span>, as all the outputs including
the warnings that might be generated for providing an out-of-range value
for a particular input parameter, are printed.

### Multiple installation of R

The user must make sure that the version of the R is 3.5.0 or above. In
case, the user has multiple installation of R, the path to the R
executable corresponding to the version specified for the QPuB should be
provided. For Windows, please click
[here](https://cran.r-project.org/bin/windows/base/rw-FAQ.html#Rcmd-is-not-found-in-my-PATH_0021)
for more details on how to include the desired R executable in the PATH
variable.

### Command line execution

Care must be taken while writing the command line to run QPuB. For
example, the following type of command will generate error:  

``` 

$ Rscript path_to_runQPuB.r -infol some_path/../examples/example_folder
```

Instead, the <span>right</span> command type is

    $ Rscript path_to_runQPuB.r -infol examples/example_folder

The same rules applies while including the titration file, for example,
assuming the name of the example folder as **P** and the corresponding
titration file name is **titration\_P**, the following command

    $ Rscript path_to_runQPuB.r -infol examples/P -titr examples/P/titration_P.csv

will generate an error, while the command

    $ Rscript path_to_runQPuB.r -infol examples/P -titr titration_P.csv

will work nicely.

### MCMC chains

It may happen that the Markov chain is not exploring the parameter
space. There might be several reasons for that. For example, the value
of `sigma_start` is so low that the chain stop exploring the parameter
space. At the same time if the value of the `GammaExponent` is higher,
due to vanishing nature of the adaptation, the covariance matrix will
remain stationary leaving the chain stuck at an undesired region.
Therefore, one solution to this is to stop the adaptation at the
beginning for sufficient exploration and then resume the adaptation
again after visually inspecting the chain. This may take several trial
runs to know the number of iterations needed before starting the
adaptation. Remember that this problem is different than the one where
the chain needs to adapt at the beginning to reach the target
covariance matrix. But, it may take over millions over millions of iterations before reaching the target region. This is
especially the case for higher dimensions . Overall, depending upon the
problem, users are recommended to experiment with the parameters
controlling the adaptation scheme and the Markov chain. In [3], the authors
provide an excellent summary of various pathological cases regarding the
adaptive scheme and optimal scaling for MCMC.

 ## Documentation
 
   A PDF version of explaining the details of the QPuB is available as a [Documentation](Documentation.pdf)
  
 ## Developers
               
  * [**Sarah Henze**](https://www.mpibpc.mpg.de/person/59990/84522)
  
  * [**Debdas Paul**](https://www.mpibpc.mpg.de/person/97709/2169)  
  
  * [**Juliane Liepe**](https://www.mpibpc.mpg.de/person/52238/15851745)
 
 ## References related to QPuB 
    
   [1] Henze, S., Paul, D., Mansurkhodzhaev, A., Henklein, P., Textoris-Taube, K., Henning, U., Mishto,
M., and Leipe, J. (2019). Quantification of in vitro peptide hydrolysis and protease-catalyzed
peptide splicing using bayesian inference. Under review.

   [2] Henze, S*., Paul, D*., Mishto, M., Liepe, J (2019). QPuB - Quantification of peptides using Bayesian inference. In preperation.  *Equal contributions
   
## Other references

   [3] Rosenthal, J. S. (2011). Optimal proposal distributions and adaptive MCMC. Handbook of Markov Chain Monte Carlo, 4(10.1201).
  
  ## LICENSE
  
  This project is licensed under the GNU General Public License v3.0 - see the [LICENSE.md](LICENSE.md) file for details
  
<!-- Markdown link & img dfn's -->
[bayestat-image]: https://en.wikipedia.org/wiki/Markov_chain_Monte_Carlo#/media/File:Bayes_icon.svg
[peptide-image]: https://en.wikipedia.org/wiki/Peptide#/media/File:Tetrapeptide_structural_formulae_v.1.png
[bayestat-url]: https://en.wikipedia.org/wiki/Bayesian_inference
[peptide-url]: https://en.wikipedia.org/wiki/Peptide
[immunoprot-url]: https://en.wikipedia.org/wiki/Immunoproteomics
