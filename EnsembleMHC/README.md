# EnsembleMHC

## Installation

### MHC-I prediction algorithms and setup

In order to run EnsembleMHC, you will need to download the following MHC-I prediction algorithms:

*  [MixMHCpred-2.0.2](https://github.com/GfellerLab/MixMHCpred/releases/tag/v2.0.2)
*  [netMHC-4.0](https://services.healthtech.dtu.dk/services/NetMHC-4.0/9-Downloads.php#)
*  [netMHCpan-4.0](https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/9-Downloads.php#)
*  [netMHCstabpan-1.0](https://services.healthtech.dtu.dk/services/NetMHCstabpan-1.0/9-Downloads.php#)
*  [pickpocket-1.1](https://services.healthtech.dtu.dk/services/PickPocket-1.1/9-Downloads.php#)
*  [mhcflurry-1.6.0](https://github.com/openvax/mhcflurry/releases/tag/1.6.0) (note: We recommend that you use the [install_mhcflurry.sh](scripts/install_mhcflurry.sh) script located in the script folder)

Be sure to follow the installation instructions of each algorithm carefully. These algorithms can be installed at any location, however, be sure to update the [ALGORITHMS_PATHS.sh](ALGORITHMS_PATHS.sh) script with the absolute path to each algorithm. 


### Other dependencies 

* [Anaconda](https://docs.anaconda.com/anaconda/install/) 

** homebrew install of anaconda
	```
	brew cask install anaconda
	```

* [R](https://www.r-project.org/)

* gnu-sed (only applicable if using macOS). This can be accomplished using homebrew
  
  Download homebrew
	```	
	/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh)"
	```
  Download gnu-sed
	```
	brew install gnu-sed
	```

  Download gnu-parallel
	```
	brew install parallel
	```

### permissions

In order to use EnsembleMHC as presented in the example, give executable permission to EnsembleMHC.

``` bash
chmod +x EnsembleMHC
```



## Usage

To see a the list of flags available to EnsembleMHC, enter `EnsembleMHC -h`

* `-p` Specify target protein/s for predictions. Proteins must be in FASTA format.
	
* `-a` Specify target HLA. A list of the supported HLAs can be seen [here.](scripts/HLA_list.txt)
	
* `-t` Specify the number of threads to accessible to EnsembleMHC. If left unassigned, the number of threads will default to the number of available cores.
	
* `-o` provide output name for the prediction folder.
    
* `-m` The parameterization summary matrix. This can be left blank when using the default algorithm and allele thresholds. The user has the option to provide a custom parameterization summary matrix. (see allele parameterization)

* `-d` add this flag to turn on debug mode. This will generate more verbose outputs.

## test 

You can test your installation of EnsembleMHC by navigating to the test directory and entering the following commands.

``` bash
cd test/

../EnsembleMHC -a HLA-A02:01 -p test.prot -o EnsembleMHC_test 

```

To check if the expected output was generated, you can enter the following command

``` bash
diff EnsembleMHC_test/HLA-A02\:01_peptideFDR_pred.csv EnsembleMHC_predict.compare | wc -l
```
numbers > 0 indicate a potential error in the install.


