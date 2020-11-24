# EnsembleMHC

## Installation

### MHC-I prediction algorithms and setup

In order to run EnsembleMHC you will need to download the following MHC-I prediction algorithms

*  [MixMHCpred-2.0.2](https://github.com/GfellerLab/MixMHCpred/releases/tag/v2.0.2)
*  [netMHC-4.0](https://services.healthtech.dtu.dk/services/NetMHC-4.0/9-Downloads.php#)
*  [netMHCpan-4.0](https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/9-Downloads.php#)
*  [netMHCstabpan-1.0](https://services.healthtech.dtu.dk/services/NetMHCstabpan-1.0/9-Downloads.php#)
*  [pickpocket-1.1](https://services.healthtech.dtu.dk/services/PickPocket-1.1/9-Downloads.php#)
*  [mhcflurry-1.6.0](https://github.com/openvax/mhcflurry/releases/tag/1.6.0) (note: We recommend that you use the [install_mhcflurry.sh](scripts/install_mhcflurry.sh) script located in the script folder)

Be sure to follow the installation instructions of each algorithmm carefully. These algorithms can be installed at any location, however, be sure to update the [ALGORITHMS_PATHS.sh](ALGORITHMS_PATHS.sh) script with the absolute path to each algorithm. 


### Other dependancies 

* [Anaconda](https://docs.anaconda.com/anaconda/install/) 

** homebrew install of anaconda
	```bash
	brew cask install anaconda
	```

* [R](https://www.r-project.org/)

* gnu-sed (only appplicable if using macOS). This can be accomplised using homebrew
  
  Download homebrew
	```bash
	/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh)"
	```
  Download gnu-sed
	```bash
	brew install gnu-sed
	```

### setup 



## usage 


