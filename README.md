# 3PG-SoNWaL_FR Calibration

This repository contains R scripts for calibrating the 3PG-SoNWaL_FR model using advanced statistical methods implemented through the `dgpsi` and `tinydancer` packages.

## Prerequisites

- R (version 4.0 or higher recommended)
- Access to the 3PG-SoNWaL_FR repository
- Internet connection for package installation

## Setup Instructions

### Step 1: Gain Access to 3PG-SoNWaL_FR Repository

To access the required 3PG-SoNWaL_FR repository, please contact **ForestResearch** for permission and access credentials.

### Step 2: Install Required Packages

#### Install dgpsi Package

The `dgpsi` package provides Deep Gaussian Process emulation functionality. Install it following the instructions at:

**https://github.com/mingdeyu/dgpsi-R**

#### Install tinydancer Package

The `tinydancer` package is required for the calibration process. Install it following the instructions at:

**https://github.com/BayesExeter/tinydancer**

### Step 3: Run Calibration

Once you have access to the 3PG-SoNWaL_FR repository and have installed the required packages, run the calibration script:

```r
source("3pg_Calibration_Tinydancer.R")
```

## Project Structure

```
.
├── 3pg_Calibration_Tinydancer.R    # Main calibration script
├── README.md                       # This file
└── [additional project files]
```

## Usage

1. Ensure you have obtained access to the 3PG-SoNWaL_FR repository from ForestResearch
2. Install the required R packages (`dgpsi` and `tinydancer`)
3. Load the calibration script in your R environment
4. Execute the calibration process

## Dependencies

- **dgpsi**: Deep Gaussian Process emulation
- **tinydancer**: Bayesian calibration framework
- **3PG-SoNWaL_FR**: Forest growth model (access required)

## Support

For issues related to:
- **3PG-SoNWaL_FR access**: Contact ForestResearch
- **dgpsi package**: Visit https://github.com/mingdeyu/dgpsi-R
- **tinydancer package**: Visit https://github.com/BayesExeter/tinydancer

## Contributing

[Add contribution guidelines if applicable]

## License

[Add license information]

## Citation

[Add citation information if applicable]

---

*This README.md file is ready to be downloaded and added to your GitHub repository.*
