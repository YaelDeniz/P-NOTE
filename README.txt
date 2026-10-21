WIP:
If PNOTE ask for ROOT_INCLUDE/MATH/interpolaor.h directory

$ find / -type f -name "Interpolator.h" 2>/dev/null
.../root_src/math/mathmore/inc/Math/Interpolator.h
../root/root/include/Math/Interpolator.h

the copy it $ cp .../root_src/math/mathmore/inc/Math/Interpolator.h ~/root-install/include/Math/

Enable Eigen from Osbrob and consruct root with MathMore

ydenizhernandez@COS-D617915:~$ cp /home/ydenizhernandez/root_src/math/mathmore/inc/Math/Interpolator.h ~/root-install/include/Math/



# Deniz et al. 2026 — Neutrino Oscillation Tomography of Earth (NOTE) code

This repository contains the CERN ROOT(c++) applications for the simulation of atmospheric neutrinos through the Earth and calculation of event distributions used in *Deniz et al. (2026)*. The code includes applications for generating   true and observed(Reconstructed) event distributions, chi2 calculations, neutrino flux and three-dimensional Earth-model visualizations.

## 1. Requirements

Before compiling the code, the paths to the required external libraries must be specified in the `Makefile`.

### CERN ROOT

This code requires CERN ROOT (follow the tutorial at https://root.cern/install/ to install it) data analysis framework. Installation flags for TGeoManger and OpenGL are required to run the program

In the `Makefile`, update according to /path/to/root:

```makefile
ROOT_LIB=/path/to/root/lib
ROOTMATH_INCLUDE=/path/to/root/include/Math
```

### OscProb

To calculate matter oscillation probabilities this framework implements calculations use the OscProb software (https://github.com/joaoabcoelho/OscProb).

Update the following line in the `Makefile`:

```makefile
OSCPROB_PATH=/path/to/OscProb
```

### Code Headers

specify the header location for the neutrino tomography simulation framework:

```makefile
NUTOMO_INCLUDE=/path/to/Neutrino-Tomography-Simulation/inc
```

Headers are located in the inc/ directory 

---

## 2. Compiling an application

After updating the paths in the `Makefile`, compile all applications with:

```
$ make all
```

This builds the following applications:

```
ObsN
ObsSenv
IntN
IntSenv
HondaFlux
Earth3D
```

To remove the compiled object files and applications:

```
$ make clean
```

Individual applications can also be compiled separately. For example:

```
$ make IntN
```

---

# 3. Applications

## 3.1 Interacting event distribution (main/IntOscGram.cpp)
```
$ ./IntN
```
Generates expected *interacting neutrino event distribution* for a range in energy and direction and bin numbers for a given Earth model. The output is a CSV file of the form E^i,theta_z^j,phi^k,Nijk (No headers)


## 3.2 Reconstructed event distribution (main/ObsOscGram.cpp)
```
$ ./ObsN
```
Generates expected *reconstructed neutrino event distribution* for a range in the reconstructed variables and bin numbers for a given Earth model. The output is a CSV file of the form E^l,theta_z^m,phi^n,Nlmn (No headers)
---

## 3.3 True Sensitivity Simulation(main/TrueSensitivity.cpp)
```
$ ./IntSenv
```
Generates the interacting sensitivity simulations by comparing the  *interacting event distribution* calculated for the reference and perturbed model, for a range of electron number density anomalies. It implements equation 15 and 17 in Deniz et al 2026. The output is a CSV file of the form \Delta n_e[%], chi2e;int , chi2mu;int, #ibins, #jbins, #kbins  (No headers)
---

## 3.4 Reconstructed Sensitivity Simulation(main/ObsSensitivity.cpp)
```
$ ./ObsSenv
```
Generates the reconstructed sensitivity simulations by comparing the  *reconstructed event distribution* calculated for the reference and perturbed model, for a range of electron number density anomalies. It implements equation 15 and 17 in Deniz et al 2026. The output is a CSV file of the form \Delta n_e[%], chi2e;reco , chi2mu;reco, #lbins, #mbins, #nbins  (No headers)


## 3.5  Atmospheric Neutrino Flux (main/HondaFlux.cpp)

```
$./HondaFlux
```

Generates visualizations of the atmospheric neutrino flux used in the simulations.

---

## 3.6 Three-Dimensional Earth Model Visualization (main/Earth3DVis.cpp)

```bash
$./Earth3D
```
"Earth3D" visualizes the three-dimensional Earth model generated with CERN ROOT TGeoManager. It implements PREM tables directly from OscProb to generate the geometry and materail of each layer



# 4. Quick Start

After cloning the repository:

### Step 1 — Configure dependencies

Edit the following variables in the `Makefile`:

```make
ROOT_LIB=/path/to/root/lib
ROOTMATH_INCLUDE=/path/to/root/include/Math
OSCPROB_PATH=/path/to/OscProb
NUTOMO_INCLUDE=/path/to/Neutrino-Tomography-Simulation/inc
```

### Step 2 — Compile

```
$ make all
```

### Step 3 — Run an application

For example:

```bash
./IntN
```

or:

```bash
./Earth3D
```

### Step 4 — Clean compiled files

```bash
make clean
```

---


