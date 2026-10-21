# P-NOTE

<!-- Short description of the project goes here -->

## Requirements

P-NOTE has been developed and tested with:

| Component | Version |
|-----------|---------|
| OS | Ubuntu 24.04 (native or WSL2), x86_64 |
| Compiler | GCC 13.3 |
| C++ standard | C++17 |
| CMake | ≥ 3.20 (tested with 3.28) |
| CERN ROOT | **6.36.00**, built with `mathmore` enabled |
| GSL | 2.7.1 |
| Python | 3.12 |

> **Important:** ROOT must be built with **MathMore** enabled, which requires GSL
> to be installed *before* building ROOT. If GSL is missing, ROOT silently
> disables MathMore and P-NOTE will fail to build.

---

## 1. Install system dependencies

```bash
sudo apt update
sudo apt install -y \
  git cmake g++ gcc binutils dpkg-dev make \
  python3 python3-dev python3-numpy \
  libx11-dev libxpm-dev libxft-dev libxext-dev \
  libssl-dev libpcre3-dev libtbb-dev libvdt-dev \
  libgsl-dev \
  libglu1-mesa-dev libglew-dev libftgl-dev libgl2ps-dev \
  libxml2-dev libsqlite3-dev libcfitsio-dev libcurl4-openssl-dev \
  libfreetype-dev libfontconfig-dev \
  libjpeg-dev libpng-dev libgif-dev libtiff-dev \
  zlib1g-dev liblzma-dev libzstd-dev liblz4-dev libxxhash-dev \
  nlohmann-json3-dev
```

Check that GSL is available:

```bash
gsl-config --version    # should print 2.7.1 or newer
```

---

## 2. Build CERN ROOT 6.36.00 from source

The build uses three folders in your home directory:
`root_src` (source), `root_build` (temporary build files), and
`root-install` (final installation).

### 2.1 Download the source

```bash
cd ~
git clone --branch v6-36-00 --depth=1 https://github.com/root-project/root.git root_src
```

### 2.2 Configure

```bash
mkdir -p ~/root_build && cd ~/root_build

cmake ../root_src \
  -DCMAKE_INSTALL_PREFIX=$HOME/root-install \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_STANDARD=17 \
  -Dfail-on-missing=ON \
  -Dmathmore=ON \
  -Dgdml=ON \
  -Dopengl=ON \
  -Dx11=ON \
  -Dasimage=ON \
  -Dimt=ON \
  -Dpyroot=ON \
  -Droofit=ON \
  -Dtmva=ON \
  -Dfitsio=ON \
  -Dsqlite=ON \
  -Dssl=ON \
  -Dxml=ON \
  -Dvdt=ON
```

Before continuing, check the line starting with `-- Enabled support for:`
near the end of the output. It **must** include `mathmore`.

### 2.3 Build and install

```bash
cmake --build . --target install -j$(nproc)
```

This can take 1–2 hours. If the build is killed due to low memory
(common on WSL or laptops), re-run it with fewer jobs, e.g. `-j4`.
It resumes where it stopped.

### 2.4 Set up the environment

```bash
echo 'source $HOME/root-install/bin/thisroot.sh' >> ~/.bashrc
source ~/.bashrc
```

### 2.5 Verify the installation

```bash
root-config --version          # 6.36.00
root-config --has-mathmore     # yes
root-config --features         # should include: mathmore geom gdml opengl
```

---

## 3. Build P-NOTE

```bash
git clone <REPOSITORY-URL> P-NOTE
cd P-NOTE
make
```

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


