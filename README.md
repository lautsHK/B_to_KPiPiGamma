# B\_to\_KPiPiGamma

## Introduction

This is a Monte-Carlo generator to generate $B \to K_{res}\gamma \to K\pi\pi\gamma$ process.
There are totally 5 $K_{res}$ included, which are :
- K(1270), namely K\_10 (10313), K\_1+ (10323) 
- K(1400), namely K'\_10 (20313), K'\_1+ (20323)
- K(1410), namely K'\*0 (100313), K'\*+ (100323)
- K(1430), namely K\_2\*0 (315), K\_2\*+ (325)
- K(1680), namely K''\*0 (30313), K''\*+ (30323)

Interference between different kaonic resonances are included. 
While for the polarisation of photon, it is controlled by the beyond standard model (BSM) dimension-6 (D6) couplings. A single parameter $\lambda$, between 1 and -1, is used here. 
When $\lambda\to 1$, the generated events are mainly right-handed polarised. While if $\lambda\to -1$, the generated events are mainly left-handed polarised.

## Installation
This package needs C++11 and cmake3. Please make sure your working environment support the mentioned version.
Once the envirnment is ready, please do :
'rm -rf build' (to clean)
'mkdir build' (to build the "build" directory)
'cd build'
'cmake3 ./..' 
'make'

The package should be compiled sucessfully.
