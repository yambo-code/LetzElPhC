![screenshot](docs/logo.png)

LetzElPhC abbreviates to _"Lëtzebuerg Electron Phonon Code"_.
_"Lëtzebuerg"_ is the Luxembourgish name for the Luxembourg Country.

LetzElPhC is distributed under the MIT license.

Please refer to the (currently evolving) [LetzElPhC documentation](https://yambo-code.github.io/LetzElPhC/) for usage.

**This code is currently under development. please use at your 
own risk**.

In case of any bugs or other issues, please open a issue.

Any contribtions to the code are always welcomed.

**Contributors note:** This project is released under the MIT License. To preserve that licensing, please avoid incorporating code from GPL-licensed software, including small snippets or translations of GPL code. Algorithmic ideas and published descriptions are fine, but implementations should be written independently or derived from permissively licensed sources (MIT, BSD, Apache-2.0, etc.) with appropriate attribution.


It should be noted that this code computes only electron–phonon matrix elements and does not evaluate any additional physical quantities. If you wish to calculate physical properties using Wannier-based interpolation techniques, you should instead use specialized tools such as [EPW](https://epw-code.org/), [Perturbo](https://perturbo-code.github.io/), [Epiq](https://the-epiq-team.gitlab.io/epiq-site/), or [Phoebe](https://phoebe-team.github.io/phoebe/).

<!-- # Citation -->
<!-- In case you wish to acknowledge the use of this code, please cite the following reference. -->
<!-- ``` -->
<!-- @phdthesis{Nalabothula2025May, -->
<!-- 	author = {Nalabothula, Muralidhar}, -->
<!-- 	title = {{Symmetries of Excitons: Implications for Exciton-Phonon Coupling and Optical Spectroscopy}}, -->
<!-- 	year = {2025}, -->
<!-- 	month = may, -->
<!-- 	address = {Luxembourg}, -->
<!-- 	school = {Unilu - Universit{\ifmmode\acute{e}\else\'{e}\fi} du Luxembourg [The Faculty of Science, Technology and Medicine]}, -->
<!-- 	url = {https://orbilu.uni.lu/handle/10993/65236} -->
<!-- } -->
<!-- ``` -->

# Authors 
Muralidhar Nalabothula\
Prof. Ludger Wirtz (Supervisor)

# Acknowledgments
Riccardo Reho (Testing and Documentation)\
Fulvio Paleari (Testing and Yambopy interface)\
University of Luxembourg (Funding)\
HPC @ Uni.lu (Computing resources)\
Henry Fried (Logo)

# TODO  
1) ~~Support XML format for dynamical matrices~~
2) ~~Work on suporting image parallelization of ph.x (preprocessor)~~
3) ~~Prepare automatic test cases~~
4) Improve OPENMP
5) Implement basic Acoustic sum rule
6) ~~Turn on different kernel options for the user~~
7) Frohlich Interpolation
8) DFT + U



## FFT Backend Licensing

This package is licensed under the MIT License. However, FFT computations are performed through a backend library, and the license terms applicable to that backend may impose additional obligations.

* **FFTW backend**: If the FFTW3 library is used, the library itself is licensed under GPLv2 or later. Programs that are distributed together with, or otherwise incorporate, the FFTW library may therefore be subject to the GPL. Merely distributing source code that depends on this package, while requiring users to obtain FFTW separately, may avoid these obligations, but users should consult the FFTW license.

* **Intel MKL backend**: If Intel oneAPI Math Kernel Library (MKL) is used via its FFTW3 compatibility interface, use of the backend is subject to Intel's licensing terms rather than the FFTW license. Users distributing software linked against MKL are responsible for complying with the applicable Intel license.

In all cases, this package itself remains licensed under MIT. Users are responsible for ensuring compliance with the license terms of the FFT backend they choose.
