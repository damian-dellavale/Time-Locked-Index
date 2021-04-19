# Time-Locked-Index
Time Locked Index (TLI) is a measure to quantify the harmonicity of noisy time series. It was developed as part of our project aimed to improve the characterization and interpretation of Cross Frequency Couplings (CFC) observed in neural recordings and multimodal time series in general.

In Velarde et al (2019), we analytically demonstrate that phase-amplitude coupling (PAC) phenomenon naturally emerges in mean-field models of biologically plausible networks, as a signature of specific bifurcation structures. The proposed analysis, based on bifurcation theory and the quantification of harmonicity using the TLI metric, allowed us to identify the mechanisms underlying oscillatory dynamics that are essentially different in the context of PAC.

In Dellavale et al (2020), the TLI allowed us to efficiently distinguish between harmonic and non harmonic phase-amplitude couplings (PAC) which coexist during the seizure dynamics of epilepsy patients observed with intracerebral iEEG recordings with macroelectrodes. Reliable identification of non harmonic PAC patterns involving High-gamma and higher frequency bands is of clinical relevance since these non harmonic PAC patterns have been previously associated with the ictal core through the paroxysmal depolarizing shifts mechanism of seizure propagation.

In Dellavale and Velarde et al (2020), we used the TLI measure to characterize different CFC types (phase-amplitude, amplitude-amplitude, phase-frequency, etc). Importantly, substantial evidence was presented supporting the conclusion that the concepts of "true" and "spurious" CFC are not intrinsically linked to low and high harmonic content in the oscillatory dynamics, respectivelly. That is, a single oscillatory dynamics characterized by a non constant oscillation period can produce "spurious" CFC with low harmonic content (i.e. non harmonic CFC). On the other hand, we found that two coupled oscillatory dynamics with independent fundamental frequencies can elicit "true" CFC with high harmonic content via the rectification mechanism (i.e. harmonic CFC).

# Description
Besides the TLI function `/matlab/function_TimeLockedIndex_v0.m`, you will find a main test script `/matlab/test_TLI.m` and all the auxiliary functions required to compute the TLI for the dynamics of the Van der Pol oscillator.

# References

- Velarde O, Urdapilleta E, Mato G, and Dellavale D (2019), Bifurcation structure determines different phase-amplitude coupling patterns in the activity of biologically plausible neural networks, NeuroImage, Vol. 202, No. 116031, DOI: [10.1016/j.neuroimage.2019.116031](https://doi.org/10.1016/j.neuroimage.2019.116031).\
Artwork selected to appear on the cover: [https://www.sciencedirect.com/journal/neuroimage/vol/202/](https://www.sciencedirect.com/journal/neuroimage/vol/202/).

- Dellavale D, Urdapilleta E, Cámpora N, Velarde O, Kochen S, Mato G (2020), Two types of ictal phase-amplitude couplings in epilepsy patients revealed by spectral harmonicity of intracerebral EEG, Clinical Neurophysiology, Vol. 131, No. 8, 1866-1885,\
DOI: [10.1016/j.clinph.2020.04.160](https://doi.org/10.1016/j.clinph.2020.04.160).\
Preprint freely available at DOI: [10.1101/2020.03.13.991299](https://doi.org/10.1101/2020.03.13.991299).

- Dellavale D, Velarde O, Mato G, Urdapilleta E (2020), Complex interplay between spectral harmonicity and different types of cross frequency couplings in non linear oscillators and biologically plausible neural network models, Phys. Rev. E 102(6), 062401,\
DOI: [10.1103/PhysRevE.102.062401](https://doi.org/10.1103/PhysRevE.102.062401).\
Preprint freely available at DOI: [10.1101/2020.10.15.341800](https://doi.org/10.1101/2020.10.15.341800).
