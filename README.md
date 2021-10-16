# SMS-EPI-Ghost-Correction
Here is the code for paper:

**Robust SENSE reconstruction of simultaneous multislice EPI with low-rank enhanced coil sensitivity calibration and slice-dependent 2D Nyquist ghost correction**

by Mengye Lyu, Markus Barth, Victor B. Xie, Yilong Liu, Xin Ma, Yanqiu Feng, and Ed X. Wu

link:  https://doi.org/10.1002/mrm.27120

## Abstract
### Purpose
To improve simultaneous multislice (SMS) EPI by robust Nyquist ghost correction in both coil sensitivity calibration and SMS reconstruction.

### Methods
To derive coil sensitivity and slice-dependent phase difference map between positive- and negative-echo images, single-band EPI reference data are fully sampled with EPI parameters matched to SMS acquisition. First, the reference data are organized into positive- and negative-echo virtual channels where missing data are estimated using low-rank-based simultaneous autocalibrating and k-space estimation (SAKE) at small matrix size. The resulting ghost-free positive- and negative-echo images are combined to generate coil sensitivity maps. Second, full-matrix positive- and negative-echo images are SENSE reconstructed from the reference data. Their phase difference or error map is then calculated. Last, SMS EPI is reconstructed using phase error correction SENSE (PEC-SENSE) that incorporates phase error map into coil sensitivity maps for negative-echo data. The proposed method was evaluated using both experimental data from 7T systems and simulations.

### Results
Virtual coil SAKE eliminated Nyquist ghosts in the single-band EPI, yielding high-quality coil sensitivity maps and phase error maps. The subsequent PEC-SENSE robustly reconstructed SMS EPI under various conditions, including presence of in-plane acceleration, with lesser artifacts and higher temporal SNR than slice-dependent 1D linear correction method.

### Conclusion
The proposed procedure of virtual coil SAKE calibration and PEC-SENSE reconstruction substantially reduces all ghost-related artifacts originating either directly from SMS EPI data or indirectly from EPI-based coil sensitivity maps. It is computationally efficient, and generally applicable to all SMS EPI-based applications.

Keywords: EPI; SMS; ghost artifact; multiband; parallel imaging; virtual coil.

![image](https://user-images.githubusercontent.com/10205514/137596076-4995fc1a-08c4-46cb-9b35-5eb0d74e5391.png)
