# Large Scale Holographic Particle 3D Imaging with Beam Propagation Method(BPM)
MatLab implementation of paper [*Large-scale holographic particle 3D imaging with the beam propagation model*](https://www.osapublishing.org/oe/fulltext.cfm?uri=oe-29-11-17159&id=451212). We provide the Beam Propagation method forward model and the reconstruction algorithm code, the simulated sample object, hologram, reconstructed particles and the sample experimental captured hologram, reconstruction result.

### Citation
If you find this project useful in your research, please consider citing our paper:

[Wang, Hao, et al. "Large-scale holographic particle 3D imaging with the beam propagation model. Opt. Express **29**(11), 17159-17172 (2021).](https://www.osapublishing.org/oe/fulltext.cfm?uri=oe-29-11-17159&id=451212)


### Abstract
We develop a novel algorithm for large-scale holographic reconstruction of 3D particle fields. Our method is based on a multiple-scattering beam propagation method (BPM)
combined with sparse regularization that enables recovering dense 3D particles of high refractive index contrast from a single hologram. We show that the BPM-computed hologram generates intensity statistics closely matching with the experimental measurements and provides up to 9× higher accuracy than the single-scattering model. To solve the inverse problem, we devise a computationally efficient algorithm, which reduces the computation time by two orders of magnitude as compared to the state-of-the-art multiple-scattering based technique. We demonstrate the superior reconstruction accuracy in both simulations and experiments under different scattering strengths. We show that the BPM reconstruction significantly outperforms the single-scattering method in particular for deep imaging depths and high particle densities.

### Overview Figure
<p align="center">
  <img src="/figure/Figure1.png">
</p>

### How to use code
Forward model: [main.m](FORWARD_CODE/main.m)<br>
Reconstruction: [main_inverse.m](INVERSE_CODE/main_inverse.m)

### Data
simulated object with dz = lambda/16: object/simulatedData/density_1.6\
ground truth object for the reconstruction object (with dz = lambda/16*100): object/simulatedDownsampledData/density_1.6.mat\
simulated sample hologram: holograms/simulatedHologram/density1.6_dn0.26.mat\
reconstructed simulated particles: results/sampleResults/simulation/density1.6_dn0.26.tiff   (need to extract files from .zip file)\
experimental captuered hologram: holograms/experimetnalHologram/density1.60.mat\
reconstruction experimental particles: results/sampleResult/experimental/density1.6.tiff   (need to extract files from .zip file)

### Data Accessibility
More data is in object.zip, holograms.zip and results.zip

