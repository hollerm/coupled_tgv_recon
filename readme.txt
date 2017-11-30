Matlab code for multi-channel regularized image reconstruction, including joint MR-PET reconstruction.

-------------------------
florian.knoll@nyumc.org
martin.holler@uni-graz.at
kristian.bredies@uni-graz.at
thomas.koesters@nyumc.org
21.8.2017
-------------------------

-------------------------
---General information---
-------------------------

This package includes the implementation and the data to reproduce some of the results presented in [1,2]. To get started, run the main control script tgv_recon_test.m. The main function that implements the algorithm from [1] is tgv_recon_joint.m. This package includes data for three test cases, which can be selected with the "testcase" switch. Examples 2 and 3, which perform simulated an in-vivo MR-PET reconstructions, require the package EMRecon [3] for evaluation of the PET-projector. This code packages includes a subset of EMRecon in the form of compiled mexfiles for windows and linux (64bit) that implement PET-projectors for the geometry of the Siemens Biograph mMR and a simple 2D example. To stay updated, and for the full functionality, we recommended to get the most version of EMRecon from the Department of Mathematics of the University of Muenster [3]. The package also uses the extremely convenient @p2DFT class to generater undersampled cartesian MR operators from Miki Lustig's SparseMRI code package: http://people.eecs.berkeley.edu/~mlustig/software/sparseMRI_v0.2.tar.gz

1) digital_brain_phantom_mr_denoising: A simple multi-contrast gaussian denoising example. The data is generated from the numerical brain phantom presented in [4], and three different MR contrasts and sequence parameters were simulated as described in [1].

2) digital_brain_phantom_mr_pet_2d: Using the same phantom as in the denoising example, this example demonstrates simulated multi-channel undersampled MR reconstruction and PET reconstruction.

3) mprage_2013_11_08: This is the in-vivo PET-MRI dataset that was used in the experiments for Figure 8 in [1]. The data (approximately 2.5GB) can be obtained from:  http://cai2r.net/resources/software. Note: This is an computationally expensive experiment. It requires roughly 30GB of memory and the runtime is several hours. The reason for this the PET projection operator, which has to be evaluated at the much smaller voxel size of the MR resolution (sidenote, since these projections could all be evaluated in parallel, this step could be substantially sped up with a GPU implementation).

-------------------------
-------Known issues------
-------------------------
EMRecon uses OpenMP, which can cause issues with Matlab mexfiles. A workaround for this is to add:

export LD_PRELOAD=/usr/lib/gcc/x86_64-linux-gnu/4.8/libgomp.so in .bashrc

or start matlab in the command line with:

LD_PRELOAD="/usr/lib/gcc/x86_64-linux-gnu/4.8/libgomp.so" matlab

-------------------------
---Relevant parameters---
-------------------------

The recommended settings for the hyperparameters for the individual experiments are given in the script tgv_recon_test.m and can also be found in [1]. They can be set in the cell-array recon_pars. The important fields are:

vec_norm: Vector norm: Options sep (separate, no coupling), frob (frobenius norm coupling), nuc (nuclear norm coupling) are available

stop_par: Number of iterations

ld: Regularization parameters for different operators, if only a single parameter is given, the data terms are weighted automatically according to their range

op_rescale: Scaling of the operator norms for faster convergence. Parameters rescale for [Gaussian, Poisson] operator. Default for real PET-MRI measurement operators is [3 10], otherwise just use 1 (no scaling).

We suggest looking at tgv_recon_joint.m for additional available parameters.

-------------------------
---New data/operators----
-------------------------

The code is designed to be easily extended to new data sets and forward operators. The data storage convention is briefly described below and we suggest looking at the provided datasets for more details.

In tgv_recon_test_script.m, the path and the filename for the new data set must be provided. The data files are expected to contain a (n_channels x 1) cell array of data structures. Each data structure must contain the name of a "parfun" member function that takes the data structure itself as input and returns a structure having at least the forward and adjoint operator and the measured data as members. A basic example how data from a given 'datafile.m' is loaded in the code is given as follows:

load('datafile.m')

for idx=1:length(data)
    
    parfun = str2func(data{idx}.parfun);
    fout = parfun(data{idx});
    
    %Set required function output    
    K{idx} = fout.K;    %Forward operator
    Kt{idx} = fout.Kt;  %Adjoint operator
    f0{idx} = fout.f0;  %Measured data
    
end

For synthetic data the structure of the required parfun can be very simple (e.g. for denoising the operators are just the identity), but for real measurement data from a PET-MR system, this is not that trivial. Essentially, the operators have to be defined such that the PET and MR images are on the same grid and spatially aligned after evaluation of the adjoint operator. We suggest looking at the provided mprage_data_parfun.m from the demo example. Essentially, it comes down to setting up the structure for the EMRecon PET-projector. Some of the parameters depend on the geometry of the used PET-MR system (the provided in-vivo example is from the Siemens Biograph mMR), some depend on the particular data (resolution, spatial offsets).

It is a requirement to guarantee convergence that the operators for the different imaging modalities have a norm of 1. This can be achieved by estimating the operator norm with a von Mises power iteration (https://en.wikipedia.org/wiki/Power_iteration), and then scaling the operators appropriately. For the demo examples we already estimated the norms and included the normalization in the parfun, for new data, this step has to be performed first. We include code to perform the power iteration (power_it.m, in the utils folder). It can be called as follows (f_pet and f_mri are PET and MRI data; K_pet, K_mri, Kt_pet, Kt_mri are corresponding forward and adjoint operators): 

iter = 10
tmp = Kt_pet(f_pet);
nrm_pet = power_it(K_pet,Kt_pet,size(tmp),iter);
display(['PET Operator norm: ',num2str(nrm_pet)])

tmp = Kt_mri(f_mri);
nrm = power_it(K_mri,Kt_mri,size(tmp),iter);
display(['MRI Operator norm: ',num2str(nrm_mri)])

If you consider this code to be useful for your research, please cite [1,2].

-------------------------
---References------------
-------------------------

[1] F Knoll, M Holler, T Koesters, R Otazo, K Bredies, DK Sodickson: Joint MR-PET reconstruction using a multi-channel image regularizer. IEEE Transactions on Medical Imaging 36: 1-16 (2017) (http://ieeexplore.ieee.org/document/7466848/)

[2] M. Holler, R. Huber, F. Knoll: Coupled regularization with multiple data discrepancies. Submitted, 2017.

[3] T. Koesters, K. Schaefers, F. Wuebbeling, "EMRECON: an expectation maximization based image reconstruction framework for emission tomography data", IEEE Nucl. Sci. Symp. Med. Imag. Conf. (NSS/MIC), pp. 4365-4368, Oct. 2011.

[4] B. Aubert-Broche, A. C. Evans, L. Collins, "A new improved version of the realistic digital brain phantom", Neuroimage, vol. 32, no. 1, pp. 138-145, Aug. 2006, [online] Available: http://dx.doi.org/10.1016/j.neuroimage.2006.03.052.








