#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun  9 19:52:14 2024

@author: shyli

"""
from platipy.imaging.registration.deformable import fast_symmetric_forces_demons_registration
from platipy.imaging.registration.utils import apply_deformable_transform
import numpy as np
import SimpleITK as sitk
import matplotlib.pyplot as plt
import os
from scipy.io import savemat
import multiprocessing


def resize_image(image, new_size, interpolator=sitk.sitkLinear):
    """
    Resize an image to a new size using the specified interpolator.
    
    Parameters:
    - image: The input SimpleITK image to be resized.
    - new_size: A tuple specifying the desired output size (width, height).
    - interpolator: The interpolator to use (default is sitk.sitkLinear).
    
    Returns:
    - The resized SimpleITK image.
    """
    original_size = image.GetSize()
    original_spacing = image.GetSpacing()
    
    new_spacing = [
        (original_size[i] * original_spacing[i]) / new_size[i] for i in range(2)
    ]
    
    resampled_image = sitk.Resample(
        image,
        new_size,
        sitk.Transform(),
        interpolator,
        image.GetOrigin(),
        new_spacing,
        image.GetDirection(),
        0.0,
        image.GetPixelID()
    )
    
    return resampled_image

# Define upper level directory
folder = '/projectnb/npbssmic/ns/Ann_Mckee_samples_55T/metrics/gsigma_1-3-5_2-3-4_3-5-7_5-7-9_7-9-11/p18/heatmaps/heatmaps_pathology_registration/'

# Define subfolders to examine
subfolders = ['AD_10382','AD_21354','AD_21424','CTE_7019','CTE_7126','NC_8095','NC_21499']

# Define the types of stains for each subfolder
stains = ['ab','pt']

# Define the filenames for the AB and p-tau stains
fnames = {
    'heatmap_file': {
        'ab': 'heatmap_A-Beta_mask_registration.tif',
        'pt': 'heatmap_p-tau_mask_registration.tif'
        },
    'ab_mask_file': {
        'ab': '_Ab_mask_registration.tif',
        'pt': '_AT8_mask_registration.tif'
        },
    'mask_analysis_file': {
        'ab': '_Ab_mask.tif',
        'pt': '_AT8_mask.tif'
        },
    'path_hm_file': {
        'ab': '_Ab_stain_heatmap_204.tif',
        'pt': '_AT8_stain_heatmap_204.tif'
        },
    'path_fig': {
        'ab': '_Ab_path_registration.png',
        'pt': '_AT8_path_registration.png',
        },
    'path_mat': {
        'ab': '_Ab_path_registered.mat',
        'pt': '_AT8_path_registered.mat'
        },
    'path_mask_mat': {
        'ab': '_Ab_path_mask_registered.mat',
        'pt': '_AT8_path_mask_registered.mat'
        },
    'mask_fig': {
        'ab': '_Ab_mask_registration.png',
        'pt': '_AT8_mask_registration.png'
        }
    }

# Iterate over each subfolder
for subfolder in subfolders[5:6]:
    print('Starting subject',subfolder,'\n')
    for stain in stains:
        niter1 = 800
        niter2 = 300
        niter3 = 300
    
        subfolder_name = os.path.basename(subfolder)
        
        # indivudual adjustment for number of iterations
        if 'CTE_6489' in subfolder_name:
            niter1 = 1200
        if 'AD_20382' in subfolder_name:
            niter1 = 2000
        
        ## Retrieve the images
        # Vasculature heatmap mask
        heatmap_file = os.path.join(subfolder, fnames['heatmap_file'][stain])
        image_1= sitk.ReadImage(heatmap_file)
        # Pathology heatmap registration mask
        ab_mask_file = [f for f in os.listdir(subfolder) if f.endswith(fnames['ab_mask_file'][stain])][0]
        image_2 = sitk.ReadImage(os.path.join(subfolder, ab_mask_file))
        # Pathology heatmap analysis mask
        analysis_mask_file = [f for f in os.listdir(subfolder) if f.endswith(fnames['mask_analysis_file'][stain])][0]
        analysis_mask = sitk.ReadImage(os.path.join(subfolder, analysis_mask_file))
        # Pathology heatmap
        path_hm_file = [f for f in os.listdir(subfolder) if f.endswith(fnames['path_hm_file'][stain])][0]
        path = sitk.ReadImage(os.path.join(subfolder, path_hm_file))        
        
        # Rotate the array by 180 degrees
        if ('NC_21499' in subfolder_name) or ('CTE_6912' in subfolder_name):
            image_array = sitk.GetArrayViewFromImage(image_2)
            rotated_array = image_array[::-1, ::-1]
            rotated_image_2 = sitk.GetImageFromArray(rotated_array)
            rotated_image_2.SetOrigin(image_2.GetOrigin())
            rotated_image_2.SetSpacing(image_2.GetSpacing())
            rotated_image_2.SetDirection(image_2.GetDirection())
            image_2 = rotated_image_2
    
        # Retrieve size of vascular heatmap mask
        hm_mask_size = image_1.GetSize()
        
        # Resample pathology images to the size and spacing of heatmap mask
        cropped_image_2 = resize_image(image_2,hm_mask_size)
        path = resize_image(path,hm_mask_size)
        analysis_mask = resize_image(analysis_mask,hm_mask_size)
        
        # Set the origin and spacing of moving (pathology) equal to fixed (OCT)
        origin = image_1.GetOrigin()
        spacing = image_1.GetSpacing()
        cropped_image_2.SetOrigin(origin)
        cropped_image_2.SetSpacing(spacing)
        path.SetOrigin(origin)
        path.SetSpacing(spacing)
        analysis_mask.SetOrigin(origin)
        analysis_mask.SetSpacing(spacing)
        
        # %% Register the pathology mask to the OCT mask
        print('Registering\n')
        image_2_deformed, tfm, dvf = fast_symmetric_forces_demons_registration(
        image_1,
        cropped_image_2, 
        resolution_staging=[8, 4, 1], 
        iteration_staging=[niter1, niter2, niter3],
        isotropic_resample=False, 
        regularisation_kernel_mm=5, 
        smoothing_sigma_factor=1,
        smoothing_sigmas=False, 
        ncores=4, 
        interp_order=2, 
        )
        
        # %% Apply the transformation to the pathology heatmap and mask
        
        # Register the pathology heatmap
        path_deformed = apply_deformable_transform(path, tfm)
        
        # Register the pathology analysis mask
        analysis_mask_deformed = apply_deformable_transform(analysis_mask, tfm)
        
        # Convert the images to numpy arrays for visualization & saving to MAT
        image_1_array = sitk.GetArrayFromImage(image_1)
        resampled_image_2_array = sitk.GetArrayFromImage(cropped_image_2)
        image_2_deformed_array = sitk.GetArrayFromImage(image_2_deformed)
        dvf_array = sitk.GetArrayFromImage(dvf)
        path_array = sitk.GetArrayFromImage(path)
        path_deformed_array = sitk.GetArrayFromImage(path_deformed)
        mask_deformed_array = sitk.GetArrayFromImage(analysis_mask_deformed)
        
        ### Before/after for pathology mask
        # Initialize before/after for pathology mask
        overlay_before = np.zeros((*image_1_array.shape, 3), dtype=np.uint8)
        overlay_after = np.zeros((*image_1_array.shape, 3), dtype=np.uint8)
        # Set red channel to image_1 and green channel to resampled_image_2 for overlay_before
        overlay_before[:,:,0] = image_1_array
        overlay_before[:,:,1] = resampled_image_2_array
        # Set red channel to image_1 and green channel to image_2_deformed for overlay_after
        overlay_after[:,:,0] = image_1_array
        overlay_after[:,:,1] = image_2_deformed_array
        
        ### Before/after for pathology heatmap
        overlay_before_path = np.zeros((*image_1_array.shape, 3), dtype=np.uint8)
        overlay_after_path = np.zeros((*image_1_array.shape, 3), dtype=np.uint8)
        # Set red channel to image_1 and green channel to resampled_image_2 for overlay_before
        overlay_before_path[:,:,0] = image_1_array
        overlay_before_path[:,:,1] = path_array
        # Set red channel to image_1 and green channel to image_2_deformed for overlay_after
        overlay_after_path[:,:,0] = image_1_array
        overlay_after_path[:,:,1] = path_deformed_array
        # Plot the figures
        plt.figure(figsize=(20, 5))
        plt.subplot(1, 2, 1)
        plt.imshow(overlay_before_path)
        plt.title('Pathology Heatmap Overlay Before Registration')
        plt.axis('off')
        plt.subplot(1, 2, 2)
        plt.imshow(overlay_after_path)
        plt.title('Pathology Heatmap Overlay After Registration')
        plt.axis('off')
        plt.tight_layout()
        plt.savefig(subfolder_name + fnames['path_fig'][stain], dpi=300)
        plt.show()
        
        ### Save pathology heatmap as .MAT
        mat_filename = subfolder_name + fnames['path_mat'][stain]
        savemat(mat_filename, {
            'path_registered': path_deformed_array,
        })
        
        ### Save pathology heatmap mask as .MAT
        mat_filename = subfolder_name + fnames['path_mask_mat'][stain]
        savemat(mat_filename, {
            'path_mask_registered': mask_deformed_array,
        })
        
        
        # Plot the images
        plt.figure(figsize=(20, 5))
        
        plt.subplot(1, 6, 1)
        plt.imshow(image_1_array, cmap='gray')
        plt.title(subfolder_name + ' Image 1')
        plt.axis('off')
        
        plt.subplot(1, 6, 2)
        plt.imshow(resampled_image_2_array, cmap='gray')
        plt.title(subfolder_name + ' Image 2')
        plt.axis('off')
        plt.subplot(1, 6, 3)
        plt.imshow(overlay_before)
        plt.title('Overlay Before Registration')
        plt.axis('off')
        
        plt.subplot(1, 6, 4)
        plt.imshow(overlay_after)
        plt.title('Overlay After Registration')
        plt.axis('off')
        
        plt.subplot(1, 6, 5)
        plt.imshow(dvf_array[:,:,0])
        plt.title('Deformation field - x')
        plt.axis('off')
        
        plt.subplot(1, 6, 6)
        plt.imshow(dvf_array[:,:,1])
        plt.title('Deformation field - y')
        plt.axis('off')
        
        plt.tight_layout()
        plt.savefig(subfolder_name + fnames['mask_fig'][stain], dpi=300)
        plt.show()