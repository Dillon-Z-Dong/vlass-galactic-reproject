#!/usr/bin/env python3
"""
Transform FITS image from RA/Dec to Galactic coordinates using astropy and reproject.
"""

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from reproject import reproject_interp
import sys
import os

def create_galactic_wcs(center_coord, image_shape, pixel_scale):
    """
    Create a WCS object for Galactic coordinates.
    
    Parameters:
    -----------
    center_coord : SkyCoord
        Center coordinate of the image
    image_shape : tuple
        (height, width) of the output image
    pixel_scale : float
        Pixel scale in degrees/pixel
    
    Returns:
    --------
    WCS object in Galactic coordinates
    """
    # Convert center to Galactic coordinates
    gal_center = center_coord.galactic
    
    # Create new WCS
    wcs_out = WCS(naxis=2)
    
    # Set coordinate system to Galactic
    wcs_out.wcs.ctype = ['GLON-TAN', 'GLAT-TAN']
    wcs_out.wcs.cunit = ['deg', 'deg']
    
    # Set reference pixel (center of image)
    wcs_out.wcs.crpix = [image_shape[1]/2, image_shape[0]/2]
    
    # Set reference coordinates (Galactic longitude and latitude)
    wcs_out.wcs.crval = [gal_center.l.degree, gal_center.b.degree]
    
    # Set pixel scale (degrees per pixel)
    wcs_out.wcs.cd = np.array([[-pixel_scale, 0], [0, pixel_scale]])
    
    return wcs_out

def transform_fits_to_galactic(input_file, output_file, output_shape, pixel_scale):
    """
    Transform FITS image from RA/Dec to Galactic coordinates.
    """
    
    print(f"Reading {input_file}...")
    
    # Read input FITS file
    with fits.open(input_file) as hdul:
        hdu = hdul[0]
        data = hdu.data
        header = hdu.header.copy()
    
    # Create input WCS
    wcs_in = WCS(header, naxis=2) 
    print(f"Input coordinate system: {wcs_in.wcs.ctype}")
    
    ny, nx = data.shape
    
    # Get image center coordinates
    center_pix = [nx/2, ny/2]
    center_coord = wcs_in.pixel_to_world(*center_pix)
    
    print(f"Image center: RA={center_coord.ra.degree:.4f}°, Dec={center_coord.dec.degree:.4f}°")
    print(f"Using pixel scale: {pixel_scale:.10f} deg/pixel ({pixel_scale*3600:.2f} arcsec/pixel)")
    
    # Create output WCS in Galactic coordinates
    wcs_out = create_galactic_wcs(center_coord, output_shape, pixel_scale)
    
    gal_center = center_coord.galactic
    print(f"Galactic center: l={gal_center.l.degree:.4f}°, b={gal_center.b.degree:.4f}°")
    
    # Perform reprojection
    print("Reprojecting using interpolation method...")
    result, footprint = reproject_interp((data, wcs_in), wcs_out, shape_out=output_shape)
    
    # Create new header with Galactic WCS
    new_header = wcs_out.to_header()
    
    # Copy all metadata from original header (except spatial WCS keywords)
    spatial_wcs_keywords = ['CTYPE1', 'CTYPE2', 'CUNIT1', 'CUNIT2', 'CRPIX1', 'CRPIX2', 
                           'CRVAL1', 'CRVAL2', 'CDELT1', 'CDELT2', 'PC1_1', 'PC1_2', 'PC2_1', 'PC2_2',
                           'LONPOLE', 'LATPOLE', 'RADESYS', 'EQUINOX']
    
    for key, value in header.items():
        if key not in spatial_wcs_keywords:
            new_header[key] = value
    
    # Add transformation info
    new_header['HISTORY'] = 'Reprojected from RA/Dec to Galactic coordinates'
    new_header['HISTORY'] = f'Original file: {input_file}'
    new_header['COMMENT'] = f'Pixel scale: {pixel_scale:.10f} deg/pixel'
    
    # Save result
    print(f"Saving to {output_file}...")
    hdu_out = fits.PrimaryHDU(data=result, header=new_header)
    hdu_out.writeto(output_file, overwrite=True)
    
    print(f"Successfully saved reprojected image to {output_file}")
    print(f"Output shape: {result.shape}")
    print(f"Data range: {np.nanmin(result):.3e} to {np.nanmax(result):.3e}")

def main():
    input_dir = 'input_images'
    output_dir = 'output_images'
    pixel_scale = 0.0002777777777778  # 1 arcsec/pixel
    output_shape = (5300, 5300)  # Larger to account for rotation

    os.makedirs(output_dir, exist_ok=True)

    fits_files = [f for f in os.listdir(input_dir) if f.lower().endswith('.fits')]

    print(f"Found {len(fits_files)} FITS files in {input_dir}.")
    print(f"Pixel scale: {pixel_scale} deg/pixel (1 arcsec/pixel)")
    print(f"Output shape: {output_shape}")
    print("-" * 50)

    for filename in fits_files:
        input_file = os.path.join(input_dir, filename)
        base, _ = os.path.splitext(filename)
        output_file = os.path.join(output_dir, f"{base}_galactic.fits")

        print(f"Transforming {input_file} -> {output_file}")
        transform_fits_to_galactic(input_file, output_file, output_shape, pixel_scale)
        print("-" * 50)

    print(f"Converted {len(fits_files)} VLASS median stack image(s) to Galactic coordinates.")

if __name__ == '__main__':
    main()