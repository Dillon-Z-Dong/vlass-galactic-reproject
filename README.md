# FITS Reprojection Tool

Transform FITS images from RA/Dec to Galactic coordinates using astropy / reproject.

## Installation

Install dependencies:

```bash
pip install -r requirements.txt
```

## Usage

1. Place FITS images to be transformed in the `input_images/` directory
2. Run the script:

```bash
python reproject_to_galactic.py
```

Output files will be saved to `output_images/` with `_galactic` suffix (e.g., `image.fits` → `image_galactic.fits`).

