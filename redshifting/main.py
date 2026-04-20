import os
import sys
import warnings

import numpy as np
import astropy.units as u
from astropy.io import fits
from astropy.wcs import WCS
from astropy.cosmology import FlatLambdaCDM
from matplotlib import pyplot as plt
from matplotlib.colors import SymLogNorm, LogNorm
from reproject import reproject_adaptive as reprojection

warnings.filterwarnings("ignore")

class RedShifting:
    def __init__(self, fitsfile='', redshift=None, spectral_index=1, gaussian_kernel=False, cosmology=None):
        """
        Parameters
        ----------
        fitsfile : Path to the FITS file to be loaded.
        redshift : Redshift of the source.
        spectral_index : float, optional
            Spectral index of the source (S ∝ ν^−α).
            Used to scale flux correctly during redshifting.
        gaussian_kernel : If True, applies a Gaussian beam convolution during reprojection.
        cosmology : Cosmological model used for distance and scaling calculations.
        """

        # Read FITS
        self.hdu = fits.open(fitsfile)[0]

        # Original redshift
        self.redshift = redshift

        # Cosmology
        if cosmology is None:
            self.cosmo = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Tcmb0=2.725 * u.K, Om0=0.3)

        # Header for simulation
        self.sim_header = self.hdu.header

        # Spectral index
        if spectral_index is not None:
            self.spectral_index = spectral_index

        # Use gaussian kernel with information from the beam to smooth image
        self.gaussian_kernel = gaussian_kernel

    @property
    def beam_area(self):
        """
        Calculate beam area by:
        Beam Area = 2 * pi * sigma^2
        FWHM = 2 * sigma * sqrt(2*ln(2))
        Number of pixels in a beam = beam area [arcsec]/(pixel length)^2
        https://www.eaobservatory.org/jcmt/faq/how-can-i-convert-from-mjybeam-to-mjy/
        """
        gfactor = 2 * np.sqrt(2*np.log(2))
        FWHM = np.sqrt(((self.hdu.header['BMAJ'] * u.deg).to(u.arcsec).value * (self.hdu.header['BMIN'] * u.deg).to(u.arcsec).value))
        sigma = FWHM/gfactor
        pix_length = (self.hdu.header['CDELT2'] * u.deg).to(u.arcsec).value
        return 2 * np.pi * sigma**2 / pix_length ** 2 / u.beam # pixels per beam

    def shift(self, dz=0, save_as=''):
        """
        Shift a source image to a higher redshift using cosmological angular
        diameter scaling and surface brightness dimming.

        Parameters
        ----------
        dz : Redshift increment (must be positive).
        save_as : Output FITS filename. If provided, the shifted image is written to disk.

        Returns
        -------
        Shifted (reprojected) image.
        """

        # ----------------------------
        # Cosmological scalings
        # ----------------------------
        pix_scaling = (self.cosmo.angular_diameter_distance(self.redshift) / self.cosmo.angular_diameter_distance(self.redshift + dz)).value
        flux_scaling = (self.cosmo.scale_factor(self.redshift) / self.cosmo.scale_factor(self.redshift + dz)) ** (3 + self.spectral_index)

        # ----------------------------
        # Extract data
        # ----------------------------
        data = np.asarray(self.hdu.data)
        original_shape = data.shape

        data_2d = np.squeeze(data)
        slice_idx = []

        while data_2d.ndim > 2:
            slice_idx.append(0)
            data_2d = data_2d[0]

        data_2d = data_2d / flux_scaling

        # ----------------------------
        # Input WCS (force 2D)
        # ----------------------------
        wcs_in = WCS(self.hdu.header).celestial
        hdu_in = fits.PrimaryHDU(data=data_2d, header=wcs_in.to_header())

        # ----------------------------
        # Output WCS
        # ----------------------------
        wcs_out = wcs_in.deepcopy()
        wcs_out.wcs.cdelt /= pix_scaling
        wcs_out.wcs.crpix *= pix_scaling

        # ----------------------------
        # Output shape
        # ----------------------------
        ny, nx = data_2d.shape
        shape_out = (int(ny * pix_scaling), int(nx * pix_scaling))

        # ----------------------------
        # Reprojection
        # ----------------------------
        if self.gaussian_kernel:
            image_2d, _ = reprojection(hdu_in, wcs_out, shape_out=shape_out, kernel='gaussian',
                                       kernel_width=np.sqrt(2 * self.beam_area.value / np.pi))
        else:
            image_2d, _ = reprojection(hdu_in, wcs_out, shape_out=shape_out)

        # ----------------------------
        # Restore to original cube shape
        # ----------------------------
        if len(original_shape) > 2:
            new_data = np.zeros(original_shape[:-2] + image_2d.shape)
            idx = tuple(slice_idx) + (slice(None), slice(None))
            new_data[idx] = image_2d
        else:
            new_data = image_2d

        # ----------------------------
        # Save output safely
        # ----------------------------
        if save_as:
            header_out = wcs_out.to_header()
            fits.PrimaryHDU(data=image_2d, header=header_out).writeto(save_as, overwrite=True)

        return new_data

    def make_image(self, dz, save_as='', video=False, same_imagescale=True):
        """
        Generate a visualisation of the source at a shifted redshift.

        Parameters
        ----------
        dz : Redshift increment used for the shift.
        save_as : Filename to save the resulting image. If not provided, the image is displayed.
        video : If True, adjusts the image dimensions for video frame compatibility.
        """
        image_data = self.shift(dz)

        # Ensure even dimensions for video compatibility
        if video:
            while image_data.shape[0] % 2 == 1:
                dz += 0.001
                image_data = self.shift(dz)

        # Choose intensity scaling
        if same_imagescale:
            vmin = 0.1 * np.nanstd(self.hdu.data)
            vmax = 20 * np.nanstd(self.hdu.data)
            norm = SymLogNorm(linthresh=np.nanstd(self.hdu.data), vmin=vmin, vmax=vmax)
        else:
            vmin = 0.001 * np.nanstd(image_data)
            vmax = 10 * np.nanstd(image_data)
            norm = LogNorm(vmin=vmin, vmax=vmax)

        # Plot image
        plt.imshow(image_data, norm=norm, cmap='CMRmap')
        plt.title(f"Redshift: {np.round(self.redshift + dz, 3)}")
        plt.tight_layout()
        plt.grid(False)
        plt.axis('off')

        if save_as:
            plt.savefig(save_as)
        else:
            plt.show()

        plt.close()

        return self

    def make_video(self, dz_max, frames=100, save_as=''):
        """
        Generate a video sequence of the source evolving with redshift using ffmpeg.

        Parameters
        ----------
        dz_max : Maximum redshift increment defining the range of evolution.
        frames: Number of frames
        save_as : Output filename for the generated video. If not provided, defaults to 'movie.mp4'.
        """
        frame_dz = np.linspace(0, dz_max, frames)
        os.system('mkdir -p video_frames')
        for n, dz in enumerate(frame_dz):
            self.make_image(dz, f'video_frames/frame_{str(n).rjust(5, "0")}.png', video=True)
        if not save_as:
            save_as='movie.mp4'
        os.system(f'ffmpeg -f image2 -r 10 -start_number 0 -i video_frames/frame_%05d.png {save_as} && rm -rf video_frames')

        return self


def move_source(input_fits='', dz=None, z=None, spectral_index=1,
                gaussian_kernel=False, make_video=False, output_fits=''):
    """
    Move a source to a higher redshift by applying cosmological scaling.

    Parameters
    ----------
    input_fits : Path to the input FITS file.
    dz : Redshift increment to apply (must be positive).
    z : Original redshift of the source.
    spectral_index : Spectral index of the source (S ∝ ν^−α), used for flux scaling.
    gaussian_kernel : If True, applies Gaussian beam smoothing during reprojection.
    output_fits : Path for the output FITS file to be written.
    """

    source = RedShifting(fitsfile=input_fits, redshift=z,
                         spectral_index=spectral_index, gaussian_kernel=gaussian_kernel)
    if dz is None:
        sys.exit('Please give dz (redshift increment).')

    if dz<=0:
        sys.exit('Redshift increment has to be larger than 0.')

    if z<=0:
        sys.exit('Original redshift has to be larger than 0')

    if make_video:
        source.make_video(dz_max=dz+z, save_as=input_fits.split('/')[-1]+'.mp4')
    elif not make_video or output_fits:
        source.shift(dz=dz, save_as=output_fits)

if __name__ == '__main__':
    sys.exit("Do not call this script directly.")
