from matplotlib import pyplot as plt
import numpy as np
import astropy.units as u
from reproject import reproject_adaptive as reprojection
from matplotlib.colors import SymLogNorm, LogNorm
import os
from astropy.cosmology import FlatLambdaCDM
from astropy.io import fits
import warnings
import pkg_resources
import sys
from past.utils import old_div

warnings.filterwarnings("ignore")

def get_noise(data):
    """
    from Cyril Tasse/kMS
    """
    maskSup = 1e-7
    m = data[np.abs(data) > maskSup]
    rmsold = np.std(m)
    diff = 1e-1
    cut = 3.
    med = np.median(m)
    for _ in range(10):
        ind = np.where(np.abs(m - med) < rmsold * cut)[0]
        rms = np.std(m[ind])
        if np.abs(old_div((rms - rmsold), rmsold)) < diff: break
        rmsold = rms
    return rms

class RedShifting:
    def __init__(self, fitsfile='', redshift=None, noise=None, keep_noise=True, noise_level=3, middle_pos=None,
                 spectral_index=1, gaussian_kernel=False, cosmology=None):
        """
        :param fits: fits file location
        :param redshift: redshift
        :param noise: give noise level in Jy/Beam
        :param middle_pos: give middle in image
        :param spectral_index: give spectral index of source. Important for pixel scaling of real emission.
        :param keep_noise: keep the noise level the same (otherwise it scales with other pixels)
                            One can decide to keep_noise=True to speed up the algorithm and if noise does not matter
        :param noise_level: factor of rms noise where above is real emmission. Example: noise_level=3 -->
                            emission above 3 times rms is real. Only works in combination with keep_noise.
        """

        self.hdu = fits.open(fitsfile)[0]

        # original redshift
        self.redshift = redshift

        # cosmology
        if cosmology is None:
            self.cosmo = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Tcmb0=2.725 * u.K, Om0=0.3)

        # use noise or calculate noise if not given
        if noise:
            self.noise = noise
        else:
            try:
                self.noise = get_noise(self.hdu.data) #noise(self.hdu.data) #flux limit is 7e-5 Jy/beam
            except ModuleNotFoundError:
                sys.exit('Give noise to RedShifting')
        if self.noise * 4 < self.hdu.data.max() / 50:
            self.noise = self.hdu.data.max() / 200

        # keep the noise level the same (True) or scale with reprojection
        self.keep_noise = keep_noise
        self.noise_level =  noise_level
        if self.noise_level!=3 and not self.keep_noise:
            print('WARNING: noise_level given but keep_noise==False, so noise_level has no effect.')

        # header for simulation
        self.sim_header = self.hdu.header

        # center position in image
        if middle_pos:
            self.hdu.header['CRVAL1'], self.hdu.header['CRVAL2'] = middle_pos
            self.hdu.header['CRPIX1'], self.hdu.header['CRPIX2'] = len(self.hdu.data) / 2, len(self.hdu.data) / 2

        # spectral index
        if spectral_index is not None:
            self.spectral_index = spectral_index

        # use gaussian kernel with information from the beam to smooth image
        self.gaussian_kernel = gaussian_kernel
        if self.gaussian_kernel and float(pkg_resources.get_distribution("reproject").version)<0.9:
            sys.exit("ERROR: Reproject 0.9 is required for gaussian kernel\nUpdate Reproject to version >0.9 or "
                     "set gaussian_kernel=False.")

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
        Shift image
        :param dz: delta redshift
        :param save_as: file name if you want to save
        :return: shifted image
        """

        # determine pixel scaling, based on redshift increment
        pix_scaling = (self.cosmo.angular_diameter_distance(self.redshift)/self.cosmo.angular_diameter_distance(self.redshift+dz)).value

        # make new header for redshift increment
        self.sim_header = self.hdu.header.copy()

        # update simulated header
        self.sim_header['NAXIS1'], self.sim_header['NAXIS2'] = int(self.hdu.data.shape[-2] * pix_scaling), int(self.hdu.data.shape[-1] * pix_scaling)
        self.sim_header['CRPIX1'], self.sim_header['CRPIX2'] = self.hdu.data.shape[-2] * pix_scaling / 2, self.hdu.data.shape[-1] * pix_scaling / 2
        self.sim_header['CDELT1'] /= pix_scaling
        self.sim_header['CDELT2'] /= pix_scaling
        self.sim_header['CDELT2_original'] = self.hdu.header['CDELT2']

        # make copy of hdu
        hdu_copy = self.hdu.copy()

        # pixel scaling, based on redshift increment
        pixel_reduction = (self.cosmo.scale_factor(self.redshift)/self.cosmo.scale_factor(self.redshift + dz)) ** (3+self.spectral_index)

        if self.keep_noise:
            # determine real emission points
            real_emission = np.argwhere(hdu_copy.data>self.noise_level*self.noise)
            for p in real_emission:
                hdu_copy.data[p[0], p[1]] /= pixel_reduction
                # if pixels reduce under noise, we replace it with a new gaussian value around the noise
                if hdu_copy.data[p[0], p[1]]<self.noise:
                    hdu_copy.data[p[0], p[1]] = abs(np.random.normal(self.noise, self.noise, 1))
        else:
            hdu_copy.data /= pixel_reduction

        #reproject
        if self.gaussian_kernel:
            image, _ = reprojection(hdu_copy, self.sim_header, kernel='gaussian', kernel_width=np.sqrt(2*self.beam_area.value/np.pi))
        else:
            image, _ = reprojection(hdu_copy, self.sim_header)

        if save_as:
            # Save extra info in header
            self.sim_header['CDELT1_original'] = self.hdu.header['CDELT1']
            self.sim_header['BMAJ_repr'] = self.hdu.header['BMAJ'] / pix_scaling
            self.sim_header['BMIN_repr'] = self.hdu.header['BMIN'] / pix_scaling
            self.sim_header['pix_scale'] = pix_scaling
            fits.writeto(save_as, image, self.sim_header, overwrite=True)

        return image

    def reduce_flux(self, factor):
        """
        Multiply flux with factor (used for RLF completeness reconstruction)
        :param factor: multiply data with this factor (to reduce flux)
        :return: shifted image
        """
        return self.hdu.data*factor


    def make_image(self, dz, save_as='', video=False, same_imagescale=True):
        """
        Make image of source
        :param dz: delta z
        :param save_as: name of image
        :param video: make image for video
        :return:
        """
        image_data = self.shift(dz)
        if video:
            while image_data.shape[0]%2==1: # we need oneven size otherwise the video will be jizzy
                dz+=0.001
                image_data = self.shift(dz)
        if same_imagescale:
            plt.imshow(image_data, norm=SymLogNorm(linthresh=np.nanstd(self.hdu.data), vmin=0.1*np.nanstd(self.hdu.data),
                                                   vmax=np.nanstd(self.hdu.data) * 20), cmap='CMRmap')
        else:
            plt.imshow(image_data, norm=LogNorm(vmin=0.001*np.nanstd(image_data),
                                                   vmax=np.nanstd(image_data) * 10), cmap='CMRmap')

        plt.title(f'Redshift: {np.round(self.redshift+dz, 3)}')
        plt.tight_layout()
        plt.grid(False)
        plt.axis('off')
        if save_as:
            plt.savefig(save_as)
        else:
            plt.show()
        plt.close()

        return self

    def make_video(self, dz_max, save_as=''):
        """
        Make a video with ffmpeg
        :param dz_max: max delta z
        :param save_as: file name
        :return:
        """
        frame_dz = np.linspace(0, dz_max, 300)
        os.system('mkdir -p video_frames')
        for n, dz in enumerate(frame_dz):
            self.make_image(dz, f'video_frames/frame_{str(n).rjust(5, "0")}.png', video=True)
        if not save_as:
            save_as='movie.mp4'
        os.system(f'ffmpeg -f image2 -r 10 -start_number 0 -i video_frames/frame_%05d.png {save_as} && rm -rf video_frames')

        return self


def move_source(input_fits='', dz=None, orig_z=None, noise=None, noise_level=3,
                spectral_index=1, gaussian_kernel=False, output_fits=''):
    """
    Move a source to a new higher redshift.

    :param input_fits: input fits name.
    :param dz: redshift increment.
    :param redshift: original redshift.
    :param noise: noise in image.
    :param noise_level: factor of rms noise where above is real emmission. Example: noise_level=3 -->
                        emission above 3 times rms is real. Only works in combination with keep_noise.
    :param spectral_index: give spectral index of source. Important for pixel scaling of real emission.
    :param gaussian_kernel: use gaussian kernel based on beam information.
    :param output_fits: output fits file name.
    """
    source = RedShifting(fitsfile=input_fits, redshift=orig_z, noise=noise, noise_level=noise_level,
                         spectral_index=spectral_index, gaussian_kernel=gaussian_kernel)
    if dz is None:
        sys.exit('Please give dz (redshift increment).')

    if dz<0:
        sys.exit('Redshift increment has to be positive.')

    if orig_z<0:
        sys.exit('Original redshift has to be positive')

    source.shift(dz=dz, save_as=output_fits)

if __name__ == '__main__':
    print("Do not call this function directly.")