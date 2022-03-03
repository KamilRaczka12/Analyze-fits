import numpy as np
import glob
import subprocess as sub
import os
from astropy.io import fits


def do_the_work(image, path_to_folder, files_to_rm):
    RA, DEC, time_data = get_image_headers(image)
    image_astrometry(image, RA, DEC)
    astrometry_image = sorted(glob.glob(image.replace('.fit', '.new')))
    new_to_sextractor(astrometry_image, image, path_to_folder, time_data)
    remove_files(path_to_folder, files_to_rm)


def get_image_headers(image):
    hdul = fits.open(image)
    ra = hdul[0].header['OBJCTRA']
    dec = hdul[0].header['OBJCTDEC']
    time_data = hdul[0].header['DATE-OBS']
    RA = ra.replace(' ', ':')
    DEC = dec.replace(' ', ':')
    return RA, DEC, time_data


def image_astrometry(image, RA, DEC):
    solve_field_command = ['solve-field', '--ra', '%s' % RA, '--dec', '%s' % DEC,
                    '--radius', '1', '--config',
                    '/home/kamilraczka/astrometry.net-0.85/etc/astrometry.cfg',
                    '--overwrite', '--no-verify', '--no-plots', '%s' % image]
    sub.Popen(solve_field_command, stdout=sub.PIPE, stderr=sub.PIPE).communicate()


def new_to_sextractor(astrometry_image, image, path_to_folder, time_data):
    image_base_name = os.path.basename(image)
    if image_base_name.endswith('.fit'):
        image_base_name = image_base_name[:-4]
    elif image_base_nam.endswith('.fits'):
        image_base_name = image_base_name[:-5]
    analyse_new_command = ['source-extractor', '-c', 
                           '/home/kamilraczka/.config/sextractor/default.sex',
                           astrometry_image[0], '-CATALOG_NAME', path_to_folder + 
                           image_base_name + time_data + '.cat']
    sub.Popen(analyse_new_command, stdout=sub.PIPE, stderr=sub.PIPE).communicate()


def remove_files(path_to_folder, files_to_rm):
    for file_ in files_to_rm:
        file_to_rm = glob.glob(path_to_folder + file_)
        os.remove(file_to_rm[0])




path_to_folder = '/home/kamilraczka/Projects/Na_zaliczenie/'
fits_images = sorted(glob.glob(path_to_folder + "*.fit*"))
# fits_images = glob.glob(path_to_folder + 'UVMon-003_R.fit')
files_to_rm = ['*.axy' , '*.corr', '*.xyls', '*.match'
                , '*.new', '*.rdls', '*.solved', '*.wcs']
sextractor_path_to_config_file = '/home/kamilraczka/.config/sextractor/default.sex'

for image in fits_images:
    do_the_work(image, path_to_folder, files_to_rm)

