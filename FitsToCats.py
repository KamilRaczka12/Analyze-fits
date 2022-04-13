import numpy as np
import glob
import logging
import subprocess as sub
import os
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits


logging.basicConfig(filename='FitsToCats.log', filemode='w', 
                    format='%(levelname)s:%(message)s',
                    encoding='utf-8', level=logging.INFO)


def main(input_dict):
    images = sorted(glob.glob(input_dict['path_to_folder'] + "*.fit*"))
    check_if_all_input_variables_filled(images, **input_dict)
    for im in images:
        im_name = get_im_name(im)
        logging.info('Program started working on image %s.' % im_name)
        hdr_values = get_im_hdrs(im, im_name, **input_dict)
        file_new = im_astr(im, im_name, **hdr_values, **input_dict)
        print(file_new)
        new_to_sextractor(im_name, file_new, **hdr_values, **input_dict)
        remove_files(**input_dict)
        logging.info('Program finished working on image %s.' % im_name)


def check_if_all_input_variables_filled(images, path_to_folder, path_to_astr_file, path_to_sex_config_file, fits_hdrs, **kwargs):
    if all([images, path_to_folder, path_to_astr_file, path_to_sex_config_file, fits_hdrs]):
        logging.info('All input variables filled.')
    else:
        logging.warning('One or more input variables empty!')


def get_im_name(image):
    pathname = os.path.splitext(image)[0]	
    image_name = pathname.split('/')[-1]
    return image_name


def get_im_hdrs(im, im_name, fits_hdrs, **kwargs):
    hdul = fits.open(im)
    dict_keys = ['RA', 'DEC', 'time_data']
    hdr_values = {}
    k = 0
    for n, key in enumerate(dict_keys):
        try:
            hdr_values[key] = hdul[0].header[fits_hdrs[n]]
        except KeyError:
            logging.warning('header for %s in image %s is incorrect!', key, im_name)
            hdr_values[key] = ''
            pass
        else:
            k += 1
    if k == len(dict_keys):
        logging.info('All headers in image %s are correct.' % im_name)
    hdr_values = change_format_ra_dec(hdr_values, dict_keys)
    return hdr_values


def change_format_ra_dec(hdr_values, dict_keys):
    date_hdrs = dict_keys[0:2]
    RA = float(hdr_values[date_hdrs[0]])
    DEC = float(hdr_values[date_hdrs[1]])
    if type(RA) == str:
        for hdr in date_hdrs:
            hdr_values[hdr] = hdr_values[hdr].replace(' ', ':')
    elif type(RA) == float:
        c = SkyCoord(ra=RA*u.degree, dec=DEC*u.degree, frame='icrs')
        a = c.to_string('hmsdms')
        a1 = a.split(' ')
        for n, i in enumerate(a1):
            i = i.replace('h', ':')
            i = i.replace('d', ':')
            i = i.replace('m', ':')
            i = i.replace('s', '')
            hdr_values[date_hdrs[n]] = i
        hdr_values['RA'] = hdr_values['RA'][0:8]
        hdr_values['DEC'] = hdr_values['DEC'][0:9]
        hdr_values['time_data'] = hdr_values['time_data'][0:22]
    else:
        pass
    return hdr_values


def im_astr(im, im_name, RA, DEC, path_to_astr_file, **kwargs):
    logging.info('Astrometry for image %s started.' % im_name)
    solve_field_command = ['solve-field', '--ra', '%s' % RA, '--dec', '%s' % DEC,
                    '--radius', '1', '--cpulimit', '30', '--config', path_to_astr_file,
                    '--overwrite', '--no-verify', '--no-plots', '%s' % im]
    sub.Popen(solve_field_command, stdout=sub.PIPE, stderr=sub.PIPE).communicate()
    file_new = im.replace('.fits', '.new')
    if os.path.exists(file_new):
        logging.info('Astrometry for image %s succesfully finished.' % im_name)
    else:
        logging.warning('Astrometry for image %s unsuccesful.' % im_name)
    return file_new


def new_to_sextractor(im_name, file_new, path_to_folder, time_data, path_to_sex_config_file, **kwargs):
    logging.info('Sextractor started working on image %s.' % im_name)
    file_cat = path_to_folder + im_name + time_data + '.cat'
    analyse_new_command = ['source-extractor', '-c', path_to_sex_config_file,
                           file_new, '-CATALOG_NAME', file_cat]
    sub.Popen(analyse_new_command, stdout=sub.PIPE, stderr=sub.PIPE).communicate()
    if os.path.exists(file_cat):
        logging.info('Sextractor finished working on image %s and returned cat file.' % im_name)
    else:
        logging.warning('Sextractor didn`t make a new file while working on image %s.' % im_name)


def remove_files(path_to_folder, files_to_rm, **kwargs):
    for file_ in files_to_rm:
        file1_ = glob.glob(path_to_folder + file_)
        os.remove(file1_[0])


input_dict = {
    'path_to_folder' : '/home/kamil/Programs/Analyze-fits/Test_files/test/',
    'path_to_astr_file' : '/home/kamil/astrometry.net-0.89/etc/astrometry.cfg',
    'path_to_sex_config_file' : '/usr/share/source-extractor/default.sex',
    # 'path_to_folder' : '/home/kamilraczka/Projects/Na_zaliczenie/'
    # 'path_to_astr_file' : '/home/kamilraczka/astrometry.net-0.85/etc/astrometry.cfg'
    # 'sextractor_path_to_config_file' : '/home/kamilraczka/.config/sextractor/default.sex'
    'files_to_rm' : ['*.axy' , '*.corr', '*.xyls', '*.match',
                   '*.rdls', '*.solved', '*.wcs'],
    # 'fits_hdrs' : ['OBJCTRA', 'OBJCTDEC', 'DATE-OBS']
    'fits_hdrs' : ['RA_OBJ', 'DEC_OBJ', 'DATE-OBS']
    }

if __name__ == '__main__':
    main(input_dict)
