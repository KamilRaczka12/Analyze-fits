from importlib.resources import path
# from inspect import _Object
from astropy.io import fits, ascii
from astropy.table import Table
from astropy.coordinates import SkyCoord, Angle
from astropy import units as u
from astroquery.jplhorizons import Horizons
from datetime import datetime
import glob
import logging
import numpy as np
import os
import pandas as pd
from photutils.datasets import make_100gaussians_image
from photutils.aperture import CircularAperture


logging.basicConfig(filename='FitsAndCatsToTable.log', filemode='w', 
                    format='%(levelname)s:%(message)s',
                    encoding='utf-8', level=logging.INFO)


def main(input_dict):
    fits_images = sorted(glob.glob(input_dict['path_to_folder'] + "*.fit*"))
    cat_files = sorted(glob.glob(input_dict['path_to_folder'] + "*.cat"))
    cs_file_exists = check_if_all_input_variables_filled(**input_dict)
    output_dict, date_ra_dec_table = get_dict_with_date_and_position(fits_images, input_dict)
    output_dict = get_dict_with_objs_info(cat_files, date_ra_dec_table, output_dict, cs_file_exists, **input_dict)
    logging.info('Succesfully created output dictionary.')
    make_output_file(output_dict, **input_dict)


def check_if_all_input_variables_filled(path_to_folder, obj_id, step, date_key_list, objects_key_list, cs_file, **kwargs):
    if not all([path_to_folder, obj_id, step, date_key_list, objects_key_list]):
        logging.critical('One or more input variables empty!')
        raise ValueError('One or more variables are missing')
    elif os.path.exists(path_to_folder + cs_file):
        logging.info('All input variables filled.')
        return True
    else:
        logging.info('All input variables filled, "cs_file" file was not found.')
        return False


def create_dict(dict_, key_list, objects):
    for adj in objects:
        for key in key_list:
            key_ = key + '_' + adj
            dict_[key_] = []
    return dict_


def get_file_name(image):
    pathname = os.path.splitext(image)[0]	
    image_name = pathname.split('/')[-1]
    return image_name


def get_dict_with_date_and_position(fits_images, input_dict):
    start_date, stop_date, images_date_list, julian_date_list = get_dates(fits_images, **input_dict)
    horizons_table = get_horizons_table(start_date, stop_date, **input_dict)
    horizon_iso_dates = get_horizon_iso_dates(horizons_table['datetime_str'])
    closest_dates_list = get_closest_dates_list(images_date_list, horizon_iso_dates)
    date_ra_dec_table = link_date_to_ra_dec(closest_dates_list, horizons_table, horizon_iso_dates)
    output_dict = get_dict_with_dates(images_date_list, julian_date_list, **input_dict)
    return output_dict, date_ra_dec_table


def get_dates(fits_images, date_hdrs_dict, **kwargs):
    images_date_list = []
    julian_date_list = []
    for image in fits_images:
        im_name = get_file_name(image)
        hdul = fits.open(image)
        try:
            date = hdul[0].header[date_hdrs_dict['date']]
        except KeyError:
            logging.critical('Couldn`t find header %s in image %s!', date_hdrs_dict['date'], im_name)
            raise KeyError('Header was not found in the targeted fits image')
        if date and len(date) == 22:
            date = date + '0'
        try:
            julian_date = hdul[0].header[date_hdrs_dict['julian_date']]
        except KeyError:
            logging.critical('Couldn`t find header %s in image %s!', date_hdrs_dict['julian_date'], im_name)
            raise KeyError('Header was not found in the targeted fits image')
        images_date_list.append(date)
        julian_date_list.append(julian_date)
    logging.info('Succesfully extracted dates from fits images.')
    images_date_list = sorted(images_date_list)
    julian_date_list = sorted(julian_date_list)
    start_date, stop_date = get_first_and_last_date(images_date_list)
    return start_date, stop_date, images_date_list, julian_date_list


def get_first_and_last_date(images_date_list):
    date1 = images_date_list[0]
    date2 = images_date_list[-1]
    start_date = date1.replace('T', ' ')
    stop_date = date2.replace('T', ' ')
    return start_date, stop_date


def get_horizons_table(start_date, stop_date, obj_id, step, **kwargs):
    obj = Horizons(id = obj_id,
               epochs={'start' : start_date,
                       'stop' : stop_date,
                       'step' : step})
    eph = obj.ephemerides(quantities=1)
    horizons_table = eph.to_pandas()
    return horizons_table


def get_dict_with_dates(images_date_list, julian_date_list, date_key_list, **kwargs):
    n_idx = np.arange(1, len(images_date_list) + 1, 1)
    list_ = [n_idx, images_date_list, julian_date_list]
    dict_ = {}
    for n, key in enumerate(date_key_list):
        dict_[key] = list_[n]
    return dict_


def change_date_format_str_to_int(date):
    date = date.replace('Jan', '01')
    date = date.replace('Feb', '02')
    date = date.replace('Mar', '03')
    date = date.replace('Apr', '04')
    date = date.replace('Mai', '05')
    date = date.replace('Jun', '06')
    date = date.replace('Jul', '07')
    date = date.replace('Aug', '08')
    date = date.replace('Sep', '09')
    date = date.replace('Oct', '10')
    date = date.replace('Nov', '11')
    date = date.replace('Dec', '12')
    date = date.replace(' ', 'T')
    if len(date) != 24:
        date = date[0:23]
    return date


def change_date_format_int_to_str(date):
    date = list(str(date))
    if date[5] == '0':
        if date[6] == '1':
            date[5:7] = 'Jan'
        elif date[6] == '2':
            date[5:7] = 'Feb'
        elif date[6] == '3':
            date[5:7] = 'Mar'
        elif date[6] == '4':
            date[5:7] = 'Apr'
        elif date[6] == '5':
            date[5:7] = 'Mai'
        elif date[6] == '6':
            date[5:7] = 'Jun'
        elif date[6] == '7':
            date[5:7] = 'Jul'
        elif date[6] == '8':
            date[5:7] = 'Aug'
        elif date[6] == '9':
            date[5:7] = 'Sep'
    elif date[5] == '1':
        if date[6] == '0':
            date[5:7] = 'Oct'
        elif date[6] == '1':
            date[5:7] = 'Nov'
        elif date[6] == '2':
            date[5:7] = 'Dec'
    date = ''.join(date)
    if len(date) != 24:
        date = date[0:24]
    date = date.replace('T', ' ')
    return date


def table_with_iso_date_format(table):
    # print(table)
    for n, date in enumerate(table['datetime_str']):
        date = change_date_format_str_to_int(date)
        print(n, date)
        table.loc[:, ('datetime_str', n)] = date
        # table['datetime_str'][n] = date
    # print(table)
    return table


def get_horizon_iso_dates(horizon_dates):
    horizon_iso_dates = []
    for date in horizon_dates:
        date = change_date_format_str_to_int(date)
        horizon_iso_dates.append(date)
    return horizon_iso_dates


def get_closest_dates_list(images_date_list, horizon_iso_dates):
    closest_dates = []
    for date in images_date_list:
        date_ = find_closest_date(date, horizon_iso_dates)
        date_ = change_date_format_str_to_int(str(date_))
        closest_dates.append(date_)
    return closest_dates


def find_closest_date(date, date_list):
    date_ = datetime.fromisoformat(date)
    date_list_ = []
    for dates in date_list:
        dates_ = datetime.fromisoformat(dates)
        date_list_.append(dates_)
    cloz_dict = { 
        abs(date_.timestamp() - date.timestamp()) : date 
        for date in date_list_}
    return cloz_dict[min(cloz_dict.keys())]


def link_date_to_ra_dec(closest_dates_list, horizons_table, horizon_iso_dates):
    dict_ = {'date_iso' : [], 'RA' : [], 'DEC' : []}
    for date in closest_dates_list:
        idx = 0
        for date_ in horizon_iso_dates:
            if date == date_:
                break
            idx += 1
        dict_['date_iso'].append(horizon_iso_dates[idx])
        dict_['RA'].append(horizons_table['RA'][idx])
        dict_['DEC'].append(horizons_table['DEC'][idx])
    return pd.DataFrame.from_dict(dict_)


def get_dict_with_objs_info(cat_files, date_ra_dec_table, output_dict, cs_file_exists, cs_file, objects_key_list, path_to_folder, **kwargs):
    objects = ['OBJ']
    if cs_file_exists:
        comp_stars_table = pd.read_csv(path_to_folder + cs_file)
        n = comp_stars_table.shape[0]
        for i in range(n):
            objects.append('COMP' + str(i+1))
    for obj in objects:
        for n, cat_file in enumerate(cat_files):
            if obj == 'OBJ':
                RA = date_ra_dec_table['RA'][n]
                DEC = date_ra_dec_table['DEC'][n]
            else:
                RA = comp_stars_table.iloc[0][2]
                DEC = comp_stars_table.iloc[0][3]
            if n == 0:
                obj_info = from_sex_get_obj(RA, DEC, cat_file, obj, objects_key_list)
            else:
                next_obj_info = from_sex_get_obj(RA, DEC, cat_file, obj, objects_key_list)
                obj_info = pd.concat([obj_info, next_obj_info], axis=0)
        output_dict = to_dict_add_obj_keys(obj, obj_info, output_dict, objects_key_list)
    return output_dict


def from_sex_get_obj(RA, DEC, cat_file, obj, key_list):
    sex_catalog = read_sex_file(cat_file)
    idx = find_closest_object(RA, DEC, sex_catalog, k=Angle('0d00m10s'))
    obj_info = sex_catalog.iloc[[idx]]
    return obj_info


def read_sex_file(file_name):
    data = ascii.read(file_name, format='sextractor')
    sex_catalog = data.to_pandas()
    return sex_catalog 


def find_closest_object(RA, DEC, sex_catalog, k):
    c = SkyCoord(ra=RA*u.degree, dec=DEC*u.degree)
    min_distance = np.inf
    # k = Angle('0d00m10s')
    for i in range(len(sex_catalog['NUMBER'])):
        ra2 = np.array([sex_catalog['ALPHA_SKY'][i]])
        dec2 = np.array([sex_catalog['DELTA_SKY'][i]])
        catalog = SkyCoord(ra=ra2*u.degree, dec=dec2*u.degree)
        _, d2d, __ = c.match_to_catalog_sky(catalog)
        if d2d < min_distance:
            min_distance = d2d
            idx_ = i
    return idx_


def to_dict_add_obj_keys(obj, obj_info, output_dict, objects_key_list):
    used_sex_catalog_keys = ['ALPHA_SKY', 'DELTA_SKY', 'X_IMAGE', 'Y_IMAGE', 'FLUX_ISO',
                             'FLUXERR_ISO', 'A_IMAGE', 'B_IMAGE', 'THETA_IMAGE']
    for n, key in enumerate(objects_key_list):
        list_ = obj_info[used_sex_catalog_keys[n]].tolist()
        to_join = [key, '_', obj]
        dict_key = ''.join(to_join)
        output_dict[dict_key] = list_
    return output_dict


def make_output_file(dict_, path_to_folder, **input_dict):
    file_ = path_to_folder + 'output_file.csv'
    s = ', '
    s = s.join(list(dict_.keys()))
    sp = '\n'
    file_headers = s + sp
    output_file = open(file_, 'w')
    output_file.write(file_headers)
    for i in range(len(dict_['IDX'])):
        file_row_list = []
        for key in list(dict_.keys()):
            file_row_list.append(str(dict_[key][i]))
        s = ', '
        s = s.join(file_row_list)
        file_line = s + sp
        output_file.write(file_line)
    output_file.close()
    if os.path.exists(file_):
        logging.info('Succesfully created output file.')
    else:
        logging.warning('Output file was not created.')





# 'path_to_folder' : '/home/kamilraczka/Projects/Fits_files/'
# 'path_to_folder' : '/home/kamilraczka/Projects/Na_zaliczenie/'
# 'path_to_folder' : '/home/kamil/Astronomia/Fits_files/'
input_dict = {'path_to_folder' : '/home/kamil/Programs/Analyze-fits/Fits_files/',
'obj_id' : 'Delia',
'step' : '1m',
# 'date_hdrs_dict': {'date': 'DATE-OBS', 'julian_date': 'JD'},
'date_hdrs_dict': {'date': 'DATE-OBS', 'julian_date': 'JD_UTC'},
'date_key_list' : ['IDX', 'DT', 'JD'],
'objects_key_list' : ['RA', 'DEC', 'X', 'Y', 'FLUX', 'FLUX_ERR', 'A_IMAGE', 'B_IMAGE', 'THETA_IMAGE'],
'cs_file' : 'comp_stars.csv'
}

if __name__ == '__main__':
    main(input_dict)
