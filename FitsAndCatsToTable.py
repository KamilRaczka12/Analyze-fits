from importlib.resources import path
# from inspect import _Object
from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord, Angle
from astropy import units as u
from astroquery.jplhorizons import Horizons
from datetime import datetime
import glob
import numpy as np
import os
import pandas as pd
from photutils.datasets import make_100gaussians_image
from photutils.aperture import CircularAperture


def do_the_work(path_to_folder, obj_id, step, comp_stars_file, key_list):
    dict_ = create_dict(key_list)
    analysed_objects = ['OBJ']
    if os.path.exists(path_to_folder + comp_stars_file):
        dict_new, n_of_comp_stars, comp_stars_table = to_dict_add_comp_star_keys(dict_, comp_stars_file, key_list)
        dict_ = dict_new
        analysed_objects = list_of_objects(analysed_objects, n_of_comp_stars)
    fits_images = sorted(glob.glob(path_to_folder + "*.fit*"))
    cat_file_list = sorted(glob.glob(path_to_folder + "*.cat"))
    images_date_list, start_date, stop_date, julian_date_list = get_dates(fits_images)
    dict_ = fill_dict_with_dates(images_date_list, julian_date_list, dict_)
    eph = give_horizons_table(obj_id, start_date, stop_date, step)
    horizons_table = eph.to_pandas()
    horizon_iso_dates = give_horizon_iso_dates(horizons_table['datetime_str'])
    closest_dates_list = give_closest_dates_list(images_date_list, horizon_iso_dates)
    date_ra_dec_table = link_date_to_ra_dec(closest_dates_list, horizons_table, horizon_iso_dates)
    dict_ = return_filled_dict(analysed_objects, cat_file_list, date_ra_dec_table, comp_stars_table, dict_)
    make_output_file(path_to_folder, dict_)
    


def create_dict(key_list):
    dict_ = {'IDX' : [], 'DT' : [], 'JD' : []}
    for key in key_list:
        key_ = key + '_OBJ'
        dict_[key_] = []
    return dict_


def to_dict_add_comp_star_keys(dict_, comp_stars_file, key_list):
    comp_stars_table = pd.read_csv(path_to_folder + comp_stars_file)
    n_of_comp_stars = len(comp_stars_table.index)
    for n in range(n_of_comp_stars):
        new_keys = []
        for key in key_list:
            new_keys.append(key + '_COMP' + str(n + 1))
        for key in new_keys:
            dict_[key] = []
    return dict_, n_of_comp_stars, comp_stars_table


def list_of_objects(analysed_objects, n_of_comp_stars):
    for i in range(n_of_comp_stars):
        analysed_objects.append('COMP' + str(i + 1))
    return analysed_objects


def get_dates(fits_images):
    images_date_list = []
    julian_date_list = []
    for image in fits_images:
        hdul = fits.open(image)
        date = hdul[0].header['DATE-OBS']
        julian_date = hdul[0].header['JD']
        if len(date) == 22:
            date = date + '0'
        images_date_list.append(date)
        julian_date_list.append(julian_date)
    date1 = images_date_list[0]
    date2 = images_date_list[-1]
    start_date = date1.replace('T', ' ')
    stop_date = date2.replace('T', ' ')
    # head1, sep1, tail1 = date1.partition('T')
    # head2, sep2, tail2 = date2.partition('T')
    # start_date = head1
    # stop_date = head2
    # if start_date == stop_date:
    #     f"{start_date} [ {tail1} ]"
    #     start_date = start_date + ' [' + tail1 + ']'
    #     stop_date = stop_date + ' [' + tail2 + ']'
    return images_date_list, start_date, stop_date, julian_date_list


def fill_dict_with_dates(images_date_list, julian_date_list, dict_):
    for n, date in enumerate(images_date_list):
        dict_['IDX'].append(str(n + 1))
        dict_['DT'].append(images_date_list[n])
        dict_['JD'].append(julian_date_list[n])
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


def give_horizons_table(obj_id, start_date, stop_date, step):
    obj = Horizons(id = obj_id,
               epochs={'start' : start_date,
                       'stop' : stop_date,
                       'step' : step})
    return obj.ephemerides(quantities=1)


def table_with_iso_date_format(table):
    # print(table)
    for n, date in enumerate(table['datetime_str']):
        date = change_date_format_str_to_int(date)
        print(n, date)
        table.loc[:, ('datetime_str', n)] = date
        # table['datetime_str'][n] = date
    # print(table)
    return table


def give_horizon_iso_dates(horizon_dates):
    horizon_iso_dates = []
    for date in horizon_dates:
        date = change_date_format_str_to_int(date)
        horizon_iso_dates.append(date)
    return horizon_iso_dates


def give_closest_dates_list(images_date_list, horizon_iso_dates):
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
            # print(date, date_)
            if date == date_:
                break
            idx += 1
        dict_['date_iso'].append(horizon_iso_dates[idx])
        # dict_['date_iso'] = table['datetime_str'][idx]
        dict_['RA'].append(horizons_table['RA'][idx])
        dict_['DEC'].append(horizons_table['DEC'][idx])
    return pd.DataFrame.from_dict(dict_)
    

def read_sextractor_files(file_name):
    sex_catalog = pd.read_fwf(file_name,
                header=None, skiprows=range(0,12),
                names=['NUMBER', 'X_IMAGE', 'Y_IMAGE', 'ALPHA_SKY', 'DELTA_SKY',
                'MAG_ISO', 'MAGERR_ISO', 'FLUX_ISO', 'FLUXERR_ISO',
                'A_IMAGE', 'B_IMAGE', 'THETA_IMAGE'])
    return sex_catalog 


def find_closest_object(RA, DEC, sex_catalog):
    c = SkyCoord(ra=RA*u.degree, dec=DEC*u.degree)
    min_distance = np.inf
    k = Angle('0d00m10s')
    for i in range(len(sex_catalog['NUMBER'])):
        ra2 = np.array([sex_catalog['ALPHA_SKY'][i]])
        dec2 = np.array([sex_catalog['DELTA_SKY'][i]])
        catalog = SkyCoord(ra=ra2*u.degree, dec=dec2*u.degree)
        _, d2d, __ = c.match_to_catalog_sky(catalog)
        if d2d < min_distance:
            min_distance = d2d
            idx_ = i
    return idx_


def from_sex_get_obj(RA, DEC, cat_file, dict_, obj):
    sex_catalog = read_sextractor_files(cat_file)
    idx = find_closest_object(RA, DEC, sex_catalog)
    obj_info = sex_catalog.iloc[[idx]]
    obj_dict = append_to_dict(obj_info, dict_, obj)
    return obj_dict


def append_to_dict(obj_info, dict_, obj):
    dict_['RA' + '_' + obj].append(obj_info.iloc[0]['ALPHA_SKY'])
    dict_['DEC' + '_' + obj].append(obj_info.iloc[0]['DELTA_SKY'])
    dict_['X' + '_' + obj].append(obj_info.iloc[0]['X_IMAGE'])
    dict_['Y' + '_' + obj].append(obj_info.iloc[0]['Y_IMAGE'])
    dict_['FLUX' + '_' + obj].append(obj_info.iloc[0]['FLUX_ISO'])
    dict_['FLUX_ERR' + '_' + obj].append(obj_info.iloc[0]['FLUXERR_ISO'])
    dict_['A_IMAGE' + '_' + obj].append(obj_info.iloc[0]['A_IMAGE'])
    dict_['B_IMAGE' + '_' + obj].append(obj_info.iloc[0]['B_IMAGE'])
    dict_['THETA_IMAGE' + '_' + obj].append(obj_info.iloc[0]['THETA_IMAGE'])
    return dict_


def return_filled_dict(analysed_objects, cat_file_list, date_ra_dec_table, comp_stars_table, dict_):
    for obj in analysed_objects:
        if obj == 'OBJ':
            for n, cat_file in enumerate(cat_file_list):
                RA = date_ra_dec_table['RA'][n]
                DEC = date_ra_dec_table['DEC'][n]
                dict_new = from_sex_get_obj(RA, DEC, cat_file, dict_, obj)
                dict_ = dict_new
        else:
            for n, cat_file in enumerate(cat_file_list):
                RA = comp_stars_table.loc[0][2]
                DEC = comp_stars_table.iloc[0][3]
                dict_new = from_sex_get_obj(RA, DEC, cat_file, dict_, obj)
                dict_ = dict_new
    return dict_



def make_output_file(path_to_folder, dict_):
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





# path_to_folder = '/home/kamilraczka/Projects/Fits_files/'
path_to_folder = '/home/kamilraczka/Projects/Na_zaliczenie/'
# path_to_folder = '/home/kamil/Astronomia/Fits_files/'
obj_id = 'Gunlod'
step = '1m'
key_list = ['RA', 'DEC', 'X', 'Y', 'FLUX', 'FLUX_ERR', 'A_IMAGE', 'B_IMAGE', 'THETA_IMAGE']
comp_stars_file = 'comp-stars.csv'
do_the_work(path_to_folder, obj_id, step, comp_stars_file, key_list)
