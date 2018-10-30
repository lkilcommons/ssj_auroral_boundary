# Copyright 2018 SEDA Group at CU Boulder
# Created by: 
# Liam Kilcommons 
# Space Environment Data Analysis Group (SEDA)
# Colorado Center for Astrodynamics Research (CCAR)
# University of Colorado, Boulder (CU Boulder)
import sys,datetime,os,requests

def test_cdf_path_and_filename():
        """ Find location of test file

        Returns
        -------
        test_data_dir : string
            Test data directory
        test_cdffn : string
            Test data filename

        """

	#Determine where this module's source file is located
	#to determine where to look for the test data
	src_file_dir = os.path.dirname(os.path.realpath(__file__))
	test_data_dir = os.path.join(src_file_dir, 'test_data')
	test_cdffn = 'dmsp-f16_ssj_precipitating-electrons-ions_20100529_v1.1.2.cdf'
	return test_data_dir, test_cdffn

def cdf_url_and_filename(dmsp_number, year, month, day):
        """ Download a DMSP file from NOAA

        Parameters
        ----------
        dmsp_number : int
            DMSP satellite number (e.g. 15)
        year : int
            Data year
        month : int
            Data month (1-12)
        day : int
            Data day of month

        Returns
        -------
        cdf_url : string
            URL to retrieve file at
        cdffn : string
            CDF filename

        """

	cdffn = ('dmsp-f%.2d' % (dmsp_number)
	 		+'_ssj_precipitating-electrons-ions_'
			+'%d%.2d%.2d_v1.1.2.cdf' % (year,month,day))

	root_url = 'https://satdat.ngdc.noaa.gov/dmsp/data/'
	one_month_ssj_data_url = 'f%.2d/ssj/%d/%.2d/' % (dmsp_number, year,
                                                         month)
	#Expected CDF file
	cdf_url = root_url + one_month_ssj_data_url + cdffn

	return cdf_url, cdffn

def download_cdf_from_noaa(cdf_url, destination_cdffn):
        """ Download CDF from NOAA

        Parameters
        ----------
        cdf_url : string
            URL (see cdf_url_and_filename) for remote file
        destination_cdffn : string
            Local name with directory destination for downloaded file

        Returns
        -------
        Void

        """

	head = requests.head(cdf_url,allow_redirects=True)
	headers = head.headers
	content_type = headers.get('content-type')
	print(headers,content_type)
	if 'html' not in content_type.lower():
		response = requests.get(cdf_url,allow_redirects=True)
		with open(destination_cdffn,'w') as f:
			f.write(response.content)
	else:
		raise RuntimeError('URL %s is not downloadable' % (cdf_url))
