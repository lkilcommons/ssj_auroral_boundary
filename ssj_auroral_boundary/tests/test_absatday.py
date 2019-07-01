import pytest
import numpy as np
from numpy import testing as nptest
import datetime,os,tempfile,shutil
from ssj_auroral_boundary.absatday import absatday
from ssj_auroral_boundary.files import test_cdf_path_and_filename

@pytest.fixture(scope='module')
def bundled_cdf_absatday(request):

    test_cdfdir,test_cdffn = test_cdf_path_and_filename()
    tmp_cdfdir = tempfile.mkdtemp()
    shutil.copyfile(
                    os.path.join(test_cdfdir,test_cdffn),
                    os.path.join(tmp_cdfdir,test_cdffn)
                    )
    cdffn = os.path.join(tmp_cdfdir,test_cdffn)
    absd = absatday(cdffn,
                    imgdir=tmp_cdfdir,
                    make_plot=False,
                    plot_failed=False,
                    csvdir=tmp_cdfdir,
                    writecsv=True,
                    csvvars=['mlat', 'mlt'])

    def remove_temporary_dir():
        shutil.rmtree(tmp_cdfdir)

    request.addfinalizer(remove_temporary_dir)

    return absd

def test_absatday_getitem_gets_date(bundled_cdf_absatday):
    #Getitem interface should allow getting the date
    #the bundled CDF is for 2010/5/29
    absd = bundled_cdf_absatday
    dts = absd['Epoch']
    dt = dts[0]
    assert dt.year==2010 and dt.month==5 and dt.day==29

def test_absatday_has_at_least_one_failed_identification(bundled_cdf_absatday):
    #We should fail at least one identification, causing a failure_reason
    #to be something other than None
    absd = bundled_cdf_absatday
    failure_reasons = [abpolarpass.failure_reason for abpolarpass in absd.polarpasses]
    assert any([reason is not None for reason in failure_reasons])
