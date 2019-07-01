import pytest
import numpy as np
from numpy import testing as nptest
import datetime,os
from ssj_auroral_boundary import files

def test_test_data_exists():
    test_data_dir, test_cdffn = files.test_cdf_path_and_filename()
    assert os.path.exists(os.path.join(test_data_dir,test_cdffn))

