# Time-stamp: <2021-12-20 09:19:28 ta32852>

"""Description: Main HMMR command

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""

# ------------------------------------
# python modules
# ------------------------------------

import os
import sys
import logging
from time import strftime

# ------------------------------------
# own python modules
# ------------------------------------
from MACS3.Utilities.Constants import *
from MACS3.Utilities.OptValidator import opt_validate_hmmratac

#from MACS3.Math.EM import HMMR_EM


# ------------------------------------
# constants
# ------------------------------------

# ------------------------------------
# Misc functions
# ------------------------------------

# ------------------------------------
# Main function
# ------------------------------------
def run( args ):
    """The HMMRATAC function/pipeline for MACS.

    """
    options = opt_validate_hmmratac( args )
    print ( options )
    return
