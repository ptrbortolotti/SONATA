#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  2 10:03:17 2019
Hides the DEBUG MEssages! 
@author: gu32kij
"""

# Core Library modules
import logging

# matplotlib logger
mpl_logger = logging.getLogger("matplotlib")
mpl_logger.setLevel(logging.WARNING)

shp_logger = logging.getLogger("shapely")
shp_logger.setLevel(logging.WARNING)
