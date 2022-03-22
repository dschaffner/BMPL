# -*- coding: utf-8 -*-
"""
Created on Sun Feb 13 14:37:12 2022

@author: dschaffner
"""

import numpy as np
import pandas as pd
sheetdirectory='C:\\Users\\dschaffner\\Dropbox\\Data\\BMPL\\BMX\\2022\\01122022\\'
spreadsheetfile = 'Separation Times by Shot.xlsx'
sheet = pd.read_excel(sheetdirectory+spreadsheetfile,header=1)
vels57 = np.array(sheet['57v'])
vels1921 = np.array(sheet['1933v'])
vels3335 = np.array(sheet['3335v'])
velsflags = np.array(sheet['flag'])



