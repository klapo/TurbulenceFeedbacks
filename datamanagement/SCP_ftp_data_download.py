####################################################################################################
# SCP_ftp_data_download 
# Karl Lapo November/2015
####################################################################################################
# Grab netcdf files from the SCP
####################################################################################################

## Import statements
# OS interaction
import sys, pickle, os
import pandas as pd

# import subplots function for plotting
import seaborn as sns
import matplotlib
from matplotlib.pyplot import subplots
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.basemap import Basemap

## Directory listing
dir_data = '/Users/karllapo/gdrive/SnowHydrology/proj/TurbulenceFeedbacks/data/SCP'

os.chdir(dir_data)
data_to_ftp = pd.read_csv('SCP_ftp_filelist.txt',sep='\t',header=None)
print(data_to_ftp)
for ftp_address in data_to_ftp.values:
	print(ftp_address[0])
	os.system('wget '+ftp_address[0])
