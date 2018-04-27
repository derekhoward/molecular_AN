import platform
import os

if platform.node() == 'RES-C02TQ1V6G8WL.local':
	DATA_DIR = "/Users/derek_howard/projects/molecular_AN/data"
	CORES = 3
else :
	raise RuntimeError('ERROR: this computer ('+platform.node()+') is not configured. Please change this in config.py')
