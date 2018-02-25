import numpy as np
import scipy # use numpy if scipy unavailable
import scipy.linalg # use numpy if scipy unavailable
import random
import math
import time
import numpy as np
import matplotlib
import matplotlib.animation as animation
matplotlib.use('TkAgg') # do this before importing pylab
import matplotlib.pyplot as plt 
import json
from pprint import pprint

def import_time_data(name = '2men.txt'):
	ids = {} #{id: timestamp}
	with open(name, "r") as f:
		for line in f:
			words = line.split(',')
			i = 0
			while i < (len(words)-1):
				ids[words[i]] = words[i+1]
				i += 2
	return ids
