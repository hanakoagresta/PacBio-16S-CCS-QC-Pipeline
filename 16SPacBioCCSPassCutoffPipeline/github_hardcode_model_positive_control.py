#! /projects/weinstocklab/comp/local/anaconda2/bin/python

from rpy2.robjects import pandas2ri
pandas2ri.activate()
import pandas as pd


from Bio import SeqIO
import os
import re
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import rpy2
import rpy2.robjects as robjects
import sys

'''
input_file=sys.argv[1]
upperlimitpassnum=sys.argv[2]
percentile=sys.argv[3]
percent_error_cutoff=sys.argv[4]
output_file_path=sys.argv[5]
original_fasta_file_path=sys.argv[6]
'''




'''Define the aggregate functions'''

def percentiledata(x):
	return np.percentile(x,[SELECT PERCENTILE OF DATA])

def SD_upper_percentile(x):
	return percentiledata(x) + (np.std(x))

def SD_lower_percentile(x):
        return percentiledata(x) - (np.std(x))


funcs = [np.mean, np.median, percentiledata, SD_upper_percentile, SD_lower_percentile]
agg_funcs = {'INS':funcs, 'DEL':funcs, 'SNP':funcs, 'TOTAL':funcs}


Error = pd.read_csv([SELECT INPUT FILE], sep='\t')


Error1 = Error
Error1 = Error1.groupby(Error1['CPN']).agg(agg_funcs)
Error1.columns = ['_'.join(x) for x in Error1.columns.get_values()]
Error1.reset_index(inplace=True)


'''This number should be modified to a custom setpoint for the upper limit of CCS passes'''
Error2 = Error1[Error1['CPN'] <= [SELECT AN UPPER LIMIT FOR PASS NUMBER]]




'''setting the x and y data as Number of CCS passes to total error rate of the 95th percentile of the data'''
xdata=np.array(Error2['CPN'])
ydata=np.array(Error2['TOTAL_percentiledata'])


'''Defining the exponential decay function with np.exp being equivalent to the function e^(-bx)'''
def func(x, a, b, c):
	return a*np.exp(-b*x) + c



'''Defining variables for the sample'''
var1, pcov1 = curve_fit(func, xdata, ydata)





'''Modeling the CCSP pass number vs error rate for sample'''


plt.plot(xdata, ydata, '-b', label='data')
plt.plot(xdata, func(xdata, *var1), 'g--',
	label='fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(var1))
plt.legend()



'''Reorganize exponential decay function to solve for x, np.log is ln'''
def revfunc(y,a,b,c):
	return np.log((y-c)/a)/(-b)




'''Determine CCS pass cutoff for sample to get 95th percentile of total error equal to or below certain % cutoff that should be modifiable
It is currently set to 1%'''
CCScutoff = revfunc([SELECT A PERCENT ERROR CUTOFF], *var1)
CCScutoff = math.ceil(CCScutoff)
print(CCScutoff)


'''Removes sequences outside upper and lower bounds of CCS pass number'''
CCSPassRemoval = open([SELECT OUTPUT FILE PATH],"w")
for sequence in SeqIO.parse([SELECT ORIGINAL FASTA FILE PATH], "fasta"):
	CCSPass = int(sequence.id.split("np:")[1])
	stringseq = str(sequence.seq)
	if CCSPass >= CCScutoff and CCSPass <= upperlimitpassnum:
		CCSPassRemoval.write('>' + sequence.id +'\n' + stringseq + '\n')
	else:
		continue
CCSPassRemoval.close()