#! /projects/weinstocklab/comp/local/anaconda2/bin/python

import pandas as pd
from Bio import SeqIO
import os
import re
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import sys
from scipy.special import expit
import argparse

'''
print("The required arguments are positially as follows: 1 infile, 2 outfile, 3 original fasta file.")
print("The optional arguments are as follows: percenterrorcutoff (abbrv -PEC) set to 1% default, percentile (abbrv -P) set to 75th percentile, , upperCCSlimit (abbrv -UL) set to 75 CCS passes, datausedformodel (abbrv -DUFM) set to the selected percentile of the data.")
'''

parser = argparse.ArgumentParser(description='Allows for the user to select the input file, upper CCS pass number limit, percentile of the data examined, percent error cutoff, output file path, and the original fasta file')
parser.add_argument('infile', type=str, help='Input file path - should be a .tsv file')
parser.add_argument('outfile', type=str, help='Output file path - should be a .fasta file')
parser.add_argument('originalfile', type=str, help='Original .fasta file that will be parsed over to remove any sequences that do not meet the user selected upper and lower CCS pass cutoffs')
parser.add_argument('-PEC','--percenterrorcutoff', type=float, default=1, help='provide an integer to be the highest percent error on the model for the nth percentile of the data (default: 1)')
parser.add_argument('-P', '--percentile', type=int, default=75, help='provide an integer for the percentile of the data (default: 75)')
parser.add_argument('-UL', '--upperCCSlimit', type=int, default=75, help='provide an integer for as the upper CCS pass limit (default: 75)')
parser.add_argument('-DUFM', '--datausedformodel', type=str, default='nthpercent', choices=['nthpercent', 'datamean'], help='select either the nth percent of the data or mean to be analyzed for the model (default: nth percent)')

arguments = parser.parse_args()




'''
# total arguments 
n = len(sys.argv) 
print("Total arguments passed:", n) 
  
# Arguments passed 
print("\nName of Python script:", sys.argv[0]) 
  
print("\nArguments passed:", end = " ") 
for i in range(1, n): 
    print(sys.argv[i], end = " ") 
'''







'''Define the aggregate functions'''

def percentiledata(x):
	return np.percentile(x,arguments.percentile)

def SD_upper_percentile(x):
	return percentiledata(x) + (np.std(x))

def SD_lower_percentile(x):
        return percentiledata(x) - (np.std(x))


funcs = [np.mean, np.median, percentiledata, SD_upper_percentile, SD_lower_percentile]
agg_funcs = {'INS':funcs, 'DEL':funcs, 'SNP':funcs, 'TOTAL':funcs}


Error = pd.read_csv(arguments.infile, sep='\t')


Error1 = Error
Error1 = Error1.groupby(Error1['CPN']).agg(agg_funcs)
Error1.columns = ['_'.join(x) for x in Error1.columns.get_values()]
Error1.reset_index(inplace=True)


'''This number should be modified to a custom setpoint for the upper limit of CCS passes'''
Error2 = Error1[Error1['CPN'] <= arguments.upperCCSlimit]


'''Defining the exponential decay function with np.exp being equivalent to the function e^(-bx)'''
def func(x, a, b, c):
        return a*expit(-b*x) + c



'''setting the x and y data as Number of CCS passes to total error rate of the nth percentile of the data or mean (default set to 75th percentile of the data)'''
xdata=np.array(Error2['CPN'])



if arguments.datausedformodel == "nthpercent":
	ydata1 = np.array(Error2['TOTAL_percentiledata'])
	var1, pcov1 = curve_fit(func, xdata, ydata1)
	
	plt.plot(xdata, ydata1, '-b', label='data')
	plt.plot(xdata, func(xdata, *var1), 'g--',
		label='fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(var1))
	plt.legend()
	
elif arguments.datausedformodel == "datamean":
	ydata2 = np.array(Error2['TOTAL_mean'])
	var1, pcov1 = curve_fit(func, xdata, ydata2)

	plt.plot(xdata, ydata2, '-b', label='data')
	plt.plot(xdata, func(xdata, *var1), 'g--',
        	label='fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(var1))
	plt.legend()


else:
	print("Please select either nthpercent or datamean as the datausedformodel argument")



'''
Defining variables for the sample


var1, pcov1 = curve_fit(func, xdata, arguments.datausedformodel)




Modeling the CCSP pass number vs error rate for sample


plt.plot(xdata, arguments.datausedformodel, '-b', label='data')
plt.plot(xdata, func(xdata, *var1), 'g--',
	label='fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(var1))
plt.legend()
'''



'''Reorganize exponential decay function to solve for x, np.log is ln'''
def revfunc(y,a,b,c):
	return np.log((y-c)/a)/(-b)




'''Determine CCS pass cutoff for sample to get mean of total error equal to or below certain % cutoff that should be modifiable
The default is 1%'''
CCScutoff = revfunc(arguments.percenterrorcutoff, *var1)
CCScutoff = math.ceil(CCScutoff)
print("The selected lower CCS Cutoff number is:", CCScutoff, "The selected upper CCS Cutoff number is:", arguments.upperCCSlimit)


'''Removes sequences outside upper and lower bounds of CCS pass number'''
CCSPassRemoval = open(arguments.outfile,"w")
for sequence in SeqIO.parse(arguments.originalfile, "fasta"):
	CCSPass = int(sequence.id.split("np:")[1])
	stringseq = str(sequence.seq)
	if CCSPass >= CCScutoff and CCSPass <= arguments.upperCCSlimit:
		CCSPassRemoval.write('>' + sequence.id +'\n' + stringseq + '\n')
	else:
		continue
CCSPassRemoval.close()


origfile = arguments.originalfile
originalcount = 0
with open(origfile, 'r') as originalfasta:
    for line in originalfasta:
        originalcount += 1


finalfile = arguments.outfile
finalcount = 0
with open(finalfile, 'r') as finalfasta:
    for line in finalfasta:
        finalcount += 1
    if finalcount > 0:
        print("The total number of sequences retained is:", int(finalcount/2))
        print("The percentage of reads retained is:", ((finalcount/2)/((originalcount/2))*100))
    else:
        print("The selected CCS pass cutoff removed all sequences. Please select a different percentile of the data or percent error cutoff.")
