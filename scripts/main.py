#!/usr/bin/python2

# Main Python Program which will be called by the controlling bash script '/start.sh'
# Inputs the PDB file name and location as a command line argument and creates another file 
# where the processed data is stored

# importing necessary modules
import sys


def main():
	filename = sys.argv[1]
	try:
		infilestream = open(filename, 'rU')		# open input pdb file under try block to tackle exceptions
	except IOError:
		print 'Error: Input PDB file was not found!'
		return 0

#	try:
#		outfilestream = open('', 'w')			# open the coordinate file under try block to handle exceptions
#	except IOError:
#		print 'Error: Coordinate file could not be opened'
#		return 0

	for line in infilestream.readlines():		#processing file line by line
		words = line.split()
		print words				#test print

# standard autocall construct for main()
if __name__ == '__main__': main()
