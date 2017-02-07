#!/usr/bin/python2

# Main Python Program which will be called by the controlling bash script '/start.sh'
# Inputs the PDB file name and location as a command line argument and creates another file 
# where the processed data is stored

# importing necessary modules
import sys
import re
import os

def coordFind(line):
	""" Returns true if the line contains atomic coordinates and false otherwise."""
	match1 = re.search('ATOM', line)
	match2 = re.search('HETATM', line)
	if match1 or match2:
		return True
	else:
		return False


def getCoords(line):
	"""Extracts and returns the coordinates of the atom / heteroatom in the line. This function assumes that the 
	coordinates exist and have been checked for before (check coordFind())"""
	words = line.split()
	coordx = words[5]
	coordy = words[6]
	coordz = words[7]
	return (coordx, coordy, coordz)


def main():
	infilename = sys.argv[1]                                 
	
	scriptPath = os.path.dirname(__file__)                 	# obtain directory name of current script in system
	scriptPath = os.path.abspath(scriptPath)				# obtain full path to the above obtained directory name

	try:
		infilestream = open(infilename, 'rU')				# open input pdb file under try block to tackle exceptions
	except IOError:
		print 'Error: Input PDB file was not found!'
		return 0

	outfilename = os.path.join(scriptPath, 'coords.tmp')	# create output coordinate file coords.tmp using absolute path
 	try:
 		outfilestream = open(outfilename, 'w')				# open the coordinate file under try block to handle exceptions
	except IOError:
		print 'Error: Coordinate file could not be opened'
		return 0

	for line in infilestream.readlines():					# processing file line by line
		if coordFind(line):									# check if line contains coordinates
			(coordx, coordy, coordz) = getCoords(line)    			# obtain cordinates of the atom in line
			outfilestream.write('{} {} {}\n'.format(coordx, coordy, coordz))			# write coordinates to output file
	infilestream.close()
	outfilestream.close()


# standard autocall construct for main()
if __name__ == '__main__': main()
