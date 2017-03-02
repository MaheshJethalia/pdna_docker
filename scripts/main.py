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
	infilename1 = sys.argv[1]                               # obtaining the pdb file names for the protein(1) and dna(2)		  
	infilename2 = sys.argv[2]
	scriptPath = os.path.dirname(__file__)                 	# obtain directory name of current script in system
	scriptPath = os.path.abspath(scriptPath)				# obtain full path to the above obtained directory name

	try:
		infilestream1 = open(infilename1, 'rU')				# open input pdb file under try block to tackle exceptions
		infilestream2 = open(infilename2, 'rU')
	except IOError:
		print 'Error: Input PDB file was not found!'
		return 0

	outfilename1 = os.path.join(scriptPath, 'coords1.tmp')	# create output coordinate file coords1.tmp and coords2.tmp using absolute path
 	outfilename2 = os.path.join(scriptPath, 'coords2.tmp')
 	try:
 		outfilestream1 = open(outfilename1, 'w')			# open the coordinate files under try block to handle exceptions
		outfilestream2 = open(outfilename2, 'w')
	except IOError:
		print 'Error: Coordinate file could not be opened'
		return 0

	for line in infilestream1.readlines():					# processing infilename1 line by line
		if coordFind(line):									# check if line contains coordinates
			(coordx, coordy, coordz) = getCoords(line)    	# obtain cordinates of the atom in line
			outfilestream1.write('{} {} {}\n'.format(coordx, coordy, coordz))			# write coordinates to output file

	for line in infilestream2.readlines():					# processing infilename2 line by line
		if coordFind(line):									# check if line contains coordinates
			(coordx, coordy, coordz) = getCoords(line)    	# obtain cordinates of the atom in line
			outfilestream2.write('{} {} {}\n'.format(coordx, coordy, coordz))			# write coordinates to output file
	
	infilestream1.close()									# close all used file objects
	infilestream2.close()
	outfilestream1.close()
	outfilestream2.close()


# standard autocall construct for main()
if __name__ == '__main__': main()