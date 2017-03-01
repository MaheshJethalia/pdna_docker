echo "Welcome to the KGP Protein Dna Docker!"

# For running 
echo "Now running the python script main.py for processing coordinates....."

#cd ./scripts
#./main.py test.pdb test.pdb
#cd ..

echo "Finished processing the coordinates."

# For c++ code compilation and running
echo "Building c++ code now..."
cd src
make
echo "Build complete! Now running the docker..."
./main
make clean
echo "Finished running the docker"

echo "End of job."
