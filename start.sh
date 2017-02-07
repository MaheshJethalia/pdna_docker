echo "Welcome to the KGP Protein Dna Docker!"
echo "Now running the python script main.py for processing coordinates....."

cd ./scripts
./main.py test.pdb
cd ..

echo "Finished processing the coordinates."
echo "End of job."
