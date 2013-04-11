#!/bin/bash 

# folder wich contains the script
folder="$(dirname $0)" || exit 0 # abort on failure
 
cd "$folder" || exit 0 # abort on failure
rm -r ./html
doxygen "ogdf-doxygen.cfg"
