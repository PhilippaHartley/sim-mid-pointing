# This script is run as root.

# Add KERN package repository, then install MeqTrees build
# dependencies, CASA and pip.

add-apt-repository -y -s ppa:kernsuite/kern-5
apt-get update
sudo apt-get build-dep -y kittens purr pyxis tigger meqtrees-timba meqtrees-cattery owlcat
sudo apt-get install -y casalite python-casacore python-pil python-pip
