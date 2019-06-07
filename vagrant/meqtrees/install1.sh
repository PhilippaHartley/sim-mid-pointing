# This script is run as root.

# Add KERN package repository, then install MeqTrees, CASA and pip
# (plus dependencies).

add-apt-repository -y -s ppa:kernsuite/kern-5
apt-get update
apt-get install -y meqtrees casalite python-pip
