# This script is run as the vagrant user.

# Clone MeqTrees repositories.

mkdir meqtrees
cd meqtrees
git clone https://github.com/ska-sa/kittens 
git clone https://github.com/ska-sa/purr 
git clone https://github.com/ska-sa/pyxis 
git clone https://github.com/ska-sa/tigger 
git clone https://github.com/ska-sa/meqtrees-timba 
git clone https://github.com/ska-sa/meqtrees-cattery 
git clone https://github.com/ska-sa/owlcat

# Compile meqtrees-timba.

cd meqtrees-timba
mkdir build
Tools/Build/bootstrap_cmake release
cd build/release 
make

# Install h5py using pip (the Ubuntu 18.04 package is too old for our
# purposes). Note this installs it in /home/vagrant/.local

pip install h5py

# Install analysis_scripts package.

cd
curl -O ftp://ftp.cv.nrao.edu/pub/casaguides/analysis_scripts.tar
tar -xf analysis_scripts.tar
rm analysis_scripts.tar

# Clone sim-mid-pointing repository.

git clone https://github.com/ska-telescope/sim-mid-pointing.git
