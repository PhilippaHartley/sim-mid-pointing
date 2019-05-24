# This script is run as the vagrant user.

# Install h5py using pip (the Ubuntu 18.04 package is too old for our
# purposes). Note this installs it in /home/vagrant/.local

pip install h5py

# Install analysis_scripts package.

curl -O ftp://ftp.cv.nrao.edu/pub/casaguides/analysis_scripts.tar
tar -xf analysis_scripts.tar
rm analysis_scripts.tar

# Clone sim-mid-pointing repository.

git clone https://github.com/ska-telescope/sim-mid-pointing.git
