# MeqTrees Vagrantfile

Vagrantfile to create an Ubuntu 18.04 virtual machine suitable for
running the SKA1-Mid simulations with MeqTrees. MeqTrees and CASA are
installed from the KERN repository. The sim-mid-pointing repository is
cloned inside the VM.

To create the VM and install the software:
```
vagrant up
```

To connect to the VM:
```
vagrant ssh
```

Inside the VM, to run the simulations:
```
. setup.sh
cd sim-mid-pointing/meqtrees_simulation
python run_meqtrees.py
```

To copy results from the VM to the host, e.g.:
```
vagrant ssh-config > ssh-config
scp -F ssh-config default:sim-mid-pointing/meqtrees_simulation/simulation1/*.png .
```

To shut down the VM:
```
vagrant halt
```
Note this preserves the state of the VM, so it can be started again
with `vagrant up`.

If you want to completely remove the VM:
```
vagrant destroy
```
