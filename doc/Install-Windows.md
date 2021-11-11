## Installing Gremlin on Windows

The below instructions were provided by Vojtech Witzany.

Install Windows Subsystem for Linux,. I will assume you have installed Ubuntu, see https://ubuntu.com/wsl

Open the WSL Ubuntu terminal, you can install all the needed packages for Gremlin by:

sudo apt-get install g++ libhdf5-dev libgsl-dev libgmp3-dev fftw3-dev make

Now go to a convenient folder, clone the Gremlin git folder, and compile using make:

git clone https://github.com/BlackHolePerturbationToolkit/GremlinEq.git
cd GremlinEq
make LIBHDF5='/usr/lib/x86_64-linux-gnu/hdf5/serial' INCLHDF5='/usr/include/hdf5/serial/' 

Notice that I have changed the LIBHDF5 and INCLHDF5 arguments for make (default values are in Makefile). This is because g++ was not able to find the hdf5 library for some reason. Similar modifications to linking may need to be made in your case.

After a succesful compilation, the binaries are in /bin and you can run them e.g. as:
./bin/Circ_Eq