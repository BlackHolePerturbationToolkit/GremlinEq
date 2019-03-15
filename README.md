# Gremlin - Teukolsky equation solver for a point particle on a bound, timelike orbit

Gremlin is a toolkit for studying solutions to the Teukolsky equation using a point-particle source. GremlinEq is a version of Gremlin that specializes the source to circular and equatorial orbits of Kerr black holes. The full Gremlin package will be released publicly via the Black Hole Perturbation Toolkit; certain proprietary libraries that were used in its development must be cleaned up and replaced with Open Source resources. GremlinEq has been so cleaned, and is hereby provided as an initial release. 

## Requirements

GremlinEq depends upon:  

- The [GNU Scientific Library (GSL)][1] version >= 2.5  
- The [GNU Multiple Precision Arithmetic Library (GMP)][2]  
- The [Fastest Fourier Transform in the West][3]  
- The [Hierarchical Data Format][4]  
	
If you are on a Mac these are can be easily installed using [Brew][5]. On Linux the package manager should be able to install them. We're not sure of the best way to install them on Windows (please let us know!)
	
To compile GremlinEq you will need:  

- a C++ compiler (we have tested it with g++)  
- Make
	
## Execution

Once you have compiled GremlinEq the binaries are placed in the bin/ folder. The main binary is 'Circ_Eq'. This will output an HDF5 file which contains information about the computed mode. Many other binaries are provided which provide additional functions -- see the documentation for details.

## Documentation

Initial documentation can be found in doc/doc.pdf  

Doxygen documentation can be found in doc/html/index.html

## Authors

Scott Hughes



[1]: https://www.gnu.org/software/gsl/
[2]: https://gmplib.org/
[3]: http://www.fftw.org/
[4]: https://www.hdfgroup.org/solutions/hdf5/
[5]: https://brew.sh/