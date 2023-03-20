1- This repo was created with the intention of implementing GMAP-TD in C++ (prior to its implementation in CUDA), with the goal of finding execution bottlenecks, etc.

2- A particularity of the implementation is that the armadillo library was used, which is a library specialized in linear algebra and is the most similar thing to Matlab that can be found in C++. The only thing that is currently essential from this library is the function to calculate eigenvalues and eigenvectors of a matrix. The rest of the things can be implemented by hand (matrix multiplication, matrix square root, etc.).

3- The file structure is as follows. On the one hand, we have the .m files, which are Matlab files whose goal is to simulate synthetic data and also to read the results from the processing and plot them. On the other hand, there are the .cpp files, which contain the implementation and testing of GMAP-TD. Below is the name of the file with its brief description.

1- genera_datos.m -> generates synthetic data. It uses the LibraryMeteo repo (addpath('../LibraryMeteo/dataGen'); addpath('../LibraryMeteo/spectrumEstimate');).

2- GMAP_TD_test.cpp -> function to test the implementation of GMAP-TD using the naive filter construction .

3- GMAP_TD_test_modificacionSebastian.cpp -> function to test the implementation of GMAP-TD using the filter construction efficiently, with the modification made by Sebastian.

4- algoritmos.cpp -> here is the implementation of the GMAP-TD algorithm. There are two possible implementations

5- common.cpp -> this is a file that contains several useful and necessary functions, which are called by the GMAP-TD algorithm.

6- GraficaResultadosMomentos.m -> Matlab script that reads the result of the moments estimated by GMAP-TD and plots them.

Way of Storing Simulated Data. Simulated data in Matlab is stored in binary form to make it easier to read the data in C++.

A version of the noise estimation is already implemented using the 1973 paper "Objective Determination of the Noise Level in Doppler Spectra."

How to Compile. To compile, you need to have the armadillo library installed. If you search http://arma.sourceforge.net/, you will find the necessary information to install it. I use Linux and compile as follows in the console: g++ -std=c++11 file_to_compile_name.cpp -o output_file_name -larmadillo
