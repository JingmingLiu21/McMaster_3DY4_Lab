#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include <cmath>
#include <chrono>

#define PI 3.14159265358979323846

// // input x is qualified as const in order not to be changed within the function
// although C++ is case sensitive it is "safer" to use Xf instead of X in order to avoid confusion
void DFT(const std::vector<float> &x, std::vector<std::complex<float>> &Xf) {
	// allocate space for the vector holding the frequency bins
	// initialize all the values in the frequency vector to complex 0
	Xf.resize(x.size(), static_cast<std::complex<float>>(0, 0));
	// "auto" keyword is used for type deduction (or inference) in C++
	for (auto m = 0; m < Xf.size(); m++) {
		for (auto k = 0; k < x.size(); k++) {
				std::complex<float> expval(0, -2*PI*(k*m) / x.size());
				Xf[m] += x[k] * std::exp(expval);
		}
	}
}

std::vector<float> inverseDiscreteFourier(std::vector<std::complex<float>> &signal) {//, std::vector<std::complex<float>> &x) {
	// We want to loop over the length of the signal we want to perform the inverse transform on
	int N = signal.size();
	std::vector<float>  x(N);

	x.resize(N, static_cast<float> (0));
	for (auto k = 0; k < N; k++) {
		for (auto m = 0; m < N; m++) {
			std::complex<float> expval(0, 2*PI*(k*m) / (float)N);
			x[k] += std::real(signal[m] * std::exp(expval));
		

		}
		x[k] = (x[k]/ (float)N);
	}
	return x;
}


// function to generate N random float values
// x is the input/output vector
// N is the number of samples
// random values are between -max and max
// precision should be capped to 3 (in the context of this type of experiments)
void generateRandomSamples(std::vector<float> &x, unsigned int N, unsigned short int max, unsigned char precision)
{
	// allocate space for the vector holding the frequency bins
	x.resize(N);
	int int_radom_max = 2*(max * static_cast<int>(pow(10,precision)));
	for (auto i = 0; i < x.size(); i++) {
		// static casting is generally safer and better for both
		// type checking as well as for documentation purposes
		x[i] = (static_cast<float>(std::rand() % int_radom_max));
		x[i] = (x[i] - (int_radom_max/2))/pow(10,precision);
	}
}

// function to print a real vector
void printRealVector(const std::vector<float> &x)
{
	std::cout << "Printing float vector of size " << x.size() << "\n";
	for (auto i = 0; i < x.size(); i++)
		std::cout << x[i] << " ";

	// starting from C++11 you can write
	// for (auto elem : x)
	//   std::cout << elem << " ";

	std::cout << "\n";
}

// function to print a complex vector
void printComplexlVector(const std::vector<std::complex<float>> &X)
{
	std::cout << "Printing complex vector of size " << X.size() << "\n";
	for (auto i = 0; i < X.size(); i++)
		std::cout << X[i] << " ";
	std::cout << "\n";
}

float average(int iterations, int N){
	float avg = 0.0;

	for(int i=0; i<iterations; i++) {
		int seed = std::time(0x0);
		// std::cout << "Starting from seed " << seed << "\n";

		// declare a vector of real values; no memory is allocated at this time
		std::vector<float> x;
		// generate 32 samples between -10 and 10; for extra flexibility
		// the last argument gives precision in terms of fraction digits
		// note: memory for x is allocated within the function called below
		generateRandomSamples(x, N, 10, 2);
		// print a vector of real numbers
		// printRealVector(x);

		// declare a vector of complex values; no memory is allocated for it
		std::vector<std::complex<float>> Xf;
		// perform DFT of x to produce Xf
		// we measure the execution time using the "chrono" class
		auto start_time = std::chrono::high_resolution_clock::now();
		DFT(x, Xf);
		auto vec = inverseDiscreteFourier(Xf);
		
		auto stop_time = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double, std::milli> DFT_run_time = stop_time-start_time;
		// print a vector of complex numbers
		// printComplexlVector(Xf);
		std::cout << "DFT + IDFT ran for " << DFT_run_time.count() << " milliseconds" << "\n";
		avg += DFT_run_time.count()/iterations; 

		

	} 
	

	return avg;



} 

// int main()
// {
// 	std::ofstream myfile;
// 	myfile.open("graphing.csv", std::ios::app );
// 	for(int i=7; i<=14; i++) {
// 		printf("%d\n", i);
// 		myfile << ("\n");
// 		int arrnum = (int)pow(2.0, (double)i);
// 		float avg = average(8, arrnum);
// 		std::cout << "Average Time ran for a DFT and IDFT pair " << avg << " milliseconds\n";
// 		std::cout << "Factor this time " << arrnum << "\n";
// 		myfile << arrnum << "," << avg << "\n"; 
// 	}
	
// 	myfile.close();

// 	// // use initialize the seed for the random number generator using current time
// 	// int seed = std::time(0x0);
// 	// std::cout << "Starting from seed " << seed << "\n";

// 	// // declare a vector of real values; no memory is allocated at this time
// 	// std::vector<float> x;
// 	// // generate 32 samples between -10 and 10; for extra flexibility
// 	// // the last argument gives precision in terms of fraction digits
// 	// // note: memory for x is allocated within the function called below
// 	// generateRandomSamples(x, 128, 10, 2);
// 	// // print a vector of real numbers
// 	// // printRealVector(x);

// 	// // declare a vector of complex values; no memory is allocated for it
// 	// std::vector<std::complex<float>> Xf;
// 	// // perform DFT of x to produce Xf
// 	// // we measure the execution time using the "chrono" class
// 	// auto start_time = std::chrono::high_resolution_clock::now();
// 	// DFT(x, Xf);
// 	// auto vec = inverseDiscreteFourier(Xf);
	
// 	// auto stop_time = std::chrono::high_resolution_clock::now();
// 	// std::chrono::duration<double, std::milli> DFT_run_time = stop_time-start_time;
// 	// // print a vector of complex numbers
// 	// // printComplexlVector(Xf);
// 	// std::cout << "DFT + IDFT ran for " << DFT_run_time.count() << " milliseconds" << "\n";




// 	// // Our idft 
// 	// auto start_time_1 = std::chrono::high_resolution_clock::now();
// 	// auto vec = inverseDiscreteFourier(Xf);
// 	// auto stop_time_1 = std::chrono::high_resolution_clock::now();
// 	// std::chrono::duration<double, std::milli> IDFT_run_time = stop_time_1-start_time_1;
// 	// // print a vector of complex numbers
// 	// // printComplexlVector(Xf);
// 	// std::cout << "IDFT ran for " << IDFT_run_time.count() << " milliseconds" << "\n";
// 	// // finished!
// 	// printRealVector(x);
// 	// printRealVector(vec);
// 	return 0;
// }
