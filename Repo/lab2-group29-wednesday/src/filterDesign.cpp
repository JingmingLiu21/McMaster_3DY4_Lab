#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <complex>
#include <cmath>
#include <chrono>
#include <numeric>
// #include "fourierTransform.cpp"
#define PI 3.14159265358979323846

// function for DFT (reused from previous experiment)
void DFT(const std::vector<float> &x, std::vector<std::complex<float>> &Xf)
{
	Xf.resize(x.size(), static_cast<std::complex<float>>(0, 0));
	for (auto m = 0; m < Xf.size(); m++) {
		for (auto k = 0; k < x.size(); k++) {
				std::complex<float> expval(0, -2*PI*(k*m) / x.size());
				Xf[m] += x[k] * std::exp(expval);
		}
	}
}

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


void large_DFT(const std::vector<float> &x, std::vector<std::complex<float>> &vec) {
	// Large DFT function to produce a very rudimentary frequency spectrum on large input
	std::vector<std::complex<float>> X_1; 
	std::vector<std::complex<float>> X_2; 
	std::vector<std::complex<float>> X_3; 
	std::vector<std::complex<float>> X_4; 
    // std::vector< std::vector <float> > element(
    //         x.size(),
    // std::vector<float>(h.size()));
	std::vector<std::complex<float>> stuff [4]= {X_1, X_2, X_3, X_4};

	int j;
	for(int i=0; i<4; i++) {
		j = 3 - i;
		auto v = (std::vector<float>(x.begin()+(i*256), x.end()-(j*256)));
		DFT(v, stuff[i]);



	}
	// after for loop runs, X_1... X_4 will be filled out
	
	vec.resize(256);
	for (int i=0; i<256; i++) {
		vec[i] = (stuff[0][i] + stuff[1][i] + stuff[2][i] + stuff[3][i])/(std::complex<float>)4;


	}
	printf("\nVEC SIZE: %d\n", vec.size());
}

// function to generate a sine with N samples per second over interval
void generateSin(std::vector<float> &t, std::vector<float> &x, float Fs, float interval, float frequency = 7.0, float amplitude = 5.0, float phase = 0.0)
{
	// we do NOT allocate memory space explicitly
	// for the time (t) vector and sample (x) vector
	t.resize(0); x.resize(0);
	float dt = 1/Fs;
	for (auto i = 0.0; i < interval; i += dt) {
		// vector size increases when pushing new elements into it
		t.push_back(i);
		x.push_back(amplitude*std::sin(2*PI*frequency*i+phase));
	}
}

// function to mix an array of sines
void mixSin(const std::vector<std::vector<float>> &sv, std::vector<float> &mixed)
{
	// assumes at least one sine passed
	// assumes all input sines are of the same size
	for (auto i = 0.0; i < sv[0].size(); i ++) {
		float mixval = 0.0;
		// note: sv.size() returns the number of sines (or rows in 2D repr)
		// sv[0].size() returns the number of samples in a sine (or cols in 2D repr)
		for (auto k = 0; k < sv.size(); k++)
			mixval += sv[k][i];
		mixed.push_back(mixval);
	}
}

// function to print a real vector (reused from previous experiment)
void printRealVector(const std::vector<float> &x)
{
	std::cout << "Printing float vector of size " << x.size() << "\n";
	for (auto i = 0; i < x.size(); i++)
		std::cout << x[i] << " ";
	std::cout << "\n";
}

// function to print a complex vector (reused from previous experiment)
void printComplexlVector(const std::vector<std::complex<float>> &X)
{
	std::cout << "Printing complex vector of size " << X.size() << "\n";
	for (auto i = 0; i < X.size(); i++)
		std::cout << X[i] << " ";
	std::cout << "\n";
} 

// function to record data in a format to be read by GNU plot
// the arguments are VERY specific to this usage in this experiment
// we have the time vector (t), a vector of sines (sv),
// input samples (x, i.e., mixed sines for this experiment),
// output samples (y, filtered input samples)
// frequency vectors (both Xf and Yf)
// note: the reference code does NOT do filtering, hence y and Yf are zeros by default
void plotMixedSinesSpectrum(const std::vector<float> &t, const std::vector<std::vector<float>> &sv, const std::vector<float> &x, const std::vector<float> &y, const std::vector<std::complex<float>> &Xf, const std::vector<std::complex<float>> &Yf)
{
	// write data in text format to be parsed by gnuplot
	const std::string filename = "../data/example.dat";
	std::fstream fd;  // file descriptor
	fd.open(filename, std::ios::out);
	fd << "#\tindex\tsine(0)\tsine(1)\tsine(2)\tdata in\tdata out\tspectrum in\tspectrum out\n";
	for (auto i = 0; i < t.size(); i++) {
		fd << "\t " << i << "\t";
		for (auto k = 0; k < sv.size(); k++)
			fd << std::fixed << std::setprecision(3) << sv[k][i] << "\t ";
		fd << x[i] << "\t "<< y[i] << "\t\t ";
		fd << std::abs(Xf[i])/Xf.size() << "\t\t " << std::abs(Yf[i])/Yf.size() <<"\n";
	}
	std::cout << "Generated " << filename << " to be used by gnuplot\n";
	fd.close();
}

// function to compute the impulse response "h" based on the sinc function
// see pseudocode from previous lab that was implemented in Python
void impulseResponseLPF(float Fs, float Fc, unsigned short int num_taps, std::vector<float> &h)
{
	// allocate memory for the impulse response
	h.resize(num_taps, 0.0);

	float norm_cutoff = Fc/(Fs/2);
	float halfway = ((num_taps - 1)/2);
	for (int i=0; i<num_taps; i++) {
		if (i == halfway) {
			h[i] = norm_cutoff;
		}

		else {
			h[i] = norm_cutoff*((sin((double)PI*norm_cutoff*halfway)) / (PI*norm_cutoff*halfway));
		}
		
		h[i] = h[i]*pow((sin((i*PI)/num_taps)), 2); // hann window
	}

}



// function to compute the filtered output "y" by doing the convolution
// of the input data "x" with the impulse response "h"; this is based on
// your Python code from the take-home exercise from the previous lab
void convolveFIR(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h)
{
	// allocate memory for the output (filtered) data
	y.resize(x.size()+h.size()-1, 0.0);
	// create 2d vector elements
    std::vector< std::vector <float> > element(
            x.size(),
    std::vector<float>(h.size()));

    int counter = 0;
    for (int i = 0; i < x.size();i++){
        for(int j = 0; j < h.size(); j++){
            element[i][j] = h[j] * x[i-counter];
            counter = counter + 1;
            if(j - counter < 0) {
                break;
            }
        }
        counter = 0;
        y[i] = std::accumulate(element[i].begin(), element[i].end(), 0);
    }




}

float average(int iterations, int N){
	float avg = 0.0;

	for(int i=0; i<iterations; i++) {
		int seed = std::time(0x0);
		// std::cout << "Starting from seed " << seed << "\n";

		// declare a vector of real values; no memory is allocated at this time


		// declare a vector of complex values; no memory is allocated for it

		// perform DFT of x to produce Xf
		// we measure the execution time using the "chrono" class
		auto start_time = std::chrono::high_resolution_clock::now();
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

int main() {


	float Fs = 1024.0;                    // samples per second
	float interval = 1.0;                 // number of seconds
	unsigned short int num_taps = 101;    // number of filter taps
	float Fc = 20.0;                      // cutoff frequency (in Hz)


	// std::vector<float> x;
	// 	// generate 32 samples between -10 and 10; for extra flexibility
	// 	// the last argument gives precision in terms of fraction digits
	// 	// note: memory for x is allocated within the function called below
	// generateRandomSamples(x, 1024, 10, 2);
	// printRealVector(x);
	// std::vector<std::complex<float>> Xf;
	// auto v = large_DFT(x, Xf);

	// generateRandomSamples(x, 1024, 10, 2);
	// std::vector<std::complex<float>> vec;
	// large_DFT(x, vec);
	// printComplexlVector(vec);

	// declare a vector of vectors for multiple sines
	std::vector<std::vector<float>> sv;
	// declare time and sine vectors
	std::vector<float> t, sine;
	// note: there is no explicit memory allocation through vector resizing
	// vector memory space will increase via the push_back method

	// generate and store the first tone
	// check the function to understand the order of arguments
	generateSin(t, sine, Fs, interval, 10.0, 5.0, 0.0);
	sv.push_back(sine);

	// generate and store the second tone
	generateSin(t, sine, Fs, interval, 40.0, 2.0, 0.0);
	sv.push_back(sine);
	// generate and store the third tone
	generateSin(t, sine, Fs, interval, 50.0, 3.0, 0.0);
	sv.push_back(sine);

	// declare the mixed sine vector and mix the three tones
	std::vector<float> x;
	mixSin(sv, x);
	// printRealVector(x);

	printf("%d\n", x.size());
	// declare a vector of complex values for DFT; no memory is allocated for it
	std::vector<std::complex<float>> Xf;
	DFT(x, Xf);
	// printComplexlVector(Xf);

	// generate the impulse response h
	// convolve it with the input data x
	// in order to produce the output data y
	std::vector<float> h;              // impulse response
	
	impulseResponseLPF(Fs, Fc, num_taps, h);
	printRealVector(h);
	std::vector<float> y;              // filter out
	convolveFIR(y, x, h);

	// compute DFT of the filtered data
	std::vector<std::complex<float>> Yf;
	std::vector<std::complex<float>> Zf;
	auto start_time_1 = std::chrono::high_resolution_clock::now();
	large_DFT(x, Yf);
	auto stop_time_1 = std::chrono::high_resolution_clock::now();
	
	std::chrono::duration<double, std::milli> dftlarge = stop_time_1-start_time_1;
	// print a vector of complex numbers
	// printComplexlVector(Xf);
	std::cout << "Large DFT ran for " << dftlarge.count() << " milliseconds" << "\n";
	printf("%d", y.size());
	
 
	// prepare the data for gnuplot
	plotMixedSinesSpectrum(t, sv, x, y, Xf, Yf);

	// naturally, you can comment the line below once you are comfortable to run gnuplot
	std::cout << "Run: gnuplot -e 'set terminal png size 1024,768' example.gnuplot > ../data/example.png\n";

	return 0;
}
