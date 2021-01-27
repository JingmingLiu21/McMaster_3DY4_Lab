# Copyright by Jingming Liu, Riley Mione and Nicola Nicolici
# Department of Electrical and Computer Engineering
# McMaster University
# Ontario, Canada
import numpy
import scipy
import cmath
from cmath import exp
import math
from fourierTransform import generateSin, plotSpectrum, plotTime, generate_square_wave
from numpy.random import randint
from math import pi
from scipy.fft import fft, ifft
def discrete_fourier_transform(x):
    """
    Discrete Fourier Transform
    """
    N = len(x)
    X = [0] * N

    for m in range(N):
        for k in range(N): 
            arg = (-2*pi*k*m*1j)/N
            X[m] = X[m] + x[k]*exp(arg)
    
    return X


def inverse_discrete_fourier_transform(signal:list):
    N = len(signal)
    x = [complex(0)] * N
    for k in range(N):
        for m in range(N):
            arg = (2*pi*k*m*1j)/N
            x[k] += (signal[m] * exp(arg))
        
        x[k] = x[k]/N

    return x


def get_square_absolute_value(signal:list):
    sum = 0
    for sample in signal:
        # Python's abs() function actually works with complex numbers too! 
        
        sum += abs(sample)**2
    return sum

def generate_random_signal(length:int):
    """
    Generates signal of specified lenght
    """
    # TODO: use randn

    return randint(-10, 10, length) 


def multi_tone_signal(num=1):
    """
    tones specifies the amount of tones to include in this multi-tone signal, with 1 being the default
    """
    tone_arr = [] 
    combined_time = numpy.array([])
    combined_signal = numpy.array([])
    for tone_number in range(num): 
        # using randint for now but this may change
        amplitude = randint(-10, 10) 
        frequency = randint(low=1, high=1e+6) # From 1Hz to 1MHz
        nyquist_frequency = 2 * frequency
        phase = randint(low=0, high= 10) 
        
        # Plots tuple with numpy
    
        # tone_arr.append((generateSin(Fs=nyquist_frequency, interval=100e-6, frequency=frequency, amplitude=amplitude, phase=phase), nyquist_frequency))
    #     tone_arr.append(generateSin(Fs=nyquist_frequency, interval=100e-6, frequency=frequency, amplitude=amplitude, phase=phase))
    # return tone_arr
        time, signal = generateSin(Fs=nyquist_frequency, interval=100e-6, frequency=frequency, amplitude=amplitude, phase=phase)
        
        if len(combined_signal) < len(signal):
            combined_signal.resize(signal.shape)
            combined_signal = combined_signal + signal

        else:
            signal.resize(combined_signal.shape)
            combined_signal = combined_signal + signal

        if len(combined_time) < len(time):
            combined_time = time
        else: 
            combined_time = combined_time




    return (combined_time, combined_signal)
if __name__ == "__main__":
    # # print(discrete_fourier_transform(x))
    # signal = generate_random_signal(1000)
    # # print(signal)

    # dft_signal = discrete_fourier_transform(signal)
    # print(f"Squared abs value of original signal: {get_square_absolute_value(signal)}")
    # print(f"Squared abs value of original signal: {get_square_absolute_value(dft_signal)}")
    
    # scipy_dft = list(fft(signal))
    
    
    # # tesing our discrete fourier transform function
  
 
    # print("Comparing dft")
    # print(numpy.allclose(dft_signal, scipy_dft))

    # idft_signal = inverse_discrete_fourier_transform(dft_signal)
    # scipy_idft = list(ifft(dft_signal))

    # # testing our inverse discrete fourier transform 
    # print("Comparing idft!")
    # print(numpy.allclose(idft_signal, scipy_idft))

    # # plot the error
    # print(multi_tone_signal(3))
    # tones = multi_tone_signal(3)

    # for tone in tones:
    #     print(len(tones))
    #     waveform, Fs = tone
    #     time, signal = waveform
    #     plotTime(signal, time)

    #     plotSpectrum(signal, Fs)

    # Testing of plotting signal from multitone signal function 
    # time, signal = multi_tone_signal(3)
    # plotTime(signal, time)
    # plotSpectrum(signal, 1e+5)
    Fs = 2e+6
    time, signal = square_wave = generate_square_wave(Fs=Fs, frequency=1000, interval=0.01, amplitude=0.5, phase=1,duty_cycle=0.5)
    plotTime(signal, time)
    plotSpectrum(signal, Fs)