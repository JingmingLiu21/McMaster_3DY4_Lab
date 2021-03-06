# Copyright by Jingming Liu, Riley Mione and Nicola Nicolici
# Department of Electrical and Computer Engineering
# McMaster University
# Ontario, Canada

import matplotlib.pyplot as plt
from scipy import signal
import numpy as np
import math
from lab1 import multi_tone_signal
# use generateSin/plotTime from the fourierTransform module
from fourierTransform import plotTime, generateSin


def freqzPlot(coeff, msg):

	# find the frequency response using freqz from SciPy:
	# https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.freqz.html
	w, h = signal.freqz(coeff)

	# plots the magnitude response where the x axis is normalized in rad/sample
	# Reminder: math.pi rad/sample is actually the Nyquist frequency
	fig, ax1 = plt.subplots()
	ax1.set_title('Digital filter frequency response (' + msg + ')')
	ax1.plot(w, 20 * np.log10(abs(h)), 'b')
	ax1.set_ylabel('Amplitude [dB]', color='b')
	ax1.set_xlabel('Frequency [rad/sample]')

# uncomment the lines below if you wish to inspect the phase response
# Note: as important as the phase response is, it is not critical
# at this stage because we expect a linear phase in the passband

# ax2 = ax1.twinx()
# angles = np.unwrap(np.angle(h))
# ax2.plot(w, angles, 'g')
# ax2.set_ylabel('Angle (radians)', color='g')

def filterSin(Fs, Fc, coeff):

	# we can control the frequency relative to the filter cutoff
	time, x = generateSin(Fs, interval = 1.0, frequency = Fc * 0.4)
	plotTime(x, time)

	# use lfilter from SciPy for FIR filtering:
	# https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.lfilter.html
	fx = signal.lfilter(coeff, 1.0, x)

	# you should cleary the effects (attenuation, delay) introduced by the filter
	plotTime(fx, time)
def band_pass_filter(Fs, N_tap, lower_bound=0, upper_bound=100):
    # Implement using firwin
    return signal.firwin(fs=Fs, numtaps=N_tap, cutoff=[lower_bound, upper_bound])
# Low_Pass
def low_pass_impluse(Fc,Fs,N_taps):
	norm_cf = Fc/(Fs/2)
	center_check = (N_taps-1)/2
	h = [0]*N_taps
	for i in range(N_taps):
		if i != center_check:
			h[i] = norm_cf * (math.sin(math.pi * norm_cf * (i - center_check))) / (math.pi * norm_cf * ((i - center_check)))
		else:
			h[i] = norm_cf
		h[i] = h[i] * math.pow((math.sin((i * math.pi) / N_taps)), 2)
	return h

def own_convolution(coe, deno_coe, sig):
	output = np.empty(shape=len(sig))
	element = np.empty(shape=(len(sig), len(coe)))
	counter = 0
	for j in range(sig.size):
		for taps in range(len(coe)):
			if deno_coe == 1:
				element[j][taps] = coe[taps] * sig[j-counter]
			else:
				break
			counter = counter + 1
			if j - counter < 0:
				break
		counter = 0
		output[j] = sum(element[j])
	return output
if __name__ == "__main__":

	Fs = 100           # sampling rate
	Fc = 15.0            # cutoff frequency
	N_taps = 40          # number of taps for the FIR
	# derive filter coefficients using firwin from Scipy:
	coe_hsci = signal.firwin(N_taps, Fc, window=('hann'), pass_zero=True, scale=True, nyq=None, fs=Fs)
	low_pass = low_pass_impluse(Fc, Fs, N_taps)
	#freqzPlot(coe_hsci, 'firwin= with ' + str(N_taps) + ' taps')
	#freqzPlot(low_pass, 'my_impulse_res = with ' + str(N_taps) + ' taps')
	# https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.firwin.html
	# second argument is the normalized cutoff frequency, i.e., the
	# cutoff frequency divided by Nyquist frequency (half of sampling rate)
	firwin_coeff = signal.firwin(N_taps, Fc/(Fs/2), window=('hann'))

	# plot the frequency response obtained through freqz
	#freqzPlot(firwin_coeff, 'firwin with ' + str(N_taps) + ' taps')

	# implement your own method for finding the coefficients for a low pass filter
	# my_own_coeff = ... provide the following arguments: Fc, Fs and N_taps
	# compare through visual inspection the frequency response against firwin
	# freqzPlot(my_own_coeff, 'my own FIR design with ' + str(N_taps) + ' taps')

	# you can confirm that a single tone has been filtered
	# filterSin(Fs, Fc, firwin_coeff)
	t, y = generateSin(Fs, interval=1.0, frequency=Fc*0.4)
	print(type(y))
	#plotTime(y,t)
	filter_result1 = signal.lfilter(low_pass,1,y)
	plotTime(filter_result1, t)
	own_filter_result = own_convolution(low_pass,1,y)
	plotTime(own_filter_result, t)
	#using impluse response filter to filter out the highest frequency
	time, sig = multi_tone_signal(3)
	filter_result = signal.lfilter(low_pass,0.00001,sig)
	#plotTime(sig, time)

	plt.show()
