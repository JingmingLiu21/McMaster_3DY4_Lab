# Copyright by Jingming Liu, Riley Mione and Nicola Nicolici
# Department of Electrical and Computer Engineering
# McMaster University
# Ontario, Canada

import matplotlib.pyplot as plt
import numpy as np
from filterDesign import low_pass_impluse
from scipy.io import wavfile
from scipy import signal
from filterDesign import low_pass_impluse, own_convolution

def filter_block_processing(audio_data, block_size, audio_Fc, audio_Fs, N_taps):

	# derive filter coefficients
	firwin_coeff = signal.firwin(N_taps, audio_Fc/audio_Fs)
	# we assume the data is stereo as in the audio test file
	filtered_data = np.empty(shape = audio_data.shape)
	# start at the first block (with relative position zero)
	position = 0

	while True:

		# filter both left and right channels
		filtered_data[position:position+block_size, 0] = \
		signal.lfilter(firwin_coeff, 1.0, \
			audio_data[position:position+block_size, 0])

		filtered_data[position:position+block_size, 1] = \
			signal.lfilter(firwin_coeff, 1.0, \
						   audio_data[position:position+block_size, 1])

		position += block_size
		if position > len(audio_data):
			break

	# to properly handle blocks you will need to use
	# the zi argument from lfilter from SciPy
	# explore SciPy, experiment, understand and learn!

	return filtered_data

def block_filtering_part3( audio_data, \
				block_size, \
				audio_Fc, \
				audio_Fs, \
				N_taps):

	# derive filter coefficients
	# we assume the data is stereo as in the audio test file
	low_pass_coeff = low_pass_impluse(Fc=audio_Fc, Fs=audio_Fs, N_taps=N_taps)
	filtered_data = np.empty(shape = audio_data.shape)
	# start at the first block (with relative position zero)
	position = 0


	# Create the initial zi value! 
	# z1 = signal.lfilter_zi(b=filtered_data[0, 0], a=1.0)
	# z2 = signal.lfilter_zi(b=filtered_data[0, 1], a=1.0)
	z1 = [0]*(len(low_pass_coeff) - 1)
	z2 = [0]*(len(low_pass_coeff) - 1)
	while True:

		# filter both left and right channels
		# filtered_data[position:position+block_size, 0], z1 = \
		# 	signal.lfilter(low_pass_coeff, 1.0, \
		# 	audio_data[position:position+block_size, 0], zi=z1)

		filtered_data[position:position+block_size, 0]= \
			own_convolution(coe=low_pass_coeff, deno_coe=1.0, \
			sig=audio_data[position:position+block_size, 0])

		# filtered_data[position:position+block_size, 1], z2 = \
		# 	signal.lfilter(low_pass_coeff, 1.0, \
		# 	audio_data[position:position+block_size, 1], zi=z2)
		filtered_data[position:position+block_size, 1]= \
			own_convolution(coe=low_pass_coeff, deno_coe=1.0, \
			sig=audio_data[position:position+block_size, 1])

		position += block_size
		if position > len(audio_data):
			break

	# to properly handle blocks you will need to use
	# the zi argument from lfilter from SciPy
	# explore SciPy, experiment, understand and learn!

	return filtered_data


def lpf_block_data_using_zi( audio_data, \
				block_size, \
				audio_Fc, \
				audio_Fs, \
				N_taps):

	# derive filter coefficients
	# we assume the data is stereo as in the audio test file
	low_pass_coeff = low_pass_impluse(Fc=audio_Fc, Fs=audio_Fs, N_taps=N_taps)
	filtered_data = np.empty(shape = audio_data.shape)
	# start at the first block (with relative position zero)
	position = 0


	# Create the initial zi value! 
	# z1 = signal.lfilter_zi(b=filtered_data[0, 0], a=1.0)
	# z2 = signal.lfilter_zi(b=filtered_data[0, 1], a=1.0)
	z1 = [0]*(len(low_pass_coeff) - 1)
	z2 = [0]*(len(low_pass_coeff) - 1)
	while True:

		# filter both left and right channels
		filtered_data[position:position+block_size, 0], z1 = \
			signal.lfilter(low_pass_coeff, 1.0, \
			audio_data[position:position+block_size, 0], zi=z1)

		filtered_data[position:position+block_size, 1], z2 = \
			signal.lfilter(low_pass_coeff, 1.0, \
			audio_data[position:position+block_size, 1], zi=z2)


		position += block_size
		if position > len(audio_data):
			break

	# to properly handle blocks you will need to use
	# the zi argument from lfilter from SciPy
	# explore SciPy, experiment, understand and learn!

	return filtered_data

def filter_single_pass(audio_data, audio_Fc, audio_Fs, N_taps):

	# derive filter coefficients
	firwin_coeff = signal.firwin(N_taps, audio_Fc/audio_Fs)
	# we assume the data is stereo as in the audio test file
	filtered_data = np.empty(shape = audio_data.shape)
	# filter left channel
	filtered_data[:,0] = signal.lfilter(firwin_coeff, 1.0, audio_data[:,0])
	# filter stereo channel
	filtered_data[:,1] = signal.lfilter(firwin_coeff, 1.0, audio_data[:,1])

	return filtered_data

# audio test file from: https://www.videvo.net/royalty-free-music/
if __name__ == "__main__":

	# use use wavfile from scipy.io for handling .wav files
	audio_Fs, audio_data = wavfile.read("../data/audio_test.wav")
	print(' Audio sample rate = {0:f} \
		\n Number of channels = {1:d} \
		\n Numbef of samples = {2:d}' \
		.format(audio_Fs, audio_data.ndim, len(audio_data)))

	# you can control the cutoff frequency and number of taps
	single_pass_data = filter_single_pass(audio_data, \
						audio_Fc = 10e3, \
						audio_Fs = audio_Fs, \
						N_taps = 51)

	# write filtered data back to a .wav file
	wavfile.write("../data/single_pass_filtered.wav", \
	 			audio_Fs, \
				single_pass_data.astype(np.int16))

	# you can control also the block size
	block_processing_data = filter_block_processing(audio_data, \
						block_size = 1000, \
						audio_Fc = 10e3, \
						audio_Fs = audio_Fs, \
						N_taps = 51)

	wavfile.write("../data/block_processing_filtered.wav", \
	 			audio_Fs, \
				block_processing_data.astype(np.int16))

	# it is suggested that you add plotting while troubleshooting
	# if you plot in the time domain, select a subset of samples,
	# from a particular channel (or both channels) e.g.,
	# audio_data[start:start+number_of_samples, 0]

	our_prosessing_data = lpf_block_data_using_zi(audio_data, \
		block_size = 1000, \
		audio_Fc = 10e3, \
		audio_Fs = audio_Fs, \
		N_taps = 51)
	
	wavfile.write("../data/our_processing_filtered.wav", \
	    audio_Fs, \
		our_prosessing_data.astype(np.int16))



	our_prosessing_data = block_filtering_part3(audio_data, \
		block_size = 1000, \
		audio_Fc = 10e3, \
		audio_Fs = audio_Fs, \
		N_taps = 51)
	
	wavfile.write("../data/part3_processing_filtered.wav", \
	    audio_Fs, \
		our_prosessing_data.astype(np.int16))

    # plt.show()
