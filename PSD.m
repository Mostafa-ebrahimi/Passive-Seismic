function [f, mag, phase] = spect(signal, sampling)
%This code is wrriten by mostfa ebrahimi
%the master student of Geophysics, in university of Tehran
%This code is about V/H in Passive Seismic Studies
%----------------------------------------------------------------------
nyquist = sampling /2;
%padding with zeros to the power of 2
N = 2.^ceil(log2(length(signal)));
f = 0:(2*nyquist/N):(nyquist-nyquist/N);
Fsignal = fft(signal, N);
%use only the first half of symetric output
Fsignal = Fsignal(1:ceil(N/2));
%magnitude
mag = abs(Fsignal)/sampling;
%energy must be the same - take into account the missing second half
%mag = 2*mag;
%scale FFT to get the real values
%mag = mag/length(signal);
%angle
phase = angle(Fsignal);
f(1)=f(2);
% phase = 2*phase;
% phase = phase/length(signal);