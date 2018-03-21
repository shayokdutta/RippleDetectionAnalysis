import scipy
import scipy.signal as signal
import numpy as np

FS=3000
#Online Filter Taps...To Use for Simulated Detection
#Bandpass FIR Filter Coeffs 150-250Hz passband
bandpassFilterTaps=signal.firwin(30, [150,250], nyq=FS/2, pass_zero=False)
#Lowpass FIR Filter Coeffs (After Absolute Value)
lowpassFilterTaps=np.asarray([0.0203770957,
        0.0108532903,
        0.0134954582,
        0.0163441640,
        0.0193546202,
        0.0224738014,
        0.0256417906,
        0.0287934511,
        0.0318603667,
        0.0347729778,
        0.0374628330,
        0.0398648671,
        0.0419196133,
        0.0435752600,
        0.0447894668,
        0.0455308624,
        0.0457801628,
        0.0455308624,
        0.0447894668,
        0.0435752600,
        0.0419196133,
        0.0398648671,
        0.0374628330,
        0.0347729778,
        0.0318603667,
        0.0287934511,
        0.0256417906,
        0.0224738014,
        0.0193546202,
        0.0163441640,
        0.0134954582,
        0.0108532903,
        0.0203770957])

def rippleBandFilterSimulated(lfp, time, FS, bpFilterTaps, lpFilterTaps):
    """
    Ripple band filter and envelope simulating real-time algorithm
    """
    #Bandpass filter into ripple band
    rippleData = signal.lfilter(bpFilterTaps,1,lfp)
    #Envelope
    rippleEnvelope = np.absolute(rippleData)
    #smooth
    smoothed_envelope = signal.lfilter(lpFilterTaps,1,rippleEnvelope)
    return smoothed_envelope, rippleData

#Ripple band filtering and envelope
def rippleBandFilter(lfp, time, FS):
    #Bandpass filter into ripple band
    b = signal.firwin(25, [150/(FS/2), 250/(FS/2)], pass_zero=False)
    rippleData = signal.filtfilt(b,1,lfp)
    #Hilbert transform
    rippleEnvelope = np.absolute(signal.hilbert(rippleData))
    #Smooth envelope with a gaussian
    EnvelopeSmoothingSD = 0.004 * FS
    smoothed_envelope = scipy.ndimage.filters.gaussian_filter1d(rippleEnvelope, EnvelopeSmoothingSD, mode = 'constant')
    return smoothed_envelope, rippleData

from itertools import groupby
from operator import itemgetter
def thresholdCrossings(data, thresh):
    aboveThreshold = np.where(data > thresh, 1, 0) #sets the value of indices (same as data) greater than thresh to 1
    eventList = []
    eventMax = []
    #loop below groups runs of zeros and ones together
    for k, v in groupby(enumerate(aboveThreshold),key=itemgetter(1)):
        if k: #k = 1 if a run of 1 is grouped together within v
            v = list(v)
            eventList.append([v[0][0],v[-1][0]]) #append indices for run of 1s
            try:
                eventMax.append(data[v[0][0]:(v[-1][0]+1)].max()) #find max value of data during run of 1s indices
            except:
                print(v, data[v[0][0]:v[-1][0]]) #print out run of 1s indices within data
    eventMax = np.asarray(eventMax)
    eventList = np.asarray(eventList)
    return eventList, eventMax