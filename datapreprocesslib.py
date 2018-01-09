import numpy as np

def PulseCounter(pulse_data,thr):
    temp2 = np.where(pulse_data > thr)
    temp2 = temp2[0]
    temp3 = np.diff(temp2)
    temp4 = np.where(temp3 > 1)
    temp4 = temp4[0]
    temp4 = temp4 + 1
    temp4 = np.append([1],temp4)

    if temp2.size == 0:
        PulseSamples = np.array([])
    else:
        PulseSamples = temp2[temp4]

    return PulseSamples

def getStimParams(M):
    M = M['stimParams']
    stimOrder  = M['stimOrder'][0,0]
    on_time_ms = M['on_time_ms'][0,0]
    on_time_ms = on_time_ms[0,0]
    off_time_ms = M['off_time_ms'][0,0]
    off_time_ms = off_time_ms[0,0]
    stim_value  = M['Value'][0,0]
    stimParams  = {'stimOrder':stimOrder,'on_time_ms':on_time_ms,'off_time_ms':off_time_ms,
                   'FSnev':30e3,'FS4':10e3,'stim_value':stim_value}

    return stimParams

def modStim4NWB2(stim_value,pulse_samples,NS5_numsamples,on_time_ms=500,off_time_ms=500,FSnev=30e3,FSns=10e3):
    """
    Formats VisaGe experiment description to NWB2.0

    modStim4NWB2(stim_value,pulse_samples,NS5_numsamples,on_time_ms=500,off_time_ms=500,FSnev=30e3)

    stim_value: matrix describing the stimulus that was run in a block, e.g. orientation
    NS5_numsamples: number of samples per channel in the NS5 file
    on_time_ms=500: stimulus presentation time in milliseconds
    off_time_ms=500: inter-stimulus-interval in milliseconds
    FSnev=30e3: sampling rate of the NS5 file in Hz
    FSns4=10e3: sampling rate of the NS4 file in Hz
    """

    # output container for stimulus informations
    stim_out = np.nan * np.ones(NS5_numsamples)

    # get rid of block end pulses
    pulse_samples = pulse_samples[::2].copy()
    print(pulse_samples.shape)
    
    # initialize containers for blank before stimulus
    stim_off_value = np.nan * np.ones((stim_value.shape[0],int(off_time_ms/1000.0*FSnev)))
    
    # repeat for each trial
    for i in range(stim_value.shape[1]):
        stim_on_value = np.repeat(stim_value[:,i],int(on_time_ms/1000.0*FSnev))
        stim_on_value = stim_on_value.reshape(stim_value.shape[0],int(on_time_ms/1000.0*FSnev))
        # catenate the blank period
        stim_on_value = np.hstack((stim_off_value,stim_on_value))
        stim_out[pulse_samples[i]:pulse_samples[i] + stim_on_value.flatten().shape[0]] = stim_on_value.flatten()

    
    return stim_out
