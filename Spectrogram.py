from obspy.fdsn import Client
from obspy import UTCDateTime

import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.mlab as ml
import scipy
from scipy.interpolate import griddata

client = Client()

def getWave(network, station, number, channel, UTC, dur):
    """
    Downloads miniseed datasets through the obspy.fdsn function.     
    """
    t = UTCDateTime(UTC)
    st = client.get_waveforms(network, station, number, channel, t, t + dur, attach_response=True)
    print st
    return st
    
def preprocess(stream):
    """Carries out simple preprocessing of trace, by first merging the stream, 
    removing instrumetn response, highpass filtering at 0.2 Hz then tapering"""
    stream.merge()
    stream.remove_response(output="vel")
    stream.filter('highpass',freq=0.02,corners=2,zerophase=True)
    stream.taper(max_percentage=0.01,type='cosine')
    return stream
    
def spec_amp(stream):
    """This produces three 1d arrays constructed after splitting the noise trace 
    into 20min intervals, with the intention that they will be used to create a 
    spectrogram. Using a while loop the routine cuts the trace, creates the time
    array using the mid point of the cut trace, carries out an fft on the cut 
    trace and creates an array of frequencies. These are then stacked onto 
    previously constructed datasets. Note, at the current time I haven't really 
    got my head around structuring np arrays, hence the fudges for reordering arrays."""
    
    data=np.array((0))
    time=np.array((0))
    freq=np.array((0))
    
    t1=stream[0].stats.starttime
    t2=stream[0].stats.endtime
    samp_rate=stream[0].stats.sampling_rate 
    
    t=t1
    
    while t < t2:
        cut=stream[0].slice(t,t+1200)
        tmid=t+600
        t=t+1200
        length = len(cut)-1
        
        time_temp=np.zeros(shape=(length,1))
        
        for a in range(0,length):
            time_temp[a]=[tmid.timestamp]

        time=np.vstack((time,time_temp))
        
        cut_data=cut.data
        cut_spec=np.fft.fft(cut_data,n=length)
    
        data_temp=np.zeros(shape=(length,1),dtype=np.complex)
        
        for b in range(0,length):
            data_temp[b]=cut_spec[b]
    
        data=np.vstack((data,data_temp))
        
        freq_temp = np.fft.fftfreq(length, d=1./samp_rate)
    
        freq_temp2=np.zeros(shape=(length,1))
        
        for c in range(0,length):
            freq_temp2[c]=freq_temp[c]   
    
    
        freq=np.vstack((freq,freq_temp2))
        
    time=np.delete(time,0,0)
    data=np.delete(data,0,0)
    freq=np.delete(freq,0,0)
    
    time=time-time.min()+600

    
    return time,data,freq
    
def Contour(time,freq,data):
    x=np.array(())
    y=np.array(())
    z=np.array(())
    
    x = np.append(x,time)
    y = np.append(y,abs(freq))
    z = np.append(z,abs(data))
    
    points = x,y
    values = z
    
    ndata = len(data)
    ny, nx = 1, 1200
    xmin, xmax = x.min(), x.max()
    ymin, ymax = y.min(), y.max()
    
    xi, yi = np.mgrid[x.min():x.max():nx, ymin:ymax:ny]
    
    zi = griddata(points, values, (xi, yi), method='nearest') 

    
    return xi,yi,zi
    
if __name__=="__main__":    
    
    st = getWave('GB','HTL','*', 'BHZ', '2010-02-15T00:00:00.000',7200)
    st_pre = preprocess(st)
    time,data,freq=spec_amp(st_pre)
    xi,yi,zi=Contour(time,freq,data)
    
    fig = plt.figure()
    plt.contourf(xi, yi, zi,1000) 
    plt.clim(0,0.000005)
    plt.show()