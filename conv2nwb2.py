# modules for NWB
from pynwb import NWBFile, get_manager
from datetime import datetime
from pynwb.form.backends.hdf5 import HDF5IO
from pynwb.ecephys import ElectricalSeries, TimeSeries
from pynwb.ecephys import ElectrodeTable, ElectrodeTableRegion
#  modules for reading Blackrock *.ns files
import sys
sys.path.append('/home/lauri/code/brPy/')
from brpylib import NsxFile
# modules for sys info
import psutil
# numpy 
import numpy as np
# some functions
import datapreprocesslib as dpl


###################################
# init NWB-file
f = NWBFile('N', 'Anesthetized marmoset recording', 'MM666', datetime.now(),
            experimenter='Dr. Lauri Nurminen',
            lab='Angelucci-lab',
            institution='University of Utah',
            experiment_description='Orientation tuning recording using Utah array in monkey V1',
            session_id='MM666')

# create device object
device = f.create_device(name='128-ch Cerebus', source="a source")

# create electrode group
electrode_name = 'Utah-array'
source = "Blackrock Microsystems"
description = "96-channel Utah array"
location = "V1"
electrode_group = f.create_electrode_group(electrode_name,
                                           source=source,
                                           description=description,
                                           location=location,
                                           device=device)



###################################

# Open file and extract headers
datafile = '/opt3/MM385/P1/20170518-223225/20170518-223225-001.ns5'
nsx_file = NsxFile(datafile)
# may cause problems with very large files
Data     = nsx_file.getdata()

# get data length and sampling rate
data_length = Data['data_headers'][0]['NumDataPoints']
sample_rate = Data['samp_per_s']
channel_num =  nsx_file.basic_header['ChannelCount']
nsx_file.close()

# set electrode table for the file
electrode_table = ElectrodeTable('V-Probe_table')
# add electrode entries
for idx in np.arange(channel_num):
    electrode_table.add_row(idx,
                    x=1.0, y=2.0, z=0.0,
                    imp=float(-idx),
                    location='V1', filtering='0.3-7500Hz',
                    description='channel %s' % idx, group=electrode_group)


# set electrode table
f.set_electrode_table(electrode_table)
# include all electrode at once
electrode_table_region = ElectrodeTableRegion(electrode_table, list(range(channel_num)), 'the first and third electrodes')
# write data for all electrodes
ephys_ts = ElectricalSeries('orientation1',
                            'an hypothetical source',
                            Data['data'],
                            electrode_table_region,
                            timestamps=np.arange(data_length)/sample_rate,
                            resolution=1.0,
                            comments="This data was randomly generated with numpy, using 1234 as the seed",
                            description="Random numbers generated with numpy.random.rand")

# stimulus_ts = TimeSeries('Stimulus timeseries',
#                          'VSG',
#                          stimulus_data,
#                          'degrees',
#                          timestamps=ephys_timestamps,
#                          resolution=0.001,
#                          comments="This data was randomly generated with numpy, using 1234 as the seed",
#                          description="Imagined orientation timeseries")
                         

filename = 'orientation-data.h5'
f.add_acquisition(ephys_ts)
f.add_acquisition(ephys_ts1)
io = HDF5IO(filename, manager=get_manager(), mode='w')
io.write(f)
io.close()
#f.add_stimulus(stimulus_ts)

