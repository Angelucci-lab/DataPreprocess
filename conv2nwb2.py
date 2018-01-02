# import modules for NWB
from pynwb import NWBFile, get_manager
from datetime import datetime
from pynwb.form.backends.hdf5 import HDF5IO
from pynwb.ecephys import ElectricalSeries, TimeSeries
from pynwb.ecephys import ElectrodeTable, ElectrodeTableRegion
# import modules for reading Blackrock *.ns files


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

# create electrode group (for machine readability only)
electrode_name = 'Utah-array'
source = "Blackrock Microsystems"
description = "96-channel Utah array"
location = "V1"
electrode_group = f.create_electrode_group(electrode_name,
                                           source=source,
                                           description=description,
                                           location=location,
                                           device=device)

# add electrode entries
for idx in [1, 2, 3, 4]:
    f.add_electrode(idx,
                    x=1.0, y=2.0, z=0.0,
                    imp=float(-idx),
                    location='V1', filtering='none',
                    description='channel %s' % idx, group=electrode_group)


electrode_table = ElectrodeTable('Utah-array_table')
electrode_table_region = ElectrodeTableRegion(electrode_table, [0, 2], 'the first and third electrodes')
###################################

ephys_ts = ElectricalSeries('test_ephys_data',
                            'an hypothetical source',
                            ephys_data,
                            electrode_table_region,
                            timestamps=ephys_timestamps,
                            # Alternatively, could specify starting_time and rate as follows
                            # starting_time=ephys_timestamps[0],
                            # rate=rate,
                            resolution=0.001,
                            comments="This data was randomly generated with numpy, using 1234 as the seed",
                            description="Random numbers generated with numpy.random.rand")



stimulus_ts = TimeSeries('Stimulus timeseries',
                         'VSG',
                         stimulus_data,
                         'degrees',
                         timestamps=ephys_timestamps,
                         resolution=0.001,
                         comments="This data was randomly generated with numpy, using 1234 as the seed",
                         description="Imagined orientation timeseries")
                         
f.add_acquisition(ephys_ts)
f.add_stimulus(stimulus_ts)

