import numpy as np
import tables as tb

animal = ('MM385','MM385','MM390')
penetr = ('P1','P2','P6')
fldrs  = ('20170519-033903','20170519-130912','20180207-230434')

results_root = '/home/lauri/projects/CorrelatedVariability/results/'
file_name    = results_root+'collated_correlation_data.h5'

# open file for read/write
data_file  = tb.open_file(file_name,'a')
# create group
data_group = data_file.create_group('/','data_group',"Data tables for each penetration") 

# loop through penetrations and animals
for p, val in enumerate(animal):
    # load processed Kilosort outputs
    spkC_L   = np.load('/opt3/'+animal[p]+'/'+penetr[p]+'/'+fldrs[p]+'/'+'Kilosorted_spkCnts_L.npy')
    spkC_NoL = np.load('/opt3/'+animal[p]+'/'+penetr[p]+'/'+fldrs[p]+'/'+'Kilosorted_spkCnts_NoL.npy')
    baseLine = np.load('/opt3/'+animal[p]+'/'+penetr[p]+'/'+fldrs[p]+'/'+'Kilosorted_baseLine.npy')
    spkR_L   = np.load('/opt3/'+animal[p]+'/'+penetr[p]+'/'+fldrs[p]+'/'+'Kilosorted_spkraster_L.npy')
    spkR_NoL = np.load('/opt3/'+animal[p]+'/'+penetr[p]+'/'+fldrs[p]+'/'+'Kilosorted_spkraster_NoL.npy')
    RF_L     = np.load('/opt3/'+animal[p]+'/'+penetr[p]+'/'+fldrs[p]+'/'+'Kilosorted_fieldParams_L.npy')
    RF_NoL   = np.load('/opt3/'+animal[p]+'/'+penetr[p]+'/'+fldrs[p]+'/'+'Kilosorted_fieldParams_NoL.npy')
    SNR      = np.load('/opt3/'+animal[p]+'/'+penetr[p]+'/'+fldrs[p]+'/'+'Kilosorted_SNR.npy')
    depth    = np.load('/opt3/'+animal[p]+'/'+penetr[p]+'/'+fldrs[p]+'/'+'Kilosorted_spikeDepths.npy')
    spkDur   = np.load('/opt3/'+animal[p]+'/'+penetr[p]+'/'+fldrs[p]+'/'+'Kilosorted_spikeWidths.npy')
    # gather stimulus info
    diams              = np.load('/opt3/'+animal[p]+'/'+penetr[p]+'/'+fldrs[p]+'/'+'diams.npy')
    contrast           = np.load(
    orientation        = np.load(
    spatial_frequency  = np.load(
    temporal_frequency = np.load(
    
    
    # prepare table
    # needs to be done in the loop to accommodate different
    # array sizes in different penetrations
    dataTable = {'spkC_L':tb.Float64Col(shape=(spkC_L.shape[1],spkC_L.shape[2])),
                 'spkC_NoL':tb.Float64Col(shape=(spkC_NoL.shape[1],spkC_NoL.shape[2])),
                 'baseLine':tb.Float64Col(shape=(baseLine.shape[1])),
                 'spkR_L':tb.Float64Col(shape=(spkR_L.shape[1],spkR_L.shape[2],spkR_L.shape[3])),
                 'spkR_NoL':tb.Float64Col(shape=(spkR_NoL.shape[1],spkR_NoL.shape[2],spkR_NoL.shape[3])),
                 'RF_L':tb.Float64Col(shape=(RF_L.shape[1])),
                 'RF_NoL':tb.Float64Col(shape=(RF_NoL.shape[1])),
                 'SNR':tb.Float64Col(1),
                 'depth':tb.Float64Col(1),
                 'layer':tb.StringCol(3),
                 'spkDur':tb.Float64Col(1),
                 'diams':tb.Float64Col(shape=(diams.shape[1]))}

    # create data table
    table = data_file.create_table(data_group, animal[p]+penetr[p], dataTable, 'Preprocessed spike-sorted data')

    # write unit data to table
    unit = table.row
    for u in range(spkC_L.shape[0]):
        unit['spkC_L']   = spkC_L[u,:,:]
        unit['spkC_NoL'] = spkC_NoL[u,:,:]
        unit['baseLine'] = baseLine[u]
        unit['spkR_L']   = spkR_L[u,:,:,:]
        unit['spkR_NoL'] = spkR_NoL[u,:,:,:]
        unit['RF_L']     = RF_L[u,:]
        unit['RF_NoL']   = RF_NoL[u,:]
        unit['SNR']      = SNR[u]
        unit['depth']    = depth[u]
        unit['spkDur']   = spkDur[u]
        unit['diams']   = diams[0,:]
        # append row
        unit.append()

    # flush data
    table.flush()

# close file
data_file.close()
