function KilosortExtractSpikeCountsRasters()
addpath(genpath('/home/lauri/matlab_functions/NPMK'));
addpath(genpath('/home/lauri/code/npy-matlab'));
addpath(genpath('/home/lauri/code/spikes'));
addpath('/home/lauri/matlab_functions/statisticals');

errb = 1;
clc; %close all

sDir = uigetdir('/opt3/');% so here you need to choose each cerebus folder
sDir = [sDir, '/'];

% array layout is the array alignment so it should be load 
load('arrayLayoutV_Probe.mat')
Mar = dir(fullfile(sDir,'stimParams*.mat'));
load([sDir,Mar(1).name])
% get pointer to dat-file
datf = dir(fullfile(sDir,'*.dat'));

% in ms defines how much before and after stimulus onset and offset
% is taken for spike-rasters
raster_trail = 100;

%% spikes 
ksf = dir(fullfile(sDir,'rez.mat'));
kilorez = load([sDir,ksf(1).name]);
NEVfreq = kilorez.rez.ops.fs;
% load spike times 
sptf = dir(fullfile(sDir,'spike_times.npy'));
spt  = readNPY([sDir,sptf(1).name]);
% load cluster identities 
cluf = dir(fullfile(sDir,'spike_clusters.npy'));
clu  = readNPY([sDir,cluf(1).name]);

% open continuous TTL channel record
a(1) = dir(fullfile(sDir,'20*.ns4'));
DataStruct  = openNSxNew([sDir,a(1).name],'read','c:1');
Fs          = DataStruct.MetaTags.SamplingFreq;
% generate pulse sample times from continuous trigger channel
pulseTimes  = PulseCounter(DataStruct.Data,20000);

%% collect stimulus info
stimOnCount    = round((stimParams.on_time_ms/1000)*Fs);
stimOffCount   = round((stimParams.off_time_ms/1000)*Fs);
stimTotalCount = stimOnCount+stimOffCount;
stimMatrix_tmp = stimParams.stimOrder;
stimMatrix     = ones(size(stimMatrix_tmp,1)/2,size(stimMatrix_tmp,2))*NaN;
% every stimulus is run twice in succession with and without the laser
for i = 1:size(stimMatrix_tmp,2)
    stimMatrix(:,i)     = unique(stimMatrix_tmp(:,i),'stable');
end
nStims    = size(stimParams.stimOrder,1);
nTrials   = size(stimParams.stimOrder,2);
% div by 2 because just every second stimulus produces a trigger
stimTimes = reshape(pulseTimes,nStims/2,nTrials); 

%% calculate number of spikes
sampFreqNEV   = double(kilorez.rez.ops.fs);
stimTimesNEV  = round(stimTimes*(sampFreqNEV/Fs));
stimLengthNEV = round(stimOnCount*(sampFreqNEV/Fs));
stimOffNEV    = round(stimOffCount*(sampFreqNEV/Fs));
% Use time after last stimulus for blank response
blnkTimes = stimTimesNEV(end,:) + stimLengthNEV + 0.1*30e3;

% conf laser info
Laser_Mat      = reshape(stimParams.laser.rnd_laser',2,(size(stimParams.laser.rnd_laser,2)*size(stimParams.laser.rnd_laser,1))/2);
Laser_Mat(2,:) = [];

N_units       = length(unique(clu));
units         = unique(clu);

spikeCnts_L   = NaN*ones(N_units, nStims/2, nTrials);
spikeCnts_NoL = NaN*ones(N_units, nStims/2, nTrials);
baseLine      = NaN*ones(N_units, nStims/2, nTrials);
spkraster_L   = NaN*ones(N_units, nStims/2, nTrials, 700);
spkraster_NoL = NaN*ones(N_units, nStims/2, nTrials, 700);

for u = 1:N_units
    spikeTimeMat = double(spt(units(u) == clu));
    
    % collate spike-counts depending on weather the laser was on or off
    for f = 1:nStims/2
        idx = find(stimMatrix==f);        
        for g = 1:length(idx)

            % laser was on during the first stimulus of the pair
            if Laser_Mat(idx(g)) == 1
                spk_tmp = (spikeTimeMat>(stimTimesNEV(idx(g)) - 0.1*30e3) & spikeTimeMat<stimTimesNEV(idx(g))+stimLengthNEV + 0.1*30e3);
                spk_tmp = (double(spikeTimeMat(spk_tmp) - stimTimesNEV(idx(g))))./30e3*1000;
                spkraster_L(u,f,g,:) = histcounts(spk_tmp,[-raster_trail:1:stimParams.on_time_ms+raster_trail]);
                spikeCnts_L(u,f,g)   = sum(spikeTimeMat>stimTimesNEV(idx(g)) & spikeTimeMat<stimTimesNEV(idx(g))+stimLengthNEV);

                spk_tmp = (spikeTimeMat>(stimTimesNEV(idx(g)) + stimLengthNEV + stimOffNEV - 0.1*30e3) & (spikeTimeMat<stimTimesNEV(idx(g))+ 2*stimLengthNEV + stimOffNEV + 0.1*30e3));
                spk_tmp = (double(spikeTimeMat(spk_tmp) - (stimTimesNEV(idx(g)) + stimLengthNEV + stimOffNEV)))./30e3*1000;
                spkraster_NoL(u,f,g,:) = histcounts(spk_tmp,[-raster_trail:1:stimParams.on_time_ms+raster_trail]);
                spikeCnts_NoL(u,f,g) = sum(spikeTimeMat>(stimTimesNEV(idx(g)) + stimLengthNEV + stimOffNEV) & spikeTimeMat<stimTimesNEV(idx(g))+ 2*stimLengthNEV + stimOffNEV);

            % laser was on during the second stimulus of the pair
            else
                spk_tmp = (spikeTimeMat>(stimTimesNEV(idx(g)) - 0.1*30e3) & (spikeTimeMat<stimTimesNEV(idx(g))+stimLengthNEV + 0.1*30e3));
                spk_tmp = (double(spikeTimeMat(spk_tmp) - stimTimesNEV(idx(g))))./30e3*1000;
                spkraster_NoL(u,f,g,:) = histcounts(spk_tmp,[-raster_trail:1:stimParams.on_time_ms+raster_trail]);
                spikeCnts_NoL(u,f,g)   = sum(spikeTimeMat>stimTimesNEV(idx(g)) & spikeTimeMat<stimTimesNEV(idx(g))+stimLengthNEV);
                
                spk_tmp = (spikeTimeMat>(stimTimesNEV(idx(g)) + stimLengthNEV + stimOffNEV - 0.1*30e3) & spikeTimeMat<stimTimesNEV(idx(g))+ 2*stimLengthNEV + stimOffNEV + 0.1*30e3);
                spk_tmp = (double(spikeTimeMat(spk_tmp) - (stimTimesNEV(idx(g)) + stimLengthNEV + stimOffNEV)))./30e3*1000;
                spkraster_L(u,f,g,:) = histcounts(spk_tmp,[-raster_trail:1:stimParams.on_time_ms+raster_trail]);
                spikeCnts_L(u,f,g)   = sum(spikeTimeMat>(stimTimesNEV(idx(g)) + stimLengthNEV + stimOffNEV) & spikeTimeMat<stimTimesNEV(idx(g))+ 2*stimLengthNEV + stimOffNEV);

            end        
            baseLine(u,f,g) = sum(spikeTimeMat > stimTimesNEV(idx(g))+stimLengthNEV+1.2*30e3 & spikeTimeMat<stimTimesNEV(idx(g))+stimLengthNEV+1.5*30e3);
        end
    end
end


% mean spike counts with and without the laser
mnSpikeCnts_L   = mean(spikeCnts_L,3)*(1000/stimParams.on_time_ms);
mnSpikeCnts_NoL = mean(spikeCnts_NoL,3)*(1000/stimParams.on_time_ms);
baseLine = mean(mean(baseLine,3),2)*(1000/300);

for u = 1:size(mnSpikeCnts_NoL,1)
    fprintf('Now fitting unit %d\n', u)
    if u == 1
        f = fit_sizetuning(mnSpikeCnts_NoL(u,:),unique(stimParams.Value)');
        fit_params     = f.final_params;
        field_params   = f.field_params;

        f = fit_sizetuning(mnSpikeCnts_L(u,:),unique(stimParams.Value)');
        fit_params_L   = f.final_params;
        field_params_L = f.field_params;
    else
        f = fit_sizetuning(mnSpikeCnts_NoL(u,:),unique(stimParams.Value)');
        fit_params     = [fit_params; f.final_params];
        field_params   = [field_params; f.field_params];
        
        f = fit_sizetuning(mnSpikeCnts_L(u,:),unique(stimParams.Value)');
        fit_params_L = [fit_params_L; f.final_params];
        field_params_L = [field_params_L; f.field_params];
    end
end

% get spike info
tempsUnW       = readNPY([sDir,'templates_unw.npy']);
spikeTemplates = readNPY([sDir,'spike_templates.npy']);
temps          = readNPY([sDir,'templates.npy']);
WInv           = readNPY([sDir,'whitening_mat_inv.npy']);
coords         = readNPY([sDir,'channel_positions.npy']);
amplitudes     = readNPY([sDir,'amplitudes.npy']);

% compute spike-widths in milliseconds using cortex-lab script
[spikeWidths, tempWidths] = computeSpikeWidths(tempsUnW, spikeTemplates);
mn_spikeWidths = NaN*ones(size(units));
for u = 1:length(units)
    mn_spikeWidths(u) = abs(mean(spikeWidths(clu == units(u))./NEVfreq*1000));
end

% get waveforms and depths
mn_spikeDepths = NaN*ones(size(units));
[spikeAmps, spikeDepths, templateDepths, tempAmps, tempsUnW, templateDuration, waveforms] = templatePositionsAmplitudes(temps, WInv, coords(:,2), spikeTemplates, amplitudes);
for u = 1:length(units)
    mn_spikeDepths(u) = mean(spikeDepths(clu == units(u)));
end

% define some parameters for SNR calculation
gwfparams.dataDir = sDir;
gwfparams.dataType = 'int16';   % Data type of .dat file (this should be BP filtered)
gwfparams.nCh      = 32;        % Number of channels that were streamed to disk in .dat file
gwfparams.wfWin    = [-40 41];  % Number of samples before and after spiketime to include in waveform
gwfparams.nWf      = 10e3;      % Max number of waveforms per unit to pull out
gwfparams.fshigh   = 300;
gwfparams.fs       = 30e3;

fileName = fullfile(sDir,datf.name);
filenamestruct = dir(fileName);
dataTypeNBytes = numel(typecast(cast(0, gwfparams.dataType), 'uint8'));
gwfparams.nSamp = filenamestruct.bytes/(gwfparams.nCh*dataTypeNBytes); 
mmf = memmapfile(fileName, 'Format', {gwfparams.dataType, [gwfparams.nCh gwfparams.nSamp], 'x'});

fprintf('Constructing high-pass filter\n')
[b1, a1] = butter(3, gwfparams.fshigh/gwfparams.fs*2, 'high');
Data = filter(b1,a1,mmf.Data.x,[],2);

fprintf('Extracting SNR\n')
SNR = NaN*ones(size(units));

for u = 1:length(units)
    gwfparams.spikeTimes = spt(clu == units(u));
    gwfparams.spikeClusters = repmat(units(u), size(gwfparams.spikeTimes));
    SNR(u) = getSNR(gwfparams,Data);
end

diams = unique(stimParams.Value)';
writeNPY(spikeCnts_L, [sDir,'Kilosorted_spkCnts_L.npy']);
writeNPY(spikeCnts_NoL, [sDir,'Kilosorted_spkCnts_NoL.npy']);
writeNPY(baseLine, [sDir,'Kilosorted_baseLine.npy']);
writeNPY(spkraster_L, [sDir,'Kilosorted_spkraster_L.npy']);
writeNPY(spkraster_NoL, [sDir,'Kilosorted_spkraster_NoL.npy']);
writeNPY(fit_params_L, [sDir,'Kilosorted_fitParams_L.npy']);
writeNPY(fit_params, [sDir,'Kilosorted_fitParams_NoL.npy']);
writeNPY(field_params_L, [sDir,'Kilosorted_fieldParams_L.npy']);
writeNPY(field_params, [sDir,'Kilosorted_fieldParams_NoL.npy']);
writeNPY(mn_spikeWidths, [sDir,'Kilosorted_spikeWidths.npy']);
writeNPY(mn_spikeDepths, [sDir,'Kilosorted_spikeDepths.npy']);
writeNPY(SNR, [sDir,'Kilosorted_SNR.npy']);
writeNPY(diams, [sDir,'diams.npy']);
writeNPY(stimParams.TF, [sDir,'TF.npy']);
writeNPY(stimParams.SF, [sDir,'SF.npy']);
writeNPY(stimParams.Contrast, [sDir,'contrast.npy']);