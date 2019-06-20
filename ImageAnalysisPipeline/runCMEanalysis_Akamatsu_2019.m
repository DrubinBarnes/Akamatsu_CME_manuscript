%% runCMEanalysis

%% Matt Akamatsu, Drubin lab, Jan 2016

% Run Danuser lab cmeAnalysis.m (Aguet 2013) with defined parameters

clear; 

% set gap length, in frames. recommended is 3 frames for a 2s interval movie. But not optimized yet. 
% Try different values and see how your data looks!! (Run manual tracking
% for a small subset of your data, and see if it looks like multiple tracks are getting stitched together.)

gapLength = 0; 

% set tracking search radius min and max, in pixels. 2.3 pixels is 250nm. 3.2 pixels is 350 nm.


trackingRadiusMin = 0; 
trackingRadiusMax = 2.3; %2; % 3.2;

curfilename = input('What do you want to call this file? ', 's');

manual = input('use manual or automatic track association? 1 for manual, 0 for automatic: ');

rfpRefChannel = input('Is RFP the reference channel? 1 for yes, 0 for no: ');

if rfpRefChannel
    disp('You chose RFP channel to be the reference channel, and GFP to be the corresponding channel.')
else
    disp('You chose GFP channel to be the reference channel, and RFP to be the corresponding channel.')

end

% These parameters usually don't change for our imaging setup.

% objectiveNA = 1.49;
objectiveNA = 1.45;

% objectiveMagnification = 60; 
objectiveMagnification = 100;

% cameraPixelSize = 6.45e-6; % in meters (sCMOS?)
cameraPixelSize = 16e-6; % in meters. EMCCD?
%cameraPixelSize = 1.09e-7; % in meters

% now cmeAnalysis will run even if you've already run it before (won't get
% caught up on the associated_tracks.mat file)

trackingData = loadConditionData('Parameters', [objectiveNA objectiveMagnification cameraPixelSize], 'IgnoreEmptyFolders', true);
for i = 1:length(trackingData)
    trackingData(i).curfilename = curfilename;
    trackingData(i).rfpRefChannel = rfpRefChannel;
end

 [res, data] = cmeAnalysis('Parameters', [objectiveNA objectiveMagnification cameraPixelSize], 'PlotAll', false, 'TrackingGapLength',...
        gapLength, 'TrackingRadius', [trackingRadiusMin trackingRadiusMax], 'data', trackingData);

numCells = length(trackingData);
numChannels = length(trackingData(1).markers)

for i = 1:numChannels
     cmeAnalysis('Parameters', [objectiveNA objectiveMagnification cameraPixelSize], 'PlotAll', false, 'TrackingGapLength',...
        gapLength, 'TrackingRadius', [trackingRadiusMin trackingRadiusMax], 'data', trackingData);
    % update the data for each cell that we have to run this on other
    % channels
    for j = 1:numCells
        trackingData(j).channels = circshift(trackingData(j).channels, 1);
        trackingData(j).markers = circshift(trackingData(j).markers, 1);
        trackingData(j).framePaths = circshift(trackingData(j).framePaths, 1);
        trackingData(j).source = trackingData(j).channels{1};  
    end
    
end

totalNbAssocPersistentRef = 0;
totalNbAssocPersistentCor = 0;
totalNbAssocComplete = 0;

%PersistentTracks creats table that provides the information of persistent tracks of all
%cells [number_ref_assoc_persistant, number_cor_assoc_persistant, complete_associated_nb]
PersistentTracks = cell(numCells,9)

for i = 1:numCells
    
    %i 
    
    for j = 1:numChannels
        extractXYcoordinates_useCellMask(trackingData(i));
        trackingData(i).channels = circshift(trackingData(i).channels, 1);
        trackingData(i).markers = circshift(trackingData(i).markers, 1);
        trackingData(i).framePaths = circshift(trackingData(i).framePaths, 1);
        trackingData(i).source = trackingData(i).channels{1};  
    end
    
    [ref_assoc_persistant_nb, cor_assoc_persistant_nb, complete_associated_nb,complete_associated_nb_woHighLow,complete_associated_nb_woPersistent,...
        end_ref_unassociated_nb,begin_ref_unassociated_nb,end_cor_unassociated_nb,begin_cor_unassociated_nb,begin_associated_nb,...
        end_associated_nb,complete_unassociated_ref,complete_unassociated_cor] = Associate_tracks_20151221(trackingData(i));
   
   
    
    ref_assoc_persistant_nb
    cor_assoc_persistant_nb
    
    complete_associated_nb
    complete_associated_nb_woHighLow
    complete_associated_nb_woPersistent
    
    end_ref_unassociated_nb 
    begin_ref_unassociated_nb
    end_cor_unassociated_nb 
    begin_cor_unassociated_nb
    
    begin_associated_nb 
    end_associated_nb
    %complete_unassociated_ref
    %complete_unassociated_cor
    
     PersistentTracks{i,1} =   ref_assoc_persistant_nb
     PersistentTracks{i,2} =  cor_assoc_persistant_nb
     PersistentTracks{i,3} =  complete_associated_nb
     PersistentTracks{i,4} =  complete_associated_nb_woHighLow
     PersistentTracks{i,5} =  complete_associated_nb_woPersistent
     PersistentTracks{i,6} = end_ref_unassociated_nb 
     PersistentTracks{i,7} = begin_ref_unassociated_nb
     PersistentTracks{i,8} = end_cor_unassociated_nb 
     PersistentTracks{i,9} =  begin_cor_unassociated_nb
     PersistentTracks{i,10} = begin_associated_nb 
     PersistentTracks{i,11} = end_associated_nb
     PersistentTracks{i,12} = complete_unassociated_ref
     PersistentTracks{i,13} = complete_unassociated_cor
    
    %totalNbAssocPersistentRef = totalNbAssocPersistentRef + number_ref_assoc_persistant;
   % totalNbAssocPersistentCor = totalNbAssocPersistentCor + number_cor_assoc_persistant;
    %totalNbAssocComplete = totalNbAssocComplete + complete_associated_nb;
    
    if manual == 1
        Manually_pick_asso_ref_cor20150312(trackingData(i));
    else
%         Clean_associated_tracks2015_03_09(trackingData(i));
%       Automatic track rejection: custom criteria (based on matching with
%       manual rejection)
        Clean_associated_tracks_MA(trackingData(i));

    end
    
     %save([curfilename '.mat']);
   
end

NumberOfTracks = cell2table(PersistentTracks,'VariableNames',{'ref_assoc_persistant_nb'...
    'cor_assoc_persistant_nb' 'complete_associated_nb' 'complete_associated_nb_woHighLow'...
    'complete_associated_nb_woPersistent' 'end_ref_unassociated_nb' 'begin_ref_unassociated_nb'...
    'end_cor_unassociated_nb' 'begin_cor_unassociated_nb' 'begin_associated_nb'...
    'end_associated_nb' 'complete_unassociated_ref' 'complete_unassociated_cor'  })

%totalNbAssocPersistentRef
%totalNbAssocPersistentCor
%totalNbAssocComplete

disp('done with all cells');



