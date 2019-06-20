function kalmanFilterInfo = kalmanReverseLinearMotion(kalmanFilterInfo,probDim)
%KALMANREVERSELINEARMOTION revese Kalman filter information in time
%
%SYNPOSIS kalmanFilterInfo = kalmanResMemLM(numFrames,numFeatures,probDim)
%
%INPUT  kalmanFilterInfo: Kalman filter information from a previous round
%                         of linking (as initialized by kalmanResMemLM).
%
%OUTPUT kalmanFilterInfo: Kalman filter information reversed in time.
%
%Khuloud Jaqaman, September 2008

%reverse time
kalmanFilterInfo = kalmanFilterInfo(end:-1:1);


%go over all frames and reverse velocity
for iFrame = length(kalmanFilterInfo) : -1 : 1
    kalmanFilterInfo(iFrame).stateVec(:,probDim+1:end) = ...
        -kalmanFilterInfo(iFrame).stateVec(:,probDim+1:end);
end
