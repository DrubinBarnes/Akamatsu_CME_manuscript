function kalmanFilterInfo = kalmanResMemLM(numFrames,numFeatures,probDim)
%KALMANRESMEMLM reserves memory for Kalman filter structure for linear motion model
%
%SYNPOSIS kalmanFilterInfo = kalmanResMemLM(numFrames,numFeatures,probDim)
%
%INPUT   numFrames  : Number of frames in movie.
%        numFeatures: An array with number of feaures in each frame.
%        probDim    : Problem dimensionality.
%
%OUTPUT   kalmanFilterInfo: Structure array with number of entries equal to 
%                           number of frames in movie. Contains the fields:
%             .stateVec        : Kalman filter state vector for each
%                                feature in frame.
%             .stateCov        : Kalman filter state covariance matrix
%                                for each feature in frame.
%             .noiseVar        : Variance of state noise for each
%                                feature in frame.
%             .stateNoise      : Estimated state noise for each feature in
%                                frame.
%             .scheme          : 1st column: propagation scheme connecting
%                                feature to previous feature. 2nd column:
%                                propagation scheme connecting feature to
%                                next feature.
%Khuloud Jaqaman, March 2007

%calculate vector size
vecSize = 2 * probDim;

%go over all frames and reserve memory
for iFrame = numFrames : -1 : 1

    kalmanFilterInfo(iFrame) = struct('stateVec',zeros(numFeatures(iFrame),vecSize),...
        'stateCov',zeros(vecSize,vecSize,numFeatures(iFrame)),...
        'noiseVar',zeros(vecSize,vecSize,numFeatures(iFrame)),...
        'stateNoise',zeros(numFeatures(iFrame),vecSize),...
        'scheme',zeros(numFeatures(iFrame),2));

end
