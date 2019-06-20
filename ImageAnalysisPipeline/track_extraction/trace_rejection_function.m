%% Function for trace rejection scheme
% Josh Feruson
% Kural Group
% The Ohio State University
% Sep 21, 2016
function accepted_allCurrentTracks = trace_rejection_function(allCurrentTracks, lifetimesCutoff, varargin)
%% Parameters

%nfr = number of frames
%rthresh = threshold r squared value
%
switch nargin
    case 2
        nfr = 4; %number of frames to consider in line fit
        rthresh = 0.5; %threshold r-squared. If any line within the trace of length
    case 3           %nfr exceeds an r-squared of rthresh it is NOT rejected
        nfr = varargin{1};
        rthresh = 0.5;
    case 4
        nfr = varargin{1};
        rthresh = varargin{2};
end
%% Main

accept = false(1,length(allCurrentTracks));
rejected = false(1,length(allCurrentTracks));
for i = 1:length(allCurrentTracks)
    int = allCurrentTracks(i).Is;
    
%     disp(int);
    
    % only consider tracks >= the number of FRAMES dictated by "lifetimescutoff". 
    
    if length(int)< lifetimesCutoff 
        continue; 
    end
    for j = 1:length(int)-nfr
%         fprintf('i %d j %d \n', i, j)

%%% I changed the "i" to "j"
%         tmpy = int(j:j+nfr);
          tmpy = int(j:j+nfr)/max(int);

%         tmpx = 1:nfr;

%%% I changed the length of tmpx so it matches tmpy, and transposed the
%%% vector
        tmpx = (1:length(tmpy))';
        
        
        
%%% I changed "tmpc" to "tmpx"
        p = polyfit(tmpx,tmpy,1); 
        rsq = 1-SSE(polyval(p,tmpx),tmpy)/SST(tmpy);
        if rsq > rthresh
            accept(i) = true;
            % add "break" so that it moves on if there's a 4-frame linear
            % element ANYWHERE in the track 
            % fprintf([num2str(i) ' accepted \n']);
            
            break;
            
        else
            rejected(i) = true;
%             fprintf([num2str(i) ' rejected \n']);

        end
    end
end

accepted_allCurrentTracks = allCurrentTracks(accept);

% display nb rejected tracks and accepted (the ones longer than the min
% length).

fprintf('%d erratic tracks rejected, ', (sum(rejected)))
fprintf('%d tracks accepted. \n', (sum(accept)))

% fprintf('%d accepted tracks and %d rejected tracks', sum(accept), length(accept)-sum(accept))

%% temporarily plot a few tracks accepted and rejected
            
%         for j = 1:length(accept)
% 
%             if accept(j)
%                 
%             figure(1); hold all;
%             plot(allCurrentTracks(j).xs, allCurrentTracks(j).ys, '.-')
% %             plot(1:length(allCurrentTracks(j).Is), allCurrentTracks(j).Is)
%             else
%                 
%             figure(2); hold all;
%             plot(allCurrentTracks(j).xs, allCurrentTracks(j).ys, '.-')
% %             plot(1:length(allCurrentTracks(j).Is), allCurrentTracks(j).Is)
% 
%             end
%             
%         end
%           figure(1); title('Accepted tracks')
%           figure(2); title('Rejected tracks')
%          

% save accepted_allCurrentTracks.mat accepted_allCurrentTracks
end

function sse = SSE(y_r,y_d)
%SSE Sum of the Squares about the Estimate
sse = 0;
ly = length(y_d);
for i = 1:ly
    sse = sse + (y_r(i)-y_d(i))^2;
end
end

function sst = SST(y_d)
%SST Sum of the Squares about the mean
sst = 0;
ly = length(y_d);
my = sum(y_d)/ly; %mean of the distribution
for i = 1:ly
    sst = sst + (y_d(i)-my)^2;
end
end