%% Make spherical coordinates 

% this will take a given radius (and breadth of sphere) and nb points and
% give you a randomly generated set of points there. 

% in the right format, exported as a txt file. 

nbPoints = 200

fractionFilled = 0.6

orientationUpsideDown = true;

curRadius = 0.045

singleToAttach = 'hip1r';

saveToText = true;

% copy curRadius to nb datapoints. 

curRadii = [];
radii = repmat(curRadius,nbPoints,1);


azimuths = [];
% azimuth is the counterclockwise angle in the x-y plane measured in radians from the positive x-axis.
elevations = [];
% elevation is the elevation angle in radians from the x-y plane.

% randomly generate an azimuth and polarAngle

u = rand(nbPoints,1);

% v = rand(nbPoints,1); % full range
v = fractionFilled*rand(nbPoints,1); % partial sphere. 

% end

azimuths = 2*pi*u;

elevations = acos(2*v-1);


% elevation = 2*pi*rand(1)

% elevation = -pi*rand(1); % just hemisphere (lower)

% elevations = minElevation + (maxElevation-minElevation)*rand(nbPoints,1); % for arbitrary

% elevations(i) = elevation;



% or set it within a range. 
xs = [];
ys = [];
zs = [];

for i = 1:length(azimuths)
    
    xs(i) = curRadius * cos(azimuths(i))*sin(elevations(i));
    ys(i) = curRadius * sin(azimuths(i))*sin(elevations(i));
    zs(i) = curRadius * cos(elevations(i));
    
end

% flip sign if you want the fraction of sphere to be positive not neegative

if orientationUpsideDown
    
    zs = -zs;
end

% matlab version
% [xs,ys,zs] = sph2cart(azimuths, elevations, radii);


%% report it in a text file for cytosim to use

if saveToText
    
%   cd('/Users/makamats/Google Drive/drubin_lab/Modeling/cytosim_current/attach_coordinate_generation')
    fileID = fopen(['beadPointArray_' num2str(fractionFilled*100) ' x' num2str(nbPoints) '.txt'],'w');
      
    fprintf(fileID, 'point0 = 0 0 0, %3f \n', curRadius);

    for i = 1:nbPoints
        
        fprintf(fileID, 'point%i = %3f %3f %3f, , %s \n', i, xs(i), ys(i), zs(i), singleToAttach);
        
    end
    
    fclose(fileID);
end
figure(1); clf;
scatter3(xs*1000,ys*1000,-zs*1000, [], -zs*1000, 'filled')
xlim([-curRadius*1000, curRadius*1000])
ylim([-curRadius*1000, curRadius*1000])
zlim([-curRadius*1000, curRadius*1000])

set(gca,'FontSize',16)
xlabel('X (nm)')
ylabel('Y (nm)')
zlabel('Z (nm)')

view([-49 13])

title({[num2str(nbPoints) 'hip1R'];[num2str(fractionFilled*100) '% coverage']})


c = colorbar;
c.Label.String = 'Z (nm)';

caxis([-45 0])
% caxislabel('Z (nm)')

saveCurGraphs([num2str(nbPoints) 'hip1r ' num2str(fractionFilled*100) 'percent'],1)
% plot3D(z,y,z)
figure(2); clf;
subplot(2,3,1)
hist(elevations)
title('elevations')

subplot(2,3,2)
hist(azimuths)
title('azimuths')

subplot(2,3,3)
hist(radii)
title('radii')

subplot(2,3,4)
hist(zs)
title('zs')

subplot(2,3,5)
hist(xs)
title('xs')

subplot(2,3,6)
hist(ys)
title('ys')

 