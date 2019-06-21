% clear;
% for all 'cell1.txt'
% allData = dir('*cell1.txt');
clear;
% update for exposure time in the title

maxX = 200000;
maxXscale = 4300;

figure(1);clf;

% choose bg subtraction method. median, JQW (concentric circles; nomenclature cell101), or VS (adjacent circles; nomenclature cell201);

medianSubtracted = false;
JQWbg = true;

onlyBGrois=false;

% select true for 6 px ROIs ("cell2.txt")

smallROIs=false;

if smallROIs
    
    footer='_cell2.txt';
else
    footer='_cell1.txt';
end

% account for small difference eGFP-tagGFP2 22mer (x1.05)

eGFPtoTagGFP2=true;
eGFPtoTagGFP2ratio = 1.05;

allI=[];

% determined empirically
binSizes=[5000 16000 16000 16000];


expectedMolecs=[12 24 60 120];

% set the file name headers, depending on which day the expt was.

% Adding the option of combining experiments. Will need some kind of
% scaling, e.g. to mean of 60mer.

% dates

% allExpts= [92718, 100918, 102118];
allExptNames= [100918, 102118];
    
    figure(2); clf;
    figure(3); clf;
    figure(4); clf;
    
% all data from multiple expts, 
% scaled by the average 60mer intensity. 
  allExpts=[];  
for p = 1:length(allExptNames)
    expt = allExptNames(p);
    
    headers=[];
    
    exposureTimes=[];
    
    switch expt
        case 92718
            cd('/Volumes/MATT4/Microscopy data/from spinning disk/20180927 24 60 120 180mer motb/Quantification/Curve Images')

            headers.Twelve=[];
            headers.MotB{1}='MotB';
            headers.TwentyFour{1}='1B';
            headers.Sixty{1}='2A';
            headers.OneHundredTwenty{1}='3A';
            
            % hard code in the exposure times since they're not in the image
            % titles here.
            
            exposureTimes.Twelve=[];
            %         exposureTimes.TwentyFour
            exposureTimes.Sixty=250;
            exposureTimes.OneHundredTwenty=150;
            
        case 100918
            cd('/Volumes/MATT4/Microscopy data/from spinning disk/20181009 12 24 60 120mer motB/quantify/Curve Images')

            headers.Twelve{1}='1C';
            headers.Twelve{2}='3B';
            headers.MotB{1}='MotB';
            headers.TwentyFour{1}='   ';
            headers.Sixty{1}='3C';
            headers.OneHundredTwenty{1}='4A';
            
        case 102118
            
            cd('/Volumes/MATT4/Microscopy data/from spinning disk/20181021 cal curve motb AP21967 titration/Quantification/Curve Images')

            headers.Twelve{1}='1A';
            headers.Twelve{2}='2A';
            headers.Twelve{3}='3A';

            headers.TwentyFour{1}='1B';
            headers.TwentyFour{2}='2B';
            headers.TwentyFour{3}='3B';
            headers.Sixty{1}='4A';
            headers.OneHundredTwenty{1}='4B';

    end
    
    % first for JQW method (called 1001.txt)
    
    % write "true" if median subtracted already
    
    calCurve = [];
    meanIntensities = [];
    medianIntensities = [];
    stdIntensities =[];
    meanGaussianIntensities = [];
    stdGaussianIntensities = [];

    
    % iterate through different proteins
    for n = 1:length(expectedMolecs)
        % for n = 1 % just  one
        disp(n)
        allIoneSecond = [];
        allI = [];
        % reset allIoneSecond each condition (e.g. 60mer)
        allIoneSecond = [];
        curSignalData=[];
        curBGData=[];
        curSignalIntensity=[];
        curBGintensity=[];
        
        % for specific subset: choose first prefix
        
        
        switch expectedMolecs(n)
            
            
            case 12
                
                allData=[];
                for i = 1:length(headers.Twelve)
                    header=headers.Twelve{i}
                    allData = [allData; dir([header '*' footer])];
                end
                
            case 24
                
                allData=[];
                for i = 1:length(headers.TwentyFour)
                    header=headers.TwentyFour{i}
                    allData = [allData; dir([header '*' footer])];
                end
            case 60
                
                allData=[];
                for i = 1:length(headers.Sixty)
                    header=headers.Sixty{i}
                    allData = [allData; dir([header '*' footer])];
                end
                
                
            case 22
                
                allData=[];
                for i = 1:length(headers.MotB)
                    header=headers.MotB{i}
                    allData = [allData; dir([header '*' footer])];
                end
                
                %             exposureTime=exposureTimes(n);
            case 120
                
                allData=[];
                for i = 1:length(headers.OneHundredTwenty)
                    header=headers.OneHundredTwenty{i}
                    allData = [allData; dir([header '*' footer])];
                end
                
                
        end
        for i = 1:length(allData)
            
            curName = allData(i).name
            
            if medianSubtracted==1
                
                if ~contains(curName,'medianFilter')
                    continue;
                end
            else
                if contains(curName,'medianFilter')
                    continue;
                end
            end
            
            
            
            curSignalData = readtable(curName,'Delimiter','tab');
            curSignalIntensity = curSignalData.RawIntDen
            curSignalArea      = curSignalData.Area;
            
            % hard code in exposure time for the two cases with no exposure
            % time in title of file.
            
            if expt == 92718
                
                if expectedMolecs(n)==60
                    exposureTime=exposureTimes.Sixty
                elseif expectedMolecs(n)==120
                    exposureTime=exposureTimes.OneHundredTwenty
                else
                    msLocation= strfind(curName,'ms');
                    exposureTime=str2double(curName(msLocation-3:msLocation-1))
                end
            else
                
                % read in the exposure time from title. NEEDS to be 3 digits.
                
                msLocation= strfind(curName,'ms');
                exposureTime=str2double(curName(msLocation-3:msLocation-1))
            end
            
            
            % read exposure time from file name if "500ms" in file name
            
            %     exposureTime = 600; %ms
            
            %         switch expectedMolecs(n)
            %             case 24
            %
            %                 if contains(curName, '500ms','IgnoreCase',true)
            %                     disp(['500ms ' num2str(i)]);
            % %                     exposureTime=500;
            %                 else
            % %                     exposureTime = 400;
            %                 end
            %         end
            
            % if median subtracted, just measure directly
            
            if medianSubtracted
                %             Comment out these two lines if you want to VS bg subtract on top of median bg subtraction:
                curBGsubtractedIntensity = curSignalIntensity;
                
            else
                
                % if not bg subtracted yet, look for bg measurement
                
                if JQWbg
                    
                    
                    curBGname = [curName(1:end-4) '01.txt']
                    curBGData=readtable(curBGname,'Delimiter','tab');
                    curBGintensity = curBGData.RawIntDen;
                    curBGarea      = curBGData.Area;
                    
                    curBGsubtractedIntensity = curSignalIntensity- (curBGintensity-curSignalIntensity).*(curSignalArea)./(curBGarea-curSignalArea);
                else
                    % VS method
                    curBGname = [curName(1:end-5) '201.txt']
                    
                    curBGData=readtable(curBGname,'Delimiter','tab');
                    curBGintensity = curBGData.RawIntDen
                    curBGarea      = curBGData.Area;
                    if curBGarea~=curSignalArea
                        error('areas are different')
                    end
                    curBGsubtractedIntensity = curSignalIntensity-curBGintensity;
                    
                end
                
                
            end
            
            % account for eGFP vs tagGFP2
            
            if eGFPtoTagGFP2
                if expectedMolecs(n)==22
                    curBGsubtractedIntensity=curBGsubtractedIntensity*eGFPtoTagGFP2ratio;
                    disp('Converted 22mer eGFP to tagGFP2 intensity')
                end
            end
            
            figure(4);subplot(1, length(allData), i);
            histogram(curBGsubtractedIntensity*1000/exposureTime);
            xlim([0 19000])
            allIoneSecond = [allIoneSecond; curBGsubtractedIntensity*1000/exposureTime];

%             allI = [allI; curBGsubtractedIntensity];
            % correct for exposure time
            
%             allIoneSecond = allI*1000/exposureTime;
            % per field
            meanI = mean(curBGsubtractedIntensity*1000/exposureTime)
            medianI = median(curBGsubtractedIntensity*1000/exposureTime)
            stdI = std(curBGsubtractedIntensity*1000/exposureTime)
        % end of "datasets" loop
        end
        
        if length(allIoneSecond)<1
            continue;
        end
        
        % save intensity and mean median in variables
        calCurve(n).Intensities = allIoneSecond;
        meanIntensities(n)=mean(allIoneSecond);
        medianIntensities(n)=median(allIoneSecond);
        stdIntensities(n) = std(allIoneSecond);
        
        % fit to Gaussian, calculate 95% CI
        
        pd = fitdist(allIoneSecond,'Normal')
        conf95 = paramci(pd)
        
        meanGaussianIntensities(n) = pd.mu;
        stdGaussianIntensities(n)  = pd.sigma;
        
        calCurve(n).pd=pd;
        calCurve(n).ci95=conf95;
        
        figure(1+p); hold on;
        subplot(length(expectedMolecs),1,n); 
        %     histogram(allIoneSecond,0:binSizes(n):maxX);
        histogram(allIoneSecond);
        title(sprintf('%2.f ± %2.f (n = %.d) ', meanI, stdI, length(allIoneSecond)))
        xlim([-5000 maxX])

    end

    figure(1); clf;
    subplot(2,1,1)
    errorbar(expectedMolecs, meanIntensities, stdIntensities, 'ok', 'LineWidth',2);
    xlim([0 150]);
    ylim([0 maxX]);
    xlabel(['Molecules per structure']);
    ylabel(['Intensity per spot (AU)']);
    title('mean ± std')
    set(gca,'FontSize',11)
    
    subplot(2,1,2)
    errorbar(expectedMolecs, medianIntensities, stdIntensities,'ok', 'LineWidth',2);
    
    xlim([0 150]);
    ylim([0 maxX]);
    xlabel(['Molecules per structure']);
    ylabel(['Intensity per spot (AU)']);
    title('median ± std')
    set(gca,'FontSize',11)
    
    figure(p+5); clf;
    errorbar(expectedMolecs, meanGaussianIntensities, stdGaussianIntensities,'ok', 'LineWidth',2);
    
    % no offset. needs curve fitting toolbox
    [rsqZero, slopeZero] = plotSlopeThruZero(expectedMolecs,meanGaussianIntensities);
    
    % with offset. uses function 'regression'
    [rsqIntercept, slopeIntercept, offset]=plotSlope(expectedMolecs, meanGaussianIntensities);
    
    
    xlim([0 150]);
    ylim([0 maxX]);
    xlabel(['Molecules per structure']);
    ylabel(['Intensity per spot (AU)']);
    title('mean (Gaussian fit) ± std')
    set(gca,'FontSize',11)
    
    % scale by mean 60mer intensity
    
    mean60merIndex = find(expectedMolecs==60);
    mean60mer = mean(calCurve(mean60merIndex).Intensities);
    
    allExpts(p).calCurve=calCurve;
    
    for n = 1:length(expectedMolecs)
        
        scaledI = calCurve(n).Intensities.*1000./mean60mer;
        calCurve(n).scaledIntensities = scaledI;
        allExpts(p).calCurve(n).scaledIntensities=scaledI;

    end

        
end   
 
% initialize "calCurveExpts" structure for coallated data
calCurveExpts=[];
for p = 1:length(allExptNames)
    
    for n = 1:length(expectedMolecs)
        calCurveExpts(n).scaledIntensities = []; 
    end
end

% coallate all the data 
for p = 1:length(allExptNames)
    for n=1:length(expectedMolecs)
        calCurveExpts(n).scaledIntensities = [calCurveExpts(n).scaledIntensities; allExpts(p).calCurve(n).scaledIntensities];
    end
end


% plot histograms and curves of all the coallated data
figure(10); clf;
figure(12); clf;
for n=1:length(calCurveExpts)
    calCurveExpts(n).scaledIntensities;
    
    pd = fitdist(calCurveExpts(n).scaledIntensities,'Normal');
    conf95 = paramci(pd)
    
    meanGaussianIntensities_allexpts(n)  = pd.mu;
    stdGaussianIntensities_allexpts(n)  = pd.sigma;
    
    calCurveExpts(n).pd=pd;
    calCurveExpts(n).ci95=conf95;
    
    figure(10); hold on;
    subplot(length(expectedMolecs),1,n);
    %     histogram(allIoneSecond,0:binSizes(n):maxX);
    histogram(calCurveExpts(n).scaledIntensities);
    title(sprintf('%2.f ± %2.f (n = %.d) ', pd.mu, pd.sigma, length(calCurveExpts(n).scaledIntensities)))
%     xlabel('Intensity (A. U.) ');
    set(gca,'FontSize',15)
    % xlim([0 max(allIoneSecond)+0.1*max(allIoneSecond)])
    xlim([-200 maxXscale])
    
    figure(12); hold on; 
    histogram(calCurveExpts(n).scaledIntensities);

end

figure(12);
    histogram(calCurveExpts(n).scaledIntensities);
%     title(sprintf('%2.f ± %2.f (n = %.d) ', pd.mu, pd.sigma, length(calCurveExpts(n).scaledIntensities)))
    set(gca,'FontSize',15)
    xlim([-200 maxXscale])
    

% plot curve of all coallated data thru zero

figure(11); clf;
errorbar(expectedMolecs, meanGaussianIntensities_allexpts, stdGaussianIntensities_allexpts,'ok', 'LineWidth',2);

% no offset. needs curve fitting toolbox
[rsqZero_allexpts, slopeZero_allexpts] = plotSlopeThruZero(expectedMolecs,meanGaussianIntensities_allexpts)

% uncomment here if you want to plot all data points

% for n = 1:length(expectedMolecs)
% xPlot = repmat(expectedMolecs(n),length(calCurveExpts(n).scaledIntensities));
% hold on;
% plot(xPlot,calCurveExpts(n).scaledIntensities,'.b');
% end


 
    xlim([0 150]);
    ylim([0 maxXscale]);
    xlabel(['Molecules per structure']);
    ylabel(['Intensity per spot (AU)']);
    title('mean (Gaussian fit) ± std')
    set(gca,'FontSize',15)

 function [rsq, zeroSlope] = plotSlopeThruZero(xs,ys)
 [slopeInfo, gof, output]=fit(xs', ys', 'a*x');
 xsThruZero=0:max(xs);
 ysThruZero=xsThruZero.*slopeInfo.a;
 rsq=gof.rsquare;
 zeroSlope=slopeInfo.a;
 
 hold on; plot(xsThruZero, ysThruZero,'k','LineWidth',2');
 end
 
 
 function [rsq,slope,offset] = plotSlope(xs,ys)
 [rsq, slope, offset]=regression(xs, ys);
 xsFromZero=0:max(xs);
 ysFromZero=xsFromZero.*slope+offset;
 
 hold on; plot(xsFromZero,ysFromZero,'b','LineWidth',2');
 end
    
    
