% Create paramters to call on AnDarksamtest
% 5/31/16
% CD
% Reads in lifetime_associated data and the field_positions to identify
% fields within the lifetime_associated data. This information can then be
% used to compare stats with an Anderson Darling test using the
% "AnDarksamtest.m" function. field_positions and lifetime_associated comes
% from Plot_stats.m function/Program.

% DONE: To do: show histograms
% DONE: To do: make it read different lifetime_associated sets
% DONE: To do: compare all combos using AD tedst
% DONE: To do: add all figures to subplots for GFP and RFP subplot(3,1,i)

close all
clear
clc

% nassoc = input('How many sets of lifetime_associated vectors would you like to read in? ');
 nassoc = input('How many experiments would you like to compare?');
        experiments = cell(nassoc);
        for n = 1:nassoc
            experiments{n} = uigetdir('Select next experiment');
        end
        if (nassoc == 1)
        end
if nassoc > 1
    stats_popups=0;
else
    stats_popups=1;
end

if nassoc == 1
    
    [field_positions, lifetime_associated] = Plot_stats20151006();
    close all
    
    field_number=size(field_positions,1);
    andarl_gfp=zeros(size(lifetime_associated));
    andarl_rfp=zeros(size(lifetime_associated));
    
    field_start=zeros(size(field_positions));%find starting positions of each field
    
    index=0;
    for i=1:field_number
        field_start(i)=index+1; %add starting position
        for j=1:field_positions(i,1)
            index=index+1;
            andarl_gfp(index,1)=lifetime_associated(index,1);
            andarl_gfp(index,2)=i;
            andarl_rfp(index,1)=lifetime_associated(index,2);
            andarl_rfp(index,2)=i;
        end
    end
    
    
    alpha = 0.05;
    
    stringa=['For all samples the Anderson Darling statistics, or chance',...
        ' that these data \n are drawn from the same distribution is: '];
    P_GFP=AnDarksamtest(andarl_gfp,alpha);
    P_RFP=AnDarksamtest(andarl_rfp,alpha);
    stringb=[stringa, sprintf('for GFP %.4f and for RFP %.4f',P_GFP,P_RFP)];
    disp(stringb);
    
    a=1:field_number; %generate unique pairs of fields to test
    b=1:field_number;
    [p,q]=meshgrid(a,b);
    mask = triu(ones(field_number), 1) > 0.5;
    pairs=[p(mask) q(mask)];
    
    
    for prz=1:size(pairs,1)
        current_pair=pairs(prz,:);
        num1=current_pair(1);
        num2=current_pair(2);
        current_gfp=[andarl_gfp(field_start(num1):field_start(num1)+field_positions(num1)-1,1) ...
            ; andarl_gfp(field_start(num2):field_start(num2)+field_positions(num2)-1,1)];
        current_gfp=[current_gfp, ones(size(current_gfp,1),1)];
        current_gfp(field_positions(num1)+1 : field_positions(num2)...
            +field_positions(num1),2)=2;
        
        current_rfp=[andarl_rfp(field_start(num1):field_start(num1)+field_positions(num1)-1,1) ...
            ; andarl_rfp(field_start(num2):field_start(num2)+field_positions(num2)-1,1)];
        current_rfp=[current_rfp, ones(size(current_rfp,1),1)];
        current_rfp(field_positions(num1)+1 : field_positions(num2)...
            +field_positions(num1),2)=2;
        P_GFP=AnDarksamtest(current_gfp,alpha);
        P_RFP=AnDarksamtest(current_rfp,alpha);
        stringc=sprintf('For fields %d and %d the AD stat is %.4f for GFP and %.4f for RFP',num1,num2,P_GFP,P_RFP);
        disp(stringc);
    end
    
    figure
    hindex=0;
    for i=1:field_number %go between fields
        h.count=0; %gfp
        h.range=[0, max(max(lifetime_associated))+1];
        h.binwidth=1; %change this to change resolution of histogram
        h2.count=0; %rfp
        h2.range=[0, max(max(lifetime_associated))+1];
        h2.binwidth=1;
        for k=1:field_positions(i,1) %advance within each field
            hindex=hindex+1;
            h=histo(h,lifetime_associated(hindex,1)); %bin in histo for gfp
            h2=histo(h2,lifetime_associated(hindex,2)); %bin in histo2 for rfp
        end
        string1=sprintf('Histogram of GFP Lifetimes in Field %d', i);
        string2=sprintf('Histogram of RFP Lifetimes in Field %d', i);
        
        subplot(field_number,1,i);
        bar(h.vals,h.hist);
        title(string1);
        xlabel('Lifetime in seconds');
        ylabel('Frequency');
    end
    
    figure
    hindex=0;
    for i=1:field_number %go between fields
        h.count=0; %gfp
        h.range=[0, max(max(lifetime_associated))+1];
        h.binwidth=1; %change this to change resolution of histogram
        h2.count=0; %rfp
        h2.range=[0, max(max(lifetime_associated))+1];
        h2.binwidth=1;
        for k=1:field_positions(i,1) %advance within each field
            hindex=hindex+1;
            h=histo(h,lifetime_associated(hindex,1)); %bin in histo for gfp
            h2=histo(h2,lifetime_associated(hindex,2)); %bin in histo2 for rfp
        end
        string1=sprintf('Histogram of GFP Lifetimes in Field %d', i);
        string2=sprintf('Histogram of RFP Lifetimes in Field %d', i);
        
        subplot(field_number,1,i);
        bar(h2.vals,h2.hist);
        title(string2);
        xlabel('Lifetime in seconds');
        ylabel('Frequency');
    end
end


%% For multiple lifetime_associated sets

if nassoc > 1
    
    for assx=1:nassoc
%         [field_positions, lifetime_associated] = Plot_stats20151006();
        
%         close all
        currentLifetimes =  open(fullfile(experiments{assx},'lifetimeData.mat'));
        s(assx).la=currentLifetimes.lifetime_associated;
    end
    
    andarl_gfp=[];
    andarl_rfp=[];
    
    sindex=0;
    la_positions=zeros(1,nassoc);
    la_start=zeros(1,nassoc);
    
    for i=1:nassoc
        la_start(i)=sindex+1;
        andarl_gfp=[andarl_gfp; s(i).la(:,1), i*ones(size(s(i).la,1),1)];
        andarl_rfp=[andarl_rfp; s(i).la(:,2), i*ones(size(s(i).la,1),1)];
        sindex=size(andarl_gfp,1);
        la_positions(i)=size(s(i).la,1);
    end
    
    alpha = 0.05;
    
    stringa=['For all samples the Anderson Darling statistics, or chance',...
        ' that these data are drawn from the same distribution is: '];
    P_GFP=AnDarksamtest(andarl_gfp,alpha);
    P_RFP=AnDarksamtest(andarl_rfp,alpha);
    stringb=[stringa, sprintf('for GFP %.4f and for RFP %.4f',P_GFP,P_RFP)];
    disp(stringb);
    
    a=1:nassoc; %generate unique pairs of fields to test
    b=1:nassoc;
    [p,q]=meshgrid(a,b);
    mask = triu(ones(nassoc), 1) > 0.5;
    pairs=[p(mask) q(mask)];
    
    
    for prz=1:size(pairs,1)
        current_pair=pairs(prz,:);
        num1=current_pair(1);
        num2=current_pair(2);
        current_gfp=[andarl_gfp(la_start(num1):la_start(num1)+la_positions(num1)-1,1) ...
            ; andarl_gfp(la_start(num2):la_start(num2)+la_positions(num2)-1,1)];
        current_gfp=[current_gfp, ones(size(current_gfp,1),1)];
        current_gfp(la_positions(num1)+1 : la_positions(num2)...
            +la_positions(num1),2)=2;
        
        current_rfp=[andarl_rfp(la_start(num1):la_start(num1)+la_positions(num1)-1,1) ...
            ; andarl_rfp(la_start(num2):la_start(num2)+la_positions(num2)-1,1)];
        current_rfp=[current_rfp, ones(size(current_rfp,1),1)];
        current_rfp(la_positions(num1)+1 : la_positions(num2)...
            +la_positions(num1),2)=2;
        P_GFP=AnDarksamtest(current_gfp,alpha);
        P_RFP=AnDarksamtest(current_rfp,alpha);
        stringc=sprintf('For fields %d and %d the AD stat is %.4f for GFP and %.4f for RFP',num1,num2,P_GFP,P_RFP);
        disp(stringc);
    end
    
    hindex=0;
    figure
    for i=1:nassoc %go between fields
        h.count=0; %gfp
        h.range=[0, max(max(max(andarl_gfp)))+1];
        h.binwidth=1; %change this to change resolution of histogram
        h2.count=0; %rfp
        h2.range=[0, max(max(max(andarl_rfp)))+1];
        h2.binwidth=1;
        for k=1:la_positions(i) %advance within each field
            hindex=hindex+1;
            h=histo(h,andarl_gfp(hindex,1)); %bin in histo for gfp
            h2=histo(h2,andarl_rfp(hindex,1)); %bin in histo2 for rfp
        end
        string1=sprintf('Histogram of GFP Lifetimes in Dataset %d', i);
        string2=sprintf('Histogram of RFP Lifetimes in Dataset %d', i);
        
        subplot(nassoc,1,i);
        bar(h.vals,h.hist);
        title(string1);
        xlabel('Lifetime in seconds');
        ylabel('Frequency');
        
    end
    
    hindex=0;
    figure
    for i=1:nassoc %go between fields
        h.count=0; %gfp
        h.range=[0, max(max(max(andarl_gfp)))+1];
        h.binwidth=1; %change this to change resolution of histogram
        h2.count=0; %rfp
        h2.range=[0, max(max(max(andarl_rfp)))+1];
        h2.binwidth=1;
        for k=1:la_positions(i) %advance within each field
            hindex=hindex+1;
            h=histo(h,andarl_gfp(hindex,1)); %bin in histo for gfp
            h2=histo(h2,andarl_rfp(hindex,1)); %bin in histo2 for rfp
        end
        string1=sprintf('Histogram of GFP Lifetimes in Dataset %d', i);
        string2=sprintf('Histogram of RFP Lifetimes in Dataset %d', i);
        
        subplot(nassoc,1,i);
        bar(h2.vals,h2.hist);
        title(string2);
        xlabel('Lifetime in seconds');
        ylabel('Frequency');
    end
    
end