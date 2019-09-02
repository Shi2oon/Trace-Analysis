function [TRACE] = SlipTrace(named)
% a Matlab code for trace analysis. 
% Select the trace on a picture
% The code decide the direction, trace incline and ‘host’ grain
% Get the ‘host’ grain mean orientation and crystal structure
% Get all possible twin and slip variation for the specific crystal 
% structure (i.e. FCC, BCC, HCP (for all systems, or each system 
% individually), tetragonal, trigonal)
% transform line formed by the slip trace into ‘host’ crystal coordinate system
% Calculate all possible Schmid factor with in ‘host’ crystal coordinate system
% Fit the race line with nearest Schmid factor vector
% Plot trace plane, direction and Schmid factor value

% Created by Abdalrhaman Mohamed Koko, as part of a DPhli/PhD in the
% University of Oxford, email abd.mohamed@stx.ox.a.cuk or abdo.aog@gmail.com
% last updated 02/09/2019

%% prerequisites
% Install MTEX
% Creat a matlab code from MTEX to import your data (you can use MTEX
% import_wizard('EBSD')
% if your crystal structure is cubic you need to add if its FCC or BCC
% add this line under CS in your MTEX function 
% CS{2}.opt.type='fcc';

% 'named' = the function full path 

% example SlipTrace('P:\Abdo\GitHub\My-Mtex-Code v2\Data\AM190406DSSB_100hr')
% named='P:\Abdo\GitHub\My-Mtex-Code v2\Data\AM190406DSSB_100hr';

% Mtex path
% addpath('C:\Users\scro3511\OneDrive - Nexus365\Documents\GitHub\mtex-5.2.beta2');  

%% the Code 
startup_mtex;   % start mtex
[file,coded,~] = fileparts(named);
addpath (file); % code path

eval(sprintf('%s;',coded));

close all; clc; warning off
%Read in ebsd file and calculate grains
if sum(CS{1}=='notIndexed')==10 
    [ebsd,CS,grains] = IndexAndFill(ebsd,CS); % fill
end

plot(ebsd,ebsd.prop.bc)
mtexColorMap white2black;       
% and superpose maybe a phase i.e. the notIndexed phase
hold on  % this keeps the present figure axis
plot(ebsd('notIndexed'),'FaceColor',[1 0.2 0],'FaceAlpha',0.6)
hold off 

%imgfile           = imgfile(1:end-4);
[filepath,~,~] = fileparts(fname);
fname = ([filepath '\' date '_Traces']);      mkdir(fname) ;

COUNTERS = 0;          answer = 'Y'; 
while answer == 'Y'
    %Display plot of grains and have user select the one that their SEM image is located in.
    COUNTERS     = COUNTERS+1;
    
    %% select twin or slip trace
    opts.Interpreter = 'tex'; % Include the desired Default answer
    opts.Default     = 'Slip';     % Use the TeX interpreter to format the question
    quest            = 'Are you tracing a Twin or a Slip?';
    Ans{COUNTERS}    = questdlg(quest,'Boundary Condition','Twin','Slip', opts);
    
    %% select it
    uiwait(msgbox('Click two points on a single slip trace.','Select Corners','modal'));
    c1{COUNTERS} = ginput(1);
    c2{COUNTERS} = ginput(1);   
    hold on; 
    h = plot([c1{COUNTERS}(1),c2{COUNTERS}(1)],[c1{COUNTERS}(2),c2{COUNTERS}(2)],...
        'Color','r','LineStyle','--','LineWidth',2); 
    
    %% decide to continue or not
    opts.Interpreter = 'tex'; % Include the desired Default answer
    opts.Default     = 'Y';     % Use the TeX interpreter to format the question
    quest            = '(Y) Chose Another Single Slip trace, (C) Remove Previous Selection, (N) Done with Trace Selection';
    answer           = questdlg(quest,'Boundary Condition','Y','N','C', opts);
    if answer=='C'
        COUNTERS=COUNTERS-1;            answer='Y';
        delete(h)
    end
end
close all;

for iv=1:COUNTERS
    line      = normalize(vector3d(c2{iv}(1)-c1{iv}(1),c2{iv}(2)-c1{iv}(2),0));
    grain     = grains(grains.findByLocation(c1{iv}));
    ori       = grain.meanOrientation;
    [CStrace] = findCS(CS,ebsd(grain)); % find corrspond CS
    if Ans{iv} == 'Twin'
        [~,~,~,sS]=decideDS(CStrace,0.3);    
    else
        [~,sS]    = decideDS(CStrace,0.3);     % find slip system    
    end
    sSlocal   = grain.meanOrientation * sS;

    %transform line formed by the slip trace into crystal coordinate system
    line    = Miller(ori*line,grain.CS);
    sigma   = stressTensor.uniaxial(vector3d.Z);
    %tau    = sS.SchmidFactor(line);
    SF      = sSlocal.SchmidFactor(sigma); %plot(SF)
    sS.CRSS = abs(SF);
    %[SFMax,~] = max(abs(SF),[],2);% take the maxium allong the rows
    %plot(grain,SFMax,'micronba r','off','linewidth',2);  hold on;    legend off
    for i=1:length(sSlocal)
        quiver(grain,sS(i).trace,'autoScaleFactor',abs(SF(i))) % 'color','k',
        %,'displayName','Burgers vector')
        %quiver(grain,sSlocal(i).b,'color','r','autoScaleFactor',1)
        %,'displayName','slip plane trace')
        index(i) = abs(sS(i).trace-line);
        hold on;
    end
    quiver(grain,line,'autoScaleFactor',max(abs(SF)),'color','k'); hold off
    axis off ;
    saveas(gcf,[fname,'\traces_' num2str(iv) '.png']);
    saveas(gcf,[fname,'\traces_' num2str(iv) '.fig']); close all

    count=0; 
    for i=1:length(sSlocal)
        if index(i) <= min(index)+(max(index)-min(index))/max(index)*0.1 %5% tolerance
            count            = count+1;
            indexed(count,1) = i;
            indexed(count,2) = abs(SF(i));
        end
    end
    [row,~]     = find(indexed==max(indexed(:,2)));
    TRACE{iv,1} = sS(indexed(row,1));
    TRACE{iv,2} = max(indexed(:,2));
    clear -vars index indexed
end

plot(ebsd,ebsd.prop.bc)
mtexColorMap white2black;   
% and superpose maybe a phase i.e. the notIndexed phase
hold on  % this keeps the present figure axis
plot(ebsd('notIndexed'),'FaceColor',[1 0.2 0],'FaceAlpha',0.6)

% plot each trace
for iv=1:COUNTERS
    hold on; 
    plot([c1{iv}(1),c2{iv}(1)],[c1{iv}(2),c2{iv}(2)],...
        'Color','r','LineStyle','-','LineWidth',2); 
    ht = text(min(c2{iv}(1),c1{iv}(1)),(c2{iv}(2)+c1{iv}(2))/2, ['[' ...
         num2str(round(TRACE{iv,1}.b.u)) num2str(round(TRACE{iv,1}.b.v)) ...
         num2str(round(TRACE{iv,1}.b.w)) '] (' ...
         num2str(round(TRACE{iv,1}.n.h)) num2str(round(TRACE{iv,1}.n.k)) ...
         num2str(round(TRACE{iv,1}.n.l)) '), ' num2str(TRACE{iv,2}) ...
         ],'Color','b');
%     iclinity = acos((c2{iv}(2) - c1{iv}(2))/sqrt((c2{iv}(2) - ...
%                c1{iv}(2))^2+(c2{iv}(1) - c1{iv}(1))^2));
%     slope = (c2{iv}(2) - c1{iv}(2))/ (c2{iv}(1) - c1{iv}(1));
%     if slope < 1;   set(ht,'Rotation',iclinity*180/pi)
%     else;           set(ht,'Rotation',iclinity*180/pi+90);     end
    set(ht,'FontSize',16)
end
hold off;  

% tight plot
% ax          = gca;
% outerpos    = ax.OuterPosition;
% ti          = ax.TightInset; 
% left        = outerpos(1) + ti(1);
% bottom      = outerpos(2) + ti(2);
% ax_width    = outerpos(3) - ti(1) - ti(3);
% ax_height   = outerpos(4) - ti(2) - ti(4);
% ax.Position = [left bottom ax_width ax_height];
saveas(gcf,[fname, '\traced.png']); % save
saveas(gcf,[fname, '\traced.fig']); close all % save