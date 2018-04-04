% Plot spikes from list of time vs cell number
% and determine quality of recall
% 
% VCUT 26-3-14

close all; clear all;

NCELL = 733;  % number of cells (neurons)
NGCELL = 100; % number of GC (output) cells
NCA3PCELL = 100; % number of CA3-PC (output) cells
NCA1PCELL = 100; % number of CA1-PC (output) cells

NPATT = 1;   % number of patterns
SPATT = 4;   % number of active cells per pattern
CPATT = 1;  % index of cue pattern

RTIME = 50+(140*3);    % run time (msecs)
STIME = 0;
ETIME = 500;

FGCPATT = '...\DGpattsN100S4P1.dat';            % GC patterns file
FCA3PCPATT = '...\CA3pattsN100S20P1.dat';     % CA3-PC patterns file
FCA1PCPATT = '...\CA1pattsN100S20P1.dat';     % CA1-PC patterns file

FSTEM = '...\Results_DG_CA3_CA1_w_inhibition_spt';   % spikes file
FSPIKE = [FSTEM '.dat'];   % spikes file

gcpatts = load(FGCPATT);   % load GC stored patterns
gccue = gcpatts(:,CPATT);   % extract GC cue pattern
ca3pcpatts = load(FCA3PCPATT);   % load CA3-PC stored patterns
ca3pccue = ca3pcpatts(:,CPATT);   % extract CA3-PC cue pattern
ca1pcpatts = load(FCA1PCPATT);   % load CA1-PC stored patterns
ca1pccue = ca1pcpatts(:,CPATT);   % extract CA1-PC cue pattern

sp = load(FSPIKE);  % load spike times
st = sp(:,1);       % extract times
cell = sp(:,2);     % extract corresponding cell indices
% extract GC spiking
gcstp = st(cell < NGCELL);
gccellp = cell(cell < NGCELL);
% extract CA3-PC spiking
ca3pcstp = st(cell >= 105 & cell <= 105+NCA3PCELL-1);
ca3pccellp = cell(cell >= 105 & cell <= 105+NCA3PCELL-1);
% extract CA1-PC spiking
ca1pcstp = st(cell >= 209 & cell <= 209+NCA1PCELL-1);
ca1pccellp = cell(cell >= 209 & cell <= 209+NCA1PCELL-1);

% Analyse spiking over time and compare with cue
DT = 1; % sliding time
% TW = 2;     % width of sliding time window
%TW = 5;    % width of sliding time window
gcTW = 5;    % width of sliding time window
ca3pcTW=5;      % width of CA3 sliding time window
ca1pcTW=5;      % width of CA1 sliding time window

gcti = 0:DT:RTIME-gcTW;
gcNW = length(gcti);            % number of time windows
ca3pcti = 0:DT:RTIME-ca3pcTW;
ca3pcNW = length(ca3pcti);      % number of time windows
ca1pcti = 0:DT:RTIME-ca1pcTW;
ca1pcNW = length(ca1pcti);      % number of time windows

gcnc = zeros(gcNW,1);
gcha = zeros(gcNW,1);
gcco = zeros(gcNW,1);
gcan = zeros(gcNW,1);
gcmc = mean(gccue);         % mean GC cue activity
ca3pcnc = zeros(ca3pcNW,1);
ca3pcha = zeros(ca3pcNW,1);
ca3pcco = zeros(ca3pcNW,1);
ca3pcan = zeros(ca3pcNW,1);
ca3pcmc = mean(ca3pccue);   % mean CA3-PC cue activity
ca1pcnc = zeros(ca1pcNW,1);
ca1pcha = zeros(ca1pcNW,1);
ca1pcco = zeros(ca1pcNW,1);
ca1pcan = zeros(ca1pcNW,1);
ca1pcmc = mean(ca1pccue);   % mean CA1-PC cue activity

% DG-GC
for i=1:gcNW
    gcrp = gccellp(gcstp>=gcti(i) & gcstp<gcti(i)+gcTW); % active cells in sliding window
    gcnc(i) = length(gcrp);    % number of active cells in window
    gcp = zeros(NGCELL,1);
    gcp(gcrp+1,1) = 1;  % recalled pattern
    gcha(i) = (sum(gcp == gccue)/NGCELL);         % hamming distance
    ca3pcha(i) = (sum(gcp == ca3pccue)/NCA3PCELL);   % hamming distance
    ca1pcha(i) = (sum(gcp == ca1pccue)/NCA1PCELL);   % hamming distance

    gcmp = mean(gcp);   % mean pattern activity
    % correlation (normalised dot product)
    if gcmp == 0
        gcco(i) = 0;
    else
        gcco(i) = dot(gcp,gccue)/sqrt(sum(gcp)*sum(gccue));
    end;    
    % angle (Graham & Willshaw 97)
    if gcmp == 0 | gcmp == 1
        gcan(i) = 0;
    else
        gcan(i) = sum((gcp-gcmp).*(gccue-gcmc))/sqrt(sum((gcp-gcmp).^2)*sum((gccue-gcmc).^2));
    end;

end;

% CA3-PC
for i=1:ca3pcNW  
    ca3pcrp = ca3pccellp(ca3pcstp>=ca3pcti(i) & ca3pcstp<ca3pcti(i)+ca3pcTW); % active cells in sliding window
    ca3pcnc(i) = length(ca3pcrp);    % number of active cells in window
    ca3pcp = zeros(NCA3PCELL,1);
    ca3pcp(ca3pcrp+1-105,1) = 1;  % recalled pattern
    ca3pcha(i) = (sum(ca3pcp == ca3pccue)/NCA3PCELL);  % hamming distance
    ca3pcmp = mean(ca3pcp);   % mean pattern activity
    % correlation (normalised dot product)
    if ca3pcmp == 0
        ca3pcco(i) = 0;
    else
        ca3pcco(i) = dot(ca3pcp,ca3pccue)/sqrt(sum(ca3pcp)*sum(ca3pccue));
    end;    
    % angle (Graham & Willshaw 97)
    if ca3pcmp == 0 | ca3pcmp == 1
        ca3pcan(i) = 0;
    else
        ca3pcan(i) = sum((ca3pcp-ca3pcmp).*(ca3pccue-ca3pcmc))/sqrt(sum((ca3pcp-ca3pcmp).^2)*sum((ca3pccue-ca3pcmc).^2));
    end;
end;

% CA1-PC
for i=1:ca1pcNW  
    ca1pcrp = ca1pccellp(ca1pcstp>=ca1pcti(i) & ca1pcstp<ca1pcti(i)+ca1pcTW); % active cells in sliding window
    ca1pcnc(i) = length(ca1pcrp);    % number of active cells in window
    ca1pcp = zeros(NCA1PCELL,1);
    ca1pcp(ca1pcrp+1-209,1) = 1;  % recalled pattern
    ca1pcha(i) = (sum(ca1pcp == ca1pccue)/NCA1PCELL);  % hamming distance
    ca1pcmp = mean(ca1pcp);   % mean pattern activity
    % correlation (normalised dot product)
    if ca1pcmp == 0
        ca1pcco(i) = 0;
    else
        ca1pcco(i) = dot(ca1pcp,ca1pccue)/sqrt(sum(ca1pcp)*sum(ca1pccue));
    end;    
    % angle (Graham & Willshaw 97)
    if ca1pcmp == 0 | ca1pcmp == 1
        ca1pcan(i) = 0;
    else
        ca1pcan(i) = sum((ca1pcp-ca1pcmp).*(ca1pccue-ca1pcmc))/sqrt(sum((ca1pcp-ca1pcmp).^2)*sum((ca1pccue-ca1pcmc).^2));
    end;
end;


% Generate figure
figure;
clf;
ms=8;
lw=2;

subplot(7,1,1);
plot(sp(:,1), sp(:,2), 'k.', 'markersize', ms);   % raster plot of Sep, EC & CA3 spiking
title('(a) Input spikes', 'fontsize', 14);
ylabel('Cell no.', 'fontsize', 14);
axis([STIME ETIME-40 NGCELL+213 NGCELL+213+420]);
subplot(7,1,2);
plot(sp(:,1), sp(:,2), 'k.', 'markersize', ms);   % raster plot of PC spiking
title('(b) Granule cell spikes', 'fontsize', 14);
ylabel('Cell no.', 'fontsize', 14);
axis([STIME ETIME-40 0 NGCELL-1]);
subplot(7,1,3);
plot(gcti, gcco, 'k-', 'LineWidth', lw); % recall quality
title('(c) GC recall quality', 'fontsize', 14);
ylabel('Quality', 'fontsize', 14);
axis([STIME ETIME-40 0 1.02]);
subplot(7,1,4);
plot(sp(:,1), sp(:,2), 'k.', 'markersize', ms);   % raster plot of PC spiking
title('(d) CA3 pyramidal cell spikes', 'fontsize', 14);
ylabel('Cell no.', 'fontsize', 14);
axis([STIME ETIME-40 105 105+NCA3PCELL-1]);
subplot(7,1,5);
plot(ca3pcti, ca3pcco, 'k-', 'LineWidth', lw); % recall quality
title('(e) CA3-PC recall quality', 'fontsize', 14);
ylabel('Quality', 'fontsize', 14);
axis([STIME ETIME-40 0 1.02]);
subplot(7,1,6);
plot(sp(:,1), sp(:,2), 'k.', 'markersize', ms);   % raster plot of PC spiking
title('(f) CA1 pyramidal cell spikes', 'fontsize', 14);
ylabel('Cell no.', 'fontsize', 14);
axis([STIME ETIME-40 209 209+NCA1PCELL-1]);
subplot(7,1,7);
plot(ca1pcti, ca1pcco, 'k-', 'LineWidth', lw); % recall quality
title('(g) CA1-PC recall quality', 'fontsize', 14);
ylabel('Quality', 'fontsize', 14);
xlabel('Time (msecs)', 'fontsize', 14);
axis([STIME ETIME-40 0 1.02]);

mean(gcco(gcco>0))
mean(ca3pcco(ca3pcco>0))
mean(ca1pcco(ca1pcco>0))

%print('-dpng', ['Images/' FSTEM]);


    