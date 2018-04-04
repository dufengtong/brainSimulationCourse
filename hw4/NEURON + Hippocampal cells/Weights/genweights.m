% Generate weight matrix through storage of random binary patterns
% with the clipped Hebbian learning rule.
% Weights stored in incoming column order for each target cell.
% BPG 19-8-08

NCELL = 100;  % number of cells (neurons)
NPATT = 1   % number of patterns
SPATT = 20;   % number of active cells per pattern
%PC = 0.5;    % percent connectivity (not normally used)

% FWGT = 'C:\Users\Vassilis\Desktop\Documents\My Research\Current\DG-CA3-CA1 microcircuit model\Code\NEURON\My model\DG-CA3-CA1-network_w_MS_recall_perf\Weights\DGwgtsN100S4P1.dat';   % weights file
% FPATT = 'C:\Users\Vassilis\Desktop\Documents\My Research\Current\DG-CA3-CA1 microcircuit model\Code\NEURON\My model\DG-CA3-CA1-network_w_MS_recall_perf\Weights\DGpattsN100S4P1.dat';   % patterns file
FWGT = 'C:\Users\Vassilis\Desktop\Documents\My Research\Current\DG-CA3-CA1 microcircuit model\Code\NEURON\My model\DG-CA3-CA1-network_w_MS_recall_perf\Weights\CA3wgtsN100S20P1.dat';   % weights file
FPATT = 'C:\Users\Vassilis\Desktop\Documents\My Research\Current\DG-CA3-CA1 microcircuit model\Code\NEURON\My model\DG-CA3-CA1-network_w_MS_recall_perf\Weights\CA3pattsN100S20P1.dat';   % patterns file
% FWGT = 'C:\Users\Vassilis\Desktop\Documents\My Research\Current\DG-CA3-CA1 microcircuit model\Code\NEURON\My model\DG-CA3-CA1-network_w_MS_recall_perf\Weights\CA1wgtsN100S20P1.dat';   % weights file
% FPATT = 'C:\Users\Vassilis\Desktop\Documents\My Research\Current\DG-CA3-CA1 microcircuit model\Code\NEURON\My model\DG-CA3-CA1-network_w_MS_recall_perf\Weights\CA1pattsN100S20P1.dat';   % patterns file


%rand('state',0);
rand('state',sum(100*clock));

%rw = rand(NCELL);
%w = ones(NCELL).*PC >= rw

w = zeros(NCELL);
p = zeros(NCELL, NPATT);

for i=1:NPATT
  % generate pattern
  pr = randperm(NCELL);
  pi = pr(1:SPATT);        % indices of active cells in pattern
  p(pi,i) = 1;
  % store in weight matrix
  w = w+(p(:,i)*p(:,i)');
end

w = w > 0;  % clip weight matrix

dlmwrite(FWGT, w, ' ');
dlmwrite(FPATT, p, ' ');

