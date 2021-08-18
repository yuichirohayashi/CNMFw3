function [filters, traces] = demoCNMFw3()

% run CNMFw3 with simulated data
% OUTPUT
% filters: spatial footprints
% traces: temporal traces

% simdata params
load('bg');
trucells = 100;
snr = [0.4];
scalef = 200;

% cell sorting method params
nPCs = 200;
d1=100;
d2=100;
d3=10;
T = 1000;

[Y, fspk] = gensim3dmovRealBG(bg, snr, scalef);
[A,C] = runCNMFw3mat(Y, nPCs);
traces = C;
filters = permute(reshape(full(A),d1,d2,d3,size(A,2)),[4 1 2 3]);





