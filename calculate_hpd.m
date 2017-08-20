

%-------------------------------------------------------------------------
%
%
%  the analysis of the simulated data for IBM,  February 14 2011
%
%
%-------------------------------------------------------------------------


clear all;
close all;
disp('please wait   .....');
disp('');
N=20000;
n=5000;

fid = fopen('C-MCMC threshold.txt', 'r');
bb = fscanf(fid, '%g %g %g  %g %g %g %g%g %g%g %g %g', [12 inf]);    % It has two rows now.
bb = bb';
fclose(fid);
a=bb(5001:end, 2:12);

size(a)

a1=mean(a)
hpdi(a, 95)
std(a)



fid = fopen('MCMC time seris normal.txt', 'r');
bb = fscanf(fid, '%g %g %g  %g %g %g ', [6 inf]);    % It has two rows now.
bb = bb';
fclose(fid);
a=bb(5001:end, 2:6);

size(a)

a1=mean(a)
hpdi(a, 95)
std(a)



