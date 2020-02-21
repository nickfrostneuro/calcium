function [include_snr] = snr_thresh(cell_sig_subtracted, thresh_snr)

%threshold signal to noise ratio

% find mean and variance of each 'cell' and throw out cells that do not
% cross 2.5x std+mean in at least 25% of 2000 frame bins

%inputs
%cell_sig_subtracted - background subtracted fluorescence signal
%cell_sig_surround - background fluorescence signal
%cell_sig - unsubtracted fluorescence signal 
%segcentroid - coordinates of segment
%thresh_snr - target threshold 

%outputs
%cell_sig_subtracted_snr - array background subtracted fluorescence signal
%cell_sig_surround_snr - background fluorescence signal
%cell_sig_snr - unsubtracted fluorescence signal 
%segcentroid_snr - coordinates of segment
%rejects - indices of rejected 'cells'



C= cell_sig_subtracted;
N = size(C,2);
win = 50;
baseline_var = 0.5;

% find the absolute change at each point
for i = 1:size(C,1)
    if mod (i,10) == 0
        i
    end
D(i,:) = C(i,2:N) - C(i,1:N-1);     %NOTE D starts at 2

    for j=1:N-win,
        x = C(i,j:j+win);
        y(j) = var(x);
        z(j) = median(abs(x - mean(x)));
    end
 
zz = tiedrank(z);
[~, quietepoch] = find (zz< (baseline_var*size(C,2))); %finds subset of epochs with lowest variance
medianchange(i) = median(z(quietepoch));   
    
end










bs = 1000; tt = []
for j = 1:floor (size(cell_sig_subtracted,2) / bs)
    for i = 1:size(cell_sig_subtracted,1)
 
        temp_diff(i,j) = max(C(i,((j-1)*bs)+1:j*bs)) - mean(C(i,((j-1)*bs)+1:j*bs));
        tt(i,j) = temp_diff(i,j) - medianchange(i)*thresh_snr;
    end
end

jj = tt>0;

include_snr = (sum(jj'>0) >15);

