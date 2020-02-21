function [cell_sig_f0 cell_sig_diff cell_sig_f_f0] = running_f_f0(cell_sig, win)

%calculates F/F0, smoothing of F0 with cubic spline function

%   INPUTS:
%       cell_sig - Cell signal to smooth
%       win      - Length of window (timepoints) over which to compute F0
%
%   OUTPUTS:
%       cell_sig_f0    - Estimate baseline for F0
%       cell_sig_diff  - F-F0
%       cell_sig_f_f0  - (F-F0)/F0
% 2016, Nicholas Frost


cell_sig_f0 = zeros(size(cell_sig));

% finds lowest point in running window, smoothed with cubic spline
        tppc = csaps([1:size(cell_sig,2)],cell_sig,0.25); 
        CC = fnval(tppc,[1:length(cell_sig)]);  

for i =   1:size(cell_sig,1)
    i
    for j = 1:size(cell_sig,2)
        if j < (win/2 + 1)
            m = min(CC(i,1:win));
        elseif j > size(cell_sig,2)-(1+win/2)
            m = min(CC(i, j - win/2 : j));
        else    
            m = min(CC(i,j - win/2: j+win/2));     
        end
        k(j) = m; 
    end
    
    cell_sig_f0(i,:) = k;
    
end

cell_sig_diff = cell_sig - cell_sig_f0;
cell_sig_f_f0 = (cell_sig_diff)./cell_sig_f0;