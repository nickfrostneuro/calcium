function [cell_sig,cell_sig_subtracted,cell_sig_surround] = CellsortApplyFilter_singleframes_subtractsurroundmask(fn, ica_segments, flims, movm, subtractsurround,backgroundmask)
% cell_sig = CellsortApplyFilter(fn, ica_segments, flims, movm, subtractmean)
%
%CellsortApplyFilter
% Read in movie data and output signals corresponding to specified spatial
% filters, Note: This works on inputs where each frame is an individual
% tiff file. 
%
% Inputs:
%     fn - file name of TIFF movie file
%     ica_segments - nIC x X matrix of ICA spatial filters
%     flims - optional two-element vector of frame limits to be read
%     movm - mean fluorescence image
%     subtractsurround - boolean specifying whether or not to subtract the
%     median of surrounding region fluorescence of each time frame
%     backgroundmask - mask consisting of all identified segments - 0 value
%     signifies background
%
% Outputs:
%     cell_sig - cellular signals
%     cell_sig_subtracted - subtract median of surrounding region
%     cell_sig_surround - signal from the surrounding region
%
%   modified from mukamel et schnitzer, 2009 by Nicholas Frost. 2016.
 
dir_use = dir('*.tif');
nt_full = length(dir_use);
nt=nt_full;
 
 
if (nargin<3)||isempty(flims)
    flims = [1,nt];
else
    nt = diff(flims)+1;
end
if nargin<5
    subtractsurround = 0;
end

[pixw,pixh] = size(imread(fn,1));
if (nargin<4)||isempty(movm)
    movm = ones(pixw,pixh);
else
    movm = double(movm);
end
k=0;

%this part defines region surrounding a segment for surround subtraction
if subtractsurround == 1
    % first identify cells bounding the IC
    ica_seg_surround = zeros(size(ica_segments,1), pixw, pixh);
    ICA_Surround_mean = zeros(size(ica_segments,1),nt);
    surround_squeeze = zeros(pixw, pixh,size(ica_segments,1));
    %Define the border surrounding each ROI, and the size of each ROI 
    %will bandpass thresholded ICA segment (nonzero) and then detect edge
    %of bandpassed image, to enlarge the region 
    for i = 1:size(ica_segments,1)
        a = ica_segments(i,:,:) > 0;
        b=bpass(squeeze(a),1,15);
        ica_seg_surround(i,:,:) = edge(squeeze(b>0));
        surround_squeeze(:,:,i) = squeeze(ica_seg_surround(i,:,:)) .* (backgroundmask == 0);
        
        %now check to make sure pixels exist in border, if not will increase radius of edge one time
        if nnz(surround_squeeze(:,:,i)) <= 1
            b=bpass(squeeze(b),1,15);
            ica_seg_surround(i,:,:) = edge(squeeze(b>0));
            surround_squeeze(:,:,i) = squeeze(ica_seg_surround(i,:,:)) .* (backgroundmask == 0);
        end
        
        %if there are still no overlapping pixels will give up and default
       %to the mean of edge pixels
        if nnz(surround_squeeze(:,:,i)) <= 1
            surround_squeeze(:,:,i) = squeeze(ica_seg_surround(i,:,:));
        end
        
        ICA_seg_size(i) = nnz(ica_segments(i,:,:));
    end
end

%preallocate large ouput files.
cell_sig = zeros(size(ica_segments,1), nt);

if subtractsurround == 1
    cell_sig_subtracted = zeros(size(ica_segments,1),nt);
    cell_sig_surround = zeros(size(ica_segments,1),nt);

    % make a copy to apply spatial filter to calculate edge median in each frame
    copy_ica_segments = ica_segments; 
    ica_seg_surround = reshape(ica_seg_surround, [], pixw*pixh);
end

ica_segments = reshape(ica_segments, [], pixw*pixh);

tic
fprintf('Loading %5g frames for %d ROIs.\n', nt, size(ica_segments,1))



while k<nt
    
    ntcurr = min(500, nt-k);
    mov = zeros(pixw, pixh, ntcurr);

    %load movie in batches of 500 frames
    for j=1:ntcurr
        temp_fname = dir_use(j+k+flims(1)-1).name;
        movcurr = imread(temp_fname);
        mov(:,:,j) = movcurr;

    %Now calculate median pixel value in edge
        if subtractsurround == 1
        for i = 1:size(ica_segments,1)
            ICA_Surround_mean(i,j+k+flims(1)-1) = mean(nonzeros(surround_squeeze(:,:,i) .* mov(:,:,j)  ));
        end 
        end
    end


    %now calculate both the cell signal and the signal in the surrounding area
    mov = reshape(mov, pixw*pixh, ntcurr);
    cell_sig(:, k+[1:ntcurr]) = (ica_segments>0)*mov;
    if subtractsurround == 1
        cell_sig_surround(:, k+[1:ntcurr]) = (ica_seg_surround>0)*mov;
    end    

    
    k=k+ntcurr;
    fprintf('Loaded %3.0f frames; ', k)
    toc
end


    if subtractsurround == 1
    %Subtracted signal is the median of the surround * number of pixels in
    %the segment ROI. Will smooth the background measurement slightly
        %tppc = csaps([1:size(ICA_Surround_mean,2)],ICA_Surround_mean,0.25); 
        %CC = fnval(tppc,[1:length(ICA_Surround_mean)]);  
        CC = ICA_Surround_mean;
        for i = 1:size(cell_sig,1)
            cell_sig_subtracted(i,:) = cell_sig(i,:) - ( CC(i,:) * ICA_seg_size(i));
           % if min(cell_sig_subtracted(i,:)') < 1
           %     cell_sig_subtracted (i,:) = cell_sig_subtracted (i,:) + (1 - min(cell_sig_subtracted(i,:)')); 
           % end
        end
    end

if subtractsurround == 0
    cell_sig_subtracted = [];
    cell_sig_surround = [];
end
