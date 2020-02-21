function [ ica_segments_edge, segcentroid_edge, exclude_edge ] = edgefilter( ica_segments, segcentroid, edge_distance, movm )
%edge filter - removes cells too close to edge as if drift correcting,
%those ICs that overlap with areas that become zeros por cell may drift off
%screen...

% Inputs:
%       ica_segments - Output from ica segmentation
%       segcentroid - Corresponding centroids for each segment
%       edge_distance - Minimum distance (in pixels) from edge

%  Outputs:
%       ica_segments_edge - Filtered ica_segments
%       segcentroid_edge - Corresponding centroids
%       exclude_edge - Indices of excluded segments

sz = size(movm);
keep = [0+edge_distance, sz(1)- edge_distance;  0+edge_distance, sz(2)- edge_distance] ;
exclude = [];

for i = 1:size(ica_segments,1)
    if segcentroid(i,1) < keep(2,1) || segcentroid(i,1) > keep(2,2) || segcentroid(i,2) < keep(1,1) || segcentroid(i,2) > keep(1,2)
        exclude = [exclude;i];
    end
end
    
exclude_edge = exclude;
segcentroid_edge = segcentroid;
ica_segments_edge = ica_segments;
segcentroid_edge(exclude,:) = [];
ica_segments_edge(exclude,:,:) = [];

