function  [dups] = clean_raster(segcentroid, raster,C, distance_threshold, corr_threshold)

%removes overlapping cells with high calculated correlation coefficients based on detected events
%inputs
%raster = event raster, n cells * x frames
%segcentroid = x/y position of identified cell (in pixels)
%distance_threshold - number of pixels below which cells are considered to
%be possibly overlapping
%corr_threshold - correlation coefficient threshold above which one cell will be discarded per pair (in cells spaced less than distance_threshold apart)

%outputs
%dups - cell id to be discarded


if nargin < 4
distance_threshold = 12;
corr_threshold = 0.25;
end
%threshold above which segments assumed to be from equivalent cells
%regardless of distance
corr_threshold2 = 0.50;

pair = []; dups = []; kl = []; cc = []; d = [];
num_cells = size(segcentroid,1);
cc = corr(full(raster)');

for i = 1:num_cells
    for j = 1:num_cells
        if i > j
        d(i,j) = ((segcentroid(i,1)-segcentroid(j,1))^2 + (segcentroid(i,2)-segcentroid(j,2))^2)^0.5;
        else
            d(i,j) = nan;
            cc(i,j) = nan;
        end  
    end      
end

close_pairs = d < distance_threshold;
corr_pairs = cc > corr_threshold;
high_corrs = cc >= corr_threshold2;
problem_cells = (close_pairs & corr_pairs) | high_corrs>0;
%indexes of cell pairs within distance threshold and above correlation cutoff
[pair(:,1) pair(:,2)] = ind2sub([num_cells num_cells],find(problem_cells(:)>0));



if isempty (pair) == 0
%identify cells overlapping with multiple other cells, if any. Otherwise
%choose based on s/n
tmp = [];
for i = 1:size(pair(:))
    tmp(i) = sum(pair(:)==pair(i));
end
tmp = tmp>1;
kl = reshape(tmp, size(pair,1),2);

%now choose one cell from each pair
tmp = [];
for i = 1:size(pair,1)
    clear tmp
    for j = 1:2
        cl = pair(i,j);
        tmp(j) = mean(C(cl,raster(cl,:)>0)');
    end
        if sum(kl(i,:)==0)
            if tmp(1) > tmp(2)
                kl(i,1) = 1;
            else kl(i,2) = 1;
            end
        end
end

dups = pair(kl);

%now remove duplicate entries
i = 1;
while i < size(dups,1)
    xx = find(dups == dups(i));
    if size(xx,1) > 1
        dups(xx(2:end)) = [];
    end
    i = i+1;
end

end
