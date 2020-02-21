function [spiketimes spt_beg spt_end] = detectspikes11(C, thresh, thresh2, thresh3, thresh_abs, spikewin, ...
                                      mindur, baseline_var)

    % C is the calcium trace to be analyzed
    %  [spiketimes spt_beg spt_end] = detectspikes5(C, thresh, thresh2, spikewin, mindur)
    % thresh is the threshold, in factors of the median
    % point-by-point change in C, for detecting candidate events
    % (usually 3)
    
    % thresh2 is threshold amplitude

    % thresh3 is the threshold, for integral under spike (usually 100-200)
    % thresh_abs is a minimum absolute value of F-F0/F0
    % spikewin is the number of points between the detection and
    % the end of an event
    
    % mindur is the minimum duration in points for an event
    
    % RETURNS:
    %               all outputs are as 1xNspikes
    %           spiketimes = the peak of each spike time
    %           spt_beg = Corresponding start time of that event
    %           spt_end = Corresponding end time
    
N = length(C)';

Dmask = ones(N-1,1);

spiketimes = [];
spt_beg = [];
spt_end = [];

% find the absolute change at each point
D = C(2:N) - C(1:N-1);     %NOTE D starts at 2nd frame

medianchange = median(abs(D));

for i=1:N-10+1,
    x = C(i:i+10-1);
    y(i) = var(x);
    z(i) = median(abs(x - mean(x)));
end
zz = tiedrank(z);
[~, quietepoch] = find (zz(1,:) < (baseline_var*length(C))); %finds subset of epochs with lowest variance
medianchange = median(z(quietepoch));



Dmask(1:spikewin) = zeros;
Dmask(end-spikewin:end) = zeros;
D = D .* Dmask;
[maxD, locmaxD] = max(D);
locmaxD = locmaxD + 1;        

while maxD && maxD > thresh*medianchange

  
    % find the beginning of the spike event
    k = 0;
    while k < spikewin-2 && locmaxD-k-1 > 5 && C(locmaxD-k) > min(C(locmaxD-k-2:locmaxD-k-1)), %look back to find location in C that has lower value than prior two bins
        k = k+1;
    end

    locbeg = locmaxD  - k;

    %locbeg is the lowest point between locmaxD and k
           
    
    % find the peak of the spike event
    %look for location where C > either of preceeding bins, 
    %then looks for max between locmaxD and locmaxD + k + 3
    k = 1;
    while locmaxD+k < (N-4) && k < spikewin && C(locmaxD+k) >= min(C(locmaxD+k-2:locmaxD+k-1)) % 
        k = k+1;
    end
    
        
    [~, rel_peak] = max(C(locmaxD:locmaxD+k)); % location of the peak
    locpeak = (rel_peak-1) + locmaxD;
    
    if locbeg+spikewin > N,            
        Dmask(locmaxD-1) = 0;
        D = C(2:N) - C(1:N-1);
        D = D .* Dmask;
        [maxD, locmaxD] = max(D);
        locmaxD = locmaxD + 1;
        continue;
    end
    
    % find the end of the spike event
    k = 1;

   while (locpeak+k < N) && (locpeak - locbeg + k) <= spikewin && C(locpeak+k)-C(locbeg) >= 0.7*(C(locpeak)-C(locbeg)) %find first point below threshold
      k = k + 1;
   end
    spikeend = locpeak + k;

    locend = spikeend;
    
        bl = min(locbeg-1, 1000);
        if locend - locbeg >= mindur...
            && C(locpeak) > mean(C(locbeg-2:locbeg)) + (thresh2*medianchange)...
            && sum(C(locbeg:locend))-(C(locbeg)*(locend-locbeg)) >= thresh3*medianchange...
            && C(locpeak) >= mean(C(locbeg-bl:locbeg-1)) + thresh_abs,
            spiketimes = [spiketimes, locpeak];
            spt_beg = [spt_beg, locbeg];
            spt_end = [spt_end, locend];
            Dmask(locbeg-1:locend-1) = 0; %locbeg + locend are indexed to C, hence n-1.
        end
    
        Dmask(locmaxD-1) = 0;
        D = C(2:N) - C(1:N-1);
        D = D .* Dmask;
        [maxD, locmaxD] = max(D);
        locmaxD = locmaxD + 1;
    
end

length(spiketimes);
