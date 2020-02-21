function res = imshowScaled(image,low,high)

% imshowScaled: displays image with low and high (fractional) values specified
% PURPOSE: read, scale, and display an image
% (well, add a choice to scale or not, later)
% low and high are proportions of the range, not intensity values

% A=imread(image);

if high<low
    temp=low;
    low=high;
    high=low;
end

minval=(max(image(:))-min(image(:)))*low+min(image(:));
maxval=(max(image(:))-min(image(:)))*high+min(image(:));

imshow(image,[minval maxval]);

