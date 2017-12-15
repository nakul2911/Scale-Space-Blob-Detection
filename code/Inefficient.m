function Inefficient(I)
%Usage: Just run the Inefficient.m script by inputting an image file I


%%Read the image
%I = imread('/Users/nakul/Documents/CSE 573/hw2/data/sunflowers.jpg');
%%Converting to grayscale
I = rgb2gray(I);

%%Convert to double
I = im2double(I);
[h,w] = size(I);


% IMPORTANT PARAMETERS TO CHANGE 
sigma = 2.5;
Threshold = 0.017;
k = 1.5;
tic
% n is number of levels in scale space
n = 10;
LogFilterSize = 2*ceil(3*sigma) + 1;
% Creating Log filter
Log = sigma^2 *  fspecial('log',LogFilterSize,sigma);

ScaleSpaceTest = zeros(h,w,n);

% Slower version: increse filter size, keep the same image 


for i = 1:n
    sigmai = sigma * k^(i-1);
    
    filter_size = 2*ceil(3*sigmai)+1;
    % Creating Log filter for each level 
    Log =  fspecial('log', filter_size, sigmai);
    Log = (sigmai^2) *Log;
    Filtered = imfilter(I,Log,'same','replicate');
    Filtered = Filtered .^ 2;
    ScaleSpaceTest(:,:,i) = Filtered; 
end
ScaleSpace = ScaleSpaceTest;
toc


% Performing Non Maximum supression between 2D scales
SuppressionSize = 3;
FullSpace = zeros(h,w,n);
for i = 1:n
    FullSpace(:,:,i) = ordfilt2(ScaleSpace(:,:,i),SuppressionSize^2,ones(SuppressionSize));
end


% Performing Non Maximum supression between thresold and scales
for i = 1:n
    FullSpace(:,:,i) = max(FullSpace(:,:,max(i-1,1):min(i+1,n)),[],3);
end
FullSpace = FullSpace .* (FullSpace == ScaleSpace);
 
    
cy = [];   
cx = [];   
BlobRadius = [];
for i=1:n
    [row,col] = find(FullSpace(:,:,i) >= Threshold);
    Blobs = length(row);
    radius =  sigma * k^(i-1) * sqrt(2); 
    radius = repmat(radius, Blobs, 1);
    cy = [cy; row];
    cx = [cx; col];
    BlobRadius = [BlobRadius; radius];
end
fprintf("Successfully Ended");
    
show_all_circles(I,cx,cy,BlobRadius);    
    
    