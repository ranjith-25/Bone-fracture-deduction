img =imread('C:\Users\HP\Documents\MATLAB\F_image\img1.jpg');

figure;
subplot(2, 2, 1);
imshow(img);
title('Original Image');

ImgBlurSigma =2; % Amount to denoise input image

MinHoughPeakDistance= 5; % Distance between peaks in Hough transform angle detection

HoughConvolutionLength = 40; % Length of line to use to detect bone regions

HoughConvolutionDilate = 2;% Amount to dilate kernel for bone detection

BreakLineTolerance= 0.25; % Tolerance for bone end detection

breakPointDilate =6; % Amount to dilate detected bone end points

img=(rgb2gray(img)); % Load image.

img=imfilter(img, fspecial('gaussian', 10, ImgBlurSigma), 'symmetric'); % Denoise

% Do edge detection to find bone edges in image

% Filter out all but the two longest lines

% This feature may need to be changed if break is not in middle of bone

boneEdges=edge(img, 'canny');

subplot(2,2,2);
imshow(boneEdges);
title('Edge Detected Image');

boneEdges=bwmorph(boneEdges, 'close');

edgeRegs=regionprops(boneEdges, 'Area', 'PixelIdxList');

AreaList=sort(vertcat(edgeRegs.Area), 'descend');

edgeRegs(~ismember(vertcat(edgeRegs.Area),AreaList(1:2))) = [];

edgeImg = zeros(size(img, 1), size(img,2));

edgeImg(vertcat(edgeRegs.PixelIdxList)) = 1;

% Do hough transform on edge image to find angles at which bone pieces are % found

% Use max value of Hough transform vs angle to find angles at which lines

% are oriented. If there is more than one major angle contribution there

% will be two peaks detected but only one peak if there is only one major

% angle contribution (ie peaks here number of located bones - Number of % breaks + 1)

[H,T,R]=hough(edgeImg, 'RhoResolution', 1, 'Theta', -90:2:89.5);

maxHough=max(H, [], 1);

HoughThresh=(max(maxHough)-min(maxHough))/2+ min(maxHough);

[~, HoughPeaks]=findpeaks(maxHough,'MINPEAKHEIGHT', HoughThresh, 'MinPeakDistance', MinHoughPeakDistance);

% Plot Hough detection results

subplot(2,2,3);

plot(T, maxHough);

hold on

plot([min(T) max(T)], [HoughThresh, HoughThresh], 'r');

plot(T(HoughPeaks), maxHough(HoughPeaks), 'rx', 'MarkerSize', 12, 'LineWidth', 2);

hold off

xlabel('Theta Value'); ylabel('Max Hough Transform');

legend({'Max Hough Transform', 'Hough Peak Threshold', 'Detected Peak'});
title('Plot Hough detection');

fprintf('HoughPeaks: %f, Theta: %f\n',HoughPeaks,T(HoughPeaks));
fl=0;
if (HoughPeaks > 30.0)
    disp('Result: Fractured');
    fl=1;
else
    disp('Result: Not Fractured');
end


% Locate site of break

if numel(HoughPeaks)>1

    BreakStack = zeros(size(img, 1), size(img, 2), numel(HoughPeaks));
% Convolute edge image with line of detected angle from hough transform for m=1:numel (Hough Peaks);
for m=1:numel(HoughPeaks)
    bonekernel = strel('line', HoughConvolutionLength, T(HoughPeaks(m)));

kern =double(bwmorph (bonekernel.getnhood(), 'dilate', HoughConvolutionDilate));
BreakStack(:,:,m) = imfilter(edgeImg, kern).*edgeImg;

end

% Take difference between convolution images. Where this crosses zero

% [within tolerance) should be where the break is. Have to filter out

% regions elsewhere where the bone simply ends.

brimg=abs(diff(BreakStack, 1, 3)) <BreakLineTolerance*max(BreakStack(:)) & edgeImg > 0;

[BpY, BpX] = find(abs(diff(BreakStack, 1, 3)) <BreakLineTolerance*max(BreakStack(:)) & edgeImg > 0);

brimg=bwmorph(brimg, 'dilate', breakPointDilate);

brReg= regionprops(brimg, 'Area', 'MajorAxisLength', 'MinorAxisLength','Orientation', 'Centroid');

brReg(vertcat(brReg.Area)~= max(vertcat(brReg.Area))) = [];

% Calculate bounding ellipse

brReg.EllipseCoords = zeros(100, 2);

t=linspace(0, 2*pi, 100);

brReg.EllipseCoords(:,1) = brReg.Centroid(1)+brReg.MajorAxisLength/2*cos(t-brReg.Orientation);

brReg.EllipseCoords(:,2)=brReg.Centroid(2)+brReg.MinorAxisLength/2*sin(t-brReg.Orientation);


else
brReg = [];

end

% Draw ellipse around break location

subplot(2,2,4);

imshow(img)
if (fl==1)
    title('Fractured location');
else
    title('Not Fractured image');
end
hold on

colormap('gray')

if ~isempty(brReg)

    plot(brReg.EllipseCoords(:,1), brReg.EllipseCoords(:,2), 'r');

end

hold off
