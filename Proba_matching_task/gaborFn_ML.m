function [im,imRGB] = gaborFn_ML(imSize, lamda, sigma, theta, phase, trim, discon)
% Simple function for making a gabor patch
%
% >>  im = gaborFn(imSize, lamda, sigma, theta, phase);
%
%     imSize = 100                           image size: n X n
%     lamda = 10                             wavelength (number of pixels per cycle)
%     theta = 15                             grating orientation
%     sigma = 10                             gaussian standard deviation in
%     pixels, if 0 then don't apply gaussian
%     phase = .25                            phase (0 -> 1)
%     [trim = .005]                          trim off edges of gaussian with values < trim
%     gnoise = .25                          add a gaussian noise
%     discon = .75                          diminish the constrast

% computed variables
freq = imSize/lamda;                    % compute frequency from wavelength
thetaRad = (theta / 360) * 2*pi;        % convert theta (orientation) to radians
% make linear ramp
X = 1:imSize;                           % X is a vector from 1 to imageSize
X0 = (X / imSize) - .5;                 % rescale X -> -.5 to .5
phaseRad = (phase * 2* pi);         % convert to radians: 0 -> 2*pi
% make 2D matrices
[Xm,Ym] = meshgrid(X0, X0);             % 2D matrices

% 4 Change orientation by adding Xm and Ym together in different proportions
Xt = Xm * cos(thetaRad);                % compute proportion of Xm for given orientation
Yt = Ym * sin(thetaRad);                % compute proportion of Ym for given orientation
XYt =  Xt + Yt;                      % sum X and Y components
XYf = XYt * freq * 2*pi;                % convert to radians and scale by frequency
grating = sin(XYf + phaseRad);                   % make 2D sinewave

im = grating;
if sigma > 0   % if sigma == 0 then don't apply a gaussian
    %Make a gaussian mask
    s = sigma / imSize;                     % gaussian width as fraction of imageSize
    gauss = exp( -(((Xm.^2)+(Ym.^2)) ./ (2* s^2)) ); % formula for 2D gaussian
    gauss(gauss < trim) = 0;                 % trim around edges (for 8-bit colour displays)
    
    % Now multply grating and gaussian to get a GABOR
    im = grating.* gauss;                % use .* dot-product
end
im = im.*discon;

max_im = ones(size(im,1),size(im,2));
im = min(im,max_im);
im = max(im,-max_im);

im = (im + 1) / 2;            % convert values -1:1 to 0:1


imRGB(:,:,1) = im;     % reassemble into 3-column array
imRGB(:,:,2) = im;     % reassemble into 3-column array
imRGB(:,:,3) = im;     % reassemble into 3-column array


