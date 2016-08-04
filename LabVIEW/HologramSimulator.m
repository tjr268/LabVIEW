
%Hologram and Reconstruction Simulation
%By: Joel Johnson
%------------------------------------------------------------------

%Simulate Object
Image = imread('moneysign.jpg') %imports image
Image = double(Image(:,:,1)) %assign image to an array
figure(1)
surf(Image) 
%colormap(gray(256)); 
title('Original Image');
[X,Y] = size(Image); %X and Y are sizes of the image
for x = 1:X;
    for y = 1:Y;
        if Image(x,y) < 10;
           Image(x,y) = 0;
           phi(x,y) = pi/2;
        else Image(x,y) > 1;
            Image(x,y) = 255
            phi(x,y) = pi;
        end
    end
end
z = 1e-6;
lambda = 633e-6;

%Phase Shift of Image
Objwave = abs(Image).*exp(1i*phi) + (0.01*lambda)*rand(64);
figure(2); 
imagesc(real(Objwave)); %displays real portion of image
colormap(gray(256)); 
title('Object Wave');

%Reference Wave
x = 1:1:X; %makes array
y = 1:1:Y; %makes array
[kx,ky] = meshgrid(x,y);
Referencewave = abs(Image).*exp(1i*(kx + ky)) + (0.01*lambda)*rand(64);

%Simulate Hologram
H = (abs(Objwave)).^2 + (abs(Referencewave)).^2 + (conj(Objwave).*Referencewave) + (conj(Referencewave).*Objwave);
Imax = max(max(H)); 
Ih = uint8(255*H/Imax);
FIh = (fftshift(fft2(Ih)));
figure(3); 
imagesc(H);
colormap(gray(256));
title('Hologram');
figure(4); 
imagesc(abs(log(FIh))); 
colormap(gray(256)); 
title('Fourier Space');

%Simulate Fourier Mask
Cursorx = 30;
Cursory = 30;
filtersize = 5;
W = zeros(X,Y);
By = (Cursorx - filtersize:Cursorx + filtersize);
Bx = (Cursory - filtersize:Cursory + filtersize);
W(round(Bx),round(By)) = 1;
FIhfilter = FIh.*W;
Bx = (X/2) - Cursorx;
By = (Y/2) - Cursory;
FIhfilter2 = circshift(FIhfilter,[round(By),round(Bx)]);
figure(5); %centers around real term (top left)
imagesc(real(FIhfilter2)); 
colormap(gray(256)); 
title('Fourier Space with Filter');

%Hologram Reconstruction using the Angular Spectrum Method
k = 2*pi/lambda;
E(x,y) = exp(1i*z*sqrt((k).^2 - (kx.*ky) - (ky.*ky)));
E = (ifft2(FIhfilter2.*E));
figure(6); 
surf(abs(E));  
title('Reconstructed Hologram Angular Spectrum');
figure(7); 
imagesc(angle(E)); 
colormap(gray(256));
title('Reconstructed Hologram Phase');
