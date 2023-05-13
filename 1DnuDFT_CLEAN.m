clear all
N=40;
L=80;

N1 = 5* N;
subplot(2,2,1)
x0 = linspace (- L / 2, L / 2, N1)';
y1=0.75*exp(-(((x0+10).^2)+(0.1^2))/2)+exp(-(((x0-10).^2)+(0.1^2))/2)+0.6*exp(-(((x0-25).^2)+(0.1^2))/2);
y1=y1';
plot(x0,y1)
title('Original Signal')

%Sinc function generation
y_sinc=sinc(0.1*x0).^2;

%Convolving the signal with the sinc function and normalising
y0=conv(y1,y_sinc, 'same');

% DFT matrix of regularly sampled data 
k = [0:floor(N/2), -floor(N/2)+1:-1];
a = exp(1i * x0 * k * pi / L);

% random sampling to get random N x N matrix from a 
randrows=randperm(size(a,1))';
randrows=sort(randrows(1:N));
a_sample = a(randrows(1:N),(1:N));

% Our corresponding y values from a
y=y0(randrows(1:N))';
x=x0(randrows(1:N));
x=sort(x);
% solve y = a * f and obtain f, i.e. nuDFT of y.
% simeq: a' * y = z * f,
% 
alpha = 0.1;
z = (a_sample' * a_sample + alpha^2 * eye(size(a_sample,1)));
f = z \ (a_sample' * y);

%iDFT
a0 = exp(1i * x0 * k* pi / L);

%y estimate using nuDFT
y_est=a0*f;
subplot(2,2,2)
y_input=real(y_est);
x0 = linspace (- L / 2, L / 2, N1)';
plot(x0,y_input);
title('Randomly Sampled in Frequency Domain')

%Initialising Max Array
yvals=zeros(1, length(x0));
yclean=zeros(1, length(x0));

%Finding the maximum y value and the corresponding x value
[y_max, index]=max(y_input);
yvals(index)=y_max;
yclean(index)=y_max;
xmax=x0(index);

%Finding the nuDFT of yvals
k = [0:floor(N/2), -floor(N/2)+1:-1];
a = exp(1i * x0 * k * pi / L);
a1=dftmtx(N1);

% solve y = a * f and obtain f, i.e. nuDFT of y.
% simeq: a' * y = (a' * a + alpha^2 * I) * f,
alpha = 0.1;
z1 = (a1' * a1 + alpha^2 * eye(size(a1,1)));
f1 = z1 \ (a1' * yvals');
%f1 is the nuDFT of yvals
% f1=f1(randrows);
f3=a0\yvals';
fnew=f-0.05*f3;
y_estnew=a0*fnew;
[y_maxnew, index]=max(y_estnew);
yclean(index)=y_maxnew;
xmax1=x0(index);

while y_maxnew>0.3*y_max
    yvals=zeros(1, length(x0));
    yvals(index)=y_maxnew;
    f3=a0\yvals';
    fnew=fnew-0.05*f3;
    y_estnew=a0*fnew;
   [y_maxnew, index]=max(y_estnew);
   if y_maxnew>yclean(index)
        yclean(index)=y_maxnew;
   end
end

%Smoothing the cleaned signal
recoveredsignal=conv(yclean, exp(-((x0).^2)/10), 'same');

%Adding back left over noise
final_image=recoveredsignal;

subplot(2,2,[3,4])
x0 = linspace (- L / 2, L / 2, N1)';

%Normalising
final_image=(1/max(final_image))*final_image;

plot(x0,final_image)
title('CLEANed Signal')

% Calculate the percentage similarity between the inpt and dirty image
avg1 = mean(y1);
avg2 = mean(y_est);
similarity = (1 - sum(abs(avg1 - avg2)) / sum(avg1 + avg2)) * 100;

fprintf('The percentage similarity between the two graphs is: %f%%\n', similarity);

% Calculate the percentage similarity between the input and clean signal
avg1 = mean(y1);
avg2 = mean(final_image);
similarity = (1 - sum(abs(avg1 - avg2)) / sum(avg1 + avg2)) * 100;

fprintf('The percentage similarity between the two graphs is: %f%%\n', similarity);


