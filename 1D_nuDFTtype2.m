%Initialising Variables
N =40 ;
L = 80;
zeta=1;
N1 = 750;
x0 = linspace (- L / 2, L / 2, N1);
y0=0.75*exp(-(((x0+10).^2)+(0.1^2))/2)+exp(-(((x0-10).^2)+(0.1^2))/2)+0.6*exp(-(((x0-25).^2)+(0.1^2))/2);
plot(x0,y0)
hold on

% random sampling
p = randperm(N1);
x = x0(p(1:N))';
y = y0(p(1:N))';

% construct DFT matrix
k = [0:floor(N/2), -floor(N/2)+1:-1];
a = exp(1i * x * k * pi / L);

% solve y = a * f and obtain f, i.e. nuDFT of y.
% exactly, to stabilize the solution, solving L2 penalized
% simeq: a' * y = (a' * a + alpha^2 * I) * f,
% and its solution minimize (y - a * f)^2 + alpha^2 * |f|^2
alpha = 0.1;
z = (a' * a + alpha^2 * eye(size(a,1)));
f = z \ (a' * y);
% iDFT
a0 = exp(1i * x0' * k* pi / L);
y_est=real(a0*f);
plot(x0, y_est, x, y, '+');
xlabel("Time (t)", 'fontsize', 18);
ylabel("Signal y", 'fontsize', 18);
title("NUDFT Type 2", "N=40", 'fontsize', 18);
legend('Original Signal','Reconstructed Signal', 'fontsize', 18)
set(gca, 'fontsize', 18, 'xticklabel', get(gca, 'xtick'), 'yticklabel', get(gca, 'ytick'));

hold off

y_diff=abs(y0'-y_est);
y_mean=(abs(y0')+abs(y_est))/2;
difference=mean((y_diff./y_mean)*100);


% Calculate the percentage similarity
avg1 = mean(y0');
avg2 = mean(y_est);
similarity = (1 - sum(abs(avg1 - avg2)) / sum(avg1 + avg2)) * 100;

fprintf('The percentage similarity between the two graphs is: %f%%\n', similarity);

 
