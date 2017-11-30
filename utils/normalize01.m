function y = normalize01(x)
% Normalizes Input data

% Only consider upper range 1
% y = abs(x)/max(abs(x(:)));

% True scaling between 0 and 1
x = abs(x);
minim = min(x(:));
maxim = max(x(:));
y = (x-minim)./(maxim-minim);
y = abs(y);

