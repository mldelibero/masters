function RLC_estimate()

% generate the data using the true model 
% here we assume the model is : R + jwL + 1/jwC 
R = 100; 
L = 0.2; 
C = 0.1; 
RLCtruth = [R, L, C]; 
freq = [0.1, 1, 10, 50].'; 
NoiseLevel = 1E-3; 
DataGeneration(freq, RLCtruth, NoiseLevel);

% estimate RLC from the data
% assume the model is known
RLC0 = [1, 1, 1];
opts = optimset('Jacobian', 'on');
[RLCnew, resnorm] = lsqnonlin(@fun_error, RLC0, [], [], opts);

disp(['Init guess RLC = ', num2str(RLC0)]); 
disp(['Estimated RLC = ', num2str(real(RLCnew))]); 
disp(['Truth RLC = ', num2str(RLCtruth)]);





% function to generate the data
function DataGeneration(freq, RLCtruth, NoiseLevel) 
R = RLCtruth(1); 
L = RLCtruth(2); 
C = RLCtruth(3);
y0 = R + j*2*pi*freq*L + 1./(j*2*pi*freq*C);
y0 = y0 + NoiseLevel*abs(y0) .* (randn(size(y0)) + j*randn(size(y0))); 
save('data.mat', 'freq', 'y0');



% optimation error function
function [err, J] = fun_error(RLC)
load('data.mat'); % load freq and y0 from data.mat
yhat = RLC(1) + j*2*pi*freq*RLC(2) + 1./(j*2*pi*freq*RLC(3)); 
err = yhat - y0; 
if nargout > 1
    J = zeros(size(freq,1), 3);
    J(:, 1) = 1;
    J(:, 2) = j*2*pi*freq;
    J(:, 3) = -1./(j*2*pi*freq) / RLC(3)^2; 
end


