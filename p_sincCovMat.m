function CMout = p_sincCovMat(P, freq, vsound)
% Produces covarience matrix for N First Order Sensors
% INPUT specifications:
%       P represents sensor positions. It is 3 x N matrix of real vectors
%           containing catresian coordinates.
%       vsound is the velocity of sound [m/s]
%       freq is the frequencies vector (= omega/c).
% OUTPUT specifications:
%       CMout is a 3 dimensional tensor (N x N x length(freq))
%               e.g., CM(:,:,f0) is the coavrience matrix for frequency f0. 

% D.Y. Levin, 2015.
% Speech and Acoustics Lab at Bar-Ilan University.

N     = size(P,2);     % Number of sensors.
F     = length(freq);  % Number of wavenumber bins
CMout = zeros(N,N,F);  % Prepare tensor for output

% Create matrix of interelement distances:
X = repmat (P(1,:),N,1);
Y = repmat (P(2,:),N,1);
Z = repmat (P(3,:),N,1);
Dist = sqrt((X - X').^2 + (Y - Y').^2 + (Z - Z').^2); % Distance matrix

% Calculate the covariance values via sinc function:
for f_index = 1:length(freq)
    CMout(:,:,f_index) = sinc(2 * freq(f_index) * Dist./vsound);
end