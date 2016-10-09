clear all;

    
% setup 
stored_atf = 0;
t60             = 0.3;
nfft            = 1024;
loading         = 1e-6;
nstart          = 20; 
M_arr           = 1:1:9;
theta_desired   = 0;
theta_desired = theta_desired*pi/180;

% DO NOT EDIT
dtheta      = 2;
radius      = 0.8;
Navg        = 100;
M           = length(M_arr);
spacing     = 0.05*(M_arr(2) - M_arr(1));

c           = 340;
f           = 1e3; FS = 16e3;
lambda      = c/f;
room        = [4,4,4];
theta_arr   = (0:dtheta:360)*pi/180; theta_arr = theta_arr(1:end-1);

hd5_filename_direct = 'orig_atf_0.h5';
% model parameters
Dx = room(1); Dy = room(2); Dz = room(3);
V = Dx*Dy*Dz;
A = 2*(Dx*Dy + Dx*Dz + Dy*Dz);

epsilon = 0.161*V/A/t60;
alpha   = (1 - epsilon)/(pi*epsilon*A);             % variance
if t60 == 0
    alpha = 0;
end
fsch = 2e3*sqrt(t60/V);


% frequency domain analysis
nfft_max    = 2^ceil(log2(0.8*FS)); 
neff        = nfft/2 + 1; 
freqs       = (0:neff - 1) / nfft * FS;   % vector of stft freqs
if (freqs(nstart) < fsch)
    error('the model doesnt hold for this nstart');
end

% noise
[mic_pos, ~]  = p_absolute_position_from_relative(M, [0,0,0], spacing, [0,0,0], 0);
phivv   = p_sincCovMat(mic_pos', freqs, c);
s       = zeros(M, M, neff);
phinn   = zeros(M, M, neff);
iphinn  = zeros(M, M, neff);

for i = nstart:neff
    phinn(:,:,i)    = (phivv(:,:,i) + loading*eye(M));
    iphinn(:,:,i)   = inv(phinn(:,:,i));
    s(:,:,i)        = chol(iphinn(:,:,i));
end

% init structures
init_zeros              = zeros(nfft, length(theta_arr));
B_expectation_theory    = init_zeros;
B_power_theory          = init_zeros;

% desired speaker 
theta_arr           = (0:dtheta:360)*pi/180; theta_arr = theta_arr(1:end-1);
theta_desired_ind	= find(theta_desired == theta_arr);


dfft = floor(nfft_max/nfft);
source_rel_position = [radius*cos(theta_desired), radius*sin(theta_desired), 0];

[mic_pos, source_position]  = p_absolute_position_from_relative(M, [0,0,0], spacing, source_rel_position, 0);
% theoretical model 
if stored_atf
    atf_direct_complete = h5read(hd5_filename_direct, '/dataset');%, [1, 1, 1, j], [M, rir.ns, 1 1]);
    atf_desired_direct = h5read(hd5_filename_direct, '/dataset', [1, 1, 1, theta_desired_ind], [inf, inf, 1 1]);
    atf_desired_direct = squeeze(atf_desired_direct);
    atf_desired_direct_fft	= fft(atf_desired_direct, nfft_max, 2);
    atf_desired_direct_fft  = atf_desired_direct_fft(M_arr, 1:dfft:end);
else
    [mic_pos, source_position]  = p_absolute_position_from_relative(M, [0,0,0], spacing, source_rel_position, 0);

    source_r = sqrt(sum((repmat(source_position, M, 1) - mic_pos).^2, 2));
    atf_desired_direct_fft = zeros(M, neff);
    for i = 1:neff
        atf_desired_direct_fft(:,i) = 1./(4*pi*source_r).*exp(1i*2*pi*source_r./c*freqs(i));
    end 
end
for j = 1:length(theta_arr)    
    source_t_rel_position = [radius*cos(theta_arr(j)), radius*sin(theta_arr(j)), 0];
    if stored_atf
        atf_direct = atf_direct_complete(:,:,:,j);
        atf_direct = squeeze(atf_direct);

        atf_direct_fft          = fft(atf_direct, nfft_max, 2);
        atf_direct_fft          = atf_direct_fft(M_arr, 1:dfft:end);
    else
        [mic_pos, source_t_position]  = p_absolute_position_from_relative(M, [0,0,0], spacing, source_t_rel_position, 0);

        source_r = sqrt(sum((repmat(source_t_position, M, 1) - mic_pos).^2, 2));
        atf_direct_fft = zeros(M, neff);
        for i = 1:neff
            atf_direct_fft(:,i) = 1./(4*pi*source_r).*exp(1i*2*pi*source_r./c*freqs(i));
        end
    end
    ddt = sqrt(sum((source_rel_position - source_t_rel_position).^2));
    gamma_dt = sinc(2*ddt*FS/(c*nfft)*[0:nfft-1].');
    
    
     
    Gamma_mm = zeros(M, M, nfft);
    for i = 1:M
        for k = 1:M
            dmm = sqrt(sum((mic_pos(i,:) - mic_pos(k,:)).^2));
            gamma_mm = sinc(2*dmm*FS/(c*nfft)*[0:nfft-1].');
            Gamma_mm(i, k, :) = reshape(gamma_mm, 1, 1, nfft);
        end
    end
    
   

    w_mvdr_direct_num = zeros(M, neff);
    T = zeros(neff,1);
    for i = nstart:neff
        w_mvdr_direct_num(:,i) = phinn(:,:,i)\atf_desired_direct_fft(:,i);
        T(i) = trace(iphinn(:,:,i)*Gamma_mm(:,:,i));
    end
    lambda_norm = (dot(w_mvdr_direct_num, atf_desired_direct_fft(:,1:neff)).' + T*alpha);
    B_expectation_theory(1:neff,j)   = (dot(w_mvdr_direct_num, atf_direct_fft(:,1:neff)).' + T.*gamma_dt(1:neff)*alpha)./lambda_norm;
    

    for kk = nstart:neff
        sGamma_mms = s(:,:,kk)*Gamma_mm(:,:,kk)*s(:,:,kk)';
        satf_direct_fft = s(:,:,kk)*atf_desired_direct_fft(:,kk);
        satf_direct_t_fft = s(:,:,kk)*atf_direct_fft(:,kk);
        xi = (norm(sGamma_mms,'fro').^2*alpha^2 + alpha*satf_direct_fft'*sGamma_mms*satf_direct_fft...
            + alpha*satf_direct_t_fft'*sGamma_mms*satf_direct_t_fft)./lambda_norm(kk)^2;
        B_power_theory(kk,j) = xi + B_expectation_theory(kk,j)*conj(B_expectation_theory(kk,j));
    end

end    
    


figure; plot(theta_arr*180/pi - 180, fftshift(squeeze(10*log10(real(mean(B_power_theory(nstart:neff,:),1))))))

% 




