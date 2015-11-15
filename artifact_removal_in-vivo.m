function data_new = artifact_removal(data_art, Fs)

% Inputs
% 1. data_art = Single channel raw in-vivo artifactual neural data as a vector
% data length must be an integer multiple of 2^10 
% 2. Fs = Sampling frequency in Hz(e.g. Fs = 40e3 if sampling freq is 40kHz
% Output
% data_new = Artifact Reduced Reconstructed Neural Data


%%  Initial Filtering and Treshold Calculation
    data_length = 2^nextpow2 (length(data_art));
    if (data_length > length(data_art))
        data_length = 2^(nextpow2 (length(data_art))-1);
    end
    sig_us = data_art(1:data_length) - mean(data_art(1:data_length));
   
    x_bpf = bandpass_filter(sig_us, Fs, 150, 300, 512);
    x_bpf_sp = bandpass_filter(sig_us, Fs, 300, 6e3, 512);

    
    sp = median(abs(x_bpf_sp))/0.6745;
    thrp = sp*sqrt(2*log10(length(x_bpf_sp)));
    
   
    s = median(abs(x_bpf))/0.6745;
    thr = s*sqrt(2*log10(length(x_bpf)));
    x_hpf = highpass_filter(sig_us, Fs, 5e3, 512);
    
    sh = median(abs(x_hpf))/0.6745;
    thrh = sh*sqrt(2*log10(length(x_hpf)));
 
%% Do Stationary wavelet decomposition/transform (SWT)
tic;
% N = floor(log2(length(sig_us)));   % No of decomposition level
N = 10;
wave_name = 'haar'; 

%%--SWT
 [A,D] = swt(sig_us,N,wave_name);
[Lo_D,Hi_D,Lo_R,Hi_R] = wfilters(wave_name);


%% Aprox Coef Thresholding
A_new = A(end,:); A_old = A(end,:); clear A;

min_ratio = min((A_new))/(median(abs(A_new))/0.6745);
max_ratio = max((A_new))/(median(abs(A_new))/0.6745);
% avg_ratio = abs((max_ratio + abs(min_ratio))/2);
avg_ratio = max(abs(A_new));
thr_ratio = 2;

    if  ( avg_ratio > thr_ratio )    
        k1 = 0.75;
    elseif  ( avg_ratio > 2*thr_ratio )
        k1 = 0.5;
    else
        k1 = 1;     
    end
% k1 = 0.8;
sigma = median(abs(A_new))/0.6745;
T = k1*sqrt(2*log10(length(A_new))*sigma^2);
%     if(T > thr_ratio)
%         T = thr_ratio;
%     end
id = find((abs(A_new)> T)==1);
tau2 = 0.8;
lamda2 = 1.1*T;
% if(abs(x_bpf(id)) < 0.75*thr) 
% %        if(abs(x_bpf(idx)) < 2.5*thr)
%             id = [];
%         else
%             id = id;
% end
A_new(id) = 0;
% A_new(id) =  T.^2./A_new(id); 
% A_new(id) =  (sign(A_new(id)).*(abs(A_new(id)-T)))./(1 + exp(-tau2*(abs(A_new(id)-lamda2))));
D_new = D;

    
for i=1:N
%% Double verification on artifacts or spikes 
    %%--Weighted Thresholding
    if (i == 3 || i == 4 || i == 5 || i == 6)     % D4, D5, D6 contains spike data, so high threshold
        k2 = 3;
    else
        k2 = 1;     % others more likely to be artifacts, so low threshold value
    end
% k2 = 1;
%%----verification ends here-------
sigma_sq = median(abs(D(i,:)))/0.6745;
%     Th = (k2 - 0.1*i)*sqrt(2*log10(length(D))*sigma_sq^2);
Th(i) = k2*sqrt(2*log10(length(D))*sigma_sq^2);

idx = find((abs(D(i,:))> Th(i))==1);
tau = 0.01;
lamda1 = 1.1*Th(i);


%% Please uncomment below to include double verification of identifying artifact
 
       if(abs(x_bpf(idx)) < 2*thr | abs(x_hpf(idx)) < 2*thrh) 
            
           if(abs(x_bpf_sp(idx)) > thrp)
                idx = [];
           end
           
       else
            idx = idx;
       end
% if (i == 4 || i == 5 || i == 6)
%     D_new(i,idx) = D(i,idx);
% else
%   D_new(i,idx) = 0;                                           % Hard

%  D_new(i,idx) =  sign(D(i,idx)).* abs(D(i,idx) - Th);  % Soft
D_new(i,idx) =     Th(i).^2./D(i,idx);               % Garrote
% D_new(i,idx) =   (sign(D(i,idx)).*(abs(D(i,idx)-Th)))./(1 + exp(-tau*(abs(D(i,idx)-lamda1))));    % SBSS

end
% end

%% Reconstruction

%%-SWT based reconstruction, ISWT
X_new = iswt(A_new, D_new, Lo_R, Hi_R); %X = ISWT(SWA(end,:),SWD,Lo_R,Hi_R)

data_new = X_new - mean(X_new);

%% Plot results

plot(data_new);hold on;plot(sig_us,'r');
xlabel('Time sample'); ylabel('Amplitude');
legend('Data New', 'Artifactual Data');




%%%%%%%%%%%%%%%%%%%%% Band Pass Filter Code %%%%%%%%%%%%%%

% Fs: sampling freq
% Fc1, Fc2: cut-off freq
% order: size of filter

function data_bpf = bandpass_filter(data, Fs, Fc1, Fc2, order)
    fraction1 = 2 * Fc1 / Fs;         
    fraction2 = 2 * Fc2 / Fs;         
    h = fir1(order, [fraction1 fraction2], 'bandpass');
    data_bpf = conv(h, data);
    to_remove = [1:order/2 length(data_bpf)-order/2+1:length(data_bpf)]; %%%% to remove side
    data_bpf(to_remove)=[];
    
%     % Adjust length
%     left_cut = round(order/2);
%     right_cut = order - left_cut;
%     data_bpf(1:left_cut) = [];
%     data_bpf(length(data_bpf) - right_cut + 1:length(data_bpf))= [];
end




%%%%%%%%%%%%%%%%%%%%% High Pass Filter Code %%%%%%%%%%%%%%

% Fs: sampling freq
% Fc: cut-off freq
% order: size of filter

function data_hpf = highpass_filter(data, Fs, Fc, order)
    fraction = 2 * Fc / Fs;         
    h = fir1(order, fraction, 'high');
    data_hpf = conv(h, data); 
     to_remove = [1:order/2 length(data_hpf)-order/2+1:length(data_hpf)]; %%%% to remove side
    data_hpf(to_remove)=[];
end


end
