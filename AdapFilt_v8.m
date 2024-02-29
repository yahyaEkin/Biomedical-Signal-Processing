fs = 360; % Sampling frequency (Hz)
t = 0:1/fs:10-1/fs; % Time vector for 10 seconds
% figure;
% plotspec(clean_ecg,1/fs,"clean ecg");

% Adding Artificial and AWGN Noises to clean ECG--------------------------
SNR_VAL = 50;
noisy_ecg = awgn(clean_ecg, SNR_VAL, "measured");
SNR_Val = snr(clean_ecg,noisy_ecg);



% Artificial Reference Signals------------------------------------------

% sinusoidal_ref = 0.25 * sin(2 * pi * freq * t);
% square_ref = 0.25 * square(2 * pi * freq * t);
% triangular_ref = 0.25 * sawtooth(2 * pi * freq * t, 0.5);


% Adaptive Filtering START ----------------------------------------------


% Initialize variables for LMS adaptive filtering
min_error = Inf; % Initialize minimum error as infinity
optimal_order = 16; % Start with the initial filter order
optimal_W = zeros(optimal_order, 1); % Initialize optimal filter coefficients
mu = 0.02; % LMS step size
amp = 0.70; freq = 1;
% Define the maximum order to test
max_order = 96;
mse_array = zeros(1,max_order-optimal_order);
artificial_reference_ecg = amp*sin(2 * pi * freq * t)
%0.25 * square(2 * pi * 50 * t); 11hz
%0.70 * sin(2 * pi * 9 * t); 11hz
%0.9*sawtooth(2 * pi * freqArray(i) * t, 0.5);

% Loop over different filter orders to find the optimal one
for test_order = optimal_order:max_order
    % Initialize filter coefficients for this test order
    W_test = zeros(test_order, 1);
    
    % Apply LMS algorithm for the current order using the artificial reference
    for n = test_order+1:length(noisy_ecg)
        x = artificial_reference_ecg(n-test_order+1:n)'; % Input signal window from the artificial reference
        d = noisy_ecg(n); % Desired output (noisy ECG)
                                                                                                                                                                                                                           
        % Calculate predicted output
        y = W_test' * x;
        
        % Error signal
        e = d - y;
        
        % Update filter coefficients using LMS algorithm
        W_test = W_test + mu * e * x;
    end
    
    % Filter the signal using the current filter coefficients
    filtered_ecg_test = filter(W_test, 1, noisy_ecg);
    
    % Calculate the mean squared error (MSE) for the current filter
    mse_test = immse(artificial_reference_ecg, filtered_ecg_test);
    mse_array(test_order)=mse_test;

    % Check if the current MSE is lower than the minimum error found so far
    if mse_test < min_error
        min_error = mse_test; % Update minimum error
        optimal_order = test_order; % Update optimal filter order
        optimal_W = W_test; % Update optimal filter coefficients
    end
    
end

% Filter the entire signal using the obtained optimal filter coefficients
filtered_ecg_optimal = filter(optimal_W, 1, noisy_ecg);
% % Calculate the group delay
% group_delay = mean(grpdelay(optimal_W, 1));
% 
% % Adjust and remove group delay
% adjusted_filtered_ecg = filtered_ecg_optimal(group_delay + 1 : end - group_delay);
% % Trim clean_ecg to match the length of adjusted_filtered_ecg
% clean_ecg_trimmed = clean_ecg(1 : length(adjusted_filtered_ecg));
% t_adjusted = t(1:length(adjusted_filtered_ecg)); % Time vector for 10 seconds
% Mean Square Error with Clean Ecg
MSEwithCleanEcg = immse(clean_ecg, filtered_ecg_optimal);

% Display the optimal filter order and coefficients
fprintf('Optimal Filter Order: %d\n', optimal_order);
disp('Optimal Filter Coefficients:');
disp(optimal_W);
disp('MSE: ');
disp(MSEwithCleanEcg);

%Plots and Figures -------------------------------------------------

% figure
% plot(t,ecg_with_combined_noise)
% figure
% plotspec(ecg_with_combined_noise,1/fs,"ecg with combined noise");

% figure(2)
% plot(mse_freqArray,"Marker","o");
% hold on;
% plot(mse_cleanEcg,"Marker","o");
% hold off;
% grid on;
% xlabel('Artificial ECG Frequency (Hz)');
% ylabel('Mean-Square Error');
% title('MSE with Artificial reference and clean-ecg');
% legend('MSE with Artificial Reference','MSE with Clean ECG');

%Plots and Figures -------------------------------------------------
figure(3)
subplot(411)
plot(t,clean_ecg);
xlabel('Time (second)');
ylabel('Amplitude (V)');
title('Clean ECG Signal');
grid on;

subplot(412)
plot(t,noisy_ecg);
xlabel('Time (second)');
ylabel('Amplitude (V)');
title(['ECG with AWGN Noise, SNR(dB): ', num2str(SNR_VAL)]);
grid on;

subplot(413)
plot(t,filtered_ecg_optimal);
xlabel('Time (second)');
ylabel('Amplitude (V)');
title(['Adaptive Filtered Signal, MSE: ', num2str(MSEwithCleanEcg), ' Filter Order:',num2str(optimal_order)]);
grid on;

subplot(414)
plot(t,artificial_reference_ecg);
xlabel('Time (second)');
ylabel('Amplitude (V)');
title(['Artificial Reference ECG, Amplitude:',num2str(amp),'V Frequency: ',num2str(freq),'Hz']);
grid on;



function plotspec(x,Ts,figureTitle)
N=length(x);                               % length of the signal x
t=Ts*(1:N);                                % define a time vector
ssf=(ceil(-N/2):ceil(N/2)-1)/(Ts*N);       % frequency vector
fx=fft(x(1:N));                            % do DFT/FFT
fxs=abs(fftshift(fx));                          % shift it for plotting
fxs=fxs*1;
fxs = 20*log10(fxs);
subplot(2,1,1), plot(t,x)                  % plot the waveform
xlabel('seconds'); ylabel('amplitude')     % label the axes
title('Amplitude of ' + figureTitle);
grid on;
subplot(2,1,2), plot(ssf,fxs)         % plot magnitude spectrum
xlabel('frequency'); ylabel('magnitude')   % label the axes
title('Frequency Response of ' + figureTitle);
grid on;
end