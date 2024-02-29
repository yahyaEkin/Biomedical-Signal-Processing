# Biomedical-Signal-Processing
Noise reduction in ECG signals with LMS method using artificial ECG signals as reference signal.
<br> As a project for Biomedical Signal Processing Lecture, me and my teammate performed LMS adaptive signal processing algorithm to reduce noise in ECG signals with using artificially generated reference signals. We took clean ECG samples from MIT-BIH Noise Stress Test Database and compare the performance of different reference signals. To show the results, we tried to write a conference paper. 

<br>There are two versions of this code: 
<br> AdapFilt_v8 : This version is used to compute the performance of artificial reference ECG signals. We tried different signals with different frequencies to see the effect of frequency and amplitude. This version also graphs the results and give MSE for each signal.
<br> AdapFilt_v10: This version is used to graph MSE performance of seperate reference signals such as sinusoidal, triangular, rectengular and combined reference signal. Combined reference signal is the combination of three references with various frequencies. The effect of noise reduction cen be seen by this version.
<br>
The clean_ecg.mat is clean ecg sanple that we use. 
