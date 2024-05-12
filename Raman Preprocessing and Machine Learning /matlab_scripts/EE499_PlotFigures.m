% EE499 Plot Figures of Pre-Processing Steps
clear all;
close all;

%% Raw Spectra
figure(1);plot(1:1024,spectra_raw,LineWidth=1.5);

hFig = figure(1);
set(hFig, 'Position', [50 50 900 700])
xlabel('Detector Index (Pixel)','Fontsize',18);
ylabel('Raman Intensity (Counts)','Fontsize',18);
set(gca,'Box','off'); 
set(gca,'Linewidth',2);
%set(gca,'YLim',[0 1.2]);
set(gca,'XLim',[1 1024]);
set(gca,'YTick',[]);
set(gca,'Fontsize',14);
set(gca,'PlotBoxAspectRatio',[1 0.6 1]);
title('Raw Spectra','Fontsize',22)
%set(gca,'PlotBoxAspectRatio',[1 1 1]);
%print('Raw_Spectra','-dpng','-r1000')

%% Step 1 - Intensity Calibrated
figure(2);plot(1:1024,spectra_intensity_calibrated,LineWidth=1.5);

hFig = figure(2);
set(hFig, 'Position', [50 50 900 700])
xlabel('Detector Index (Pixel)','Fontsize',18);
ylabel('Raman Intensity (A.U.)','Fontsize',18);
set(gca,'Box','off'); 
set(gca,'Linewidth',2);
%set(gca,'YLim',[0 1.2]);
set(gca,'XLim',[1 1024]);
set(gca,'YTick',[]);
set(gca,'Fontsize',14);
set(gca,'PlotBoxAspectRatio',[1 0.6 1]);
title('Intensity Calibrated Spectra','Fontsize',22)
%set(gca,'PlotBoxAspectRatio',[1 1 1]);
%print('Raw_Spectra','-dpng','-r1000')

%% Step 2 - Wavenumber Calibration Spectra

figure(3);plot(wavenumber_axis,spectra_calibrated,LineWidth=1.5);

hFig = figure(3);
set(hFig, 'Position', [50 50 900 700])
xlabel('Raman Shift (cm^-^1)','Fontsize',18);
ylabel('Raman Intensity (A.U.)','Fontsize',18);
set(gca,'Box','off'); 
set(gca,'Linewidth',2);
%set(gca,'YLim',[0 1.2]);
set(gca,'XLim',[200 1800]);
set(gca,'YTick',[]);
set(gca,'Fontsize',14);
set(gca,'PlotBoxAspectRatio',[1 0.6 1]);
title('Wavenumber and Peak Calibrated Spectra','Fontsize',22)
%set(gca,'PlotBoxAspectRatio',[1 1 1]);
%print('Calibrated_Spectra','-dpng','-r1000')

%% Step 3 - Cosmic Ray Removal Spectra
figure(4);plot(wavenumber_axis,spectra_cosmic_ray_removed,LineWidth=1.5);

hFig = figure(4);
set(hFig, 'Position', [50 50 900 700])
xlabel('Raman Shift (cm^-^1)','Fontsize',18);
ylabel('Raman Intensity (A.U.)','Fontsize',18);
set(gca,'Box','off'); 
set(gca,'Linewidth',2);
%set(gca,'YLim',[0 1.2]);
set(gca,'XLim',[200 1800]);
set(gca,'YTick',[]);
set(gca,'Fontsize',14);
set(gca,'PlotBoxAspectRatio',[1 0.6 1]);
title('Cosmic Ray Removed Spectra','Fontsize',22)
%set(gca,'PlotBoxAspectRatio',[1 1 1]);
%print('spectra_cosmic_ray_removed','-dpng','-r1000')

%% Step 3 - De-noisied Spectra
figure(5);plot(wavenumber_axis,spectra_denoised,LineWidth=1.5);

hFig = figure(5);
set(hFig, 'Position', [50 50 900 700])
xlabel('Raman Shift (cm^-^1)','Fontsize',18);
ylabel('Raman Intensity (A.U.)','Fontsize',18);
set(gca,'Box','off'); 
set(gca,'Linewidth',2);
%set(gca,'YLim',[0 1.2]);
set(gca,'XLim',[200 1800]);
set(gca,'YTick',[]);
set(gca,'Fontsize',14);
set(gca,'PlotBoxAspectRatio',[1 0.6 1]);
title('Denoised Spectra','Fontsize',22)
%set(gca,'PlotBoxAspectRatio',[1 1 1]);
%print('spectra_denoised','-dpng','-r1000')

X = min(min(spectra_baseline_removed));
Y = max(max(spectra_baseline_removed-X));

%% Step 4 - Baseline Subtraction Spectra
figure(6);plot(wavenumber_axis,spectra_baseline_removed,LineWidth=1.5);

hFig = figure(6);
set(hFig, 'Position', [50 50 900 700])
xlabel('Raman Shift (cm^-^1)','Fontsize',18);
ylabel('Raman Intensity (A.U.)','Fontsize',18);
set(gca,'Box','off'); 
set(gca,'Linewidth',2);
%set(gca,'YLim',[0 1.2]);
set(gca,'XLim',[200 1800]);
set(gca,'YTick',[]);
set(gca,'Fontsize',14);
set(gca,'PlotBoxAspectRatio',[1 0.6 1]);
title('Baseline Removed Spectra','Fontsize',22)
%set(gca,'PlotBoxAspectRatio',[1 1 1]);
%print('spectra_baseline_removed','-dpng','-r1000')


%% Step 5 - Normalisation Spectra

X = min(min(spectra_baseline_removed));
Y = max(max(spectra_baseline_removed-X));

normalised_spectra = (spectra_baseline_removed-X)/Y;

figure(6);plot(wavenumber_axis,normalised_spectra,LineWidth=1.5);

hFig = figure(6);
set(hFig, 'Position', [50 50 900 700])
xlabel('Raman Shift (cm^-^1)','Fontsize',18);
ylabel('Raman Intensity (A.U.)','Fontsize',18);
set(gca,'Box','off'); 
set(gca,'Linewidth',2);
set(gca,'YLim',[0 1.2]);
set(gca,'XLim',[200 1800]);
set(gca,'YTick',[]);
set(gca,'Fontsize',14);
set(gca,'PlotBoxAspectRatio',[1 0.6 1]);
title('Normalised Spectra','Fontsize',22)
%set(gca,'PlotBoxAspectRatio',[1 1 1]);
%print('spectra_baseline_removed','-dpng','-r1000')
clear X Y

%% Subplot of all Spectra in Pre-processing steps
subplot(4,2,[1 2]);
plot(1:1024,spectra_raw(:,1),'Linewidth',1.5);
set(gca,'XLim',[1 1024]);
set(gca,'YTick',[]);
xlabel('Detector Index (Pixel)');
ylabel('Intensity (Counts)');
title('Raw Spectrum');

subplot(4,2,3);
plot(1:1024,spectra_intensity_calibrated(:,1),'Linewidth',1.5);
set(gca,'XLim',[1 1024]);
set(gca,'YTick',[]);
xlabel('Detector Index (Pixel)');
ylabel('Raman Intensity (A.U.)');
title('Intensity Calibrated Spectrum');

subplot(4,2,4);
plot(wavenumber_axis,spectra_calibrated(:,1),'Linewidth',1.5);
set(gca,'XLim',[200 1800]);
set(gca,'YTick',[]);
xlabel('Raman Shift (cm^-^1)');
ylabel('Raman Intensity (A.U.)');
title('Wavenumber Calibrated Spectrum');


subplot(4,2,5);
plot(wavenumber_axis,spectra_cosmic_ray_removed(:,1),'Linewidth',1.5);
set(gca,'XLim',[200 1800]);
set(gca,'YTick',[]);
xlabel('Raman Shift (cm^-^1)');
ylabel('Raman Intensity (A.U.)');
title('Removed Cosmic Ray Spectrum');

subplot(4,2,6);
plot(wavenumber_axis,spectra_denoised(:,1),'Linewidth',1.5);
set(gca,'XLim',[200 1800]);
set(gca,'YTick',[]);
xlabel('Raman Shift (cm^-^1)');
ylabel('Raman Intensity (A.U.)');
title('De-noised Spectrum');

subplot(4,2,7);
plot(wavenumber_axis,spectra_baseline_removed(:,1),'Linewidth',1.5);
set(gca,'XLim',[200 1800]);
set(gca,'YTick',[]);
xlabel('Raman Shift (cm^-^1)');
ylabel('Raman Intensity (A.U.)');
title('Baseline Subtraction Spectrum');

subplot(4,2,8);
X = min(spectra_baseline_removed(:,1));
Y = max(spectra_baseline_removed(:,1))-X;
plot(wavenumber_axis,(spectra_baseline_removed(:,1)-X)/Y,'Linewidth',1.5);
set(gca,'XLim',[200 1800]);
set(gca,'YTick',[]);
xlabel('Raman Shift (cm^-^1)');
ylabel('Intensity (AU)');
title('Normalised Spectrum');

sgtitle('Pre-processing Steps on a Nylon Raman Spectrum')

%% Plot of Recorded Pre-Processed Datasets
clear all
clc
PVC = load('PVC_matfiles/PVC_Dataset.mat');
PVC_dataset = PVC.dataset;
PVC_wn = PVC.wavenumber_axis;


PS = load('PS_matfiles/PS_Dataset.mat');
PS_dataset = PS.dataset;
PS_wn = PS.wavenumber_axis;

PMMA = load('PMMA_matfiles/PMMA_Dataset.mat');
PMMA_dataset = PMMA.dataset;
PMMA_wn = PMMA.wavenumber_axis;

figure;
plot(PVC_wn,PVC_dataset);
set(gca,'XLim',[200 1800]);
set(gca,'YTick',[]);
xlabel('Raman Shift (cm^-^1)');
ylabel('Intensity (AU)');
title('Raman Spectra of Polyvinyl Chloride (PVC)');
set(gca,'Fontsize',14);

figure;
plot(PMMA_wn,PMMA_dataset);
set(gca,'XLim',[200 1800]);
set(gca,'YTick',[]);
xlabel('Raman Shift (cm^-^1)');
ylabel('Intensity (AU)');
title('Raman Spectra of Polymethyl Methacrylate (PMMA)');
set(gca,'Fontsize',14);

figure;
plot(PS_wn,PS_dataset);
set(gca,'XLim',[200 1800]);
set(gca,'YTick',[]);
xlabel('Raman Shift (cm^-^1)');
ylabel('Intensity (AU)');
title('Raman Spectra of Polystyrene (PS)');
set(gca,'Fontsize',14);


%% Plot of dataset with mean and standard deviation
clear all
clc

data = load('PVC_matfiles/PVC_Calibrated_Dataset.mat');
dataset = data.dataset;
wn = data.wavenumber_axis;

% Calculate mean and standard deviation
spectra_mean = mean(dataset,2);
spectra_std_dev = std(dataset,0,2);

% Calculate upper and lower bounds for standard deviation
std_pos = spectra_mean + spectra_std_dev;
std_neg = spectra_mean  - spectra_std_dev;


% Plot mean and fill between upper and lower bounds
figure;
plot(wn, spectra_mean, LineWidth=1.5);
hold on;
fill([wn, fliplr(wn)],[std_pos', fliplr(std_neg')], 'b', 'FaceAlpha', 0.3);


% Set axis labels and remove y-axis ticks
xlabel('Raman Shift (cm^-^1)');
ylabel('Raman Intensity (A.U.)');
title('Raman Spectrum of Polyvinyl chloride (PVC)')
set(gca,'XLim',[200 1800]);
set(gca, 'ytick', []);
legend('$\bar{x}$','$\bar{x} \pm \sigma$','Interpreter','Latex');
set(gca,'Fontsize',14);

