close all;

X = min(min(spectra_raw));
Y = max(max(spectra_raw-X));

figure(1);plot(wavenumber_axis,(spectra_raw-X)/Y);

hFig = figure(1);
set(hFig, 'Position', [50 50 900 700])
xlabel('Wavenumber cm^-^1','Fontsize',18);
ylabel('AU','Fontsize',18);
set(gca,'Box','off'); 
set(gca,'Linewidth',2);
set(gca,'YLim',[0 1.2]);
set(gca,'XLim',[200 1800]);
set(gca,'YTick',[]);
set(gca,'Fontsize',14);
set(gca,'PlotBoxAspectRatio',[1 0.6 1]);
title('Raw Spectra')
%set(gca,'PlotBoxAspectRatio',[1 1 1]);
%print('Raw_Spectra','-dpng','-r1000')

X = min(min(spectra_calibrated));
Y = max(max(spectra_calibrated-X));

figure(2);plot(wavenumber_axis,(spectra_calibrated-X)/Y);

hFig = figure(2);
set(hFig, 'Position', [50 50 900 700])
xlabel('Wavenumber cm^-^1','Fontsize',18);
ylabel('AU','Fontsize',18);
set(gca,'Box','off'); 
set(gca,'Linewidth',2);
set(gca,'YLim',[0 1.2]);
set(gca,'XLim',[200 1800]);
set(gca,'YTick',[]);
set(gca,'Fontsize',14);
set(gca,'PlotBoxAspectRatio',[1 0.6 1]);
title('Intensity Calibrated Spectra')
%set(gca,'PlotBoxAspectRatio',[1 1 1]);
%print('Calibrated_Spectra','-dpng','-r1000')

figure(3);plot(wavenumber_axis,(spectra_cosmic_ray_removed-X)/Y);

hFig = figure(3);
set(hFig, 'Position', [50 50 900 700])
xlabel('Wavenumber cm^-^1','Fontsize',18);
ylabel('AU','Fontsize',18);
set(gca,'Box','off'); 
set(gca,'Linewidth',2);
set(gca,'YLim',[0 1.2]);
set(gca,'XLim',[200 1800]);
set(gca,'YTick',[]);
set(gca,'Fontsize',14);
set(gca,'PlotBoxAspectRatio',[1 0.6 1]);
title('Cosmic Rays Removed')
%set(gca,'PlotBoxAspectRatio',[1 1 1]);
%print('spectra_cosmic_ray_removed','-dpng','-r1000')

figure(4);plot(wavenumber_axis,(spectra_denoised-X)/Y);

hFig = figure(4);
set(hFig, 'Position', [50 50 900 700])
xlabel('Wavenumber cm^-^1','Fontsize',18);
ylabel('AU','Fontsize',18);
set(gca,'Box','off'); 
set(gca,'Linewidth',2);
set(gca,'YLim',[0 1.2]);
set(gca,'XLim',[200 1800]);
set(gca,'YTick',[]);
set(gca,'Fontsize',14);
set(gca,'PlotBoxAspectRatio',[1 0.6 1]);
title('Denoised')
%set(gca,'PlotBoxAspectRatio',[1 1 1]);
%print('spectra_denoised','-dpng','-r1000')

X = min(min(spectra_baseline_removed));
Y = max(max(spectra_baseline_removed-X));

figure(5);plot(wavenumber_axis,(spectra_baseline_removed-X)/Y);

hFig = figure(5);
set(hFig, 'Position', [50 50 900 700])
xlabel('Wavenumber cm^-^1','Fontsize',18);
ylabel('AU','Fontsize',18);
set(gca,'Box','off'); 
set(gca,'Linewidth',2);
set(gca,'YLim',[0 1.2]);
set(gca,'XLim',[200 1800]);
set(gca,'YTick',[]);
set(gca,'Fontsize',14);
set(gca,'PlotBoxAspectRatio',[1 0.6 1]);
title('Baseline Removed')
%set(gca,'PlotBoxAspectRatio',[1 1 1]);
%print('spectra_baseline_removed','-dpng','-r1000')

clear X Y