function [fitresult, gof] = createFit_polynomial(noise_seg_adret, mean_slope_seg_adret)
%CREATEFIT1(NOISE_SEG_ADRET,MEAN_SLOPE_SEG_ADRET)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : noise_seg_adret
%      Y Output: mean_slope_seg_adret
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  另请参阅 FIT, CFIT, SFIT.

%  由 MATLAB 于 14-May-2020 16:42:35 自动生成


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( noise_seg_adret, mean_slope_seg_adret );

% Set up fittype and options.
ft = fittype( 'poly3' );

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft );

% Plot fit with data.
% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult, xData, yData );
% legend( h, 'mean_slope_seg_adret vs. noise_seg_adret', 'untitled fit 1', 'Location', 'NorthEast' );
% % Label axes
% xlabel noise_seg_adret
% ylabel mean_slope_seg_adret
% grid on


