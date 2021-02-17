function [fitresult, gof] = createFit(slope_seg_ubac, mean_noise_seg_ubac)
%CREATEFIT(SLOPE_SEG_UBAC,MEAN_NOISE_SEG_UBAC)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : slope_seg_ubac
%      Y Output: mean_noise_seg_ubac
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  另请参阅 FIT, CFIT, SFIT.

%  由 MATLAB 于 13-May-2020 16:37:18 自动生成


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( slope_seg_ubac, mean_noise_seg_ubac );

% Set up fittype and options.
ft = fittype( 'a*cosd(b*x+c)+d', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.849129305868777 0.933993247757551 0.678735154857773 0.757740130578333];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult, xData, yData );
% legend( h, 'mean_noise_seg_ubac vs. slope_seg_ubac', 'untitled fit 1', 'Location', 'NorthEast' );
% % Label axes
% xlabel slope_seg_ubac
% ylabel mean_noise_seg_ubac
% grid on


