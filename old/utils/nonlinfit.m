function [fitresult, gof] = nonlinfit(zts, Xdat, per1guess, per2guess)
%  Nonlinear regression for a two period harmonic regression model with
%  unknown periods
%
%  Auto-generated by MATLAB on 26-Jul-2022 17:33:11


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( zts, Xdat );

% Set up fittype and options.
ft = fittype( 'a0+a1*sin(2*pi*x/per1)+a2*cos(2*pi*x/per1)+a3*sin(2*pi*x/per2)+a4*cos(2*pi*x/per2)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Algorithm = 'Trust-Region';
opts.Display = 'Off';
opts.Lower = [0 -10 -10 -10 -10 .1 .1];
opts.StartPoint = [0 0 0 0 0  per1guess per2guess];
opts.Upper = [1 10 10 10 10 1 1];
opts.TolFun=1e-2;
opts.TolX=1e-2;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% % Plot fit with data.
% %figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult, xData, yData );
% legend( h, 'Xdat_unif_example vs. zts_unif', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
% % Label axes
% xlabel( 'zts_unif', 'Interpreter', 'none' );
% ylabel( 'Xdat_unif_example', 'Interpreter', 'none' );
% grid on

