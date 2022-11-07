function fitobj = nonlinfit(zts, Y, model, lin_result)
switch model
    case 'cosinorOneFreq'
        per1guess = lin_result.per1guess;
        a1guess   = lin_result.a1guess;
        a2guess   = lin_result.a2guess;
        [xData, yData] = prepareCurveData( zts, Y );
        % Set up fittype and options.
 %       ft = fittype( 'a0*cos(2*pi*x/a2-a1)', ...
%            'independent', 'x', 'dependent', 'y' );
        ft = fittype( 'a1*sin(2*pi*x/per1)+a2*cos(2*pi*x/per1)', ...
            'independent', 'x', 'dependent', 'y' );
        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        opts.Algorithm = 'Trust-Region';
        opts.Display = 'Off';
        %opts.Lower =      [0 0 0];
        %opts.StartPoint = [1 lin_result pi];
        %opts.Upper =      [100 1 2*pi];
        opts.Lower =      [-10 -10 0];
        opts.StartPoint = [a1guess a2guess per1guess];
        opts.Upper =      [10 10 1];
        opts.TolFun=1e-2;
        opts.TolX=1e-2;
        
        % Fit model to data.
        fitobj = fit( xData, yData, ft, opts );
end

end

