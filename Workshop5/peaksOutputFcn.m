function varargout = peaksOutputFcn(varargin)
%   Copyright 2015 The MathWorks, Inc.
% Output function that plots the iterates of the optimization algorithm.

% pattern search
% [stop,options,optchanged] = psoutputhistory(optimvalues,options,flag,interval)
%
% simulated annealing
% [stop,options,optchanged] = myfun(options,optimvalues,flag)
%
% genetic algorithms
% [state,options,optchanged] = myfun(options,state,flag,interval)
%
% optimization toolbox
% stop = outfun(x, optimValues, state)

stop = false;
varargout{1} = stop;
switch nargin
    case 3 % which solver type?
        if isnumeric(varargin{1}) % optimization toolbox
            % stop = outfun(x, optimValues, state)
            x = varargin{1};
            f = varargin{2}.fval;
            state = varargin{3};
        elseif isfield(varargin{1},'x') % pattern search
            % [stop,options,optchanged] =
            % psoutputhistory(optimvalues,options,flag,interval)
            optimValues = varargin{1};
            state = varargin{3};
            x = optimValues.x;
            f = peaksObj(x);
            varargout{3} = false;
        elseif isfield(varargin{2},'Population') % genetic algorithms
            % [state,options,optchanged] = myfun(options,state,flag,interval)
            x = varargin{2}.Population;
            state = varargin{3};
            f = peaksObj(x);
            varargout{1} = varargin{2};
            varargout{2} = varargin{1};
            varargout{3} = false;
        else
            % simulated annealing
            % [stop,options,optchanged] = myfun(options,optimvalues,flag)
            optimValues = varargin{2};
            x = optimValues.x;
            state = varargin{3};
            varargout{3} = false;
            f = peaksObj(x);
        end
    case 4 % pattern search or genetic algorithms
        if isstruct(varargin{2}) % pattern search
            % [stop,options,optchanged] =
            % psoutputhistory(optimvalues,options,flag,interval)
            x = varargin{2}.Population;
            state = varargin{3};
            f = peaksObj(x);
            varargout{1} = varargin{2};
            varargout{2} = varargin{1};
            varargout{3} = false;
        else % genetic algorithms
            % [state,options,optchanged] = myfun(options,state,flag,interval)
            x = varargin{2}.Population;
            state = varargin{3};
            f = peaksObj(x);
            varargout{3} = false;
        end
    otherwise
        error('peaksOutputFcn:undefinedInputs',...
            'Input syntax is invalid.  Type ''help peaksOutputFcn'' for correct calling syntax');
end

switch state
    case 'init'
        % Plot objective function surface
        PlotSurface(x,f,state);
        pause
    case 'iter'
        % Update surface plot to show current solution
        PlotUpdate(x,f,state);
    case 'done'
        % After optimization, display solution in plot title
        PlotUpdate(x,f,state);
        DisplayTitle(x,f,state)
        try
            rmpref('peaksNonlinearPlot');
        catch
            % do nothing
        end
end

end

% -------------------------------------------------------------------------
% helper function PlotSurface
% -------------------------------------------------------------------------
function PlotSurface(x,z,state)

% Check to see if figure exists, replace if it is
h = findobj('Tag','PlotIterates');
if ~isempty(h)
    figure(h), hold off
end

% Plot the Peaks function
ezsurfc(@(x,y) peaksObj([x y]),[-3 3],30)
set(gcf,'Tag','PlotIterates')
%colormap 'autumn'
xlabel('X1')
ylabel('X2')
view([-33 22])
hold on

% check to see if nonlinear constraint should be plotted
try
doPlotCon = getpref('peaksNonlinearPlot','doplot');
if ~isempty(doPlotCon)
    if doPlotCon
        % plot contraint once
        axScale = axis;
        [xc,yc] = meshgrid(axScale(1):0.04:axScale(2),...
            axScale(3):0.04:axScale(4));
        
        zc = peaksCon([xc(:),yc(:)]);
        % convert constraint to matris for plotting
        [r,c] = size(xc);
        zc = reshape(zc,r,c);
        % only want to know where positive values reside
        %indx = find(zc<=0);
        zc(zc<=0) = NaN;
        zc(zc>0)= min(get(gca,'ZLim'));
        surf(xc,yc,zc,'EdgeColor','none','FaceAlpha',0.5);
    end
    % do nothing, purposely set to false so don't plot again
end
catch %do nothing
end
hold off
PlotUpdate(x,z,state)
DisplayTitle(x,z,state)
end

% -------------------------------------------------------------------------
% helper function PlotUpdate
% -------------------------------------------------------------------------
function PlotUpdate(x,z,state)

% Check to see if figure exists, if not, create
h = findobj('Tag','PlotIterates');
if isempty(h)
    PlotSurface(x,z,state)
    set(gcf,'Tag','PlotIterates')
    h = findobj('Tag','PlotIterates');
end

% Update Plot with New Point
figure(h)
hold on

switch state
    case 'init'
        ptColor = [0 0 0]; % black
        mrkSize = 14;
    case 'iter'
        ptColor = 'm'; % magenta
        mrkSize = 4;
    case 'done'
        ptColor = [1 0 0]; 
        mrkSize = 14;
    otherwise
        ptColor = [0 0 1]; % blue
        mrkSize = 14;
end
plot3(x(:,1),x(:,2),z,'MarkerFaceColor',ptColor,'MarkerSize',mrkSize,...
    'Marker','diamond',...
    'LineStyle','none',...
    'Color',ptColor);
ax1 = findobj('Tag','LowerContour');
minz = min(get(gca,'ZLim'));
minz = repmat(minz,size(x(:,1)));
if isempty(ax1) || length(x)>2
    if length(x) > 2
        marktype = 'k.';
    else
        marktype = 'k.-';
    end
    ax1 = plot3(x(:,1),x(:,2),minz,marktype,'MarkerSize',16);
    set(ax1,'Tag','LowerContour');
else
    set(ax1, 'XData', [get(ax1,'XData'),x(1)]);
    set(ax1, 'YData', [get(ax1,'YData'),x(2)]);
    set(ax1, 'ZData', [get(ax1,'ZData'),minz]);       
end

switch state
    case 'init'
        hold on
        plot3(x(:,1),x(:,2),minz,'k.','MarkerSize',30); % init point
        hold off
    case 'done'
        hold on
        plot3(x(:,1),x(:,2),minz,'r.','MarkerSize',30); % final point
        hold off
end

alpha 0.5

h_surface = findobj('EdgeColor', [0 0 0]);
if(double(h_surface))
    h_surface.EdgeColor = [0.5020    0.5020    0.5020];
end

DisplayTitle(x,z,state)
end

% -------------------------------------------------------------------------
% helper function DisplayTitle
% -------------------------------------------------------------------------
function DisplayTitle(x,z,state)

% Check to see if figure exists, if not, create
h = findobj('Tag','PlotIterates');
if isempty(h)
    PlotUpdate(x,z,state)
end

h = findobj('Tag','PlotIterates');
% Update Plot Title
figure(h), hold on
z0 = min(get(gca,'ZLim'));

switch state
    case 'iter'
        str = 'Current';
    case 'init'
        str = 'Iniital';
        text(x(1),x(2),z0*1.1,'Start')
    case 'done'
        str = 'Final';
        text(x(1),x(2),z0*1.1,'End')
end

str = sprintf('%s x = [%6.4f %6.4f], Objective = %6.4f',str, x(1),x(2),z);
title(str)

end