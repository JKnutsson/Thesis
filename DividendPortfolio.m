clc, clear all, format long g
%% Dividend portfolio

%% Load data set from excel.
% Data.Properties.VariableNames(2);

Data         = readtable('Data.xlsx','Sheet','Pivot 2008-2022');
DivData      = readtable('Dividend.xlsx','Sheet','Pivot 2008-2022');

OMX     = readtable('Data.xlsx','Sheet','Pivot OMX');
RF      = readtable('RfSE10Y.xlsx');

clc
%% Sets of data.

% Mean risk free rate of swedish 10 year bonds
RF              = mean(table2array(RF(197:372,2)))/100;
% Clean up of OMXS data to fit the time period of 2008-2022, This also removes the "OMXS's" that have no data during our selected time period.

OMX(:,7)        = [];
OMX(:,3:5)      = [];
OMX(1:4265, : ) = [];
OMX(end, : )    = [];
OMX             = rmmissing(OMX);



dateEntry    = Data.RowLabels(148:end,:);        % Contain all dates from the first period of data
priceEntry   = Data( 148:end ,2:end);            % Contain all prices from the first period of data

DivDate      = DivData.RowLabels(148:end,:);    
DivYieldData = DivData( 148:end ,2:end);


% define two sets of column variable names
set1 = DivYieldData.Properties.VariableNames;
set2 = priceEntry.Properties.VariableNames;

% find the variable names that are in set1 but not in set2 and the other
% way around.
missingVars  = setdiff(set1, set2);
missingVars2 = setdiff(set2, set1);

% remove the columns in priceEntry and DivYieldData  that correspond to the missing variable names
priceEntry   = removevars(priceEntry, missingVars2);
DivYieldData = removevars(DivYieldData, missingVars);

% Converts prices and div yield to the "double" format
priceData    = table2array( priceEntry );  
DivYield     = table2array(DivYieldData);

DivYield(1:1403,17) = 0;


clear set1 set2 missingVars missingVars2

%%

PortSize = [10,20,30,40,50,100,150,170];
%%

for i = 1:width(PortSize)
    [time,PortValueWeighted{i,1}, PortValueUniform{i,1}, PortValueOMX] = DividedYield(priceData, DivYield, dateEntry, DivDate, OMX, PortSize(i));
end

%%

OMXS30     = PortValueOMX(:,1);
OMXSPI     = PortValueOMX(:,2);
%%

PortValueWeighted = cell2mat(PortValueWeighted);
PortValueUniform  = cell2mat(PortValueUniform);

%% Annual returns monthly

for i = 1:width(PortSize)

    % Calculate the total return over the entire investment period for the Wwighted portfolio.
    ReturnWM(i) = (PortValueWeighted(i,end)-PortValueWeighted(i,1))/PortValueWeighted(i,1);

    % Calculate the total return over the entire investment period for the Uniform portfolio.
    ReturnUM(i) = (PortValueUniform(i,end)-PortValueUniform(i,1))/PortValueUniform(i,1);

    % Calculate the total return over the entire investment period for the OMXS30 index.
    ReturnOMXS30 = (OMXS30(end)-OMXS30(1))/OMXS30(1);

    % Calculate the total return over the entire investment period for the OMXSPI index.
    ReturnUMOMXSPI = (OMXSPI(end)-OMXSPI(1))/OMXSPI(1);

    % Calculate the annual return for the Wwighted portfolio, assuming that there are 14 years.
    AnnualReturnW(i) = ((1 + ReturnWM(i))^(1 / 14) - 1 )*100;

    % Calculate the annual return for the Uniform portfolio, assuming that there are 14 years.
    AnnualReturnU(i) = ((1 + ReturnUM(i))^(1 / 14) - 1 )*100;

    % Calculate the annual return for the OMXS30 index, assuming that there are 14  years.
    AnnualReturnOMXS30 = ((1 + ReturnOMXS30)^(1 / 14) - 1 )*100;

    % Calculate the annual return for the OMXSPI index, assuming that there are 14 years.
    AnnualReturnOMXSPI = ((1 + ReturnUMOMXSPI)^(1 / 14) - 1 )*100;


    % The following lines initialize a counter variable to 1, and then loop through
    % the rows of the PortValueW matrix in increments of 12 (representing one year),
    % stopping at the second-to-last year.

    j=1;
    for l = 1:12:height(PortValueWeighted')-12

           % Calculate the return for the W portfolio over the next 12 months.
           ReturnsWeighted(i,j) = (PortValueWeighted(i,l+12) - PortValueWeighted(i,l))/ PortValueWeighted(i,l);

           % Calculate the return for the U portfolio over the next 12 months.
           ReturnsUniform(i,j)  = (PortValueUniform(i,l+12) - PortValueUniform(i,l))/ PortValueUniform(i,l);

           % Calculate the return for the OMXS30 index over the next 12 months.
           ReturnsOMXS30(j)   = (OMXS30(l+12) - OMXS30(l))/ OMXS30(l);

           % Calculate the return for the OMXSPI index over the next 12 months.
           ReturnsOMXSPI(j)   = (OMXSPI(l+12) - OMXSPI(l))/ OMXSPI(l);

           % Increment the counter variable by 1.
           j = j + 1;

    end


    % Calculate the annualized standard deviation (volatility) for each portfolio,
    % which is the standard deviation of monthly returns multiplied by the square root of 12.
    stdWeightedM(i) = std(ReturnsWeighted(i,:))*sqrt(12);
    stdUniformM(i)  = std(ReturnsUniform(i,:))*sqrt(12);
    stdOMX30        = std(ReturnsOMXS30)*sqrt(12);
    stdOMXSPI       = std(ReturnsOMXSPI)*sqrt(12);

    % Calculate the Sharpe ratio as a measure of risk-adjusted return
    % for each portfolio, where AnnualReturnW/U/OMXS30/OMXSPI is the annualized return 
    % for each portfolio and RF is the risk-free rate.
    SharpeWeightedMonthly(i) = ((AnnualReturnW(i)/100) - RF)/stdWeightedM(i);
    SharpeUniformMonthly(i)  = ((AnnualReturnU(i)/100) - RF)/stdUniformM(i);
    SharpeOMXS30             = (AnnualReturnOMXS30/100 - RF)/stdOMX30;
    SharpeOMXSPI             = (AnnualReturnOMXSPI/100 - RF)/stdOMXSPI;

end
%% Individual Plots
color_order = [0.8500 0.3250 0.0980; % Orange
               0.0 0.1 0.0; % Dark green
               0.1176 0.9 1.0000; % Blue
               0.5412 0.1686 0.8863; % Purple
               0.0 0.0 1; % Turquoise
               1.0000 0.0 0.0000; % Red
               0.8431 0.0784 1; % Magenta
               0.1176 0.5176 0.2627; % Green
               1, 0, 1]; % Magenta


%% Uniform Normal

% Set custom color order as default
set(groot, 'defaultAxesColorOrder', color_order);


fig = gcf;
set(fig, 'Units', 'normalized', 'Position', [0.15, 0.15, 0.9, 0.9]); % adjust the numbers to scale down the figure
ax = gca;
set(ax, 'Units', 'normalized', 'Position', [0.15, 0.15, 0.8, 0.75]); % adjust the numbers to scale down the axis

plot(time, PortValueUniform, 'LineWidth', 1);
hold on
plot(time, OMXS30, 'LineWidth', 1)
plot(time, OMXSPI, 'LineWidth', 1)
hold off
xticks(time([1:12:end]));
set(gca, 'XTick', time([1:12:end]),'FontSize', 10);

datetick('x', 'dd-mm-yyyy', 'keepticks');
xtickangle(45);
ax = gca;
ax.Position(2) = ax.Position(2) + 0.05; % adjust the number to create enough space
grid on;

legend('Uniform model (10 stocks)', 'Uniform model (20 stocks)', 'Uniform model (30 stocks)', 'Uniform model (40 stocks)', 'Uniform model (50 stocks)', 'Uniform model (100 stocks)', 'Uniform model (150 stocks)', 'OMXS30', 'OMXSPI')
legend('Location', 'South', 'Orientation', 'Horizontal','NumColumns', 5, 'Units', 'normalized', 'Position', [0.48 0 0.1 0.1], 'Box', 'off', 'FontSize', 16)


%title('Rebalance monthly - Normal scale', 'FontSize', 16);
ytickformat('%,.0f SEK')
ylim([50000 2600000])
set(gca, 'yticklabel', sprintfc('%gM SEK', get(gca, 'ytick')/1e6), 'FontSize', 14 )
% Add a y-axis label
ylabel('Portfolio value (SEK)')
ylabel_handle = get(gca, 'ylabel');
current_position = get(ylabel_handle, 'position');
set(ylabel_handle, 'position', [current_position(1)-300, current_position(2), current_position(3)]);


print -depsc2 Dividend_Uniform_Normal

%% Uniform Log


% Set custom color order as default
set(groot, 'defaultAxesColorOrder', color_order);


fig = gcf;
set(fig, 'Units', 'normalized', 'Position', [0.15, 0.15, 0.9, 0.9]); % adjust the numbers to scale down the figure
ax = gca;
set(ax, 'Units', 'normalized', 'Position', [0.15, 0.15, 0.8, 0.75]); % adjust the numbers to scale down the axis


semilogy(time, PortValueUniform, 'LineWidth',1);
hold on
semilogy(time, OMXS30, 'LineWidth',1)
semilogy(time, OMXSPI, 'LineWidth',1)
hold off
xticks(time([1:12:end]));
set(gca, 'XTick', time([1:12:end]),'FontSize', 10);
datetick('x', 'dd-mm-yyyy', 'keepticks');
xtickangle(45);
ax = gca;
ax.Position(2) = ax.Position(2) + 0.05; % adjust the number to create enough space
grid on;

legend('Uniform model (10 stocks)', 'Uniform model (20 stocks)', 'Uniform model (30 stocks)', 'Uniform model (40 stocks)', 'Uniform model (50 stocks)', 'Uniform model (100 stocks)', 'Uniform model (150 stocks)', 'OMXS30', 'OMXSPI')
legend('Location', 'South', 'Orientation', 'Horizontal','NumColumns', 5, 'Units', 'normalized', 'Position', [0.48 0 0.1 0.1], 'Box', 'off', 'FontSize', 16)

set(gca, 'yticklabel', sprintfc('%gM SEK', get(gca, 'ytick')/1e6), 'FontSize', 14 )
% Add a y-axis label
yticks = logspace(log10(min(PortValueUniform(6,:))), log10(max(PortValueUniform(6,:))), 10);
yticklabels = arrayfun(@(x) sprintf('%.0f%%', x), ...
                        (yticks - min(PortValueUniform(6,:))) ./ min(PortValueUniform(6,:)) * 100, ...
                        'UniformOutput', false);
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels);
ylabel('Portfolio Value Increase (%)', 'FontSize', 14)


print -depsc2 Dividend_Uniform_Log

%% Weighted Normal


% Set custom color order as default
set(groot, 'defaultAxesColorOrder', color_order);


fig = gcf;
set(fig, 'Units', 'normalized', 'Position', [0.15, 0.15, 0.9, 0.9]); % adjust the numbers to scale down the figure
ax = gca;
set(ax, 'Units', 'normalized', 'Position', [0.15, 0.15, 0.8, 0.75]); % adjust the numbers to scale down the axis


plot(time, PortValueWeighted,'LineWidth',1);
hold on
plot(time, OMXS30, 'LineWidth',1)
plot(time, OMXSPI, 'LineWidth',1)
hold off
xticks(time([1:12:end]));
set(gca, 'XTick', time([1:12:end]),'FontSize', 10);

datetick('x', 'dd-mm-yyyy', 'keepticks');
xtickangle(45);
ax = gca;
ax.Position(2) = ax.Position(2) + 0.05; % adjust the number to create enough space
grid on;

legend('Weighted model (10 stocks)', 'Weighted model (20 stocks)', 'Weighted model (30 stocks)', 'Weighted model (40 stocks)', 'Weighted model (50 stocks)', 'Weighted model (100 stocks)', 'Weighted model (150 stocks)', 'OMXS30', 'OMXSPI')
legend('Location', 'South', 'Orientation', 'Horizontal','NumColumns', 5, 'Units', 'normalized', 'Position', [0.48 0 0.1 0.1], 'Box', 'off', 'FontSize', 16)

ytickformat('%,.0f SEK')
ylim([50000 2500000])
set(gca, 'yticklabel', sprintfc('%gM SEK', get(gca, 'ytick')/1e6), 'FontSize', 14 )
% Add a y-axis label
ylabel('Portfolio value (SEK)')
ylabel_handle = get(gca, 'ylabel');
current_position = get(ylabel_handle, 'position');
set(ylabel_handle, 'position', [current_position(1)-300, current_position(2), current_position(3)]);


print -depsc2 Dividend_Weighted_Normal


%% Weighted Log


% Set custom color order as default
set(groot, 'defaultAxesColorOrder', color_order);


fig = gcf;
set(fig, 'Units', 'normalized', 'Position', [0.15, 0.15, 0.9, 0.9]); % adjust the numbers to scale down the figure
ax = gca;
set(ax, 'Units', 'normalized', 'Position', [0.15, 0.15, 0.8, 0.75]); % adjust the numbers to scale down the axis


semilogy(time, PortValueWeighted, 'LineWidth', 1);
hold on
semilogy(time, OMXS30,'LineWidth', 1)
semilogy(time, OMXSPI,'LineWidth', 1)
hold off
xticks(time([1:12:end]));
set(gca, 'XTick', time([1:12:end]),'FontSize', 10);
datetick('x', 'dd-mm-yyyy', 'keepticks');
xtickangle(45);
ax = gca;
ax.Position(2) = ax.Position(2) + 0.05; % adjust the number to create enough space
grid on;

legend('Weighted model (10 stocks)', 'Weighted model (20 stocks)', 'Weighted model (30 stocks)', 'Weighted model (40 stocks)', 'Weighted model (50 stocks)', 'Weighted model (100 stocks)', 'Weighted model (150 stocks)', 'OMXS30', 'OMXSPI')
legend('Location', 'South', 'Orientation', 'Horizontal','NumColumns', 5, 'Units', 'normalized', 'Position', [0.48 0 0.1 0.1], 'Box', 'off', 'FontSize', 16)

set(gca, 'yticklabel', sprintfc('%gM SEK', get(gca, 'ytick')/1e6), 'FontSize', 14 )
% Add a y-axis label
yticks = logspace(log10(min(PortValueWeighted(1,:))), log10(max(PortValueWeighted(1,:))), 10);
yticklabels = arrayfun(@(x) sprintf('%.0f%%', x), ...
                        (yticks - min(PortValueWeighted(1,:))) ./ min(PortValueWeighted(1,:)) * 100, ...
                        'UniformOutput', false);
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels);
ylabel('Portfolio Value Increase (%)', 'FontSize', 14)


print -depsc2 Dividend_Weighted_Log

%%


