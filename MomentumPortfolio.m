clc, clear, format long g
%% Load data set from excel.
% Data.Properties.VariableNames(2);

Data1   = readtable('Data.xlsx','Sheet','Pivot 1991-2012');
Data2   = readtable('Data.xlsx','Sheet','Pivot 2012-2022');
OMX     = readtable('Data.xlsx','Sheet','Pivot OMX');
RF      = readtable('RfSE10Y.xlsx');


%% Set variables.

RF          = mean(table2array(RF(1:372,2)))/100;
OMX(:,7)    = [];
OMX(:,3:5)  = [];
OMX         = rmmissing(OMX);

%%

dateEntry1   = Data1.RowLabels;             % Contain all dates from the first period of data
priceEntry1  = Data1( : ,2:end);            % Contain all prices from the first period of data
priceData1   = table2array( priceEntry1 );  % Converts prices to the "double" format

% Since some of the data from Eikon is incomplete and have varying "periods" of missing data for either/both price or ev/ebitda and are removed manually
% if they can't "pass" the while loop in the Backtesting section. The while loop checks if the price at a specific date is 0 
% if it is then it checks if there is a availabe price in +2/-2 steps from the current position that can be used.
% The following code removes the stocks or periods which could not pass the while loop.

priceData1(2116:2753, 133)  =   0;          % Gap in the data so removed the prices until the set it complete.

% Keep prices or remove?
 priceData1(4401:4411, 200)  =   59.25;
% priceData1(1:4411, 200)  =   0;  


dateEntry2   = Data2.RowLabels;             % Contain all dates from the second period of data
priceEntry2  = Data2( : ,2:end);            % Contain all prices from the first period of data
priceData2   = table2array( priceEntry2 );  % Converts prices to the "double" format

clc
%%
PortSize = [10,20,30,40,50,100,150];

%% Momentum portfolio with rebalancing based on seasonality (USA)

%%
for i = 1:width(PortSize)
    
    [time,PortfolioValueW{i,1}, PortfolioValueU{i,1}, MonthlySeasonalityWeighted{i,1}, MonthlySeasonalityUniform{i,1}] = Momentum(priceData1, priceData2, dateEntry1 , dateEntry2, priceEntry1, priceEntry2,PortSize(i));
    
end

%% Momentum portfolio with rebalancing at the end of each month.
for i = 1:width(PortSize)
    
    [PortfolioValueMonthlyW{i,1}, PortfolioValueMonthlyU{i,1}, PortValueOMX ] = MomentumMonthly(priceData1, priceData2, dateEntry1 , dateEntry2, OMX,PortSize(i));
end

%%
OMXS30     = PortValueOMX(:,1);
OMXSPI     = PortValueOMX(:,2);

%clear data1 data2 dateEntry1 dateEntry2 priceData1 priceData2 priceEntry1 priceEntry2

%%
PortfolioValueW = cell2mat(PortfolioValueW);

PortfolioValueU = cell2mat(PortfolioValueU);

MonthlySeasonalityWeighted = cell2mat(MonthlySeasonalityWeighted);

MonthlySeasonalityUniform = cell2mat(MonthlySeasonalityUniform);

PortfolioValueMonthlyW = cell2mat(PortfolioValueMonthlyW);

PortfolioValueMonthlyU = cell2mat(PortfolioValueMonthlyU);
%% Annual returns seasonality
for i = 1:width(PortSize)
    
    ReturnWS(i) = (MonthlySeasonalityWeighted(i,end)- MonthlySeasonalityWeighted(i,1)) / MonthlySeasonalityWeighted(i,1);
    ReturnUS(i) = (MonthlySeasonalityUniform(i,end)- MonthlySeasonalityUniform(i,1)) / MonthlySeasonalityUniform(i,1);

    AnnualReturnSeasonalityWeighted(i) = ((1 + ReturnWS(i))^(1 / 31) - 1 )*100;
    AnnualReturnSeasonalityUniform(i)  = ((1 + ReturnUS(i))^(1 / 31) - 1 )*100;

    for j = 1:height(MonthlySeasonalityUniform')-1


           ReturnsWeightedSeasonality(i,j) = (MonthlySeasonalityWeighted(i,j+1) - MonthlySeasonalityWeighted(i,j))/ MonthlySeasonalityWeighted(i,j);
           ReturnsUniformSeasonality(i,j) = (MonthlySeasonalityUniform(i,j'+1) - MonthlySeasonalityUniform(i,j))/ MonthlySeasonalityUniform(i,j);

    end


    % Calculate the annualized standard deviation of each portfolio by 
    % multiplying the monthly standard deviation with the square root of 12.
    stdSeasonalityWeightedM(i) = std(ReturnsWeightedSeasonality(i,:))*sqrt(12);
    stdSeasonalityUniformM(i)  = std(ReturnsUniformSeasonality(i,:))*sqrt(12);

    % Calculate the Sharpe ratio of each portfolio, which measures the excess return 
    % of the portfolio above the risk-free rate relative to its risk.
    % It is calculated as the difference between the portfolio's annualized return (in decimal)
    % and the risk-free rate divided by the annualized standard deviation.
    SharpeSeasonalityWeightedMonthly(i) = ((AnnualReturnSeasonalityWeighted(i)/100) - RF)/stdSeasonalityWeightedM(i);
    SharpeSeasonalityUniformMonthly(i)  = ((AnnualReturnSeasonalityUniform(i)/100) - RF)/stdSeasonalityUniformM(i);
end
%% Monthly returns for the monthly portfolio

for i = 1:width(PortSize)
    
    % Calculate the monthly returns for each portfolio,
    % which is the difference between the end and start values of the portfolio
    % divided by the start value of the portfolio.
    ReturnWM(i) = (PortfolioValueMonthlyW(i,end)-PortfolioValueMonthlyW(i,1))/PortfolioValueMonthlyW(i,1);
    ReturnUM(i) = (PortfolioValueMonthlyU(i,end)-PortfolioValueMonthlyU(i,1))/PortfolioValueMonthlyU(i,1);
   

    % Calculate the annualized return for each portfolio,
    % where ReturnWM/U/OMXS30/OMXSPI is the monthly return
    % for each portfolio, raised to the power of 1/31 (the number of months in a year),
    % minus 1, then multiplied by 100 to convert to percentage.
    AnnualReturnMonthlyWeighted(i) = ((1 + ReturnWM(i))^(1 / 31) - 1 )*100;
    AnnualReturnMonthlyUniform(i)  = ((1 + ReturnUM(i))^(1 / 31) - 1 )*100;

    % Iterate through the rows of the transpose of PortfolioValueMonthlyW,
    % from the first row to the second to last row.
    % (height(PortfolioValueMonthlyW')-1 is the number of rows - 1, 
    % since we need to skip the last row to avoid an "out of bounds" error)
    for j = 1:height(PortfolioValueMonthlyW')-1

           % Calculate the monthly return for weighted portfolio,
           % which is the difference between the value in the next row and current row
           % divided by the value in the current row.
           ReturnsWeightedMonthly(i,j) = (PortfolioValueMonthlyW(i,j+1) - PortfolioValueMonthlyW(i,j))/ PortfolioValueMonthlyW(i,j);

           % Calculate the monthly return for uniform portfolio,
           % which is the difference between the value in the next row and current row
           % divided by the value in the current row.
           ReturnsUniformMonthly(i,j)  = (PortfolioValueMonthlyU(i,j+1) - PortfolioValueMonthlyU(i,j))/ PortfolioValueMonthlyU(i,j);

    end


    % Calculate the annualized standard deviation of each portfolio by 
    % multiplying the monthly standard deviation with the square root of 12.
    stdWeightedM(i) = std(ReturnsWeightedMonthly(i,:))*sqrt(12);
    stdUniformM(i)  = std(ReturnsUniformMonthly(i,:))*sqrt(12);


    % Calculate the Sharpe ratio of each portfolio, which measures the excess return 
    % of the portfolio above the risk-free rate relative to its risk.
    % It is calculated as the difference between the portfolio's annualized return (in decimal)
    % and the risk-free rate divided by the annualized standard deviation.
    SharpeWeightedMonthly(i) = ((AnnualReturnMonthlyWeighted(i)/100) - RF)/stdWeightedM(i);
    SharpeUniformMonthly(i)  = ((AnnualReturnMonthlyUniform(i)/100) - RF)/stdUniformM(i);
end
%% OMX

    ReturnOMXS30 = (OMXS30(end-1)-OMXS30(1))/OMXS30(1);
    ReturnUMOMXSPI = (OMXSPI(end-1)-OMXSPI(1))/OMXSPI(1);
    
    for j = 1:height(PortfolioValueMonthlyW')-1
        
           % Calculate the monthly return for OMXS30 index,
           % which is the difference between the value in the next row and current row
           % divided by the value in the current row.
           ReturnsOMXS30(j)  = (OMXS30(j+1) - OMXS30(j))/ OMXS30(j);

           % Calculate the monthly return for OMXSPI index,
           % which is the difference between the value in the next row and current row
           % divided by the value in the current row.
           ReturnsOMXSPI(j)  = (OMXSPI(j+1) - OMXSPI(j))/ OMXSPI(j);
    end
    AnnualReturnOMXS30 = ((1 + ReturnOMXS30)^(1 / 31) - 1 )*100;
    AnnualReturnOMXSPI = ((1 + ReturnUMOMXSPI)^(1 / 31) - 1 )*100;
    
    stdOMXS30    = std(ReturnsOMXS30)*sqrt(12);
    stdOMXSPI    = std(ReturnsOMXSPI)*sqrt(12);
    
    SharpeOMXS30 = (AnnualReturnOMXS30/100 - RF)/stdOMXS30;
    SharpeOMXSPI = (AnnualReturnOMXSPI/100 - RF)/stdOMXSPI;



%% Figures for LaTeX (eps format)

color_order = [0.8500 0.3250 0.0980; % Orange
               0.0 0.1 0.0; % Dark green
               0.1176 0.9 1.0000; % Blue
               0.5412 0.1686 0.8863; % Purple
               0.0 0.0 1; % Turquoise
               1.0000 0.0 0.0000; % Red
               0.8431 0.0784 1; % Magenta
               0.1176 0.5176 0.2627; % Green
               1, 0, 1]; % Magenta


%% LOG SCALE
%% Uniform


% Set custom color order as default
set(groot, 'defaultAxesColorOrder', color_order);

fig = gcf;
set(fig, 'Units', 'normalized', 'Position', [0.15, 0.15, 0.9, 0.9]); % adjust the numbers to scale down the figure
ax = gca;
set(ax, 'Units', 'normalized', 'Position', [0.15, 0.15, 0.8, 0.75]); % adjust the numbers to scale down the axis


semilogy(time, PortfolioValueMonthlyU,'LineWidth',1);
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

yticks = logspace(log10(min(PortfolioValueMonthlyU(1,:))), log10(max(PortfolioValueMonthlyU(1,:))), 10);
yticklabels = arrayfun(@(x) sprintf('%.0f%%', x), ...
                        (yticks - min(PortfolioValueMonthlyU(1,:))) ./ min(PortfolioValueMonthlyU(1,:)) * 100, ...
                        'UniformOutput', false);
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels);
ylabel('Portfolio Value Increase (%)', 'FontSize', 14)



print -depsc2 MomentumMonthlyUniformLogScale

%% WEIGHTED

% Set custom color order as default
set(groot, 'defaultAxesColorOrder', color_order);

fig = gcf;
set(fig, 'Units', 'normalized', 'Position', [0.15, 0.15, 0.9, 0.9]); % adjust the numbers to scale down the figure
ax = gca;
set(ax, 'Units', 'normalized', 'Position', [0.15, 0.15, 0.8, 0.75]); % adjust the numbers to scale down the axis


semilogy(time, PortfolioValueMonthlyW,'LineWidth',1);
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
legend('Weighted model (10 stocks)', 'Weighted model (20 stocks)', 'Weighted model (30 stocks)', 'Weighted model (40 stocks)', 'Weighted model (50 stocks)', 'Weighted model (100 stocks)', 'Weighted model (150 stocks)', 'OMXS30', 'OMXSPI')
legend('Location', 'South', 'Orientation', 'Horizontal','NumColumns', 5, 'Units', 'normalized', 'Position', [0.48 0 0.1 0.1], 'Box', 'off', 'FontSize', 16)

yticks = logspace(log10(min(PortfolioValueMonthlyW(1,:))), log10(max(PortfolioValueMonthlyW(1,:))), 10);
yticklabels = arrayfun(@(x) sprintf('%.0f%%', x), ...
                        (yticks - min(PortfolioValueMonthlyW(1,:))) ./ min(PortfolioValueMonthlyW(1,:)) * 100, ...
                        'UniformOutput', false);
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels);
ylabel('Portfolio Value Increase (%)', 'FontSize', 14)



print -depsc2 MomentumMonthlyWeightedLogScale
%% Normal Scale
%% UNIFORM 

% Set custom color order as default
set(groot, 'defaultAxesColorOrder', color_order);

fig = gcf;
set(fig, 'Units', 'normalized', 'Position', [0.15, 0.15, 0.9, 0.9]); % adjust the numbers to scale down the figure
ax = gca;
set(ax, 'Units', 'normalized', 'Position', [0.15, 0.15, 0.8, 0.75]); % adjust the numbers to scale down the axis


plot(time, PortfolioValueMonthlyU,'LineWidth',1);
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

legend('Uniform model (10 stocks)', 'Uniform model (20 stocks)', 'Uniform model (30 stocks)', 'Uniform model (40 stocks)', 'Uniform model (50 stocks)', 'Uniform model (100 stocks)', 'Uniform model (150 stocks)', 'OMXS30', 'OMXSPI')
legend('Location', 'South', 'Orientation', 'Horizontal','NumColumns', 5, 'Units', 'normalized', 'Position', [0.48 0 0.1 0.1], 'Box', 'off', 'FontSize', 16)

ylim([50000 110000000])
set(gca, 'yticklabel', sprintfc('%gM SEK', get(gca, 'ytick')/1e6), 'FontSize', 14 )
% Add a y-axis label
ylabel('Portfolio value (SEK)', 'FontSize', 14)
ylabel_handle = get(gca, 'ylabel');
current_position = get(ylabel_handle, 'position');
set(ylabel_handle, 'position', [current_position(1)-300, current_position(2), current_position(3)]);

print -depsc2 MomentumMonthlyUniform


%% WEIGHTED

%% UNIFORM
% Set custom color order as default
set(groot, 'defaultAxesColorOrder', color_order);

fig = gcf;
set(fig, 'Units', 'normalized', 'Position', [0.15, 0.15, 0.9, 0.9]); % adjust the numbers to scale down the figure
ax = gca;
set(ax, 'Units', 'normalized', 'Position', [0.15, 0.15, 0.8, 0.75]); % adjust the numbers to scale down the axis


plot(time, PortfolioValueMonthlyW,'LineWidth',1);
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
ylim([50000 100000000])
set(gca, 'yticklabel', sprintfc('%gM SEK', get(gca, 'ytick')/1e6), 'FontSize', 14 )
% Add a y-axis label
ylabel('Portfolio value (SEK)', 'FontSize', 14)
ylabel_handle = get(gca, 'ylabel');
current_position = get(ylabel_handle, 'position');
set(ylabel_handle, 'position', [current_position(1)-300, current_position(2), current_position(3)]);

print -depsc2 MomentumMonthlyWeighted

%% Seasonality

%% Normal Scale

% Set custom color order as default
set(groot, 'defaultAxesColorOrder', color_order);

fig = gcf;
set(fig, 'Units', 'normalized', 'Position', [0.15, 0.15, 0.9, 0.9]); % adjust the numbers to scale down the figure
ax = gca;
set(ax, 'Units', 'normalized', 'Position', [0.15, 0.15, 0.8, 0.75]); % adjust the numbers to scale down the axis


plot(time, MonthlySeasonalityWeighted,'LineWidth',1);
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
ylim([50000 27500000])
set(gca, 'yticklabel', sprintfc('%gM SEK', get(gca, 'ytick')/1e6), 'FontSize', 14 )
% Add a y-axis label
ylabel('Portfolio value (SEK)')
ylabel_handle = get(gca, 'ylabel');
current_position = get(ylabel_handle, 'position');
set(ylabel_handle, 'position', [current_position(1)-300, current_position(2), current_position(3)]);

print -depsc2 MomentumSeasonalityWeighted

%% UNIFORM

% Set custom color order as default
set(groot, 'defaultAxesColorOrder', color_order);

fig = gcf;
set(fig, 'Units', 'normalized', 'Position', [0.15, 0.15, 0.9, 0.9]); % adjust the numbers to scale down the figure
ax = gca;
set(ax, 'Units', 'normalized', 'Position', [0.15, 0.15, 0.8, 0.75]); % adjust the numbers to scale down the axis


plot(time, MonthlySeasonalityUniform,'LineWidth',1);
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

legend('Uniform model (10 stocks)', 'Uniform model (20 stocks)', 'Uniform model (30 stocks)', 'Uniform model (40 stocks)', 'Uniform model (50 stocks)', 'Uniform model (100 stocks)', 'Uniform model (150 stocks)', 'OMXS30', 'OMXSPI')
legend('Location', 'South', 'Orientation', 'Horizontal','NumColumns', 5, 'Units', 'normalized', 'Position', [0.48 0 0.1 0.1], 'Box', 'off', 'FontSize', 16)

ytickformat('%,.0f SEK')
ylim([50000 18000000])
set(gca, 'yticklabel', sprintfc('%gM SEK', get(gca, 'ytick')/1e6), 'FontSize', 14 )
% Add a y-axis label
ylabel('Portfolio value (SEK)')
ylabel_handle = get(gca, 'ylabel');
current_position = get(ylabel_handle, 'position');
set(ylabel_handle, 'position', [current_position(1)-300, current_position(2), current_position(3)]);

print -depsc2 MomentumSeasonalityUniform

%% Log Scale

%% Weighted

% Set custom color order as default
set(groot, 'defaultAxesColorOrder', color_order);

fig = gcf;
set(fig, 'Units', 'normalized', 'Position', [0.15, 0.15, 0.9, 0.9]); % adjust the numbers to scale down the figure
ax = gca;
set(ax, 'Units', 'normalized', 'Position', [0.15, 0.15, 0.8, 0.75]); % adjust the numbers to scale down the axis


semilogy(time, MonthlySeasonalityWeighted,'LineWidth',1);
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

legend('Weighted model (10 stocks)', 'Weighted model (20 stocks)', 'Weighted model (30 stocks)', 'Weighted model (40 stocks)', 'Weighted model (50 stocks)', 'Weighted model (100 stocks)', 'Weighted model (150 stocks)', 'OMXS30', 'OMXSPI')
legend('Location', 'South', 'Orientation', 'Horizontal','NumColumns', 5, 'Units', 'normalized', 'Position', [0.48 0 0.1 0.1], 'Box', 'off', 'FontSize', 16)

set(gca, 'yticklabel', sprintfc('%gM SEK', get(gca, 'ytick')/1e6), 'FontSize', 14 )
% Add a y-axis label
yticks = logspace(log10(min(MonthlySeasonalityWeighted(1,:))), log10(max(MonthlySeasonalityWeighted(1,:))), 10);
yticklabels = arrayfun(@(x) sprintf('%.0f%%', x), ...
                        (yticks - min(MonthlySeasonalityWeighted(1,:))) ./ min(MonthlySeasonalityWeighted(1,:)) * 100, ...
                        'UniformOutput', false);
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels);
ylabel('Portfolio Value Increase (%)', 'FontSize', 14)

print -depsc2 MomentumSeasonalityWeightedLogScale


%% Uniform

% Set custom color order as default
set(groot, 'defaultAxesColorOrder', color_order);

fig = gcf;
set(fig, 'Units', 'normalized', 'Position', [0.15, 0.15, 0.9, 0.9]); % adjust the numbers to scale down the figure
ax = gca;
set(ax, 'Units', 'normalized', 'Position', [0.15, 0.15, 0.8, 0.75]); % adjust the numbers to scale down the axis


semilogy(time, MonthlySeasonalityUniform,'LineWidth',1);
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
yticks = logspace(log10(min(PortfolioValueMonthlyU(1,:))), log10(max(PortfolioValueMonthlyU(1,:))), 10);
yticklabels = arrayfun(@(x) sprintf('%.0f%%', x), ...
                        (yticks - min(PortfolioValueMonthlyU(1,:))) ./ min(PortfolioValueMonthlyU(1,:)) * 100, ...
                        'UniformOutput', false);
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels);
ylabel('Portfolio Value Increase (%)', 'FontSize', 14)

print -depsc2 MomentumSeasonalityUniformLogScale
%%

