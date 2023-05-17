clc, clear all, format long g
%% EV/EBIDTA portfolio

%% Load data set from excel.
% Data.Properties.VariableNames(2);

Data           = readtable('Data2.xlsx','Sheet','2001-2022');
EVEBITDASheet  = readtable('EVEBITDAWEEKLY.xlsx','Sheet','PivotEVTest');


OMX     = readtable('Data.xlsx','Sheet','Pivot OMX');
RF      = readtable('RfSE10Y.xlsx');

clc
%% Sets of data.

Data(5528:end, :) = [];

% Mean risk free rate of swedish 10 year bonds
RF              = mean(table2array(RF(109:372,2)))/100;
% Clean up of OMXS data to fit the time period of 2001-2022, This also removes the "OMXS's" that have no data during our selected time period.
OMX(:,7)        = [];
OMX(:,3:5)      = [];
OMX(1:2252, : ) = [];
OMX(end, : )    = [];
OMX             = rmmissing(OMX);

dateEntryOMX   = OMX.RowLabels;
dat3           = datevec(dateEntryOMX);
[~,idxOMX]     = unique(dat3(:,1:2),'rows','last'); % idx2 indexes the row coresponding to the end of each month in priceData 
DateIndexOMX   = dateEntryOMX(idxOMX);
OMX            = table2array(OMX( : , 2:end));

clear dat3
%%

dateEntry    = Data.RowLabels(:,:);        % Contain all dates from the first period of data
priceEntry   = Data(:,2:end);            % Contain all prices from the first period of data

EVEBITDADate = EVEBITDASheet.RowLabels(:,:);
EVEBITDAData = EVEBITDASheet( : ,2:end);


% define two sets of column variable names
set1 = EVEBITDAData.Properties.VariableNames;
set2 = priceEntry.Properties.VariableNames;

% find the variable names(i.e assets) that are in set1 but not in set2 and the other
% way around.
missingVars  = setdiff(set1, set2);
missingVars2 = setdiff(set2, set1);

% remove the columns in priceEntry that correspond to the missing variable names
priceEntry   = removevars(priceEntry, missingVars2);
EVEBITDAData = removevars(EVEBITDAData, missingVars);

% Converts prices to the "double" format
priceData       = table2array( priceEntry );  
EVEBITDAData    = table2array(EVEBITDAData);

% Since some of the data from Eikon is incomplete and have varying "periods" of missing data for either/both price or ev/ebitda and are removed manually
% if they can't "pass" the while loop in the Backtesting section. The while loop checks if the price at a specific date is 0 
% if it is then it checks if there is a availabe price in +2/-2 steps from the current position that can be used.
% The following code removes the stocks or periods which could not pass the while loop.
%%

EVEBITDAData(101:end, 17)   = 0;

EVEBITDAData(403, 20)   = 0;

EVEBITDAData(20, 23)   = 0;
EVEBITDAData(120:121, 23)   = 0;

EVEBITDAData(174:176, 28)   = 0;

EVEBITDAData(80:82, 49)   = 0;
EVEBITDAData(105:106, 49)   = 0;

EVEBITDAData(149:end, 98)   = 0;

EVEBITDAData(153:154, 107)    = 0;

EVEBITDAData(115, 115)  = 0;
EVEBITDAData(192:193, 115)  = 0;
EVEBITDAData(232:250, 115)  = 0;

EVEBITDAData(177:180, 138)  = 0;

EVEBITDAData(251, 251)  = 0;
EVEBITDAData(306:320, 251)  = 0;

EVEBITDAData(9, 242)  = 0;
EVEBITDAData(44, 242)  = 0;
EVEBITDAData(70:75, 242)  = 0;
EVEBITDAData(80:86, 242)  = 0;

EVEBITDAData(86:88, 256)  = 0;

clear set1 set2 missingVars missingVars2

%%

PortSize = [10,20,30,40,50,100,150];

%%

for i = 1:width(PortSize)
   [time,PortValueWeighted{i,1}, PortValueUniform{i,1}] = EVEBITDA(priceData, EVEBITDAData, dateEntry, EVEBITDADate, PortSize(i));
end

%%

PortValueWeighted = cell2mat(PortValueWeighted);
PortValueUniform  = cell2mat(PortValueUniform);

%% OMX "Portfolio"
% Using the same logic as the backtesting, a portfolio for the OMXS is created. This allows the value of the portfolio
% to update at the same intervals as the backtesting so that it can be plotted on the same graph later.


for j = 1:width(OMX)
    
    numberOfOMX(1,j)  = floor( 100000 / OMX( idxOMX(1) , j) );

    PortValueOMX(1,j) = numberOfOMX(1,j) * OMX( idxOMX(1) , j);

    CashOMX(1,j)      = 100000 - PortValueOMX(1,j);
    
end
%%
for i = 1:height(idxOMX)-1

    for j = 1:width(OMX)
        
       PortValueOMX(i+1,j) = numberOfOMX(1,j) * OMX( idxOMX(i) , j);
       
    end
end
%%
for j = 1:2

    PortValueOMX(end,j) = PortValueOMX(end,j)+ CashOMX(1,j)  ;

end

PortValueOMX(isinf(PortValueOMX)) = 0;
PortValueOMX(isnan(PortValueOMX)) = 0;

OMXS30     = PortValueOMX(:,1);
OMXSPI     = PortValueOMX(:,2);
%%
clear PortValueOMX PortValueOMX numberOfOMX CashOMX
clc

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
    AnnualReturnW(i) = ((1 + ReturnWM(i))^(1 / 22) - 1 )*100;

    % Calculate the annual return for the Uniform portfolio, assuming that there are 14 years.
    AnnualReturnU(i) = ((1 + ReturnUM(i))^(1 / 22) - 1 )*100;

    % Calculate the annual return for the OMXS30 index, assuming that there are 14  years.
    AnnualReturnOMXS30 = ((1 + ReturnOMXS30)^(1 / 22) - 1 )*100;

    % Calculate the annual return for the OMXSPI index, assuming that there are 14 years.
    AnnualReturnOMXSPI = ((1 + ReturnUMOMXSPI)^(1 / 22) - 1 )*100;


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


plot(time, PortValueUniform,'LineWidth',1);
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
ylim([50000 34000000000])
set(gca, 'yticklabel', sprintfc('%gM SEK', get(gca, 'ytick')/1e6), 'FontSize', 14 )
% Add a y-axis label
ylabel('Portfolio value (SEK)')
ylabel_handle = get(gca, 'ylabel');
current_position = get(ylabel_handle, 'position');
set(ylabel_handle, 'position', [current_position(1)-300, current_position(2), current_position(3)]);


print -depsc2 EVEBITDA_Uniform_Normal
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
yticks = logspace(log10(min(PortValueUniform(1,:))), log10(max(PortValueUniform(1,:))), 10);
yticklabels = arrayfun(@(x) sprintf('%.0f%%', x), ...
                        (yticks - min(PortValueUniform(1,:))) ./ min(PortValueUniform(1,:)) * 100, ...
                        'UniformOutput', false);
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels);
ylabel('Portfolio Value Increase (%)', 'FontSize', 14)


print -depsc2 EVEBITDA_Uniform_Log

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
ylim([50000 360000000000])
set(gca, 'yticklabel', sprintfc('%gM SEK', get(gca, 'ytick')/1e6), 'FontSize', 14 )
% Add a y-axis label
ylabel('Portfolio value (SEK)')
ylabel_handle = get(gca, 'ylabel');
current_position = get(ylabel_handle, 'position');
set(ylabel_handle, 'position', [current_position(1)-300, current_position(2), current_position(3)]);


print -depsc2 EVEBITDA_Weighted_Normal


%% Weighted Log


% Set custom color order as default
set(groot, 'defaultAxesColorOrder', color_order);

fig = gcf;
set(fig, 'Units', 'normalized', 'Position', [0.15, 0.15, 0.9, 0.9]); % adjust the numbers to scale down the figure
ax = gca;
set(ax, 'Units', 'normalized', 'Position', [0.15, 0.15, 0.8, 0.75]); % adjust the numbers to scale down the axis
semilogy(time, PortValueWeighted, 'LineWidth', 1);
hold on
semilogy(time, OMXS30, 'LineWidth', 1)
semilogy(time, OMXSPI, 'LineWidth', 1)
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
                        (yticks - min(PortValueWeighted(i,:))) ./ min(PortValueWeighted(i,:)) * 100, ...
                        'UniformOutput', false);
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels);
ylabel('Portfolio Value Increase (%)', 'FontSize', 14)


print -depsc2 EVEBITDA_Weighted_Log



