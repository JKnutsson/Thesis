function [time,PortValueWeighted, PortValueUniform, PortValueOMX] = DividedYield(priceData, DivYield, dateEntry, DivDate, OMX, PortSize)


%% Index the end of each month in the price data and EV/EBITDA variable.
% We are buying/selling at the end of each month so instead of looping through all prices until we get to X date
% we index the position of the end of each month in a vector.

dat1            = datevec(dateEntry);              
[~,idxPrice]    = unique(dat1(:,1:2),'rows','last');   % idxPrice indexes the row coresponding to the end of each month in priceData
DateIndex       = dateEntry(idxPrice);                 % idxPrice in 'dd-mm-yyyy' format


dat2           = datevec(DivDate);              
[~,idxDiv]     = unique(dat2(:,1:2),'rows','last'); % idxDiv indexes the row coresponding to the end of each month in DivData
DateIndexDiv   = DivDate(idxDiv);                   % idxDiv in 'dd-mm-yyyy' format

dateEntryOMX   = OMX.RowLabels;
dat3           = datevec(dateEntryOMX);
[~,idxOMX]     = unique(dat3(:,1:2),'rows','last'); % idx2 indexes the row coresponding to the end of each month in priceData 
DateIndexOMX   = dateEntryOMX(idxOMX);
OMX            = table2array(OMX( : , 2:end));
[~, width]     = size(OMX);


clear dat1 dat2 dat3

%% Find the X largest dividend yield stocks and their position


% For each row in the table "idxDiv", find the X largest dividend yield stocks 
% and their positions within the array "DivYield".
for i = 1:height(idxDiv)   
    % largestDivYield{i} contains the X largest dividend yields for row i, 
    % while Position{i} contains the corresponding positions of these yields
    % within the row.
    [largestDivYield{i},Position{i}] = maxk(DivYield(idxDiv(i), : ),PortSize);
end

% Convert the cell arrays into matrices for ease of use.
largestDivYield = largestDivYield.';
largestDivYield = cell2mat(largestDivYield);

Position = Position.';
Position = cell2mat(Position);


%% Weights

% For each row in "largestDivYield", sum up all the dividend yields.
for i = 1:height(largestDivYield)
    SumDiv(i,1) = sum(largestDivYield(i,:));
end

% For each row in "largestDivYield", calculate the weight of each stock by 
% dividing its dividend yield by the sum of all the dividend yields.
for i = 1:height(largestDivYield) 
    
    for k = 1:PortSize
        Weights(i,k) = largestDivYield(i,k) / SumDiv(i,1);
    end
end




%% Portfolio-backtesting

% Initialise portfolio values for each simulation.
PortValueWeighted(1) = 100000;      % Portfolio value - weighted strategy
PortValueUniform(1) = 100000;      % Portfolio value - uniform strategy
m  = 0;                     % Variable used to adjust transaction time-stamp
ok = 0;                     % Flag variable

% Loop through all rows in "largestDivYield" (except the last one).
for i = 1:height(largestDivYield)-1
    
    % Loop through each stock in the portfolio.
    for k = 1:PortSize
        
            if Weights(i,k) > 0
        
            % While there are no price data available, adjust the time-stamp 
            % until price data is found. If too many adjustments have been
            % made, flag an error.
            while priceData( idxPrice(i,1)+m , Position(i,k) )  == 0

                if m < 3 
                    m = m+1;
                    Error = 87
                elseif ok == 1
                    m = -3;
                else
                    m = -1;
                    ok = 1;
                    Error = 93
                 end
            end

            % Calculate the number of stocks to buy based on the weighting 
            % for this particular stock.
            NumStocksWeighted(i,k) = floor( ( PortValueWeighted(i) * Weights(i,k) ) / priceData( idxPrice(i)+m, Position(i,k) ) );

            % Calculate the number of stocks to buy under a uniform investment 
            % strategy.
            NumStocksUniform(i,k)  = floor( ( PortValueUniform(i) * (1/PortSize) ) / priceData( idxPrice(i)+m, Position(i,k) ) );

            % Calculate the cost of buying the weighted number of stocks.
            BuyWeighted(i,k)       = NumStocksWeighted(i,k) * priceData( idxPrice(i), Position(i,k) );

            % Calculate the cost of buying an equal number of each stock.
            BuyUniform(i,k)        = NumStocksUniform(i,k) * priceData( idxPrice(i), Position(i,k) );

            m = 0;     % Reset time-stamp adjustment variable.
            end
    end

     
    
        % Calculate remaining cash after buying stocks in weighted and uniform
        % strategies.
        CashWeighted(i) = PortValueWeighted(i) - sum(BuyWeighted(i, : ));

        CashUniform(i)  = PortValueUniform(i) - sum(BuyUniform(i, : ));
        
    
        % Loop through each stock in the portfolio again (on the next day).
        for k = 1:PortSize     
             if Weights(i,k) > 0 
            % While there are no price data available, adjust the time-stamp 
            % until price data is found. If too many adjustments have been made,
            % flag an error.
            while priceData( idxPrice(i+1,1)+m , Position(i,k) )  == 0

                if m < 3 
                    m = m+1;
                    Error = 130
                elseif ok == 1
                    m = -3;
                else
                    m = -1;
                    ok = 1;
                    Error = 136;
                 end
            end

            % Calculate the revenue from selling the weighted number of stocks.
            SellWeighted(i,k) = NumStocksWeighted(i,k) * priceData( idxPrice(i+1)+m, Position(i,k) );

            % Calculate the revenue from selling an equal number of each stock.
            SellUniform(i,k)  = NumStocksUniform(i,k) * priceData( idxPrice(i+1)+m, Position(i,k) );

            m = 0;     % Reset time-stamp adjustment variable.
             end
         end  

        % Calculate the total portfolio value after selling stocks.
        PortValueWeighted(i+1) = CashWeighted(i) + sum(SellWeighted(i, :)) ;
        PortValueUniform(i+1) = CashUniform(i) + sum(SellUniform(i, :)) ;

    
end

%%
clear ok m


% The following loop will go through each stock in the OMX index and
% initialize a portfolio for it

for j = 1:width
    
    % Calculate the number of shares that can be bought with an initial 
    % investment of 100000 for the j-th stock, rounded down to the nearest integer.
    numberOfOMX(1,j)  = floor( 100000 / OMX( idxOMX(1) , j) );
    
    % Calculate the initial value of the portfolio by multiplying the 
    % number of shares purchased and the initial price of the share for
    % the j-th stock
    PortValueOMX(1,j) = numberOfOMX(1,j) * OMX( idxOMX(1) , j);
    
    % Calculate the remaining cash in the portfolio after buying the
    % j-th stock. This is the starting cash of 100000 minus the initial
    % value of the portfolio for the j-th stock.
    CashOMX(1,j)      = 100000 - PortValueOMX(1,j);
    
end

%%
% The following loop will go through each row in the idxOMX matrix,
% which contains the dates for which we have data on the OMX index. For
% each date, the loop will calculate the value of the portfolio.

for i = 1:height(idxOMX)-1

    % Nested loop to go through each stock in the OMX index and update
    % its value in the portfolio.
    for j = 1:width
        
       % Calculate the new value of the j-th stock in the portfolio by
       % multiplying the number of shares held and its current price.
       % Update the PortValueOMX matrix with this new value.
       PortValueOMX(i+1,j) = numberOfOMX(1,j) * OMX( idxOMX(i) , j);
       
    end
end

%%
% This loop updates the final row of the PortValueOMX matrix with the
% total value of the portfolio at the end of the investment period.
% It does this by adding the remaining cash in the portfolio (stored in
% the CashOMX matrix) to the value of the stocks in the final row of the
% PortValueOMX matrix.

for j = 1:width

    PortValueOMX(end,j) = PortValueOMX(end,j)+ CashOMX(1,j);

end


% This line sets any infinite values in the PortValueOMX matrix to zero.
PortValueOMX(isinf(PortValueOMX)) = 0;

% This line sets any NaN (Not a Number) values in the PortValueOMX matrix to zero.
PortValueOMX(isnan(PortValueOMX)) = 0;

% The following lines extract the OMXS30 and OMXSPI portfolios from the
% PortValueOMX matrix.



% Time is the X-axis and since it's based on the DateIndex which is the end
% last date of the month that has a price during the period I add 2 to it
% so that the X-axis start at 01-01-xxxx in the graphs
time = datenum(DateIndex);
time = time + 2;
time(49) = time(49)+2;
time(109) = time(109)+1;
time(121) = time(121)+2;
end

