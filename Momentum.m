function [time,PortfolioValueW, PortfolioValueU, MonthlySeasonalityWeighted, MonthlySeasonalityUniform ] = Momentum(priceData1, priceData2, dateEntry1 , dateEntry2, priceEntry1, priceEntry2,PortSize)


% Convert date strings to date vectors
dat1        = datevec(dateEntry1);  

% Find the unique rows corresponding to each month end in dat1 and 
% get the indexes of these rows
[~,idx1]    = unique(dat1(:,1:2),'rows','last'); % idx1 indexes the row coresponding to the end of each month in priceData
DateIndex1  = dateEntry1(idx1);

% Get the number of columns in priceData1
[~, width1] = size(priceData1);

% Repeat the same process for the second set of date entries
dat2        = datevec(dateEntry2);
[~,idx2]    = unique(dat2(:,1:2),'rows','last'); % idx2 indexes the row coresponding to the end of each month in priceData                 
[~, width2] = size(priceData2);
DateIndex2  = dateEntry2(idx2);


%% Find the momentum of each stock at the end of each month which is then stored in the "Momentum 1 or 2" matrices

% Initialize a zero matrix for Momentum1
Momentum1 = zeros(height(idx1), width1);

% Loop through each month in idx1 starting from 12 (skipping first year)
for l = 12:height(idx1)     % l is the index number for each month.
    
    % Loop through each column (stock) in priceData1
    for k = 1:width1        % k keeps track of the colums (stocks).         
    
        % Check if starting and ending day price is not equal to 0.
        if priceData1( idx1(l) ,k) ~= 0 && priceData1( idx1(l-11) ,k) ~= 0 
            
            % Calculate momentum using a 1-year window excluding the last month.
            Momentum1(l,k) = ((priceData1( idx1(l-11) ,k) - priceData1(1,k)) / priceData1(idx1(l-11),k))+1;
    
            % Loop through the previous 10 months using i as the index
            for i = 10:-1:1            
                
                m = 0; 
                
                % Check if data point is missing
                if priceData1( idx1(l-i) ,k) == 0 
                    % Change row position by 1 incase one data point in the set is missing.
                    m = 1;              
                end
                
                % Calculate cumulative returns for momentum calculation
                Momentum1(l,k) = Momentum1(l,k) * (((priceData1( idx1(l-i)+m ,k) - priceData1( idx1(l-(i+1))+m ,k)) / priceData1( idx1(l-i)+m,k))+1);
                                
            end 
            
            % Subtract 1 from cumulative returns to calculate generic momentum
            Momentum1(l,k) = Momentum1(l,k)-1; % Cumulative returns (generic momentum) for each stock
        end  
    end
end


%% Same as above just for the second time period.

Momentum2 = zeros(132 , 281); 

for l = 12:height(idx2)           % l is the index number for each month.
    
    for k = 1:width2        % k keeps track of the colums (stocks).         
    
        if priceData2( idx2(l) ,k) ~= 0 && priceData2( idx2(l-11) ,k) ~= 0  % Filter periods where the starting and ending day price is 0.
        
                                Momentum2(l,k) = ( (priceData2( idx2(l-11) ,k) - priceData2( 1 ,k) ) / priceData2( idx2(l-11) ,k))+1;
    
            for i = 10:-1:1            % Calculate 1 years momentum skipping the last month.
                
                                m = 0;  
                                
                    if priceData2( idx2(l-i) ,k) == 0 % Change row position by 1 incase one data point in the set is missing.
                    
                                m = 1;              
                   
                    end
                    
                                Momentum2(l,k) = Momentum2(l,k) * (((priceData2( idx2(l-i)+m ,k) - priceData2( idx2(l-(i+1))+m ,k)) / priceData2( idx2(l-i)+m ,k))+1);
                                
            end 
                                Momentum2(l,k) = Momentum2(l,k)-1; % Cumulative returns (generic momentum) for each stocks
        end  
    end
end
%% Clean up inf / NaN's.


% Replace all infinities in Momentum1 with 0
Momentum1(isinf(Momentum1)) = 0; % 663

% Replace all NaNs in Momentum1 with 0
Momentum1(isnan(Momentum1)) = 0; % 445


% Replace all infinities in Momentum2 with 0
Momentum2(isinf(Momentum2)) = 0; % 131

% Replace all NaNs in Momentum2 with 0
Momentum2(isnan(Momentum2)) = 0; % 32





%%
positive1    = zeros(height(idx1) , width1);
negative1    = zeros(height(idx1) , width1);
positive2    = zeros(height(idx1) , width1);
negative2    = zeros(height(idx1) , width1);
ID1          = zeros(height(idx1) , width1); % Information discreteness (Frog-in-the-pan theory)

positive3    = zeros(height(idx2) , width2);
negative3    = zeros(height(idx2) , width2);
positive4    = zeros(height(idx2) , width2);
negative4    = zeros(height(idx2) , width2);
ID2          = zeros(height(idx2) , width2); % Information discreteness (Frog-in-the-pan theory)


% This loop iterates over each month (l is the index number for each month).
for l = 12:height(idx1)
    
    % This loop iterates over each stock (k keeps track of the columns/stocks).
    for k = 1:width1              
        
        % This loop calculates the number of positive and negative trading days 1 trading year (251 days) back.
        for p = 251:-1:2
                                 
            % If we're in December, set m equal to 4, otherwise it's zero.
            if l == 12
                 m = 4;
            end
            
            % Check if momentum isn't zero for both previous and current trading days.
            if priceData1( idx1(l)-p+m , k) ~= 0 && priceData1( idx1(l)-p+1+m , k) ~= 0
                
                % If the price on the first day is greater than the second day, increment the number of positive days for that stock.
                if priceData1( idx1(l)-p+m , k) > priceData1( idx1(l)-p+1+m , k)
                    positive1(l,k) = positive1(l,k) + 1;
                      
                % If the price on the first day is less than the second day, increment the number of negative days for that stock.
                elseif priceData1( idx1(l)-p+m , k) < priceData1( idx1(l)-p+1+m , k)
                    negative1(l,k) = negative1(l,k) + 1;
                end                               
            end
            
            % Set m back to zero for the next iteration.
            m = 0;
        end        
    end
end

% Same as above but for the next time period.

for l = 12:height(idx2)                  % l is the index number for each month.
    
    for k = 1:width2               % k keeps track of the colums (stocks).
        
        for p = 251:-1:2        % Cacluates the number of positive and negative trading days in a 1 trading year 
                                % starting in the end of each month 
             if l == 12
                 m = 2;
             end
            if priceData2( idx2(l)-p+m , k) ~= 0 && priceData2( idx2(l)-p+1+m , k) ~= 0 % ADD IF THE MOMENTUM ISN'T 0???
                
                        
                      if priceData2( idx2(l)-p+m , k) > priceData2( idx2(l)-p+1+m , k) % Adds up the number of positive days

                                positive3(l,k) = positive3(l,k) + 1;
                      end
                      
                      if priceData2( idx2(l)-p +m, k) < priceData2( idx2(l)-p+1+m , k) % Adds up the number of negative days

                                negative3(l,k) = negative3(l,k) + 1;      
                      end
                               
            end
            m = 0;
        end
        
           
    end
end

  % This loop iterates over each month (l is the index number for each month).
for l = 12:height(idx1)
    
    % This loop iterates over each stock (k keeps track of the columns/stocks).
    for k = 1:width1 
        
        % Calculate the percentage of positive days for each stock.
        positive2(l,k) = positive1(l,k) / (positive1(l,k) + negative1(l,k));
        
        % Calculate the percentage of negative days for each stock.
        negative2(l,k) = negative1(l,k) / (positive1(l,k) + negative1(l,k));

        % Create a filter based on the "Frog-in-the-pan" theory.
        ID1(l,k) = sign(Momentum1(l,k)) * (negative2(l,k) - positive2(l,k));            
    end
end

% This loop iterates over each month (l is the index number for each month).
for l = 12:height(idx2)
    
    % This loop iterates over each stock (k keeps track of the columns/stocks).
    for k = 1:width2 
        
        % Calculate the percentage of positive days for each stock.
        positive4(l,k) = positive3(l,k) / (positive3(l,k) + negative3(l,k));
        
        % Calculate the percentage of negative days for each stock.
        negative4(l,k) = negative3(l,k) / (positive3(l,k) + negative3(l,k));
        
        % Create a filter based on the "Frog-in-the-pan" theory.
        ID2(l,k) = sign(Momentum2(l,k)) * (negative4(l,k) - positive4(l,k));            
    end
end

 

%% Replaces inf / NaN's with zeros.

ID1(isinf(ID1)) = 0;
ID1(isnan(ID1)) = 0;

ID2(isinf(ID2)) = 0;
ID2(isnan(ID2)) = 0;


% This line creates an empty matrix of zeros with dimensions 265x217.
momentumQualityFiltered1 = zeros(265, 217);

% This line creates an empty matrix of zeros with dimensions 132x281.
momentumQualityFiltered2 = zeros(132, 281);

% This loop iterates over each month (l is the index number for each month).
for l = 12:height(idx1)
    
    % This loop iterates over each stock (k keeps track of the columns/stocks).
    for k = 1:width1 
        
        % Check if the stock has positive momentum and high quality.
        if Momentum1(l,k) > 0 && ID1(l,k) < 0  
                           
            % If it does, add its Momentum value to the filtered matrix.
            momentumQualityFiltered1(l,k) = Momentum1(l,k);
                       
        end  
    end
end

% This loop iterates over each month (l is the index number for each month).
for l = 12:height(idx2)
    
    % This loop iterates over each stock (k keeps track of the columns/stocks).
    for k = 1:width2 
        
        % Check if the stock has positive momentum and high quality.
        if Momentum2(l,k) > 0 && ID2(l,k) < 0  
                           
            % If it does, add its Momentum value to the filtered matrix.
            momentumQualityFiltered2(l,k) = Momentum2(l,k);
                       
        end  
    end
end

%%
% Create a counter variable x and a matrix to store row indices for each rebalancing.
x = 1;
monthRowIndex1 = zeros(45,1);

% Find the PortSize(i) stocks with the highest momentum during the first rebalancing period (in February).
largestMomentum1{x} = maxk(momentumQualityFiltered1(12,: ), PortSize);
monthRowIndex1(x) = 12;
x = x + 1;
%%
% Loop over the remaining 21 years of data, and for each year, find the PortSize(i) stocks with the highest momentum during each rebalancing period (in May, August, and November).
for y = 1:21
    largestMomentum1{x} = maxk(momentumQualityFiltered1(2+12*y,: ), PortSize);
    monthRowIndex1(x) = 2 + 12*y;
    x = x + 1;
    
    largestMomentum1{x} = maxk(momentumQualityFiltered1(5+12*y,: ), PortSize);
    monthRowIndex1(x) = 5 + 12*y;
    x = x + 1;   
    
    largestMomentum1{x} = maxk(momentumQualityFiltered1(8+12*y,: ), PortSize);
    monthRowIndex1(x) = 8 + 12*y;
    x = x + 1;
    
    largestMomentum1{x} = maxk(momentumQualityFiltered1(11+12*y,: ), PortSize);
    monthRowIndex1(x) = 11 + 12*y;
    x = x + 1;
end

% Transpose the arrays and convert from cell to matrix.
largestMomentum1 = largestMomentum1.';
largestMomentum1 = cell2mat(largestMomentum1);
%%
% Create a counter variable x and a matrix to store row indices for each rebalancing.
x = 1;
monthRowIndex2 = zeros(41,1);

% Find the PortSize(i) stocks with the highest momentum during the first rebalancing period (in February).
largestMomentum2{x} = maxk(momentumQualityFiltered2(12,: ), PortSize);
monthRowIndex2(x) = 12;
x = x + 1;

% Loop over the remaining 10 years of data, and for each year, find the PortSize(i) stocks with the highest momentum during each rebalancing period (in May, August, and November).
for y = 1:10
    largestMomentum2{x} = maxk(momentumQualityFiltered2(2+12*y,: ), PortSize);
    monthRowIndex2(x) = 2 + 12*y;
    x = x + 1;
    
    largestMomentum2{x} = maxk(momentumQualityFiltered2(5+12*y,: ), PortSize);
    monthRowIndex2(x) = 5 + 12*y;
    x = x + 1;   
    
    largestMomentum2{x} = maxk(momentumQualityFiltered2(8+12*y,: ), PortSize);
    monthRowIndex2(x) = 8 + 12*y;
    x = x + 1;
    
    largestMomentum2{x} = maxk(momentumQualityFiltered2(11+12*y,: ), PortSize);
    monthRowIndex2(x) = 11 + 12*y;
    x = x + 1;
end

% Transpose the arrays and convert from cell to matrix.
monthRowIndex2 = monthRowIndex2.';
largestMomentum2 = largestMomentum2.';
largestMomentum2 = cell2mat(largestMomentum2);

%%

% Loop over each row of the largestMomentum1 matrix to compute the sum of momentums for each stock.
for i = 1:85
    sumMomentum1(i,1) = sum(largestMomentum1(i,:));
end

% Loop over each row of the largestMomentum2 matrix to compute the sum of momentums for each stock.
for i = 1:41
    sumMomentum2(i,1) = sum(largestMomentum2(i,:));
end



% Loop over each row i of the largestMomentum1 matrix and compute the position of the k-th largest momentum stock.
for i = 1:85
    for k = 1:PortSize

        % Check if the k-th largest momentum stock has positive momentum and record its position.
        if largestMomentum1(i,k) > 0
            [~, c1(i,k)] = find(momentumQualityFiltered1(monthRowIndex1(i):monthRowIndex1(i),:) == largestMomentum1(i,k));
        end
        
    end
end

% Convert the month and position data into a table and then to an array.
positionMatrix1 = table(monthRowIndex1, c1);
positionMatrix1 = table2array(positionMatrix1);

% Repeat the above code but for largestMomentum2 data.
for i = 1:41
    for k = 1:PortSize

        if largestMomentum2(i,k) > 0
            [~, c2(i,k)] = find(momentumQualityFiltered2(12:132,:) == largestMomentum2(i,k));
        end
        
    end
end

% Convert the month and position data into a table and then to an array.
monthRowIndex2 = monthRowIndex2';
positionMatrix2 = table(monthRowIndex2, c2);
positionMatrix2 = table2array(positionMatrix2);



% Iterate over each row i of largestMomentum1 matrix and calculate the weight for each stock.
% Also calculate the number of non-zero weights in each row of the matrix.

UniformWeights1 = zeros(85,1);
UniformWeights2 = zeros(41,1);
Weights1        = zeros(85, PortSize);
Weights2        = zeros(41, PortSize);

for i = 1:85

    for j = 1:PortSize
         
        % Check if the j-th largest momentum stock has positive momentum and calculate its weight.
        if largestMomentum1(i,j) > 0 

            Weights1(i,j) = largestMomentum1(i,j) / sumMomentum1(i);
            UniformWeights1(i) = UniformWeights1(i) + 1;

        end
    end
end

% Repeat the above loop but for largestMomentum2 data.
for i = 1:41

    for j = 1:PortSize
         
        if largestMomentum2(i,j) > 0 

            Weights2(i,j) = largestMomentum2(i,j) / sumMomentum2(i);
            UniformWeights2(i) = UniformWeights2(i) + 1;

        end
    end
end

% delete the first row of the copied Weights and Position Matrices.
Weights2(1,:) = [];
positionMatrix2(1,:) = [];


%% Portfolio Backtesting

% Smart Rebalance: This rebalanced portfolio is based on seasonality and quality momentum stocks.
% The portfolio is rebalanced on the close of trading in February, May, August, and November.
% The PositionMatrix contains information about the rows and columns of the largest quality momentum stocks.
% The IDX variable is used to find the end of each month in the priceData matrix, specifically for February, May, August, and November months.


% This code initializes the portfolio value for two portfolios (PortfolioValueW and PortfolioValueU) in the beginning, based on an initial investment of $100,000.
% It then implements a smart rebalance approach to buying and selling stocks in each portfolio. 
% The code uses a loop that runs i times (i.e., once for each time period), and within each iteration, it does the following:
% For buying stocks:
%   - Loops through PortSize(i) different stock weights assuming PortSize(i) stocks passed our filter in the blocks above otherviwse it loops PortSize(i) times but "skips" the "stocks" with 0 momentum.
%   - Calculates the number of stocks to buy based on momentum and quality criteria, as well as the number of stocks to buy with uniform weights.
%   - Checks if the stock data is missing or partial by searching for the next available market day. If the data is missing or partial, it increments a counter variable to try again on the next day. Once the counter is reset to 0, the code calculates the total amount spent to buy stocks using both approaches.
%   - Finally, it calculates the cash available after buying the stocks for each portfolio.
% For selling stocks:
%   - Loops through the same PortSize(i) different stock weights.
%   - Calculates the amount received for selling the stocks using both Smart and Uniform Rebalance approaches. The counter variable is incremented if the stock data is missing or partial.
%   - The new portfolio value is calculated using the cash available after selling the stocks.
% The final result is the updated PortfolioValueW and PortfolioValueU matrices, containing the updated value of each portfolio at each time period.

% Initial Portfolio Value

PortfolioValueW(1) = 100000;
PortfolioValueU(1) = 100000;
ok = 0;

z = 1;

% Initialize a counter variable
m = 0;

% Loop for Buying Stocks using Smart Rebalance approach
for i = 1:85
    
    % Buy section
    for j = 1:PortSize

        % Check if the stock weight is greater than 0 for the given row
        if Weights1(i,j) > 0
            
            % Loop to check if the specific stock data is missing or partial
            m = 0;
        while priceData1(idx1(positionMatrix1(i,1))+m , positionMatrix1(i,j+1)) == 0
            % Loop through price data for that stock until valid price data is encountered
            if m < 3 
                % Adjust time stamp by one day (24 hours)
                m = m + 1;
            elseif ok == 1
                % If more than 3 adjustments have been made, set m value to -3 to flag an error
                m = -3;
            else
                % If only three adjustments were made, adjust timestamp by a larger amount
                % to find price data
                m = -1;
                ok = 1;
            end
        end

            % Calculate the number of stocks to buy based on momentum and quality
            numberOfStocks1W(i,j) = floor( ( Weights1(i,j) * PortfolioValueW(i) ) / priceData1( idx1(positionMatrix1(i,1))+m , positionMatrix1(i,j+1) ) );

            % Calculate the total amount spent to buy stocks
            BuySum1W(i,j)         = numberOfStocks1W(i,j) * priceData1( idx1(positionMatrix1(i,1))+m , positionMatrix1(i,j+1) );

            % Calculate the number of stocks to buy with uniform weights
            numberOfStocks1U(i,j) = floor( ( (1/UniformWeights1(i)) * PortfolioValueU(i) ) / priceData1( idx1(positionMatrix1(i,1))+m , positionMatrix1(i,j+1) ) );

            % Calculate the total amount spent to buy stocks with uniform weights
            BuySum1U(i,j)         = numberOfStocks1U(i,j) * priceData1( idx1(positionMatrix1(i,1))+m , positionMatrix1(i,j+1) );

            % Reset the counter variable to 0 for next iteration
            m = 0;

        end
    end

    % Calculate the cash available after buying the stocks in Weitghed and Uniform Rebalance
    Cash1W(i) = PortfolioValueW(i) - sum(BuySum1W(i,:));
    Cash1U(i) = PortfolioValueU(i) - sum(BuySum1U(i,:));
    
        
       % The following for-loop is used to keep track of the monthly price
        % changes of the portfolio.
    
        MonthlySeasonalityWeighted(1,z) = PortfolioValueW(i);    
        MonthlySeasonalityUniform(1,z)  = PortfolioValueU(i);
        
        z = z +1;
        
        if i <   85 && i > 1
            
            MonthUpdate = 2;
        else 
            MonthUpdate = 1;
        end
    for q = 1:MonthUpdate
        
        for j = 1:PortSize

            % Check if the stock weight is greater than 0 for the given row
            if Weights1(i,j) > 0

                % Loop to check if the specific stock data is missing or partial
               m = 0;
        while priceData1(idx1(positionMatrix1(i,1)+q)+m , positionMatrix1(i,j+1)) == 0
            % Loop through price data for that stock until valid price data is encountered
            if m < 3 
                % Adjust time stamp by one day (24 hours)
                m = m + 1;
            elseif ok == 1
                % If more than 3 adjustments have been made, set m value to -3 to flag an error
                m = -3;
            else
                % If only three adjustments were made, adjust timestamp by a larger amount
                % to find price data
                m = -1;
                ok = 1;
            end
        end
        ok = 0;
                            
              MonthlyValueWeighted(q,j) = numberOfStocks1W(i,j) * priceData1( idx1(positionMatrix1(i,1)+q)+m , positionMatrix1(i,j+1) ) ;
                
                
              MonthlyValueUniform(q,j)  = numberOfStocks1U(i,j) * priceData1( idx1(positionMatrix1(i,1)+q)+m , positionMatrix1(i,j+1) ) ;

            end
        end
        
        MonthlySeasonalityWeighted(1,z) = sum(MonthlyValueWeighted(q,:)) + Cash1W(i);
        
        MonthlySeasonalityUniform(1,z)  = sum(MonthlyValueUniform(q,:)) + Cash1U(i);

        z = z +1;
        
        MonthlyValueWeighted = [];
        MonthlyValueUniform  = []; 
    end
        
    
    

    % Sell section
    % Reset the counter variable to 0 for next iteration
    m = 0;
    
    % Check if it is not the last iteration
    if i < 85
         
        % Loop for selling stocks in Smart and Uniform Rebalance
        for j = 1:PortSize

            % Check if the stock weight is greater than 0 for the given row
            if Weights1(i,j) > 0
                 
                m = 0;
        while priceData1(idx1(positionMatrix1(i+1,1))+m , positionMatrix1(i,j+1)) == 0
            % Loop through price data for that stock until valid price data is encountered
            if m < 3 
                % Adjust time stamp by one day (24 hours)
                m = m + 1;
            elseif ok == 1
                % If more than 3 adjustments have been made, set m value to -3 to flag an error
                m = -3;
            else
                % If only three adjustments were made, adjust timestamp by a larger amount
                % to find price data
                m = -1;
                ok = 1;
            end
        end
            ok = 0;
                % Calculate the amount received for selling the stocks in Smart Rebalance
                CashSell1W(i,j) = numberOfStocks1W(i,j) * priceData1( idx1(positionMatrix1(i,1)+3)+m , positionMatrix1(i,j+1) ) ;

                % Calculate the amount received for selling the stocks in Uniform Rebalance
                CashSell1U(i,j) = numberOfStocks1U(i,j) * priceData1( idx1(positionMatrix1(i,1)+3)+m , positionMatrix1(i,j+1) ) ;

                % Reset the counter variable to 0 for next iteration
                m = 0;

            end
        end
    
        % Calculate the new Portfolio Value after selling the stocks in Smart and Uniform Rebalance
        PortfolioValueW(i+1) = Cash1W(i) + sum(CashSell1W(i,:));
        PortfolioValueU(i+1) = Cash1U(i) + sum(CashSell1U(i,:));

    end   
end

%% Bridge period 1 to 2
% This code updates the portfolio values to ensure that the second time period starts with the correct values.
% It achieves this by checking the positions of our stocks from Set 1 in Set 2.
% Since we buy our last portfolio in November and need it to "run" for 3 months,
% it is important to adjust the values correctly.
% The code ensures that the updated portfolio values accurately reflect the positions of the stocks in each set at the beginning of the second time period.

% Create tables to hold the stock symbols for the last portfolio

last_Row = positionMatrix1(end,2:end)';
non_Zero_Pos = last_Row(last_Row ~= 0);

set1 = table(priceEntry1.Properties.VariableNames(non_Zero_Pos)');
set2 = table(priceEntry2.Properties.VariableNames()');

numSmallElements = height(set1);
[~,Bridge] = ismember(set1,set2);





    % The following for-loop is used to keep track of the monthly price
    % changes of the portfolio.
   for q = 1:2     
    for j = 1:height(Bridge)

        % Check if the stock weight is greater than 0 for the given row
        if Weights1(i,j) > 0

            % Loop to check if the specific stock data is missing or partial
            while priceData2(idx2(positionMatrix2(q,1))+m , Bridge(j,1)) == 0
    % Loop through price data for that stock until valid price data is encountered
            if m < 3 
                % Adjust time stamp by one day (24 hours)
                m = m + 1;
            elseif ok == 1
                % If more than 3 adjustments have been made, set m value to -3 to flag an error
                m = -3;
            else
                % If only three adjustments were made, adjust timestamp by a larger amount
                % to find price data
                m = -1;
                ok = 1;
            end
            end

          MonthlyValueWeighted(1,j) = numberOfStocks1W(i,j) * priceData2( idx2(positionMatrix2(q,1))+m , Bridge(j,1) ) ;


          MonthlyValueUniform(1,j)  = numberOfStocks1U(i,j) * priceData2( idx2(positionMatrix2(q,1))+m , Bridge(j,1) ) ;

        end
    end

    MonthlySeasonalityWeighted(1,z) = sum(MonthlyValueWeighted(1,:)) + Cash1W(i);

    MonthlySeasonalityUniform(1,z)  = sum(MonthlyValueUniform(1,:)) + Cash1U(i);

    z = z +1;

    MonthlyValueWeighted = [];
    MonthlyValueUniform  = []; 
   end


% Loop through the columns of Weights1 matrix from 1 to PortSize(i)
for j = 1:height(Bridge)
    
    % Check if the element in Weights1 matrix is greater than zero
    if Weights1(i,j) > 0
        
        % Initialize a variable m to 0
        m = 0;
        
        % Loop until a non-zero stock price is found in priceData2 matrix
        while priceData2( idx2(positionMatrix2(1,1))+m , Bridge(j,1) ) == 0
            
              % Loop through price data for that stock until valid price data is encountered
            if m < 3 
                % Adjust time stamp by one day (24 hours)
                m = m + 1;
            elseif ok == 1
                % If more than 3 adjustments have been made, set m value to -3 to flag an error
                m = -3;
            else
                % If only three adjustments were made, adjust timestamp by a larger amount
                % to find price data
                m = -1;
                ok = 1;
            end
        end
        
        % Calculate the total value of the stocks sold for portfolio 1W
        CashSell1W(i,j) = numberOfStocks1W(i,j) * priceData2( idx2(positionMatrix2(2,1))+m , Bridge(j,1) );
        
        % Calculate the total value of the stocks sold for portfolio 1U
        CashSell1U(i,j) = numberOfStocks1U(i,j) * priceData2( idx2(positionMatrix2(2,1))+m , Bridge(j,1) );
        
        % Set m back to 0
        m = 0;
    end
end

% Calculate the portfolio value for portfolio 1W for the next period
PortfolioValueW(i+1) = Cash1W(i) + sum(CashSell1W(i,:));

% Calculate the portfolio value for portfolio 1U for the next period
PortfolioValueU(i+1) = Cash1U(i) + sum(CashSell1U(i,:));


%% Period 2 

% Loop through the rows of Weights2 matrix from 1 to 40

for i = 1:40
    ok = 0;
    % Buy stocks for each column in Weights2 matrix that has a value greater than zero
    for j = 1:PortSize

        % Check if the element in Weights2 matrix is greater than zero
        if Weights2(i,j) > 0
            
            % Initialize a variable m to 0
            m = 0;
            
            % Loop until a non-zero stock price is found in priceData2 matrix
            while priceData2( idx2(positionMatrix2(i,1))+m , positionMatrix2(i,j+1) ) == 0
                
                % If m is less than 3, increment m by 1 and continue looping
                if m < 3 
                    m = m+1;
                
                % If m is greater than or equal to 3, set m to -1 and test
                % if the price is non 0 at the previous date.
                else 
                    m = -1;
                end
            end
            
            % Calculate the number of stocks to buy for portfolio 2W
            numberOfStocks2W(i,j) = floor( ( Weights2(i,j) * PortfolioValueW(i+85) ) / priceData2( idx2(positionMatrix2(i,1))+m , positionMatrix2(i,j+1) ) );
            
            % Calculate the total cost of buying the stocks for portfolio 2W
            BuySum2W(i,j) = numberOfStocks2W(i,j) * priceData2( idx2(positionMatrix2(i,1))+m , positionMatrix2(i,j+1) );

            % Calculate the number of stocks to buy for portfolio 2U
            numberOfStocks2U(i,j) = floor( ( (1/UniformWeights2(i)) * PortfolioValueU(i+85) ) / priceData2( idx2(positionMatrix2(i,1))+m , positionMatrix2(i,j+1) ) );

            % Calculate the total cost of buying the stocks for portfolio 2U
            BuySum2U(i,j) = numberOfStocks2U(i,j) * priceData2( idx2(positionMatrix2(i,1))+m , positionMatrix2(i,j+1) );
            
            % Set m back to 0
            m = 0;
        end
    end
    
    % Calculate the cash amount remaining after buying stocks for portfolio 2W
    Cash2W(i) = PortfolioValueW(85+i) - sum(BuySum2W(i,:));
    
    % Calculate the cash amount remaining after buying stocks for portfolio 2U
    Cash2U(i) = PortfolioValueU(85+i) - sum(BuySum2U(i,:));
    
    
        % The following for-loop is used to keep track of the monthly price
        % changes of the portfolio.
    
        MonthlySeasonalityWeighted(1,z) = PortfolioValueW(85+i);    
        MonthlySeasonalityUniform(1,z)  = PortfolioValueU(85+i);
        
        z = z +1;
        
        if i < 40 && i > 1
            
            MonthUpdate = 2;
        else 
            MonthUpdate = 1;
        end
    for q = 1:MonthUpdate
        
        for j = 1:PortSize

            % Check if the stock weight is greater than 0 for the given row
            if Weights2(i,j) > 0

                % Loop to check if the specific stock data is missing or partial
                while priceData2(idx2(positionMatrix2(i,1)+q)+m , positionMatrix2(i,j+1)) == 0

                  % If m is less than 3, increment m by 1 and continue looping
                if m < 3 
                    m = m+1;
                
                % If m is greater than or equal to 3, set m to -1 and test
                % if the price is non 0 at the previous date.
                else 
                    m = -1;
                end
                end

              MonthlyValueWeighted(q,j) = numberOfStocks2W(i,j) * priceData2( idx2(positionMatrix2(i,1)+q)+m , positionMatrix2(i,j+1) ) ;
                
                
              MonthlyValueUniform(q,j)  = numberOfStocks2U(i,j) * priceData2( idx2(positionMatrix2(i,1)+q)+m , positionMatrix2(i,j+1) ) ;

            end
        end
        
        MonthlySeasonalityWeighted(1,z) = sum(MonthlyValueWeighted(q,:)) + Cash2W(i);
        
        MonthlySeasonalityUniform(1,z)  = sum(MonthlyValueUniform(q,:)) + Cash2U(i);

        z = z +1;
        
        MonthlyValueWeighted = [];
        MonthlyValueUniform  = []; 
    end

% Sell stocks for current time period if it is not the last time period (i.e. i < 40)
if i < 40
    
    % Loop through each column of Weights2 matrix
    for j = 1:PortSize
        
        % Check if the element in Weights2 matrix is greater than zero
        if Weights2(i,j) > 0
            
            % Initialize a variable m to 0
            m = 0;
            
            % Loop until a non-zero stock price is found in priceData2 matrix
            while priceData2( idx2(positionMatrix2(i,1)+3)+m , positionMatrix2(i,j+1) ) == 0
                
                % If m is less than 3, increment m by 1 and continue looping
                if m < 3 
                    m = m+1;
                
                % If m is greater than or equal to 3, set m to -1 and exit the loop
                else 
                    m = -1;
                end
                   
            end
            
            % Calculate the cash amount earned from selling stocks for portfolio 2W
            CashSellW2(i,j) = numberOfStocks2W(i,j) * priceData2( idx2(positionMatrix2(i,1)+3)+m , positionMatrix2(i,j+1) ) ;

            % Calculate the cash amount earned from selling stocks for portfolio 2U
            CashSellU2(i,j) = numberOfStocks2U(i,j) * priceData2( idx2(positionMatrix2(i,1)+3)+m , positionMatrix2(i,j+1) ) ;
            
            % Set m back to 0
            m = 0;
        end
    end
    
    % Update the portfolio value for portfolio 2W for the next time period
    PortfolioValueW(85+i+1) = Cash2W(i) + sum(CashSellW2(i,:));
    
    % Update the portfolio value for portfolio 2U for the next time period
    PortfolioValueU(85+i+1) = Cash2U(i) + sum(CashSellU2(i,:));
 end
end

%% 
% Initialize variable m to 0
m = 0;

% Loop through each column of Weights2 matrix
for j = 1:PortSize
    
    % Check if the element in Weights2 matrix is greater than zero
    if Weights2(i,j) > 0
        
        % Loop until a non-zero stock price is found in priceData2 matrix,
        % starting from the end of the matrix and decrementing m by 1
         while priceData2( idx2(positionMatrix2(i,1))+m , positionMatrix2(i,j+1) ) == 0
                    if m < 3 
                    m = m+1;
                    elseif ok == 1
                        m = -2;
                    else
                        m = -1;
                        ok = 1;
                    end
         end
        
        % Calculate the cash amount earned from selling stocks for portfolio 2W using the last non-zero stock price
        CashSellW2(i,j) = numberOfStocks2W(i,j) * priceData2( end-m , positionMatrix2(i,j+1) ) ;

        % Calculate the cash amount earned from selling stocks for portfolio 2U using the last non-zero stock price
        CashSellU2(i,j) = numberOfStocks2U(i,j) * priceData2( end-m , positionMatrix2(i,j+1) ) ;
        
        % Reset m to 0
        m = 0;
    end
end

% Update the portfolio value for portfolio 2W for the last sell order
PortfolioValueW(i+85+1) = Cash2W(i) + sum(CashSellW2(i,:));

% Update the portfolio value for portfolio 2U for the last sell order
PortfolioValueU(i+85+1) = Cash2U(i) + sum(CashSellU2(i,:));
%%
time = datenum([dateEntry1(idx1(12:end)) ; dateEntry2(idx2(13:end))])';
time = time+2;
time(49)  = time(49)+1;
time(121) = time(121)+2;
time(109) = time(109)+1;
time(181) = time(181)+1;
time(193) = time(193)+2;
time(253) = time(253)+2;
time(313) = time(313)+1;
time(325) = time(325)+2;
end

