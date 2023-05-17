function [PortfolioValueMonthlyW, PortfolioValueMonthlyU, PortValueOMX ] = MomentumMonthly(priceData1, priceData2, dateEntry1 , dateEntry2, OMX,PortSize)

%% Momentum portfolio with rebalancing at the end of each month.


% Convert dateEntry1 into a matrix format (year, month, day) using datevec.
dat1        = datevec(dateEntry1);

% Find the unique rows of 'dat1' that correspond to the end of each month
% and store the index corresponding to these rows in 'idx1'.
% The 'last' option specifies that the latest occurrence of each row should be used.
[~,idx1]    = unique(dat1(:,1:2),'rows','last');


% Extract the number of columns in priceData1 and store it in width1.
[~, width1] = size(priceData1);

% Repeat the above steps for 'dat2' and 'priceData2'.
dat2        = datevec(dateEntry2);
[~,idx2]    = unique(dat2(:,1:2),'rows','last');
%idx2(end)   = [];
[~, width2] = size(priceData2);

% Extract the date entries in the OMX table and convert them to matrix format.
dateEntryOMX   = OMX.RowLabels;
dat3           = datevec(dateEntryOMX);

% Find the unique rows of 'dat3' that correspond to the end of each month
% and store the index corresponding to these rows in 'idxOMX'. 
% The 'last' option specifies that the latest occurrence of each row should be used.
[~,idxOMX]     = unique(dat3(:,1:2),'rows','last');

% Convert the OMX table to a numeric array, excluding the date column.
OMX            = table2array(OMX( : , 2:end));

% Extract the number of columns in OMX and store it in width3.
[~, width3]    = size(OMX);
%%

% Initialize a new matrix 'Momentum1' with the same number of rows as 'idx1'
% and the same number of columns as 'priceData1', to store momentum data.
Momentum1 = zeros(height(idx1) , width1); 

Start_date = 12;

% Loop through each month starting from the 12th month in idx1.
for l = Start_date:height(idx1)
    
    % Loop through each stock in priceData1.
    for k = 1:width1       
    
        % Check if end-of-month prices exist for both the current month and
        % the month 12 months ago, and skip if they don't.
        if priceData1( idx1(l) ,k) ~= 0 && priceData1( idx1(l-11) ,k) ~= 0  
        
            % Calculate the momentum for the current stock.
            if l > 12
            
            Momentum1(l,k) = ( (priceData1( idx1(l-11) ,k) - priceData1(idx1(l-12)+1 ,k) ) / priceData1( idx1(l-12)+1 ,k)) +1;
            
            
            else
            
            Momentum1(l,k) = ( (priceData1( idx1(l-11) ,k) - priceData1( 1 ,k) ) / priceData1( 1 ,k)) +1;
            
            end
            % Loop through the past 11 months, skipping the most recent month.
            for i = 10:-1:1
                
                m = 0;  
                
                % If one data point is missing in the set, change row position by 1.
                if priceData1( idx1(l-i) ,k) == 0 
                    m = 1;              
                end
                
                % Calculate the momentum for the current stock based on the past i months.
                Momentum1(l,k) = Momentum1(l,k) * (((priceData1( idx1(l-i)+m ,k) - priceData1( idx1(l-(i+1))+m ,k)) / priceData1( idx1(l-(i+1))+m ,k))+1);
                
            end 
            % Compute the cumulative returns (generic momentum) for each stock.
            Momentum1(l,k) = Momentum1(l,k)-1;
        end  
    end
end

% Repeat the above steps for priceData2 to calculate 'Momentum2' for the second time-period.
Momentum2 = zeros(height(idx2) , width2); 

for l = Start_date:height(idx2)        
    
    for k = 1:width2            
        
        if priceData2( idx2(l) ,k) ~= 0 && priceData2( idx2(l-11) ,k) ~= 0  
            
            if l > 12
            
            Momentum2(l,k) = ( (priceData2( idx2(l-11) ,k) - priceData2( idx2(l-12)+1 ,k) ) / priceData2( idx2(l-12)+1 ,k))+1;   
            
            else
                
                Momentum2(l,k) = ( (priceData2( idx2(l-11) ,k) - priceData2( 1 ,k) ) / priceData2( 1 ,k))+1;
            
            end
            
            for i = 10:-1:1               
                m = 0;                  
                if priceData2( idx2(l-i) ,k) == 0                  
                    m = 1;              
                end                
                
                Momentum2(l,k) = Momentum2(l,k) * (((priceData2( idx2(l-i)+m ,k) - priceData2( idx2(l-(i+1))+m ,k)) / priceData2( idx2(l-(i+1))+m ,k))+1); 
                
            end 
            
            Momentum2(l,k) = Momentum2(l,k)-1; 
        end  
    end
end

%% Clean up inf / NaN's.


% Replace infinite values in 'Momentum1' with 0.
Momentum1(isinf(Momentum1)) = 0;

% Replace NaN values in 'Momentum1' with 0.
Momentum1(isnan(Momentum1)) = 0;

% Compute the total error for Momentum1, assuming a total of 663+445 data points were affected.
% The error is calculated as the sum of all the replaced values divided by the total number of data points.
% In this case, the error is 0.001408.
 
% Repeat for Momentum2.
Momentum2(isinf(Momentum2)) = 0;
Momentum2(isnan(Momentum2)) = 0;
% Compute the total error for Momentum2, assuming a total of 131+32 data points were affected.
% The error is calculated as the sum of all the replaced values divided by the total number of data points.
% In this case, the error is 0.0044.

% Initialize new matrices to store positive and negative momentum values for each month and stock, as well as
% information discreteness measures ID1 and ID2, which are not yet calculated in this section of the code.
positive1    = zeros(height(idx1) , width1);
negative1    = zeros(height(idx1) , width1);
positive2    = zeros(height(idx1) , width1);
negative2    = zeros(height(idx1) , width1);
ID1          = zeros(height(idx1) , width1);

positive3    = zeros(height(idx2) , width2);
negative3    = zeros(height(idx2) , width2);
positive4    = zeros(height(idx2) , width2);
negative4    = zeros(height(idx2) , width2);
ID2          = zeros(height(idx2) , width2);



% Loop over months, stocks and the past year's trading days (251 to 2).
% Calculate the number of positive and negative trading days for each month and stock
% based on the comparison between the stock price p days ago and the stock price (p+1) days ago.
% If a stock price is zero or NaN at either point, that data point will be skipped.
% Store the count of positive and negative days in separate matrices, 'positive1' and 'negative1'.
% The initial value of m changes when l is 12 because the first year's momentum is
% calculated using data from the end of the fourth month of trading rather than the first.
% Repeat the same process for another dataset, storing the count of positive and negative days in 
% matrices named 'positive3' and 'negative3'.
for l = Start_date:height(idx1)                  
    for k = 1:width1               
        for p = 251:-1:2        
            if l == Start_date
                 m = 4;
            end
            if priceData1( idx1(l)-p+m , k) ~= 0 && priceData1( idx1(l)-p+1+m , k) ~= 0 
                if priceData1( idx1(l)-p+m , k) > priceData1( idx1(l)-p+1+m , k)
                    positive1(l,k) = positive1(l,k) + 1;
                end  
                if priceData1( idx1(l)-p +m, k) < priceData1( idx1(l)-p+1+m , k) 
                    negative1(l,k) = negative1(l,k) + 1;      
                end            
            end
            m = 0;
        end             
    end
end

% Similar to the above loop, calculate the number of positive and negative trading days, but for another dataset,
% storing the count in matrices named 'positive3' and 'negative3'.
% The initial value of m changes when l is 12 because the first year's momentum
% is calculated using data from the end of the second month of trading rather than the first.
for l = Start_date:height(idx2)                  
    for k = 1:width2               
        for p = 251:-1:2        
            if l == Start_date
                 m = 2;
            end
            if priceData2( idx2(l)-p+m , k) ~= 0 && priceData2( idx2(l)-p+1+m , k) ~= 0 
                if priceData2( idx2(l)-p+m , k) > priceData2( idx2(l)-p+1+m , k) 
                    positive3(l,k) = positive3(l,k) + 1;
                end  
                if priceData2( idx2(l)-p +m, k) < priceData2( idx2(l)-p+1+m , k) 
                    negative3(l,k) = negative3(l,k) + 1;      
                end            
            end
            m = 0;
        end             
    end
end


%%            Convert the number of positive / negative days to a '%'.
            
% Loop over months and stocks to calculate the percentage of positive and negative trading days
% for each month and stock, based on the count of those days calculated in the previous loop.
% Store the percentage of positive and negative days in separate matrices, 'positive2' and 'negative2'.
% Then create a filter, ID1, that captures how momentum changes with respect to positive and negative days.
for l = Start_date:height(idx1) 
    for k = 1:width1               
            positive2(l,k)          = positive1(l,k) / (positive1(l,k) + negative1(l,k));
            negative2(l,k)          = negative1(l,k) / (positive1(l,k) + negative1(l,k));
            % Create a filter based on the "Frog-in-the-pan"
            % theory i.e. if the momentum is based on many small 
            % changes in price or fewer large changes.
            ID1(l,k) = sign(Momentum1(l,k)) * (negative2(l,k) - positive2(l,k));       
    end
end

% Similar to the above loop, calculate the percentage of positive and negative trading days,
% store them in separate matrices named 'positive4' and 'negative4', and create a filter, ID2.
for l = Start_date:height(idx2) 
    for k = 1:width2               
            positive4(l,k)          = positive3(l,k) / (positive3(l,k) + negative3(l,k));
            negative4(l,k)          = negative3(l,k) / (positive3(l,k) + negative3(l,k));
            % Create a filter based on the "Frog-in-the-pan"
            % theory i.e. if the momentum is based on many small 
            % changes in price or fewer large changes.
            ID2(l,k) = sign(Momentum2(l,k)) * (negative4(l,k) - positive4(l,k));     
    end
end

 

% Replace infinite values in matrix ID1 with zero.
ID1(isinf(ID1)) = 0;

% Replace NaN (Not a Number) values in matrix ID1 with zero.
ID1(isnan(ID1)) = 0;

% Replace infinite values in matrix ID2 with zero.
ID2(isinf(ID2)) = 0;

% Replace NaN (Not a Number) values in matrix ID2 with zero.
ID2(isnan(ID2)) = 0;



% Create matrices to filter stocks based on positive momentum and high quality.
momentumQualityFiltered1 = zeros(height(idx1) , width1);
momentumQualityFiltered2 = zeros(height(idx2) , width2);

% Loop over the monthly data from the 12th month to the end of the data.
for l = Start_date:height(idx1)           % l is the index number for each month.
    
    % Loop over each stock in the dataset.
    for k = 1:width1        % k keeps track of the colums (stocks).         
    
        % Check if the current stock has positive momentum and high quality.
        if Momentum1( l ,k) > 0 && ID1( l ,k) < 0  
                           
            % Index the current stock in the filtered matrix.
            momentumQualityFiltered1(l,k) = Momentum1(l,k);  
                       
        end  
    end
end

% Repeat the process for the second dataset.
for l = Start_date:height(idx2)           % l is the index number for each month.
    
    for k = 1:width2        % k keeps track of the colums (stocks).         
    
        if Momentum2( l ,k) > 0  && ID2( l ,k) < 0  
                           
            momentumQualityFiltered2(l,k) = Momentum2(l,k);  
                       
        end  
    end
end


% These lines of code are used to filter stocks based on positive momentum and high quality.
% Two matrices, `momentumQualityFiltered1` and `momentumQualityFiltered2`,
% are created with the same size as their corresponding datasets to store the filtered results.
% Within each dataset, this code loops through each month and each stock to identify which stocks
% meet the filtering criteria. For a stock to be included in the filtered result,

%% FIND LARGEST PortSize(i) MOMENTUM STOCKS AND THEIR POSITIONS.

% Portfolio is rebalanced at the end of each month.
largestMomentum1{1} = maxk(momentumQualityFiltered1(Start_date, : ),PortSize);

% Loop through each year of data and identify the 30 stocks with the largest momentum for each month.
for y = 1:height(idx1)-12 % number of years in the set
    
    largestMomentum1{y+1} = maxk(momentumQualityFiltered1(Start_date+y, : ),PortSize);
  
end
%%
% Convert to a matrix and transpose the rows and columns for easier calculation.
largestMomentum1 = largestMomentum1.';
largestMomentum1 = cell2mat(largestMomentum1);

% Repeat the process using the second dataset.
largestMomentum2{1}   = maxk(momentumQualityFiltered2(Start_date, : ),PortSize);

for y = 1:height(idx2)-12 % number of years in the set
    
    largestMomentum2{y+1} = maxk(momentumQualityFiltered2(Start_date+y, : ),PortSize);
  
end

largestMomentum2 = largestMomentum2.';
largestMomentum2 = cell2mat(largestMomentum2);

% Calculate the sum of the momentum values for each row (month) of the two datasets.
for i = 1:height(largestMomentum1)
    
    sumMomentum1(i,1) = sum(largestMomentum1(i, : ));
    
end

for i = 1:height(largestMomentum2)
    
    sumMomentum2(i,1) = sum(largestMomentum2(i, : ));
    
end


%% Create matrices to store the positions of the top PortSize(i) momentum stocks for each month.

for i = 1:height(largestMomentum1)
    count = zeros(1, PortSize);
    for k = 1:PortSize
        % Check if the current momentum value is greater than 0.
        % If so, find the corresponding column index in the filtered dataset and store it.
        if largestMomentum1(i,k) > 0
            idx = i+Start_date-1 : i+Start_date-1;
            if count(k) == 0
                [ ~ , c1(i,k)] = find(momentumQualityFiltered1(idx, :) == largestMomentum1(i,k), 1);
                count(k) = 1;
            else
                [ ~ , c1(i,k)] = find(momentumQualityFiltered1(idx, :) == largestMomentum1(i,k), 1, 'last');
            end
%             positionMatrix1(i,k) = c1(i,k);
        end
    end
end

clear count idx 

% Convert the matrices to a table with the corresponding month-year indexes.
positionMatrix1 = table(idx1(Start_date:end),c1);
positionMatrix1 = table2array(positionMatrix1);


for i = 1:height(largestMomentum2)
    count = zeros(1, PortSize);
    for k = 1:PortSize
        % Check if the current momentum value is greater than 0.
        % If so, find the corresponding column index in the filtered dataset and store it.
        if largestMomentum2(i,k) > 0
            idx = i+Start_date-1 : i+Start_date-1;
            if count(k) == 0
                [ ~ , c2(i,k)] = find(momentumQualityFiltered2(idx, :) == largestMomentum2(i,k), 1);
                count(k) = 1;
            else
                [ ~ , c2(i,k)] = find(momentumQualityFiltered2(idx, :) == largestMomentum2(i,k), 1, 'last');
            end
%             positionMatrix1(i,k) = c1(i,k);
        end
    end
end


% Convert the matrices to a table with the corresponding month-year indexes.
positionMatrix2 = table(idx2(Start_date:end),c2);
positionMatrix2 = table2array(positionMatrix2);


% This code calculates the weights for a portfolio of stocks, based on their momentum values.
% The `largestMomentum` matrices identify the 30 stocks with the highest momentum value each month,
% for each dataset. The `sumMomentum` matrices calculate the sum of the momentum values
% for each month for each dataset. The `positionMatrix` matrices store the positions of the top 30 momentum stocks
% for each month along with their respective month-year



%% Calculate the wheights in each month.

% Create matrices to store uniform weights and momentum-based weights for each dataset.
UniformWeights1 = zeros(height(largestMomentum1),1);
UniformWeights2 = zeros(height(largestMomentum2),1);
Weights1        = zeros(height(largestMomentum1), PortSize);
Weights2        = zeros(height(largestMomentum2), PortSize);

% Loop through each row (month), and calculate the momentum-based weight and the uniform weight for each stock.
for i = 1:height(largestMomentum1)

    for j = 1:PortSize
         
        % Check if the current momentum value is greater than 0.
        % If so, calculate the momentum-based weight and add 1 to the uniform weight.
        if largestMomentum1(i,j) > 0 

            Weights1(i,j)      = largestMomentum1(i,j) / sumMomentum1(i);
            UniformWeights1(i) = UniformWeights1(i) + 1;

        end
    end
end

% Repeat the process using the second dataset.
for i = 1:height(largestMomentum2)

    for j = 1:PortSize
         
        if largestMomentum2(i,j) > 0 

            Weights2(i,j)      = largestMomentum2(i,j) / sumMomentum2(i);
            UniformWeights2(i) = UniformWeights2(i) + 1;

        end
    end
end



%% Portfolio Backtesting


PortfolioValueMonthlyW(1) = 100000;
PortfolioValueMonthlyU(1) = 100000;


BuySum1W    = zeros(height(largestMomentum1), PortSize);
BuySum1U    = zeros(height(largestMomentum1), PortSize);
CashSell1W  = zeros(height(largestMomentum1), PortSize);
CashSell1U  = zeros(height(largestMomentum1), PortSize);

m  = 0;

for i = 1:height(largestMomentum1)-1
    ok = 0;
      
            
               % Buy
% This section of code is responsible for calculating the buying related parameters such as 
% number of stocks to buy, and buying sums for different instruments.

for j = 1:PortSize
    % Loop through all columns of Weights1 matrix
    if Weights1(i,j) > 0
        % If weight value is greater than zero for a particular stock
        % While there are no price data available, adjust the time-stamp 
        % until price data is found. If too many adjustments have been made,
        % flag an error.
        m = 0;
        ok = 0;
        while priceData1( positionMatrix1(i,1)+m , positionMatrix1(i,j+1) )  == 0
            % Loop through price data for that stock until valid price data is encountered
            if m < 3 
                % Adjust time stamp by one day (24 hours)
                m = m+1;
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
        
        % Calculate number of stocks to buy based on respective weights and portfolio values
        numberOfStocks1W(i,j) = floor( (Weights1(i,j) * PortfolioValueMonthlyW(i)) / priceData1( positionMatrix1(i,1)+m , positionMatrix1(i,j+1) ) );
        
        % Calculate the total buy sum for each instrument
        BuySum1W(i,j) = numberOfStocks1W(i,j) * priceData1( positionMatrix1(i,1)+m , positionMatrix1(i,j+1) );
        
        % Calculate number of stocks to buy based on uniform weights and portfolio values
        numberOfStocks1U(i,j) = floor( ((1/UniformWeights1(i)) * PortfolioValueMonthlyU(i)) / priceData1( positionMatrix1(i,1)+m , positionMatrix1(i,j+1) ) );
        
        % Calculate the total buy sum for each instrument
        BuySum1U(i,j) = numberOfStocks1U(i,j) * priceData1( positionMatrix1(i,1)+m , positionMatrix1(i,j+1) );
                    
    end
end


    % Calculate the cash remaining after buying
        Cash1W(i) = PortfolioValueMonthlyW(i) - sum(BuySum1W(i,:));
        Cash1U(i) = PortfolioValueMonthlyU(i) - sum(BuySum1U(i,:));


                                                  % Sell
% This section of code is responsible for calculating the selling related parameters such as 
% cash sell amounts for different instruments.

for j = 1:PortSize
    % Loop through all columns of Weights1 matrix
    ok = 0;
    if Weights1(i,j) > 0
        % If weight value is greater than zero for a particular stock
        % While there are no price data available, adjust the time-stamp 
        % until price data is found. If too many adjustments have been made,
        % flag an error.
        m = 0;
        while priceData1(positionMatrix1(i+1,1)+m , positionMatrix1(i,j+1)) == 0
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
        
        % Calculate the cash sell amount for each instrument
        CashSell1W(i,j) = numberOfStocks1W(i,j) * priceData1(positionMatrix1(i+1,1)+m , positionMatrix1(i,j+1));
        
        % Calculate the cash sell amount for each uniform weight instrument
        CashSell1U(i,j) = numberOfStocks1U(i,j) * priceData1(positionMatrix1(i+1,1)+m , positionMatrix1(i,j+1));
        
        % Reset timestamp to zero for each instrument
        m = 0;
                    
    end
end

      
            
    % Calculate the new portfolio value after selling
        PortfolioValueMonthlyW(i+1) = Cash1W(i) + sum(CashSell1W(i,:));
        PortfolioValueMonthlyU(i+1) = Cash1U(i) + sum(CashSell1U(i,:));
        
      
end





%% The following code works the same way just for the seconds set of data.


BuySum2W    = zeros(height(largestMomentum2), PortSize);
BuySum2U    = zeros(height(largestMomentum2), PortSize);
CashSell2W  = zeros(height(largestMomentum2), PortSize);
CashSell2U  = zeros(height(largestMomentum2), PortSize);

for i = 1:height(largestMomentum2)-1
    ok = 0;
                                                    % Buy
    
        for j = 1:PortSize

            if Weights2(i,j) > 0
                
                while priceData2( positionMatrix2(i,1)+m , positionMatrix2(i,j+1) )  == 0
                    
                     if m < 3 
                    m = m+1;
                    else 
                        m = -1;
                    end
                end

                    numberOfStocks2W(i,j) = floor( ( Weights2(i,j) * PortfolioValueMonthlyW(height(largestMomentum1)-1+i) ) / priceData2( positionMatrix2(i,1)+m , positionMatrix2(i,j+1) ) );

                    BuySum2W(i,j) = numberOfStocks2W(i,j) * priceData2( positionMatrix2(i,1)+m , positionMatrix2(i,j+1) );

                    numberOfStocks2U(i,j) = floor( ( (1/UniformWeights2(i)) * PortfolioValueMonthlyU(height(largestMomentum1)-1+i) ) / priceData2( positionMatrix2(i,1)+m , positionMatrix2(i,j+1) ) );

                    BuySum2U(i,j) = numberOfStocks2U(i,j) * priceData2( positionMatrix2(i,1)+m , positionMatrix2(i,j+1) );
                    
                    m = 0;

            end
        end

            Cash2W(i) = PortfolioValueMonthlyW(height(largestMomentum1)+i-1) - sum(BuySum2W(i,:));
            
            Cash2U(i) = PortfolioValueMonthlyU(height(largestMomentum1)+i-1) - sum(BuySum2U(i,:));
             

                                                    % Sell
 

        for j = 1:PortSize
            ok = 0;
            if Weights2(i,j) > 0
                
                while priceData2( positionMatrix2(i+1,1)+m , positionMatrix2(i,j+1) ) == 0
                    if m < 3 
                    m = m+1;
                    elseif ok == 1
                        m = -2;
                    else
                        m = -1;
                        ok = 1;
                    end
                end

                    CashSell2W(i,j) = numberOfStocks2W(i,j) * priceData2( positionMatrix2(i+1,1)+m , positionMatrix2(i,j+1) ) ;

                    CashSell2U(i,j) = numberOfStocks2U(i,j) * priceData2( positionMatrix2(i+1,1)+m , positionMatrix2(i,j+1) ) ;

                        m = 0;

            end
        end
            
     
        
        
                     PortfolioValueMonthlyW(height(largestMomentum1)+i) = Cash2W(i) + sum(CashSell2W(i,:));
                     
                     PortfolioValueMonthlyU(height(largestMomentum1)+i) = Cash2U(i) + sum(CashSell2U(i,:));
           
end




% Calculates the number of OMX shares that can be bought with 100,000 SEK 
% and initializes portfolio values as well as cash amounts
for j = 1:width3
    
    numberOfOMX(1,j)  = floor( 100000 / OMX( positionMatrix1(1) , j) );

    PortValueOMX(1,j) = numberOfOMX(1,j) * OMX( positionMatrix1(1) , j);

    CashOMX(1,j)      = 100000 - PortValueOMX(1,j);
    
end

% Loop through all timestamps in the idxOMX array and calculate the 
% latest portfolio values based on the adjusted timestamp values
for i = 1:height(idxOMX)
    
    % Loop through all OMX instruments
    for j = 1:width3
        
        ok = 0;
                
                % If no price data is found at a particular time-stamp, 
                % adjust the timestamp until valid data is found.
                while OMX( idxOMX(i)+m , j ) == 0
                    if m < 3 
                        % Adjust timestamp by one day (24 hours)
                        m = m+1;
                    elseif ok == 1
                        % If more than three adjustments have been made, set m value to -2
                        m = -2;
                    else
                        % If only three adjustments were made, adjust timestamp by a larger amount
                        % to find valid data
                        m = -1;
                        ok = 1;

                    end
                end
        
       % Calculate new portfolio value from the latest price data
       PortValueOMX(i+1,j) = numberOfOMX(1,j) * OMX( idxOMX(i)+m , j);
       
    end
end

% Add the remaining cash amounts to the end of the PortValueOMX matrix
for j = 1:width3

    PortValueOMX(end,j) = PortValueOMX(end,j)+ CashOMX(1,j)  ;

end

% Handle cases where portfolio value is infinity or NaN
PortValueOMX(isinf(PortValueOMX)) = 0;
PortValueOMX(isnan(PortValueOMX)) = 0;



