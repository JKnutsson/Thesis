function [time,PortValueWeighted, PortValueUniform] = EVEBITDA(priceData, EVEBITDAData, dateEntry, EVEBITDADate, PortSize)

    %% Index the end of each month in the price data and EV/EBITDA variable.
    % We are buying/selling at the end of each month so instead of looping through all prices until we get to X date
    % we index the position of the end of each month in a vector and add the begining date of our data so that the 
    % portfolio starts buying in the begining of the first year.

    dat1            = datevec(dateEntry);              
    [~,idxPrice]    = unique(dat1(:,1:2),'rows','last');   % idxPrice indexes the row coresponding to the end of each month in priceData
    idxPrice        = vertcat(4, idxPrice);                % This line adds the first date in january that has both price and EV/EBITDA data
    DateIndex       = dateEntry(idxPrice);                 % idxPrice in 'dd-mm-yyyy' format

    dat2           = datevec(EVEBITDADate);              
    [~,idxEV]      = unique(dat2(:,1:2),'rows','last'); % idxEV indexes the row coresponding to the end of each month in DivData
    idxEV          = vertcat(1, idxEV);
    DateIndexEV    = EVEBITDADate(idxEV);               % idxEV in 'dd-mm-yyyy' format



    clear dat1 dat2 dat3
    %% Find the smallest EV/EBITDA ratios and positions.
    % In order to find the lowest ev/ebitda values we first need to "mask" the non-positive values with "NaN".
    % Then we can use "mink()" to find the lowest values and their position in the price matrix at the end of each month
    % using the index variable idxEV.

    mask    = EVEBITDAData <= 0;
    idxMask = find(mask);

    EVEBITDAData(mask) = NaN;
%%
    for i = 1:height(idxEV)   

        [minimumEV{i},Position{i}]   = mink(EVEBITDAData(idxEV(i), : ),PortSize);

    end
    
    % Change data type to "double" otherwise we won't be able to do operations on the matrix.
    minimumEV = minimumEV.';
    minimumEV = cell2mat(minimumEV);
    minimumEV(isnan(minimumEV)) = 0;

    Position = Position.';
    Position = cell2mat(Position);

    clear idxMask mask

    %% Weightes based on lowest EV/EBITDA

    total_EV_EBITDA = zeros(height(minimumEV),1);
    for i = 1:height(minimumEV)  
       for k = 1:width(minimumEV)
           if minimumEV(i,k) > 0
                total_EV_EBITDA(i) = total_EV_EBITDA(i) + 1/minimumEV(i,k);
           end
       end  

       for k = 1:width(minimumEV)
           if minimumEV(i,k) > 0
                Weights(i,k) = (1/minimumEV(i,k)) / total_EV_EBITDA(i);
           end
       end
    end

%%
    clear total_EV_EBITDA

    %% Portfolio-backtesting
    % This code simulates a portfolio backtesting scenario using two methods (weighted and uniform) to allocate funds among a set of stocks in a given portfolio,
    % and then tracks the resulting portfolio values over a specific time period.

    % Two initial values, PortValueW(1) and PortValueU(1), are set to the value 100,000,
    % which represent the starting portfolio value for the weighted and uniform methods, respectively.

    % Next, two loops are executed that iterate through each row of the minimumEV matrix up to the second-last row.
    % The inner loop iterates through each column of the weight matrix (which has dimensions of number of rows in minimumEV X PortSize).

    % The code calculates the maximum number of shares (NumStocksWeighted(i,k) and NumStocksUniform(i,k))
    % of each stock that can be bought based on the current available fund values (PortValueW(i) and PortValueU(i),
    % respectively) and the current stock prices.
    % It then records the amount spent buying each stock under each method (BuyWeighted(i,k) and BuyUniform(i,k)),
    % as well as the remaining cash amounts after all purchases (CashWeighted(i) and CashUniform(i)).

    % After these calculations are done for all columns,
    % the code then enters another loop that iterates again through each column but now for the next row in the minimumEV matrix. 
    % Using the same logic as before, it determines the proceeds that would be realized by selling all shares of each stock under each method,
    % and records these results as SellWeighted(i,k) and SellUniform(i,k).

    % Finally, it calculates the new total portfolio values for both methods based on the sums of cash plus sales proceeds minus purchase costs (PortValueW(i+1) and PortValueU(i+1)),
    % and the loops start again for the next iteration. 

    % Overall, this code simulates a backtesting investment strategy that alternates between buying and selling stocks
    % based on price movements and specific fund allocation strategies.


    PortValueWeighted(1) = 100000;
    PortValueUniform(1) = 100000;
    m  = 0;
    ok = 0;


    for i = 1:height(minimumEV)-1

        for k = 1:PortSize

                    % This while loop checks the prices of stock k in row i. If it's 0, then it checks the prices 
                    % +2/-2 steps from row i. If it gets stuck in an infinite loop, it's used as an indication 
                    % for me to check the data for errors or decide to exclude a specific period or remove the stock entirely.

                if Weights(i,k) > 0
                   while priceData( idxPrice(i,1)+m , Position(i,k) )  == 0

                    if m < 3 
                        m = m+1;
                        elseif ok == 1
                            m = -3;
                         Error = 128
                        else
                            m = -1;
                            ok = 1;
                            Error = 132

                     end
                   end
                    
                   
                   NumStocksWeighted(i,k) = floor( (PortValueWeighted(i) * Weights(i,k)) / priceData( idxPrice(i)+m, Position(i,k) ) );
                   
                   NumStocksUniform(i,k)  = floor( (PortValueUniform(i) * (1/PortSize)) / priceData( idxPrice(i)+m, Position(i,k) ) );

                   BuyWeighted(i,k)       = NumStocksWeighted(i,k) * priceData( idxPrice(i), Position(i,k) );

                   BuyUniform(i,k)        = NumStocksUniform(i,k) * priceData( idxPrice(i), Position(i,k) );

                   m = 0;
                end
        end


                   CashWeighted(i) = PortValueWeighted(i) - sum(BuyWeighted(i, : ));
                   

                   CashUniform(i)  = PortValueUniform(i) - sum(BuyUniform(i, : ));

              

              for k = 1:PortSize 
                    if Weights(i,k) > 0
                    % This while loop checks the prices of stock k in row i. If it's 0, then it checks the prices 
                    % +2/-2 steps from row i. If it gets stuck in an infinite loop, it's used as an indication 
                    % for me to check the data for errors or decide to exclude a specific period or remove the stock entirely.

                   while priceData( idxPrice(i+1,1)+m , Position(i,k) )  == 0

                        if m < 3 
                            m = m+1;
                            elseif ok == 1
                                m = -3;
                                Error = 167
                            else
                                m = -1;
                                ok = 1;
                                Error = 171
                         end
                   end

                   
                   SellWeighted(i,k) = NumStocksWeighted(i,k) * priceData( idxPrice(i+1)+m, Position(i,k) );

                   SellUniform(i,k)  = NumStocksUniform(i,k) * priceData( idxPrice(i+1)+m, Position(i,k) );

                   m = 0;
                   
                   end
              end  
                   PortValueWeighted(i+1) = CashWeighted(i) + sum(SellWeighted(i, :)) ;
                   PortValueUniform(i+1) = CashUniform(i) + sum(SellUniform(i, :)) ;
                   

    end


    %%
    clear ok m


    % Time is the X-axis and since it's based on the DateIndex which is the end
    % last date of the month that has a price during the period I add 2 to it
    % so that the X-axis start at 01-01-xxxx in the graphs
    time = datenum(DateIndex);
    time = time + 2; 
    time(1) = time(1)-2;
    time(13) = time(13)+2;
    time(73) = time(73)+1;
    time(85) = time(85)+2;
    time(145) = time(145)+2;
    time(205) = time(205)+1;
    time(217) = time(217)+2;
end

