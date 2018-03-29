%% Loadings
clear all;
%loading files
load('SIZE.mat');
%load('IndustryPortfolios.mat');
load('monthlyReturns.mat');
%formatting NaN 
Size(Size==-99.9900)=NaN;
monthlyReturns1(monthlyReturns1==-99.9900)=NaN;
monthlyReturns = monthlyReturns1/100; 
clear monthlyReturns1
%IndustryPortfolios(IndustryPortfolios==-99.9900)=NaN;
%monthlyReturns = IndustryPortfolios/100;
 
%% Estimating Returns of Size strategy 
numMonths = size(Size,1);
ewR=nan(numMonths - 1,1);
vwR=nan(numMonths - 1,1);
for iInvestmentDate = 1:numMonths-1
    % Sort the neccessary criterion
    % What is the sortCriterion : Size
    sortCriterion = Size(iInvestmentDate,:);
    % go to net period returns 
    sortCriterion(isnan(monthlyReturns(iInvestmentDate+1,:)))=nan;
    
    % Sort the sortCriterion
    [sortedSize,indexSorted] = sort(sortCriterion ,2);
    
    % number of active assets and number of assets in decile 
    numOfAssetsOfMonth=sum(~isnan(sortedSize),2);
    numInDecileThisIteration=round(numOfAssetsOfMonth/10);
    
    % Find the long and the short assets 
    ixShortThisIteration=indexSorted(numOfAssetsOfMonth-numInDecileThisIteration+1:numOfAssetsOfMonth);
    ixLongThisIteration=indexSorted(1:numInDecileThisIteration);
    
    %realized returns of equally weighted strategy 
    ewR(iInvestmentDate)...
        =mean(monthlyReturns(iInvestmentDate+1,ixLongThisIteration))-...
         mean(monthlyReturns(iInvestmentDate+1,ixShortThisIteration)); 
    
    weightLong= Size(iInvestmentDate,ixLongThisIteration)/sum(Size(iInvestmentDate,ixLongThisIteration));
    weightShort= Size(iInvestmentDate,ixShortThisIteration)/sum(Size(iInvestmentDate,ixShortThisIteration));
 
    vwR(iInvestmentDate)...
        =sum(monthlyReturns(iInvestmentDate+1,ixLongThisIteration).* weightLong)-...
         sum(monthlyReturns(iInvestmentDate+1,ixShortThisIteration).*weightShort);      
end

MeanEwr = mean(ewR,1);
MeanVwr = mean(vwR,1);

%% CAPM Regression of Size Strategy  
load('FFResearchDataFactors.mat');
FFResearchDataFactors = FFResearchDataFactors( 2 : end , :);
MarketPremium = FFResearchDataFactors(:,1) - FFResearchDataFactors(:,4);
[EstCov, NwstdEqCapm, CapmEqBetas ] =  hac(MarketPremium, ewR );
tstatEqCapm = CapmEqBetas./NwstdEqCapm;
[EstCov, NwstdVwCapm, CapmVwBetas ] =  hac(MarketPremium, vwR );
tstatVwCapm = CapmVwBetas./NwstdEqCapm;
% We get no significant a. The returns of size factor in industry
% portfolios are explained by market.

%% Figure 1: plot equally vs value weighted portfolio returns based on size strategy
dm=load('DATESIZE');
figure 
plot(datenum(num2str(dm.Date2(2:1099)),'yyyymm'),cumprod(1+ewR)-1);hold on
plot(datenum(num2str(dm.Date2(2:1099)),'yyyymm'),cumprod(1+vwR)-1)
legend('equally weighted','value weighted');
datetick

%% Figure 2: Compounded  Return and Cumulative log Excess Return
figure 
cumReturns = cumprod( ewR +1  ) -1 ;
cumLogReturns = cumsum( log( 1+ ewR) ) ;
yyaxis left
plot(datenum(num2str(dm.Date2(2:1099)),'yyyymm'),cumReturns); hold on
ylabel('Cumulative Size')
yyaxis right
plot(datenum(num2str(dm.Date2(2:1099)),'yyyymm'),cumLogReturns); 
legend( {'Compounded  Return' , 'Cumulative log Excess Return'   } );
datetick
xlabel('Date')
hold off

%% Dependent Sort 
load('booktomarket.mat');
booktomarket(booktomarket==-99.9900)=NaN;
%format BM matrix 
t=0;
for i=1:92
    for j=t:t+11
        BM(i+j,:) = booktomarket(i,:);    
    end
    t=j;
end
clear i j t booktomarket;
% delete five last months, we want matrix to end 
BM=BM(1:1099,:); 
FinalSize = nan(1099,49);
newSize = nan(1099,49);
newMonthlyReturns = nan(1099,49);
FinalReturns = nan(1099,49);
for iDate= 1:1099-1
    criterion = BM(iDate,:);
    criterion(isnan(Size(iDate,:)))=nan; 
    criterion(isnan(monthlyReturns(iDate+1,:)))=nan;
    
    [sortedBM,indexBM] = sort(criterion(1,:));
    numOfAssetsMonth = sum(~isnan(BM),2);
    numInQuantile = round(numOfAssetsMonth/5);
    for j=1:numOfAssetsMonth 
        newSize(iDate,j) = Size(iDate,indexBM(1,j));
        newMonthlyReturns(iDate+1,j) = monthlyReturns(iDate+1,indexBM(1,j));
    end
    stepQuantile = 1 ;
    
for iQuantile=1:5
    if iQuantile==5
    [FinalSize(iDate,stepQuantile:end),indexfinal] = sort(newSize(iDate, ...
        stepQuantile:end),2);
         
          FinalReturns(iDate+1,stepQuantile:end) = ...
          newMonthlyReturns(iDate+1,indexfinal); 
    else    
    [FinalSize(iDate,stepQuantile:stepQuantile+numInQuantile-1),indexfinal] = sort(newSize(iDate, ...
        stepQuantile:stepQuantile+numInQuantile-1),2);
      
          FinalReturns(iDate+1,stepQuantile:stepQuantile+numInQuantile-1) = ...
          newMonthlyReturns(iDate+1,indexfinal);   
      
    stepQuantile = stepQuantile+numInQuantile;
    end
end
            
end
  
FinalReturns(isnan(newMonthlyReturns))=nan;
DoubleSortedReturns = FinalReturns;
DoubleSortedSize = FinalSize; 
clear iQuantile iDate newmonthlyReturns j newSize FinalReturns FinalSize;
 

%eqW

DoubleSortedReturns = DoubleSortedReturns(2:1099,:);
for i=1:1098 
 numOfAssetsMonth = sum(~isnan(DoubleSortedReturns),2);
 numInQuantile = round(numOfAssetsMonth/5);
    EqRetQuantilesDouble(i,1) = mean(DoubleSortedReturns(i,1:round(0.04*numOfAssetsMonth)));
    t=0.04;
    l=0.08;
    for j=2:24
        EqRetQuantilesDouble(i,j) = nanmean(DoubleSortedReturns(i,round(t*numOfAssetsMonth)+1 : round(l*numOfAssetsMonth)));
        t=t+0.04; 
        l=l+0.04;    
    end
    EqRetQuantilesDouble(i,25) = nanmean(DoubleSortedReturns(i,round(0.96*numOfAssetsMonth)+1 : numOfAssetsMonth)); 
end

stili = 1;
for j=0:5:20
  DiferenceWithinQuantileEq(:,stili) = EqRetQuantilesDouble(:,j+1)-EqRetQuantilesDouble(:,j+4);
  stili = stili + 1;
end
clear stili 
%means for the tables 
meanDifferenceEq = mean(DiferenceWithinQuantileEq,1);
meanInQuantilesEq = mean(EqRetQuantilesDouble,1);
%regression
%give the same format 
MarketPremium = MarketPremium/100; 

[EstCov, NwstdEqCapm1, CapmEqBetas1 ] =  hac(MarketPremium, DiferenceWithinQuantileEq(:,1) );
tstatEqCapm1 = CapmEqBetas1./NwstdEqCapm1;

[EstCov, NwstdEqCapm2, CapmEqBetas2 ] =  hac(MarketPremium, DiferenceWithinQuantileEq(:,2) );
tstatEqCapm2 = CapmEqBetas2./NwstdEqCapm2;

[EstCov, NwstdEqCapm3, CapmEqBetas3 ] =  hac(MarketPremium, DiferenceWithinQuantileEq(:,3) );
tstatEqCapm3 = CapmEqBetas3./NwstdEqCapm3;

[EstCov, NwstdEqCapm4, CapmEqBetas4 ] =  hac(MarketPremium, DiferenceWithinQuantileEq(:,4) );
tstatEqCapm4 = CapmEqBetas4./NwstdEqCapm4;

[EstCov, NwstdEqCapm5, CapmEqBetas5 ] =  hac(MarketPremium, DiferenceWithinQuantileEq(:,5) );
tstatEqCapm5 = CapmEqBetas5./NwstdEqCapm5;

%% vw 
DoubleSortedSize = DoubleSortedSize(1:1098,:);
WeightMatrix = DoubleSortedReturns .* DoubleSortedSize ;  
for i=1:1098
    numOfAssetsMonth = sum(~isnan(DoubleSortedReturns),2);
    numInQuantile = round(numOfAssetsMonth/5);
    VwRetQuantilesDouble(i,1) = sum(WeightMatrix(i,1:round(0.04*numOfAssetsMonth))) / sum(DoubleSortedSize(i,1:round(0.04*numOfAssetsMonth)));   
    t=0.04;
    l=0.08;
    for j=2:24
        VwRetQuantilesDouble(i,j) = sum(WeightMatrix(i,round(t*numOfAssetsMonth)+1 : round(l*numOfAssetsMonth))) / sum(DoubleSortedSize(i,round(t*numOfAssetsMonth)+1 : round(l*numOfAssetsMonth)));
        t=t+0.04; 
        l=l+0.04;    
    end
        VwRetQuantilesDouble(i,25) = sum(WeightMatrix(i,round(0.96*numOfAssetsMonth)+1 : numOfAssetsMonth)) / sum(DoubleSortedSize(i,round(0.96*numOfAssetsMonth)+1 : numOfAssetsMonth));
end

stili = 1;
for j=0:5:20
  DiferenceWithinQuantileVw(:,stili) = VwRetQuantilesDouble(:,j+1)-VwRetQuantilesDouble(:,j+4);
  stili = stili + 1;
end
clear stili 

%means for the tables 
meanDifferenceVw = mean(DiferenceWithinQuantileVw,1);
meanInQuantilesVw = mean(VwRetQuantilesDouble,1);


%regression Vw 
[EstCov, NwstdVwCapm1, CapmVwBetas1 ] =  hac(MarketPremium, DiferenceWithinQuantileVw(:,1) );
tstatVwCapm1 = CapmVwBetas1./NwstdVwCapm1;

[EstCov, NwstdVwCapm2, CapmVwBetas2 ] =  hac(MarketPremium, DiferenceWithinQuantileVw(:,2) );
tstatVwCapm2 = CapmVwBetas2./NwstdVwCapm2;

[EstCov, NwstdVwCapm3, CapmVwBetas3 ] =  hac(MarketPremium, DiferenceWithinQuantileVw(:,3) );
tstatVwCapm3 = CapmVwBetas3./NwstdVwCapm3;

[EstCov, NwstdVwCapm4, CapmVwBetas4 ] =  hac(MarketPremium, DiferenceWithinQuantileVw(:,4) );
tstatVwCapm4 = CapmVwBetas4./NwstdVwCapm4;

[EstCov, NwstdVwCapm5, CapmVwBetas5 ] =  hac(MarketPremium, DiferenceWithinQuantileVw(:,5) );
tstatVwCapm5 = CapmVwBetas5./NwstdVwCapm5;


%%
for i=1:1098
    MeanRetInQuantile(i,1) =100* mean(DoubleSortedReturns(i,1:round(0.2*numOfAssetsMonth)));
    MeanRetInQuantile(i,2) =100* mean(DoubleSortedReturns(i,1:round(0.2*numOfAssetsMonth)+1 : round(0.4*numOfAssetsMonth)));
    MeanRetInQuantile(i,3) =100* mean(DoubleSortedReturns(i,1:round(0.4*numOfAssetsMonth)+1 : round(0.6*numOfAssetsMonth)));
    MeanRetInQuantile(i,4) =100* mean(DoubleSortedReturns(i,1:round(0.6*numOfAssetsMonth)+1 : round(0.8*numOfAssetsMonth)));
    MeanRetInQuantile(i,5) =100* mean(DoubleSortedReturns(i,1:round(0.8*numOfAssetsMonth)+1 :end));
end

average=mean(MeanRetInQuantile,1); 

