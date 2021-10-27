% set 1, prob 8 

x = rand(100000,1); 

y = zeros(100000,1); 
for i = 1:length(y)
    y(i) = mean(rand(100,1)); 
end 

figure
subplot(2,1,1) 
    histogram(x, 'Normalization', 'probability') 
    title('samples from uniform distribution') 
    ylabel('pdf units') 
subplot(2,1,2) 
    histogram(y, 'Normalization', 'probability') 
    title('average of 100 uniform random variables')
    ylabel('pdf units') 
    xlabel('value of random variable') 
sgtitle('100,000 random variables') 

%% 

x = rand(10000,1);

figure
histogram(x,'Normalization','probability'); 