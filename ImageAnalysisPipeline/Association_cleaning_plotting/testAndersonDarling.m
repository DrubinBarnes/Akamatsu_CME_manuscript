% Create paramters to call on AnDarksamtest

close all
clear
clc

sizal=10;
sizbet=50;
alpha = 0.1;
X = zeros(sizal*sizbet,2);

adindex=0;
for al=1:sizal
    for bet=1:sizbet
        adindex=adindex+1;
        X(adindex,1)=exp(rand);
        X(adindex,2)=al;
        
    end
end

for bet=1:sizbet
    X(bet+sizbet,1)=X(bet+sizbet,1)+1;
end

AnDarksamtest(X,alpha)