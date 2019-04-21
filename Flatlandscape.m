
function [Probabilities, ExpectedValue, Variance] = Flatlandscape(H)
%Function to calculate the Probability, Expected Value, and Variance of a one dimensional lightning leader with a height of H
Probabilities = ones(1,H);
%This initializes the vector associated with the probability values.
E = ones(1,H);
%This initializes a vector to make computing the expected value easier.
E2 = ones(1,H);
%This initializes a vector to make computing the expected value easier.
for i = 0:(H-1)
    P = nchoosek(H-1,i)*(1/2)^(H-1);
    %This calculates the probability at the specific point
    Probabilities(1,i+1) = P;
    %Puts the calculated Probability into the vector that will be output
    E(H-i)= (H-(i))*P;
    %Calculates the expected value
    E2(H-i) = ((H-(i))^2)*P;
    %Used in the calculation for the variance
end
ExpectedValue = sum(E);
%Final calculation for the expected value
Variance = sum(E2) - (ExpectedValue)^2;
%Final Calculation for the Variance
