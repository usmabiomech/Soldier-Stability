function [shortenedMatrix] = zeroKillerChicken(largeMatrix)
i=1;
shortenedMatrix=[];
while i<=length(largeMatrix(:,1))
    if largeMatrix(i,1)~=0
        i=i+1000;
    elseif largeMatrix(i+1)==0
                start=i-1000;
                break
    end
end
i=start;
while i<=length(largeMatrix(:,1))
    if largeMatrix(i,1)~=0
        i=i+1;
    elseif largeMatrix(i,1)==0
                shortenedMatrix=largeMatrix(1:i-1,:);
                break
    end
end
end