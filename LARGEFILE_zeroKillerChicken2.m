function [shortenedMatrix] = zeroKillerChicken2(largeMatrix)
i=1;
shortenedMatrix=[];
while i<=length(largeMatrix)
    if largeMatrix(1,i)~=0
        i=i+100;
    elseif largeMatrix(1,i+1)==0
                start=i-100;
                break
    end
end
i=start;
while i<=length(largeMatrix)
    if largeMatrix(1,i)~=0
        i=i+1;
    elseif largeMatrix(1,i)==0
                shortenedMatrix=largeMatrix(1,1:i-1);
                break
    end
end
end