function [shortenedMatrix] = zeroKillerChicken(largeMatrix)
ii=1;
shortenedMatrix=[];
start=1;
while ii<=length(largeMatrix(:,1))
    if largeMatrix(ii,1)~=0
        ii=ii+100;
    elseif largeMatrix(ii+1)==0
                start=ii-100;
                if start<=0
                    start=1;
                end
                break
    end
end
ii=start;
while ii<=length(largeMatrix(:,1))
    if largeMatrix(ii,1)~=0
        ii=ii+1;
    elseif largeMatrix(ii,1)==0
                shortenedMatrix=largeMatrix(1:ii-1,:);
                break
    end
end
end