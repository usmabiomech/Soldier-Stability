function [shortenedMatrix] = zeroKillerChicken2(largeMatrix)
ii=1;
shortenedMatrix=[];
start=1;
while ii<=length(largeMatrix)
    if largeMatrix(1,ii)~=0
        ii=ii+10;
    elseif largeMatrix(1,ii+1)==0
                start=ii-10;
                if start<=0
                    start=1;
                end
                break
    end
end
ii=start;
while ii<=length(largeMatrix)
    if largeMatrix(1,ii)~=0
        ii=ii+1;
    elseif largeMatrix(1,ii)==0
                shortenedMatrix=largeMatrix(1,1:ii-1);
                break
    end
end
end