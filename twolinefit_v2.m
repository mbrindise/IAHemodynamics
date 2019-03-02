function [bestind,lineinfo] = twolinefit_v2(y)

% Get signal info and initialze output
if size(y,2) == 1
    y = y';
end
numy = length(y);
ersum = zeros(numy,1);
r2 = ersum;

% Compute best two line fit
for q = 1:1:numy-2
    if q == 86
        w = 1;
    end
    if q == 1
        x1 = q+1:numy;
        p1 = polyfit(x1,y(x1),1);
        yfit1 = polyval(p1,x1);
        ersum(q) = sum(abs(yfit1-y(x1)));
        r2(q,1) = 1-sum((y(x1)-yfit1).^2)/sum((y(x1)-mean(y(x1))).^2);
    elseif q == numy
        x1 = 1:q-1;
        y1 = y(1:q-1);
        
        p1 = polyfit(x1,y(x1),1);
        yfit1 = polyval(p1,x1);
        ersum(q) = sum(abs(yfit1-y(x1)));
        r2(q,1) = 1-sum((y(x1)-yfit1).^2)/sum((y(x1)-mean(y(x1))).^2);
    else
        x1 = 1:q;
        x2 = q+1:numy;
        
        p1 = polyfit(x1,y(x1),1);
        yfit1 = polyval(p1,x1);
        er1(q,1) = sum(abs(yfit1-y(x1)));
        r2_1(q,1) = 1-sum((y(x1)-yfit1).^2)/sum((y(x1)-mean(y(x1))).^2);
        
        p2 = polyfit(x2,y(x2),1);
        yfit2 = polyval(p2,x2);
        er2(q,1) = sum(abs(yfit2-y(x2)));
        r2_2(q,1) = 1-sum((y(x2)-yfit2).^2)/sum((y(x2)-mean(y(x2))).^2);
        
        ersum(q) = er1(q,1) + er2(q,1);
        r2(q,1) = (r2_1(q,1)+r2_2(q,1))/2;
        
        % For use with RANSAC
        %y1 = y(x1);
        %y2 = y(x2);
        %[p1R_s1,p1R_i1] = ransac([x1;y1],2,500,10,0.90);
        %[p2R_s1,p2R_i1] = ransac([x2;y2],2,500,10,0.90);
        %yfit1_1 = polyval([p1R_s1,p1R_i1],x1);
        %er1_1(q,1) = sum(abs(yfit1_1-y(x1)));
        %r2_1_1(q,1) = 1-sum((y(x1)-yfit1_1).^2)/sum((y(x1)-mean(y(x1))).^2);
        %yfit2_1 = polyval([p2R_s1,p2R_i1],x2);
        %er2_1(q,1) = sum(abs(yfit2_1-y(x2)));
        %r2_2_1(q,1) = 1-sum((y(x2)-yfit2_1).^2)/sum((y(x2)-mean(y(x2))).^2);
        %ersum_1(q) = er1_1(q,1) + er2_1(q,1);
        %r2p1(q,1) = (r2_1_1(q,1)+r2_2_1(q,1))/2;
    end
end

% Determine index with minimal error
[minerr,bestind] = max(r2./ersum);
%[minerr,bestind] = min(ersum);
%[minerr,bestind] = max(r2);

%Determine quality of best ind
[pks,locs] = findpeaks(r2./ersum,'MinPeakDistance',round(length(y)/25));
if length(pks) > 1
    srtpks = sort(pks,'descend');
    met = r2./ersum; met(99) = 0;
    pkrat = srtpks(1)/srtpks(2);
    %pkrat = minerr/mean(met);
    if pkrat < 1.3
        bestind = find(er2(2:end) < 0.1,1,'first') + 1;
        if isempty(bestind)
            bestind = length(y);
        end
        met = r2./ersum;
        minerr = met(bestind);

    end
end



% Obtain Info for Best Line
x1 = bestind+1:numy;
p1 = polyfit(x1,y(x1),1);
yfit1 = polyval(p1,x1);
x2 = 1:bestind;
p2 = polyfit(x2,y(x2),1);
yfit2 = polyval(p2,x2);

% Save Line Info for Output
lineinfo.x1 = x1;
lineinfo.x2 = x2;
lineinfo.fit1 = yfit1;
lineinfo.fit2 = yfit2;
lineinfo.l1 = p1;
lineinfo.l2 = p2;
lineinfo.error = minerr;


end