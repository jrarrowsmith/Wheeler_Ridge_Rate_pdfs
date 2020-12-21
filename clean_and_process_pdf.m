function [most_common_rate, max_rate, min_rate] = clean_and_process_pdf(h, sigmas)
%funcation taks input historgram of pdf of rates and normalizes it and then
%clips it by the amount in sigmas
bincenters = h.BinEdges(1):h.BinWidth:h.BinEdges(h.NumBins);

wbc=h.BinCounts;
area_of_each_bar=wbc.*h.BinWidth;
orginal_area_under_curve=sum(area_of_each_bar);
normalized_area_of_each_bar = area_of_each_bar./orginal_area_under_curve;
plot(bincenters, normalized_area_of_each_bar,'k-')
na = sum(normalized_area_of_each_bar);

while na>sigmas
    tf=find(wbc>0);
    temp=wbc(tf)-1;
    wbc(tf)=temp;
    area_of_each_bar=wbc.*h.BinWidth;
    area_under_curve=sum(area_of_each_bar);
    normalized_area_of_each_bar = area_of_each_bar./orginal_area_under_curve;
    na=sum(normalized_area_of_each_bar);
end
hold on
fill(bincenters, normalized_area_of_each_bar,[192/255 192/255 192/255])
max_prob=max(normalized_area_of_each_bar);
loc=find(normalized_area_of_each_bar==max_prob);
most_common_rate=bincenters(loc);
plot(most_common_rate,max_prob, 'ro')
tf=normalized_area_of_each_bar>0;
min_rate = min(bincenters(tf));
loc=find(bincenters==min_rate)-1;%this is the first non-zero but we want the last zero
min_rate = bincenters(loc);
plot(min_rate,normalized_area_of_each_bar(loc), 'ro')
max_rate=max(bincenters(tf));
loc=find(bincenters==max_rate)+1;%this is the las non-zero but we want the first zero on the high side
if loc>length(bincenters)
    loc=length(bincenters)
end
max_rate=bincenters(loc);
plot(max_rate,normalized_area_of_each_bar(loc), 'ro')

atxt = sprintf('Most common rate = %0.1f\nsig range = %0.4f\nmin = %0.1f to max = %0.1f',most_common_rate,sigmas,min_rate,max_rate);
gtext(atxt)
ylabel('probability')

end

