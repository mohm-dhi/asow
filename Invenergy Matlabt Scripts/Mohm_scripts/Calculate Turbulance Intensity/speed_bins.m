function [indices,lim] = speed_bins(data,bins)
%jbej 18032020
%Divide speed in bins


%data = 0 + (20-0).*rand(10000,1); 
%bins=[0,6.8,10,20]


nbins=length(diff(bins));

 
for i = 1:nbins+1
    
    if i == nbins+1
        indices{nbins+1}=transpose([1:1:length(data)]);
        lim{nbins+1}=strcat('all');
    else
        indices{i} = find(bins(i) < data & data < bins(i+1) );
        lim{i}=strcat(num2str(bins(i)),'-',num2str(bins(i+1)),'{ }','m/s');
    end
   
%     histogram(mod(data(indices{i}),360),[0:0.1:20])
%     xticks([0:1:20])
%     title(lim{i})
%     i
end


end