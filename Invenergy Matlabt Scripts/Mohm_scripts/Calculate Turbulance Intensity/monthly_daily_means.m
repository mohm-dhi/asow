function [weight_total,weight_months,weight_days] = monthly_daily_means(time_vector,data_per_hour)


%MONTHLY MEANS AND DAILY MEANS ADD WEIGTHING 
%JABE 7/12/2020

data_per_day=data_per_hour*24;

days_in_months=[[31 28 31 30 31 30 31 31 30 31 30 31] ; [31 29 31 30 31 30 31 31 30 31 30 31]];

time=time_vector;
len=length(time);

weight_total=zeros(len,1);
weight_days=zeros(len,1);
weight_months=zeros(len,1);



years=unique(year(time));

ly=leapyears(years);

monthly_num=zeros(1,12);
for iy=1:length(years)
    monthly_num=monthly_num+days_in_months(ly(iy)+1,:)/length(years);
end
    
%loop over months
for im=1:12
    ind_m=find(month(time)==im);
    len_ind=length(ind_m);
    data_per_month=data_per_day*monthly_num(im);
    weight_m(im)=data_per_month/len_ind;
    weight_months(ind_m)=weight_m(im);
    
end
weight_months=weight_months/sum(weight_months)*len;

%loop over hours
for ih=1:24
    ind_h=find(hour(time)==ih-1);
    len_ind=length(ind_h);
    weight_h(ih)=data_per_hour/len_ind;
    weight_days(ind_h)=weight_h(ih);
        
end
weight_days=weight_days/sum(weight_days)*len;


weight_total=weight_months.*weight_days;

end



function out = leapyears(YEARS)
    list=zeros(1,length(YEARS));
    x = rem(YEARS,[4,100,400]);
    true = find(x(:,1)==0 & x(:,2) | x(:,3) == 0);
    list(true)=1;
    out = list;
end
