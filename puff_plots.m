% clear all
% close all
% TESTING THE CUDA PUFF MODEL  

Calcium = csvread("Calcium.csv");
Time = csvread("Time.csv");
B = csvread("B.csv");

trials = 5; % change when you change the number of trials

BT = 20;

F = BT - B    ;

%%%%%%%%%%%%% DEFINE PUFF CUTOFF %%%%%%%%%%%%%

cutoff = 8; %(F-F0)/F0 = 3


for i = 1:length(F)
    for j = 1:trials
        if F(i,j) < cutoff
            F2(i,j) = 0;
        else
            F2(i,j) = F(i,j);
        end
    end
end

%%%%%%%%%%%%% PUFF DURATION %%%%%%%%%%%%%

%Find where it starts to go above 2 and below 2
%If puff peak indices is between those points count it
%Take difference between start and stop time

dur_start_i = zeros(500, trials);
dur_start_t = zeros(500,trials);
dur_end_i = zeros(500,trials);
dur_end_t = zeros(500,trials);


for x = 1:trials
    count = 1;
    count2 = 1;
    for y = 2:length(F2)
        if F2((y-1), x) < 2 && F2(y,x) > 2
            dur_start_i(count,x) = y;
            dur_start_t(count,x) = Time(y,x);
            count = count + 1;
        end
        
        if F2(y,x) < 2 && F2((y-1),x) > 2
            dur_end_i(count2,x) = y;
            dur_end_t(count2,x) = Time(y,x);
            count2 = count2 + 1;
        end
    end
end


if length(dur_start_i) ~= length(dur_end_i)
    
    if length(dur_start_i) > length(dur_end_i)
        dur_start_i = dur_start_i(:,1:length(dur_end_i));
        dur_start_t = dur_start_t(:,1:length(dur_end_t));
    elseif length(dur_end_i) > length(dur_start_i)
        dur_end_i = dur_end_i(:,2:length(dur_start_i)+1);
        dur_end_t = dur_end_t(:,1:length(dur_start_t));
    end

end


%Calculate puff durations

puffdur = zeros(500,trials);

for p = 1:trials
    for m = 1:length(dur_end_t) %go thru all start/end indices
        if dur_end_t(m,p) ~= 0 && dur_start_t(m,p) ~= 0
        puffdur(m,p) = dur_end_t(m,p) - dur_start_t(m,p);
        end
    end 
end

pufftimes = zeros(500,trials);
puffamps  = zeros(500,trials);
for i = 1:trials
   for jj = 1:length(dur_end_i)
       if dur_end_i(jj,i) ~= 0 && dur_start_i(jj,i) ~= 0
           A = F2(dur_start_i(jj,i):dur_end_i(jj,i),i);
           [M, Ind] = max(A);
           disp(M)
           puffamps(jj,i) = M;
           loc = dur_start_i(jj,i) + Ind;
           pufftimes(jj,i) = Time(loc,i);
       end
   end
end


%%%%%%%%%%%% MEAN INTER-PUFF INTERVAL %%%%%%%%%%%%%

IPI = zeros(500,trials);
for j = 1:trials
    for k = 2:length(pufftimes)
        IPI(k-1,j) = pufftimes(k,j)-pufftimes(k-1,j);
        if IPI(k-1,j) <= 0
            IPI(k-1,j) = NaN;
        end
    end
end
meanIPI = nanmean(IPI);


%%%%%%%%%%%%%%%%%% PLOT %%%%%%%%%%%%%%%%%%%%


% Cleaving pre-allocated zeros that will affect plots

for q = 2:length(Time)
    for r = 1:trials
        if F(q,r) == 20
            F(q,r) = NaN;
            Time(q,r) = NaN;
        end
    end
end

for s = 2:length(puffdur)
    for t = 1:trials
        if puffdur(s,t) == 0
            puffdur(s,t) = NaN;
        end
    end
end

for u = 2:length(puffamps)
    for v = 1:trials
        if puffamps(u,v) == 0
            puffamps(u,v) = NaN;
        end
    end
end

%  Fluorescence Time Course (subplots for each trial)
figure
hold on
for x = 1:trials
    subplot(2,3,x); % change dimensions based on how many trials are run in parallel
    plot(Time(:,x), F(:,x));
    title(['Trial: ',num2str(x)]);
    xlabel("Time in ms");
    ylabel("F/F0");
end


hold off


% %Plot Fluorescence Time Course
figure
hold on
for c = 1:trials
    plot(Time(:,c), F(:,c));
end
plot([0 100],[cutoff cutoff],'r-')
plot([0 100],[2 2],'b-')
ylabel('F/F0');
xlabel('Time');
title('F/F0 over Time');
hold off

mn = mean(meanIPI);
%Plot IPI histogram
figure
hold on
histogram(IPI);
xlabel('IPI')
ylabel('Occurences')
title(['Mean Inter-Puff Interval = ',num2str(mn)])

%Plot peak amplitude histogram
figure
hold on
histogram(puffamps,10);
xlabel('Puff Amplitude')
ylabel('Occurences')

%Plot histogram of puff durations
figure
hold on
histogram(puffdur);
xlabel('Puff Duration')
ylabel('Occurences')
title("Duration");

% 
% 
