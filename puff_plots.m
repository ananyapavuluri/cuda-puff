% clear all
% close all
% TESTING THE CUDA PUFF MODEL  

Calcium = csvread("Calcium.csv");
Time = csvread("Time.csv");
B = csvread("B.csv");

trials = 5; % change when you change the number of trials

%%%%%%%% DEANALISA's CODE (modified for 2D data matrices)%%%%%%%%%%%

BT = 20;

F = BT - B    ;

%%%%%%%%%%%%% DEFINE PUFF CUTOFF %%%%%%%%%%%%%

cutoff = 8; %(F-F0)/F0 = 3
%F2 = F';

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

% for l = 2:length(F)
%     
%     if F(l-1) < 2 && F(l) > 2   
%         dur_start_i = [dur_start_i l];  
%         dur_start_t = [dur_start_t totaltime(l)];
%     end
%         
%     if F(l) < 2 && F(l-1) > 2
%         dur_end_i = [dur_end_i l];
%         dur_end_t = [dur_end_t totaltime(l)];
%     end
% 
% end

%Make sure they are the same length

% if dur_end_t(y,x) < dur_start_t(y,x) %The first start should be less than the first end
%     dur_end_t = dur_end_t(2:length(dur_end_t),x); %Cleave first end value   
% end

if length(dur_start_i) ~= length(dur_end_i)
    
    if length(dur_start_i) > length(dur_end_i)
        dur_start_i = dur_start_i(:,1:length(dur_end_i));
        dur_start_t = dur_start_t(:,1:length(dur_end_t));
    elseif length(dur_end_i) > length(dur_start_i)
        dur_end_i = dur_end_i(:,2:length(dur_start_i)+1);
        dur_end_t = dur_end_t(:,1:length(dur_start_t));
    end

end








%Calculate puff duration

puffdur = zeros(500,trials);
% [pks2,locs2] = findpeaks(pks);
% [pks3,locs3] = findpeaks(pks2);
% [pks4,locs4] = findpeaks(pks3);
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

%%%%%%%%%%%%%%%%% CLEAVE %%%%%%%%%%%%%%%%%%

%Only want to characterize puffs that we have all info for
%The limiting variable is puff duration because a part of a puff may be
%simulated i.e. there is a start but no end to calculate duration

% pufftime = pufftime(:,1:length(puffdur));
% pks = pks(1:length(puffdur),:);

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

% Finding the time taken by each trial
% 
% trialtimes = zeros(trials);
% for 

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

%  Calcium Time Course
figure
hold on
for x = 1:trials
    subplot(2,3,x);
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
% plot(totaltime,F(:,1:length(totaltime)),'k')
% plot(Time,F2)
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
% pd = fitdist(IPI','Normal');
% x_values = 0:0.1:30;
% y = pdf(pd,x_values);
%plot(x_values,y,'LineWidth',3)
%disp(IPI)
xlabel('IPI')
ylabel('Occurences')
title(['Mean Inter-Puff Interval = ',num2str(mn)])

%Plot peak amplitude histogram
figure
hold on
histogram(puffamps,10);
% pd = fitdist(puffamps,'Normal');
% x_values = 0:0.1:18;
% y = pdf(pd,x_values);
% plot(x_values,y,'LineWidth',3)
% meanpeakamp = mean(pks);
% stdpeakamp = std(pks);
%title(['Mean Puff Amplitude = ',num2str(meanpeakamp),' +/- ',num2str(stdpeakamp)])
xlabel('Puff Amplitude')
ylabel('Occurences')

%Plot histogram of puff durations

figure
hold on
histogram(puffdur);
% pd = fitdist(puffdur','Normal');
% x_values = 0:0.01:1.5;
% y = pdf(pd,x_values);
% plot(x_values,y,'LineWidth',3)
xlabel('Puff Duration')
ylabel('Occurences')
title("Duration");

% 
% 
