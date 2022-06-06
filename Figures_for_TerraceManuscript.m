%% Code for replicating figures in 
% A multi-proxy assessment of terrace formation in the lower Trinity River 
% Valley, Texas 
% This file contains the code to load and plot plot the data for 
% Figure 4 - 5, 7 - 12 

%% Figure 4 
% Median terrace elevation (A) and detrended elevation (B) with 
%interquartile range versus median terrace latitude (UTM) and 
%interquartile range for the 52 terraces along the N-S trending valley. 

% load elevation and detrrended elevation data for terraces 
% columns include:
% FID, Elevation, TerraceID, GarvinClass, ValleyPlane, DetrendedElevation, 
% XCoordinate, YCoordinate
data_elev = csvread('Terrace_5m_elevation_detrended_elevation.csv',1,0);
% convert -9999 elevations to NaN
data_elev(data_elev(:,2) < -100, 2) = NaN; 
% find the terrace number for each entry
num_terrace = unique(data_elev(:,3)); 

% load valley floor plane 
% columns inlcude: FID,pointid,Elevation,POINT_X,POINT_Y
fid = fopen('Terrace_10m_valley_plane_elevation.csv');
valley_data = textscan(fid,'%*f%*f%f%f%f','Delimiter',',','Headerlines',1);
valley_x = valley_data{:,2};
valley_y = valley_data{:,3};
valley_elev = valley_data{:,1};
fclose(fid);

% load data previously run fitting results 
% See plane fitting code 
% Model_Run	Terrace_Level	Param_A	Param_B	Param_C	Param_A_Error	
% Param_B_Error	Param_C_Error	Stirke_deg	Strike_Error_deg	
% Dip_deg	Dip_Error_deg	RMSE	Num_Terraces_in_Group
data_fitting = csvread('Trinity_Terraces_Plane_Fitting_Output.csv',1,0);

low_terrace = data_fitting(1,:);
intermediate_terrace = data_fitting(2,:);
high_terrace = data_fitting(3,:);

% find elevation statistics for each terrace 
for i = 1:length(num_terrace)

    terrace_data = data_elev(data_elev(:,3)==num_terrace(i),:);
    % Find the unique id for each terrace
    terrace_id(i) = terrace_data(1,3);
    % Find the terrace classification based on Garvin 2008
    gavin_class_terrace(i) = terrace_data(1,4);
    % Find the median elevation for each terrace
    elevation_median(i) = quantile(terrace_data(:,2),.5);
    % Find the 25th percent quantile elevation for each terrace
    elevation_1st(i) = quantile(terrace_data(:,2),.25);
    % Find the 75th percent quantile elevation for each terrace
    elevation_3rd(i) = quantile(terrace_data(:,2),.75);
    % Find the median detrended elevation for each terrace
    detrendelevation_median(i) = quantile(terrace_data(:,6),.5);
    % Find the 25th percent quantile detrended elevation for each terrace
    detrendelevation_1st(i) =quantile(terrace_data(:,6),.25);
    % Find the 75th percent quantile detrended elevation for each terrace
    detrendelevation_3rd(i) = quantile(terrace_data(:,6),.75);
    % Find the median x coordinate for each terrace
    X_median(i) = quantile(terrace_data(:,7),.5);
    % Find the 25th percent quantile X coordinate for each terrace
    X_1st(i) = quantile(terrace_data(:,7),.25);
    % Find the 75th percent quantile X coordinate for each terrace
    X_3rd(i) = quantile(terrace_data(:,7),.75);
    % Find the median Y coordinate for each terrace
    Y_median(i) = quantile(terrace_data(:,8),.5);
    % Find the 25th percent quantile Y coordinate for each terrace
    Y_1st(i) =  quantile(terrace_data(:,8),.25);
    % Find the 75th percent quantile Y coordinate for each terrace
    Y_3rd(i) =  quantile(terrace_data(:,8),.75);
end


% define color scheme used 
color_terrace =[178, 178, 178;
                224, 243, 219;
                168, 221, 118;
                67, 162, 202;
                ]./255;
            
color_valley = [42,181,115]./255;

% set up figure size characteristics
fig= figure(1);
fig.Position = [100, 100, 1120,  420];

% plot the terrace median Y coordinate vs elevation
subplot(1,4,1:2)

hold on

% valley slope infromation from best-fit
valley_abc = [3.15E-05	0.000294437	-983.462052];

% define the range of locations of valley lidar points for which
% the best-fit will be plotted
range_y = [min(valley_y),max(valley_y)];
range_x = [median(valley_x)-1000,median(valley_x)+1000];
[range_X,range_Y]=meshgrid(range_x, range_y);

%add valley surface to plot 
surf(range_X,range_Y, [range_X*valley_abc(1)+range_Y*valley_abc(2)+ ...
     valley_abc(3)],'FaceColor',color_valley(1,:),'EdgeColor','none');

% add slopes of high, intermediate, low Deweyville with Z = aX + bY + c.
surf(range_X,range_Y, [range_X*low_terrace(3)+range_Y*low_terrace(4)+ ...
    low_terrace(5)],'FaceColor',color_terrace(2,:),'EdgeColor','none');
surf(range_X,range_Y, [range_X*intermediate_terrace(3)+...
    range_Y*intermediate_terrace(4)+ ...
    intermediate_terrace(5)],'FaceColor',color_terrace(3,:),...
    'EdgeColor','none');
surf(range_X,range_Y,[range_X*high_terrace(3)+range_Y*high_terrace(4)+...
    high_terrace(5)],'FaceColor',color_terrace(4,:),'EdgeColor','none');

% now draw the vertical errorbar for each point

% loop through all terrace median points 
for i=1:length(X_median)
        
        % create a vector that repeats the median twice
        xV = [X_median(i); X_median(i)];
        yV = [Y_median(i); Y_median(i)];
        zV = [elevation_median(i); elevation_median(i)];
        
        % create a vector that ranges around 
        yV_cap = [Y_median(i)-500; Y_median(i)+500];
        xV_cap = [X_median(i)-500; X_median(i)+500];
        zV_cap = [elevation_median(i)-0.25; elevation_median(i)+0.25];

        xMin = X_1st(i);
        xMax = X_3rd(i);
        yMin = Y_1st(i);
        yMax = Y_3rd(i);
        zMin = elevation_3rd(i);
        zMax = elevation_1st(i);

        xB = [xMin, xMax];
        yB = [yMin, yMax];
        zB = [zMin, zMax];

        % draw error bars
        h=plot3(xV, yV, zB, '-k','Color',[0.55,0.55,0.55]);
        set(h, 'LineWidth', 1);
        %h=plot3(xB, yV, zV, '-k','Color',[0.55,0.55,0.55]);
        %set(h, 'LineWidth', 1);
        h=plot3(xV, yB, zV, '-k','Color',[0.55,0.55,0.55]);
        set(h, 'LineWidth', 1);
        
        % caps on error bars 
        h=plot3(xV, yV_cap, [zMin,zMin], '-r','Color',[0.55,0.55,0.55]);
        set(h, 'LineWidth', 1);
        h=plot3(xV, yV_cap, [zMax,zMax], '-r','Color',[0.55,0.55,0.55]);
        set(h, 'LineWidth', 1);
        h=plot3(xV, [yMin, yMin], zV_cap, '-r','Color',[0.55,0.55,0.55]);
        set(h, 'LineWidth', 1);
        h=plot3(xV, [yMax, yMax], zV_cap, '-r','Color',[0.55,0.55,0.55]);
        set(h, 'LineWidth', 1);
end

%plot only median values for each terrace 
scatter3(X_median,Y_median,elevation_median,20,...
    color_terrace(gavin_class_terrace+1,:),'filled')

%adjust view to along Y coordinates 
view(90,0)

% add x and z lables 
zlabel('elevation (m)')
xlabel('northing (m)')

%plot the terrace Y vs detrend elevation
subplot(1,4,3:4)

hold on
% add error bars that represent the 25th and 75th quantile detrended 
% elevation and Y coordinates
errorbar(Y_median,detrendelevation_median,...
    detrendelevation_median-detrendelevation_1st,...
    detrendelevation_3rd-detrendelevation_median,...
    Y_median-Y_1st,Y_3rd-Y_median,'.','Color',[0.55,0.55,0.55])
%plot the terrace Y vs detrend elevation
scatter(Y_median,detrendelevation_median,20,...
    color_terrace(gavin_class_terrace+1,:),'filled')

% add x and y lables 
xlabel('northing (m)')
ylabel('detrended elevation (m)')


%% Figure 5 
% Distributions of detrended elevations based on the Garvin(2008)
% classification

% Low intermediate, high, unclassified color scheme
color_terrace =[224, 243, 219;
                168, 221, 118;
                67, 162, 202;
                178, 178, 178]./255;

% scaling based on 
% https://clauswilke.com/dataviz/histograms-density-plots.html

% set up figure 
figure(2)
pd_detrendedelev = fitdist(data_elev(:,6), 'Kernel','BandWidth',0.2);
x_detrendedelev = -10:0.1:20;
y_detrendedelev = pdf(pd_detrendedelev,x_detrendedelev);
num_all = length(data_elev(:,6));
ax5 =area(x_detrendedelev,y_detrendedelev,'FaceColor','k',...
    'FaceAlpha',0.6);

% loop through each Deweyville terrace and unclassified data points
for i = [3,2,1,0]
    pd_detrendedelev = fitdist(data_elev(data_elev(:,4)==i,6),...
        'Kernel','BandWidth',0.2);
    %median of each classification
    %median(data_elev(data_elev(:,4)==i,6))
    x_detrendedelev = -10:0.1:20;
    y_detrendedelev = pdf(pd_detrendedelev,x_detrendedelev);
    hold on
    num_comp(i+1) = length(data_elev(data_elev(:,4)==i,6));
    if i == 0 % unclassified 
        ax1=area(x_detrendedelev,...
            y_detrendedelev.*(num_comp(i+1)/num_all),...
            'FaceColor',color_terrace(4,:),'FaceAlpha',.8);
    elseif i == 1 % low Deweyville
        ax2 =area(x_detrendedelev,...
            y_detrendedelev.*(num_comp(i+1)/num_all),...
            'FaceColor',color_terrace(i,:),'FaceAlpha',.8);
    elseif i ==2 % intermediate Deweyville 
        ax3 =area(x_detrendedelev,...
            y_detrendedelev.*(num_comp(i+1)/num_all),...
            'FaceColor',color_terrace(i,:),'FaceAlpha',.6);
    elseif i  == 3 % high Deweyville 
        ax4 =area(x_detrendedelev,...
            y_detrendedelev.*(num_comp(i+1)/num_all),...
            'FaceColor',color_terrace(i,:),'FaceAlpha',.6);
    end
end

% set elevation limits on x axis
xlim([-4,14])
xticks([-4:2:14])

legend([ax5,ax1,ax2,ax3,ax4],...
    'all elevation','unclassified','low','intermediate','high')
xlabel('detrended elevation (m)')
ylabel('scaled density')

%% Figure 7
% Root mean square error (RMSE) of a plane fitted to elevation points of 
% terraces previously classified as high Deweyville, intermediate 
% Deweyville and low Deweyville in the Trinity River valley compared to a 
% distribution of RMSE from 150,000 randomly grouped terraces.

% See plane fitting code 
% load data previously run fitting results 
% Model_Run	Terrace_Level	Param_A	Param_B	Param_C	Param_A_Error	
% Param_B_Error	Param_C_Error	Stirke_deg	Strike_Error_deg	
% Dip_deg	Dip_Error_deg	RMSE	Num_Terraces_in_Group
data_fitting = csvread('Trinity_Terraces_Plane_Fitting_Output.csv',1,0);

% find probability of RMSE falling within autogenic distribution 
% low
index_8 = find(data_fitting(:,14)== 8);
% middle 
index_19 = find(data_fitting(:,14) == 19);
% high 
index_22 = find(data_fitting(:,14) == 22);

figure(10) %create a figure to capture output of histogram
% low Deweyville
h_low = histogram(data_fitting(index_8,13),'Normalization','cdf',...
    'BinMethod','auto');
% find first bin edge where value is greater 
terrace_low_bin=find(h_low.BinEdges>= data_fitting(1,13),1); 
% probability of occurance in percent 
terrace_low_prob = h_low.Values(terrace_low_bin)*100;
% 21.4256% that an RMSE of the lower deweyville terraces or less would
% occure

% middle Deweyville
h_middle = histogram(data_fitting(index_19,13),'Normalization','cdf',...
'BinMethod','auto');
% find first bin edge where value is greater 
terrace_middle_bin=find(h_middle.BinEdges>= data_fitting(2,13),1);
% probability of occurance in percent 
terrace_middle_prob = h_middle.Values(terrace_middle_bin)*100; 
% 0.5140% that an RMSE of the middle dewevylille terraces or less would
% occure

% high Deweyville
h_high = histogram(data_fitting(index_22,13),'Normalization','cdf',...
'BinMethod','auto');
% find first bin edge where value is greater 
terrace_high_bin=find(h_high.BinEdges>= data_fitting(3,13),1); 
% probability of occurance in percent 
terrace_high_prob =  h_high.Values(terrace_high_bin)*100; 
%   0.0080 % that an RMSE of the high dewevylille terrace or less would
%  occure


% inter quartile range 
%8 terraces 
[quantile(data_fitting(index_8,13),0.25),...
    quantile(data_fitting(index_8,13),0.5),...
    quantile(data_fitting(index_8,13),0.75)];
%19 terraces
[quantile(data_fitting(index_19,13),0.25),...
    quantile(data_fitting(index_19,13),0.5),...
    quantile(data_fitting(index_19,13),0.75)];
%22 terraces
[quantile(data_fitting(index_22,13),0.25),...
    quantile(data_fitting(index_22,13),0.5),...
    quantile(data_fitting(index_22,13),0.75)];


% plot histogram of RMSE 

figure(3);
set(gcf,'Position',[0,0,1100,350]);

% High Deweyville Terrace
subplot(1,3,3)
% 22 terraces
v = get(gca,'Position');
set(gca,'Position',[v(1) v(2)*1.8 v(3:4)])
histogram(data_fitting(index_22,13),25,'Normalization','count',...
    'FaceColor','k');
hold on
ax = gca;
ax3 = plot([data_fitting(3,13) data_fitting(3,13)],...
    [ax.YLim(1),ax.YLim(2)],'Color',color_terrace(3,:),'LineWidth',1.5); 
axis([0.5,4,0,6000])
xticks([0:0.5:4])
ax.FontSize = 12;
ax.LineWidth = 1;
ax.FontName = 'Arial';
ax.Box = 'off';
text(2.3,5500,{'number of terraces ','in Monte Carlo = 22'})
legend(ax3,'High Deweyville','Location','northoutside')

%Intermediate Deweyville Terrace 
subplot(1,3,2)
v = get(gca,'Position');
set(gca,'Position',[v(1) v(2)*1.8 v(3:4)])
hold on
% 19 terraces
histogram(data_fitting(index_19,13),25,'Normalization','count',...
    'FaceColor','k')
ax=gca;
ax1 = plot([data_fitting(2,13)' data_fitting(2,13)'],[0,6000],...
    'Color',color_terrace(2,:),'LineWidth',1.5); % intermediate
axis([0.5,4,0,6000])
xticks([0:0.5:4])
xlabel({'Root Mean Square Error of terrace points fit to plane (m)',''})
ax.FontSize = 12;
ax.LineWidth = 1;
ax.FontName = 'Arial';
ax.Box = 'off';
text(2.25,5500,{'number of terraces ','in Monte Carlo = 19'})
legend(ax1,'Intermediate Deweyville','Location','northoutside')

%Low Deweyville Terrace
subplot(1,3,1)
v = get(gca,'Position');
set(gca,'Position',[v(1) v(2)*1.8 v(3:4)])
hold on
% 8 terraces
histogram(data_fitting(index_8,13),25,'Normalization','count',...
    'FaceColor','k');
ax=gca;
ax2 = plot([data_fitting(1,13)' data_fitting(1,13)'],[0,6000],...
    'Color',color_terrace(1,:),'LineWidth',1.5); % low
axis([0.5,4,0,6000])
xticks([0:0.5:4])
ax.FontSize = 12;
ax.LineWidth = 1;
ax.FontName = 'Arial';
legend(ax2, 'Low Deweyville','Location','northoutside')
ylabel('number of Monte Carlo runs')
text(2.3,5500,{'number of terraces ','in Monte Carlo = 8'})
ax.Box = 'off';

%% Figure 8A 
% Paleochannel width extraction.

% Go to folder with transect files 
foldername = 'PaleochannelTransectsforFigure7';
cd(foldername)

% find names of all transects 
filenames = dir('**/*.txt');

% plot transects
fig =figure(4);
fig.Position = [100,100,650,500];

% loop thru each transect
for i = 1:length(filenames)
    %open transect file
    hold on
    fid =fopen(fullfile(foldername,filenames(i).name));
    % each file contains distance along transect and elevation
    data_transect = textscan(fid,'%f%f','Headerlines',1,'Delimiter',',');
    fclose(fid);
    % plot transect in an new subplot 
    subplot(12,1,2*i-1:2*i)
    % plot the distance with the center of transect centered around 300m
    % assuming that the largest transect is 600m
    plot(data_transect{:,1}+(600-max(data_transect{:,1}))./2,...
        data_transect{:,2},'k')
    % set the x-axis limits to 600
    xlim([0,600])
    % only label the x-axis ticks for the last subplot
    if i ~= 6
        xticklabels('')
        % only labe the middle subplot with y-axis label
        if i == 3
            ylabel('elevation (m)')
        end
    end
end
xlabel('distance from centerline (m)')
% assign x tick labels of transects such that 300m is 0m from the center of
% the transects
xticklabels([-300,-200,-100,0,100,200,300])

%return to original folder
cd ..
%% Figure 8E and 9

% load data 
fid = fopen('Paleochannel_summary.csv');
%Paleochannel ID	mean width (m)	std width (m)	channel length (m)	
%Radius (m)	Terrace Number	Terrace level	Gavin Classification	
%grain size (Gavin Upper Bar High)(mm)	grain size(Gavin Upper Bar Low)(mm)
%grain size (Gavin Lower Bar High)(mm)	grain size(Gavin Lower Bar Low)(mm)
% number of bends %bends preserved as cut off 
data_paleo = textscan(fid,'%f%f%f%f%f%f%s%f%f%f%f%f%f%f%s',...
    'Delimiter',',',...
    'Headerlines',1);
fclose(fid);

% assign variables to data columns
paleochannel_num = data_paleo{1};
terrace_num = data_paleo{6};
length_terrace_num = length(terrace_num);
paleochannel_width_mean = data_paleo{2};%m
paleochannel_width_std = data_paleo{3};%m
paleochannel_length = data_paleo{4};%m
Gavin_class_channel = data_paleo{8}; 
d50_upbar_urange = data_paleo{9}./1000; %m
d50_upbar_lrange = data_paleo{10}./1000; %m
d50_lbar_urange = data_paleo{11}./1000; %m
d50_lbar_lrange = data_paleo{12}./1000; %m
num_bends = data_paleo{13};
cutoff = data_paleo{14}; %likelihood of cutoff occurance 

% variables 
nu =10^-6;  %m^2/s at 20 degree
p = 0.998;% g/cm^3 water density at 20 degrees
ps = 2.65; % g/cm^3 sediment density
R = ps/p - 1; %reynolds number non dimensional 
g = 9.81; % m/s^2

d50 =[d50_upbar_urange,d50_upbar_lrange,d50_lbar_urange,d50_lbar_lrange];%m

% Wilkerson and Parker 2011
beta = 0.00398; 
n_B = 0.494;
n_B_SE = 0.14;
m_B = 0.269;
m_B_SE =0.031;

% Monte Carlo Variables

num_iter = 50000; %for error analysis 

% find the elevation info correlated with the terrace that the paleochannel
% is sitting on
for i = 1:length_terrace_num
    indx = find(terrace_id == terrace_num(i));
    if ~isempty(indx)
        terrace_delev(i) = detrendelevation_median(indx); % m 
        terrace_delev_1st(i) = detrendelevation_1st(indx); % m 
        terrace_delev_3rd(i) = detrendelevation_3rd(indx); % m 
    else 
        terrace_delev(i) = NaN; % m 
        terrace_delev_1st(i) = NaN; % m 
        terrace_delev_3rd(i) = NaN; % m 
    end
end

% paleodischarge equation from dimensional bankfull equation in Wilkerson 
% and Parker (2011) 22a which is missing a zero. Bbf is bankful width
% 
% Bbf = 0.00398*(sqrt(R)/nu)^0.494*g^(-0.0875)*Qbf^0.669*D_50^0.0685;
% rearranging for Qbf =
% ((g^(1/5)*Bbf*(1/0.00398)*((sqrt(R*g*D_50)*D_50)/nu)^(-n_B))
% Bbf^(1/0.669).*...
% (0.00398*(sqrt(R)/nu)^0.494*g^(-0.0875)*D_50^0.0685)^(-1/0.669)

% paleodischarge calculation 

% find bankful width from random normal distribution
rng('shuffle')    
Bbf_rand=paleochannel_width_std.'.*randn(num_iter,length_terrace_num)+...
    paleochannel_width_mean.';

% remove Bbf_rand that are less than 0
for i = 1:length_terrace_num
    while sum(Bbf_rand(:,i)<=0)
        Bbf_zero = find(Bbf_rand(:,i)<=0);
        Bbf_rand(Bbf_zero,i)=paleochannel_width_std(i).'.*...
            randn(length(Bbf_zero),1)+paleochannel_width_mean(i).';
    end
end

% find n_B from random normal distribution
rng('shuffle')    
n_B_rand=repmat(n_B_SE,1,length_terrace_num).*...
    randn(num_iter,length_terrace_num)+repmat(n_B,1,length_terrace_num);

% find m_B from random normal distribution
rng('shuffle')    
m_B_rand=repmat(m_B_SE,1,length_terrace_num).*...
    randn(num_iter,length_terrace_num)+repmat(m_B,1,length_terrace_num);

% find D_50 from uniform distibution
rng('shuffle')
D_50_rand=min(d50,[],2).'+range(d50.').*rand(num_iter,length_terrace_num);

% calculate Qbf based on random values 
% bankful discharge for for each terrace and each iteration
Qbf = zeros(num_iter,length_terrace_num); 
Qbf = ...
    ((Bbf_rand.*g^(1/5))./0.00398.*...
    (sqrt(R*g.*D_50_rand).*D_50_rand./nu).^(-n_B_rand).*...
    (sqrt(g.*D_50_rand).*D_50_rand.^2).^m_B_rand).^(5./(5.*m_B_rand+2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 7E plot discharge versus width 
fig2 = figure(5);
fig2.Position = [0,0,300,300];%[0,0,800,600];
ax5=errorbar(paleochannel_width_mean,median(Qbf,1),median(Qbf,1)-...
    quantile(Qbf,0.25),quantile(Qbf,0.75)-median(Qbf,1),'o',...
    'Color',[0.75,0.75,0.75],'MarkerSize',7,...
    'MarkerEdgeColor',[0.55 .55 .55],'MarkerFaceColor',[67,162,202]./250); 

hold on 

xlabel('width (m)')
ylabel('paleodischarge (m^3/s)')

ax = gca;
ax.FontSize = 8;
ax.FontName = 'Arial';
ax.Position = [ax.Position(1:2),0.8,0.8];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 9 plot paleodischarge vs detrended elevation
x = terrace_delev;
y = median(Qbf,1);
err_neg = terrace_delev_1st-terrace_delev;
err_pos = terrace_delev_3rd-terrace_delev;

%color scheme 
opacity = [0:range(paleochannel_length)]./range(paleochannel_length);
cmap1 =[flipud([235, 245, 250]./255.*opacity.'+...
    (1-opacity.').*[67,162,202]./250);...
    [26, 79, 101]./255.*opacity.'+(1-opacity.').*[67,162,202]./250];

fig =figure(6);
fig.Position = [0,0,800,600];
hold on
colormap(cmap1);
ax = gca;
errorbar(x,y,y-quantile(Qbf,0.25),quantile(Qbf,0.75)-y,...
    err_neg,err_pos,'o','Color',[0.75,0.75,0.75],'MarkerSize',1,...
    'MarkerEdgeColor','blue','MarkerFaceColor','blue');

scatter(x,y,30,paleochannel_length./1000,'filled',...
'MarkerEdgeColor',[0.55 .55 .55]);

xlabel('elevation above modern valley floor (m)')
ylabel('paleo-discharge (m^3/s)')
ax.FontSize = 16;
ax.FontName = 'Arial';
set(ax,'YScale','log')

colormap(cmap1);
c = colorbar;
ylabel(c,'paleo-channel length (km)')

%% Figure 10
% Distributions using elevation (A) and paleo-channels (B) that support 
% interpretations of an allogenic forcing in terrace abandonment due to 
% decrease in discharge (B) and two individual distributions of detrended 
% elevations (A). 

% Evaluate distributions of elevation

% Variables considered
k = 1:9;
nK = numel(k);
Sigma = {'diagonal','full'};
nSigma = numel(Sigma);
SharedCovariance = {true,false};
SCtext = {'true','false'};
nSC = numel(SharedCovariance);
RegularizationValue = 0.01;
options = statset('MaxIter',1000);

% Preallocation of memory
gm_elev = cell(nK,nSigma,nSC);
aic_elev = zeros(nK,nSigma,nSC);
bic_elev = zeros(nK,nSigma,nSC);
converged_elev = false(nK,nSigma,nSC);

% Detrend elevation data 
% Fit all models. This can take up to 300s.
%tic
for i = 1:nK
    for m = 1:nSC
        for j = 1:nSigma    
            gm_elev{i,j,m} = fitgmdist(data_elev(:,6),k(i),...
                'CovarianceType',Sigma{j},...
                'SharedCovariance',SharedCovariance{m},...
                'RegularizationValue',RegularizationValue,...
                'Options',options);
            aic_elev(i,j,m) = gm_elev{i,j,m}.AIC;
            bic_elev(i,j,m) = gm_elev{i,j,m}.BIC;
            converged_elev(i,j,m) = gm_elev{i,j,m}.Converged;
        end
    end
end
%toc
allConverge_elev = (sum(converged_elev(:)) == nK*nSigma*nSC);

% select the best model
BestModel_elev = gm_elev{2,2,2};

% Paleochannel data
Qbf_median_notnan_temp = [nanmedian(Qbf,1).'];
Qbf_median_notnan_temp = ...
    Qbf_median_notnan_temp(~isnan(Qbf_median_notnan_temp));
Qbf_median_notnan = (Qbf_median_notnan_temp - ...
    mean(Qbf_median_notnan_temp))/std(Qbf_median_notnan_temp);

% Preallocation
gm_channel = cell(nK,nSigma,nSC);
aic_channel = zeros(nK,nSigma,nSC);
bic_channel = zeros(nK,nSigma,nSC);
converged_channel = false(nK,nSigma,nSC);

%Paleodischarge data
% Fit all models
%tic
for i = 1:nK
    for m = 1:nSC
        for j = 1:nSigma    
            gm_channel{i,j,m} = fitgmdist(Qbf_median_notnan,k(i),...
                'CovarianceType',Sigma{j},...
                'SharedCovariance',SharedCovariance{m},...
                'RegularizationValue',RegularizationValue,...
                'Options',options);
            aic_channel(i,j,m) = gm_channel{i,j,m}.AIC;
            bic_channel(i,j,m) = gm_channel{i,j,m}.BIC;
            converged_channel(i,j,m) = gm_channel{i,j,m}.Converged;
        end
    end
end
%toc

allConverge_channel = (sum(converged_channel(:)) == nK*nSigma*nSC);

% find best model
BestModel_channel = gm_channel{3,2,2};

% plot data 
fig7 =figure(7);
fig7.Position=[265   269   994   768];

% A Detrended elevation data Gaussian Mix Model fit
subplot(4,4,[1:2,5:6])
hold on 
ax3 =histogram(data_elev(:,6),'Normalization','pdf','FaceColor','k');
ax3.EdgeColor = 'none';
ax = gca;
ax1 =plot(ax.XLim(1):0.01:ax.XLim(2),...
    pdf(BestModel_elev,[ax.XLim(1):0.01:ax.XLim(2)].'),...
    'Color',[0.2,0.5,0.5],'Linewidth',2);

% loop through reach fit and plot it 
for i = 1:length(BestModel_elev.mu)
       ax2 =plot(ax.XLim(1):0.01:ax.XLim(2),...
           normpdf(ax.XLim(1):0.01:ax.XLim(2),BestModel_elev.mu(i),...
           BestModel_elev.Sigma(i)),'Color',[0.5,0.5,0.5]);
end
xlabel('detrended elevation (m)')
ylabel('pdf')
legend([ax1,ax2,ax3],'Gaussian Mix Model','individual componets',...
    'histogram of data')

%B Paleochannel data gaussian mix model
subplot(4,4,[3:4,7:8])
hold on 
ax3 =histogram((Qbf_median_notnan.'*std(Qbf_median_notnan_temp))+...
    mean(Qbf_median_notnan_temp),'Normalization','pdf',...
    'FaceColor',[67,162,202]./255,'NumBins',50);
ax1 =plot([min(Qbf_median_notnan):0.01:max(Qbf_median_notnan)]*...
    std(Qbf_median_notnan_temp)+mean(Qbf_median_notnan_temp),...
    pdf(BestModel_channel,...
    [min(Qbf_median_notnan):0.01:max(Qbf_median_notnan)].')./1000,...
    'Color',[0.2,0.5,0.5],'Linewidth',2);

% loop through each fit and plot it
for i = 1:3
    A = trapz([min(Qbf_median_notnan):0.01:max(Qbf_median_notnan)]*...
        std(Qbf_median_notnan_temp)+mean(Qbf_median_notnan_temp),...
        normpdf(min(Qbf_median_notnan):0.01:max(Qbf_median_notnan)));
    
    ax2 =plot([min(Qbf_median_notnan):0.01:max(Qbf_median_notnan)]*...
        std(Qbf_median_notnan_temp)+mean(Qbf_median_notnan_temp),...
        normpdf(min(Qbf_median_notnan):0.01:max(Qbf_median_notnan),...
        BestModel_channel.mu(i),BestModel_channel.Sigma(i))./A,...
        'Color',[0.5,0.5,0.5]);
end
xlabel('calculated discharge (m^3/s)')
ylabel('pdf')
legend([ax1,ax2,ax3],'Gaussian Mix Model','individual componets',...
    'histogram of data','Location','northwest')
ylim([0,0.006])

%C AIC results for number of groups in mixing model for elevation
subplot(4,4,[9:10,13:14])
plot(reshape(aic_elev,nK,nSigma*nSC));xlabel('$k$','Interpreter','Latex');
xlabel('$k$','Interpreter','Latex');
ylabel('AIC');
legend({'Diagonal-shared','Full-shared','Diagonal-unshared',...
    'Full-unshared'});
xlim([1,9])

%D AIC results for number of groups in mixing model for paleochannel
subplot(4,4,[11:12,15:16])
plot(reshape(aic_channel,nK,nSigma*nSC));
xlabel('$k$','Interpreter','Latex');
ylabel('AIC');
legend({'Diagonal-shared','Full-shared','Diagonal-unshared',...
    'Full-unshared'});
xlim([1,9])

%% Figure 11

% load terrace length data
%Terrace ID 	GarMS_Clas	median elevation (m)	
% 1st quartile elevation(m)	3rd quartile elevation(m)	
%median detrended elevation (m)	1st quartile detrended elevation(m)	
%3rd quartile detrended elevation(m)	median UTM easting (m)	
%1st quartile UTM easting (m)	3rd quartile UTM easting (m)	
%median UTM northing (m)	1st quartile UTM northing (m)	
%3rd quartile UTM northing (m)	Terrace ID of adjacent higher terraces	
%Terrace ID of lower adjacent terraces	
%minimum bounding box width (m)	minimum bounding box length (m)
fid = fopen('TrinityTerraces_characteristics.csv');
data_terrace_info = textscan(fid,'%f%f%f%f%f%f%f%f%f%f%f%f%f%f%s%s%f%f',...
    'Delimiter',',','Headerlines',1);
fclose(fid);
length_NSenvel = data_terrace_info{:,18};
width_NSenvel = data_terrace_info{:,17};

x = terrace_delev;
err_neg = terrace_delev_1st-terrace_delev;
err_pos = terrace_delev_3rd-terrace_delev;

fig11 = figure(8);
fig11.Position = [0,0,400,900];

% A) plot terrace length vs detrended elevation
subplot(6,1,1:2);
ax=gca;
hold on

%plot only median values for each terrace 
errorbar(detrendelevation_median,length_NSenvel./1000,...
     detrendelevation_median-detrendelevation_1st,...
 detrendelevation_3rd-detrendelevation_median,'horizontal',...
 '.','Color',[0.55,0.55,0.55])

scatter(detrendelevation_median,length_NSenvel./1000,20,'k','filled')

ylabel('terrace length (km)')
ax.FontSize = 10;
ax.FontName = 'Arial';
xlim([-2,16])
ylim([0,14])

% B) plot paleochannel length vs detrended elevation
subplot(6,1,3:4);
ax = gca;
hold on

%plot only median values for each terrace 
errorbar(x,paleochannel_length./1000,...
     err_neg,err_pos,'horizontal','.','Color',[0.55,0.55,0.55])

scatter(x,paleochannel_length./1000,20,'k','filled',...
    'MarkerEdgeColor',[0.55 .55 .55],...
    'MarkerFaceColor',[67,162,202]./250);

ylabel('paleochannel length (km)')
ax.FontSize = 10;
ax.FontName = 'Arial';
xlim([-2,16])

% C) plot paleochannel width vs detrended elevation
subplot(6,1,5:6);
y = paleochannel_width_mean;

ax = gca;
hold on
 errorbar(x,paleochannel_width_mean,...
     paleochannel_width_std,...
     paleochannel_width_std,err_neg,err_pos,'.','Color',[0.55,0.55,0.55])

scatter(x,paleochannel_width_mean,20,'k','filled',...
    'MarkerEdgeColor',[0.55 .55 .55],...
    'MarkerFaceColor',[67,162,202]./250);

xlabel('elevation above modern valley floor (m)')
ylabel('paleochannel width (m)')
ax.FontSize = 10;
ax.FontName = 'Arial';
xlim([-2,16])

%% Figure 12
% Proxy measurements to assess the importance of meander bend-cutoff for
% terrace formation. 

fig9 =figure(9);
fig9.Position = [280   278   643   241];

%A: elevation difference adjacent terrace 

% adjacent terrace dataset if part of the data_terrace_info dataset

% load data from plane fitting
low_terrace = data_fitting(1,:);
low_terrace_slope = tan(data_fitting(1,11))*pi()/180;
low_terrace_slope_error  = tan(data_fitting(1,12))*pi()/180;

intermediate_terrace = data_fitting(2,:);
intermediate_terrace_slope = tan(data_fitting(2,11))*pi()/180;
intermediate_terrace_slope_error  = tan(data_fitting(2,12))*pi()/180;

high_terrace = data_fitting(3,:);
high_terrace_slope = tan(data_fitting(3,11))*pi()/180;
high_terrace_slope_error  = tan(data_fitting(3,12))*pi()/180;

% loop through terraces and find elevation difference to adjacent terrace
for i = 1:length(num_terrace)
    indx =find(data_terrace_info{:,1} ==i-1);
    if ~isempty(indx)
        
        if ~strcmp(data_terrace_info{1,16}(indx), 'NaN')
            lower_terraces = cellfun(@(s) str2double(s),...
                (strsplit(cell2mat(data_terrace_info{1,16}(indx)),';')));
            for k = 1:length(lower_terraces)
                elevation_median_diff_lower(i,k) = elevation_median(i) -...
                    elevation_median(lower_terraces(k)+1);
            end
        end
        if ~strcmp(data_terrace_info{1,15}(indx), 'NaN')
            upper_terraces = cellfun(@(s) str2double(s),...
                (strsplit(cell2mat(data_terrace_info{1,15}(indx)),';')));
            for k = 1:length(upper_terraces)
                elevation_median_diff_upper(i,k) = elevation_median(i)-...
                    elevation_median(upper_terraces(k)+1);
            end
        end
   end
end
% convert zeros to nan
elevation_median_diff_lower(elevation_median_diff_lower == 0) = NaN;
elevation_median_diff_upper(elevation_median_diff_upper == 0) = NaN;

% paleochannel length for high Dewevyille Terraces
high_channel_length =nanmean(paleochannel_length...
    (Gavin_class_channel==3&cutoff==1));
high_channel_length_std = nanstd(paleochannel_length...
    (Gavin_class_channel==3&cutoff==1));
high_channel_length_num = length(paleochannel_length...
    (Gavin_class_channel==3&cutoff==1));

% paleochannel length for intermediate Dewevyille Terraces
intermediate_channel_length = nanmean(paleochannel_length...
    (Gavin_class_channel==2&cutoff==1));
intermediate_channel_length_std =  nanstd(paleochannel_length...
    (Gavin_class_channel==2&cutoff==1));
intermediate_channel_length_num =  length(paleochannel_length(...
    Gavin_class_channel==2&cutoff==1));

% paleochannel length for low Dewevyille Terraces
low_channel_length = nanmean(paleochannel_length...
    (Gavin_class_channel==1&cutoff==1));
low_channel_length_std =  nanstd(paleochannel_length...
    (Gavin_class_channel==1&cutoff==1));
low_channel_length_num =  length(paleochannel_length...
    (Gavin_class_channel==1&cutoff==1));

% elevation drop due to cutoff high Dewevyille Terraces
high_height = high_channel_length*high_terrace_slope;
high_height_error_min = (high_channel_length-high_channel_length_std)*...
    (high_terrace_slope-high_terrace_slope_error);
high_height_error_max = (high_channel_length+high_channel_length_std)*...
    (high_terrace_slope+high_terrace_slope_error);

% elevation drop due to cutoff intermediate Dewevyille Terraces
intermediate_height = intermediate_channel_length*...
    intermediate_terrace_slope;
intermediate_height_error_min = (intermediate_channel_length-...
    intermediate_channel_length_std)*(intermediate_terrace_slope-...
    intermediate_terrace_slope_error);
intermediate_height_error_max = (intermediate_channel_length+...
    intermediate_channel_length_std)*(intermediate_terrace_slope+...
    intermediate_terrace_slope_error);

% elevation drop due to cutoff low Dewevyille Terraces
low_height = low_channel_length*low_terrace_slope;
low_height_error_min = (low_channel_length-low_channel_length_std)*...
    (low_terrace_slope-low_terrace_slope_error);
low_height_error_max = (low_channel_length+low_channel_length_std)*...
    (low_terrace_slope+low_terrace_slope_error);

% plot data
subplot(1,2,1)
histogram(abs([elevation_median_diff_lower]),'FaceColor','k')
ax = gca;
hold on

ax1=area(...
    [low_height_error_min,low_height_error_max,low_height_error_min],...
    [ax.YLim(1),ax.YLim(2),ax.YLim(2)],...
    'FaceColor',color_terrace(1,:),'FaceAlpha',0.5,'EdgeColor','none') ;
ax2=area(...
    [intermediate_height_error_min,intermediate_height_error_max,...
    intermediate_height_error_min],[ax.YLim(1),ax.YLim(2),ax.YLim(2)],...
    'FaceColor',color_terrace(2,:),'FaceAlpha',0.3,'EdgeColor','none') ;
ax3=area(...
    [high_height_error_min,high_height_error_max,high_height_error_min],...
    [ax.YLim(1),ax.YLim(2),ax.YLim(2)],...
    'FaceColor',color_terrace(3,:),'FaceAlpha',0.3,'EdgeColor','none') ;

plot([low_height,low_height],[ax.YLim(1),ax.YLim(2)],...
    'Color',color_terrace(1,:),'LineWidth',1.5);
plot([intermediate_height,intermediate_height],[ax.YLim(1),ax.YLim(2)],...
    'Color',color_terrace(2,:),'LineWidth',1.5);
plot([high_height,high_height],[ax.YLim(1),ax.YLim(2)],...
    'Color',color_terrace(3,:),'LineWidth',1.5) ;

xlabel('absolute elevation difference to adjacent terrace (m)')
ylabel('number of terraces')
yticks([0:2:16])
xticks([0:2:12])
ylim([0,16])
legend([ax1,ax2,ax3],'low','intermediate','high')

%B: number of bend per terrace 
% load data 
fid = fopen('Paleochannel_bends_on_terraces.csv');
% Terrace ID  Gavin classification   Sum of bends on terrace	
% Max number of bends from one generation of channel 	
% Max number of bends from one channel segment 

data_terrace_channel = textscan(fid,'%f%f%f%f%f','Delimiter',',',...
    'Headerlines',1);
fclose(fid);

bend_count_terrace = data_terrace_channel{:,5};

%plot number of bends for longest paleochannel on each terraces 
subplot(1,2,2)
edges = 0:1:14;
h1 = histcounts(bend_count_terrace(gavin_class_terrace.'==1),edges);
h2 = histcounts(bend_count_terrace(gavin_class_terrace.'==2),edges);
h3 = histcounts(bend_count_terrace(gavin_class_terrace.'==3),edges);
h4 = histcounts(bend_count_terrace(gavin_class_terrace.'==0),edges);

b = bar(edges(1:end-1),[h1;h2;h3;h4]','hist');
b(1).FaceColor = color_terrace(1,:);
b(2).FaceColor = color_terrace(2,:);
b(3).FaceColor = color_terrace(3,:);
b(4).FaceColor = color_terrace(4,:);
xlim([-0.5,6.5])
legend('low','intermediate','high','unclassified')
xlabel('number of bends per terrace')
ylabel('number of terraces')