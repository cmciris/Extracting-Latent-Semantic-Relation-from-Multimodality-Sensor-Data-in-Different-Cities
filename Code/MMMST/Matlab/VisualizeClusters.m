%***************Visualize the clustering result****************************
%Created by Ying Wei
%2014/12/27
function [T, values,instances] = VisualizeClusters(clusteringFilename)
%load the data in 
inData = load(clusteringFilename);
data = inData(:,1:end-3);
cluster = inData(:,end-2);
maxLabel = inData(:,end-1);
purity = inData(:,end);
T = table(cluster,maxLabel,purity);
T = unique(T);
values = unique(cluster);
instances = histc(cluster,values);
T = [T table(instances)];
%pca on the data so that the data can be visualized
% [coeff, data] = pca(data);
%start to visualize
% marker = {'o','+','d','*','v','s'};
% color = colormap(lines);
% for i = 1:1%size(inData,1)
% if(maxLabel(i)~=-1)
% markers = marker{maxLabel(i)};
% else
% markers = '^';
% end
% colors = color(cluster(i) * 2 + 1,:);
% scatter(data(i,1),data(i,5),[],markers); %colors,
% hold on;
% end
end