function [avg,st] = CalculateMeanStd()
    city = {'Baoding','Beijing','Shanghai','Shenzhen','Tianjin'};
    all = [];

    for c = 1:1:5
        for ts = 0:1:23
            dir = strcat('D:\Dropbox\Dropbox\Project\MSRA\Result\Features\',city(c),'\labeledSubset2\',num2str(ts,'%.2d'));
            road = load(strcat(dir{1},'\road.txt'));
            poi = load(strcat(dir{1},'\poi.txt'));
            met = load(strcat(dir{1},'\meterology.txt'));
            if(c == 2)
                mob = load(strcat(dir{1},'\mobility.txt'));
            else
                mob = zeros(size(road,1),11);
            end
            all = [all; [road poi met mob]];
        end
    end
    avg = mean(all);
    st = std(all,0,1);
end