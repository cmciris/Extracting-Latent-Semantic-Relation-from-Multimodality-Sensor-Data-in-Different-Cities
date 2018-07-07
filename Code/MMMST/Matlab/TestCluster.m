function TestCluster(clusterDirectory,timeslot)
    dict1Filename = strcat(clusterDirectory,num2str(timeslot,'%.2d'),'\00.txt');
	dict2Filename = strcat(clusterDirectory,num2str(timeslot,'%.2d'),'\01.txt');
	dict3Filename = strcat(clusterDirectory,num2str(timeslot,'%.2d'),'\02.txt');
	dict4Filename = strcat(clusterDirectory,num2str(timeslot,'%.2d'),'\03.txt');

	D1 = load(dict1Filename); 
	D2 = load(dict2Filename); 
	D3 = load(dict3Filename); 
	D4 = load(dict4Filename); 
    
    figure(1)
    X = D1(:,1:end-3);
    Y = D1(:,end-1)+1;
    stem(IndFeat(X,Y))
    grid on
    
    figure(2)
    X = D2(:,1:end-3);
    Y = D2(:,end-1)+1;
    stem(IndFeat(X,Y))
    grid on
    
    figure(3)
    X = D3(:,1:end-3);
    Y = D3(:,end-1)+1;
    stem(IndFeat(X,Y))
    grid on
    
    figure(4)
    X = D4(:,1:end-3);
    Y = D4(:,end-1)+1;
    stem(IndFeat(X,Y))
    grid on
    
    T1 = VisualizeClusters(strcat(clusterDirectory,num2str(timeslot,'%.2d'),'\00.txt'));
    T2 = VisualizeClusters(strcat(clusterDirectory,num2str(timeslot,'%.2d'),'\01.txt'));
    T3 = VisualizeClusters(strcat(clusterDirectory,num2str(timeslot,'%.2d'),'\02.txt'));
    T4 = VisualizeClusters(strcat(clusterDirectory,num2str(timeslot,'%.2d'),'\03.txt'));
    disp(size(T1,1));
    disp(size(T3,1));
end