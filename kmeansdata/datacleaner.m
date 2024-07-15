function [A,true_labels] = datacleaner(problem,N,k)

if ( ...
    strcmp(problem, 'wdbc')    == 1  || ...
    strcmp(problem, 'cloud')   == 1  ||...
    strcmp(problem, 'ecoli')   == 1  ||...
    strcmp(problem, 'iris')    == 1  ||...
    strcmp(problem, 'pima')    == 1  ||...
    strcmp(problem, 'raisin')  == 1  || ...
    strcmp(problem, 'seeds')   == 1  || ...
    strcmp(problem, 'SPECTF')  == 1  || ...
    strcmp(problem, 'thyroid') == 1  ||...
    strcmp(problem, 'wine')    == 1 ... 
    )
    treat_data = true;
elseif ( ...
    strcmp(problem, 'ionosphere')  == 1 || ...
    strcmp(problem, 'parkinsons')  == 1 || ...
    strcmp(problem, 'transfusion') == 1 ||...
    strcmp(problem, 'synthetic') == 1 ...
     )
    treat_data = false;
end

% ==========================================  

if ( strcmp(problem, 'wdbc') == 1 )
    
    formatSpec = '%d %s';
    formatSpec = strcat(formatSpec,repmat('%f', 1, 30));
    
    fileID = fopen('wdbc.data', 'r');

    A = textscan(fileID, formatSpec, 'Delimiter', ',');
    fclose(fileID);

    true_labels = A{:,2};

    A = cell2mat(A(:,3:end));
end

% ==========================================  

if ( strcmp(problem, 'cloud') == 1 )
    
    formatSpec = repmat('%f', 1, 10);
    
    fileID = fopen('cloud.data', 'r');

    A = textscan(fileID, formatSpec, 'Delimiter', ' ','MultipleDelimsAsOne', 1);
    fclose(fileID);

    true_labels = repmat({'A'}, 1024, 1);
    true_labels = [true_labels; repmat({'B'}, 1024, 1)];

    A = cell2mat(A);
 end

% ==========================================  

if ( strcmp(problem, 'ecoli') == 1 )
    
    formatSpec = '%s';
    formatSpec = strcat(formatSpec,repmat('%f', 1, 7),'%s');
    
    fileID = fopen('ecoli.data', 'r');

    A = textscan(fileID, formatSpec, 'Delimiter', ' ','MultipleDelimsAsOne', 1);
    fclose(fileID);

    true_labels = A{:,end};

    A = cell2mat(A(:,2:8));
end

% ==========================================    

if ( strcmp(problem, 'ionosphere') == 1 )

    formatSpec = repmat('%f', 1, 34);
    formatSpec = strcat(formatSpec, '%s');

    fileID = fopen('ionosphere.data', 'r');

    A = textscan(fileID, formatSpec, 'Delimiter', ',');
    fclose(fileID);

    true_labels = A{:,end};

    A = cell2mat(A(:,1:34));
end

% ==========================================  

if ( strcmp(problem, 'iris') == 1 )

    formatSpec = repmat('%f', 1, 4);
    formatSpec = strcat(formatSpec, '%s');

    fileID = fopen('iris.data', 'r');

    A = textscan(fileID, formatSpec, 'Delimiter', ',');
    fclose(fileID);

    true_labels = A{:,5};

    A = cell2mat(A(:,1:4));
end

% ==========================================  

if ( strcmp(problem, 'parkinsons') == 1 )
    
    formatSpec = strcat('%s',repmat('%f', 1, 16),'%s',repmat('%f', 1, 6));

    fileID = fopen('parkinsons.data', 'r');

    A = textscan(fileID, formatSpec, 'Delimiter', ',','HeaderLines', 1);
    fclose(fileID);

    true_labels = A{:,18};

    A{:,18} = [];

    A = cell2mat(A(:,2:end));
end

% ==========================================  

if ( strcmp(problem, 'pima') == 1 )

    formatSpec = repmat('%f', 1, 8);
    formatSpec = strcat(formatSpec, '%s');

    fileID = fopen('pimadiabetes.csv', 'r');

    A = textscan(fileID, formatSpec, 'Delimiter', ',', 'HeaderLines', 1);
    fclose(fileID);

    true_labels = A{:,end};

    A = cell2mat(A(:,1:8));
end

% ==========================================  

if ( strcmp(problem, 'raisin') == 1 )

    formatSpec = '%f %f %f %f %f %f %f %s';

    fileID = fopen('Raisin_Dataset.arff', 'r');

    A = textscan(fileID, formatSpec, 'Delimiter', ',');
    fclose(fileID);

    true_labels = A{:,end};

    A = cell2mat(A(:,1:7));
end

% ==========================================  

if ( strcmp(problem, 'seeds') == 1 )

    formatSpec = strcat(repmat('%f', 1, 7),'%s');

    fileID = fopen('seeds_dataset.txt', 'r');

    A = textscan(fileID, formatSpec, 'Delimiter', ' ', 'MultipleDelimsAsOne', 1);
    fclose(fileID);

    true_labels = A{:,end};

    A = cell2mat(A(:,1:7));

    csvwrite('segmentation.csv', A);
 end

 % ==========================================    

 if ( strcmp(problem, 'SPECTF') == 1 )

    formatSpec = strcat('%s',repmat('%f', 1, 44));

    fileID = fopen('SPECTF.train', 'r');

    A = textscan(fileID, formatSpec, 'Delimiter', ',');
    fclose(fileID);

    true_labels = A{:,1};

    A = cell2mat(A(:,2:end));
 end

 % ==========================================    

 if ( strcmp(problem, 'thyroid') == 1 )

    formatSpec = strcat('%s',repmat('%f', 1, 5));

    fileID = fopen('new-thyroid.data', 'r');

    A = textscan(fileID, formatSpec, 'Delimiter', ',');
    fclose(fileID);

    true_labels = A{:,1};

    A = cell2mat(A(:,2:end));
 end

% ==========================================    

 if ( strcmp(problem, 'transfusion') == 1 )

    formatSpec = strcat(repmat('%f', 1, 4), '%s');

    fileID = fopen('transfusion.data', 'r');

    A = textscan(fileID, formatSpec, 'Delimiter', ',','HeaderLines', 1);
    fclose(fileID);

    true_labels = A{:,end};

    A = cell2mat(A(:,1:4));
 end    

% ==========================================    

if ( strcmp(problem, 'wine') == 1 )

    formatSpec = '%s';
    formatSpec = strcat(formatSpec,repmat('%f', 1, 13));

    fileID = fopen('wine.data', 'r');

    A = textscan(fileID, formatSpec, 'Delimiter', ',');
    fclose(fileID);

    true_labels = A{:,1};

    A{:,1} = [];

    A = cell2mat(A); 
end

% ==========================================    

if ( strcmp(problem, 'synthetic') == 1 )

    rng(2028)

    % Generate data points
    
    a = 1;
    b = 4;
    
    A = [];

    % References points
    centroids = [3, 3; -3, -3; 6, -6];

    points_per_cluster = floor(N / k);
    
    for i = 1:k
        dispersion_factor = ( b - a ) * rand(points_per_cluster,1) + a;
        A = [A; centroids(i, :) + dispersion_factor .* randn(points_per_cluster, 2)];
    end
    
    % Add extra points if N is not exactly divisible by k

    remaining_points = N - points_per_cluster * k;
    if remaining_points > 0
        dispersion_factor = ( b - a ) * rand(remaining_points,1) + a;
        A = [A; centroids(1:remaining_points, :) + dispersion_factor .* randn(remaining_points, 2)];
    end

    true_labels = [];
end

% ==========================================    

if ( treat_data )
    [nrow, ~] = size(A);
    stdA = std(A);
    meanA = mean(A);
    A = (A-repmat(meanA, [nrow,1]));
    A = A./repmat(stdA,[nrow,1]);

    if ( strcmp(problem, 'ionosphere') == 1 )
        A(:,2) = 0;
    end
end    