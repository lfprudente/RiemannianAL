function plot_kmeans(x,N,k,data)

if ( size(x,2) == 1 )
    Y = reshape(x, [N, k]);
else
    Y = x;
end

for i = 1:N
    [~,j] = max(Y(i,:));
    Y(i,:) = 0;
    Y(i,j) = 1;
end

centroids = (Y' * data) ./ sum(Y', 2);

figure;
hold on;

colors = {'k','r','b'};
marker = {'o','s','^'};

axis off

for j = 1:k
    cluster_points = data(Y(:, j) == 1, :);
    scatter(cluster_points(:, 1), cluster_points(:, 2), 36, colors{j},marker{j},'filled');
    %scatter(centroids(j, 1), centroids(j, 2), 300, 'y', 'pentagram', 'filled');
end