function [accuracy] = analyses(x,N,k,true_labels)

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

unique_labels = unique(true_labels);
num_labels = length(unique_labels);

% Initialize an array to hold the count of correct predictions for each label-column pair
counts = zeros(num_labels, num_labels);

% Count the correct predictions for each label-column pair
for i = 1:num_labels
    for j = 1:num_labels
        counts(i, j) = sum(Y(:, j) == 1 & strcmp(true_labels, unique_labels{i}));
    end
end

% Initialize a mapping from columns to labels
col_to_label = cell(num_labels, 1);
used_labels = false(num_labels, 1);  % Track which labels have been used
used_columns = false(num_labels, 1); % Track which columns have been used

% Flatten the counts matrix and sort the indices by descending order of counts
[~, sorted_indices] = sort(counts(:), 'descend');

% Determine the best matching columns for each label
for idx = 1:length(sorted_indices)
    [i, j] = ind2sub(size(counts), sorted_indices(idx));
    if ~used_labels(i) && ~used_columns(j)
        col_to_label{j} = unique_labels{i};
        used_labels(i) = true;
        used_columns(j) = true;
        % Set the counts for this row and column to zero to avoid repeated assignments
        counts(i, :) = 0;
        counts(:, j) = 0;
    end
end

% Initialize the converted predicted labels
N = length(true_labels);
predicted_labels = cell(N, 1);

% Convert prediction matrix to class labels using the determined column indices
for i = 1:N
    [~, max_index] = max(Y(i, :));
    predicted_labels{i} = col_to_label{max_index};
end

% Calculate the number of correct predictions
num_correct = sum(strcmp(predicted_labels, true_labels));

% Calculate the accuracy percentage
accuracy = (num_correct / N) * 100;