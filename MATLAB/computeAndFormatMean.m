function formatted_str = computeAndFormatMean(array)
    % Compute the mean and standard error of the mean (SEM)
    mean_value = mean(array);
    error_value = std(array) / sqrt(length(array));

    % Format the output string
    formatted_str = fmtMeanUnc(mean_value, error_value);
end
