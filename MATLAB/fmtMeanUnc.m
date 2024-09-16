function formatted_str = fmtMeanUnc(mean_value, error_value)

    error_value = ceil(error_value * 10^(ceil(-log10(error_value)))) / (10^(ceil(-log10(error_value))));

    % Extract the first significant digit of the uncertainty
    first_digit_error = floor(error_value * 10^(ceil(-log10(error_value))));

    % Determine the number of significant digits for the mean
    num_significant_digits = ceil(-log10(error_value));

    % Round the mean to the appropriate number of significant digits
    rounded_mean = round(mean_value, num_significant_digits);

    % Format the output string
    formatted_str = sprintf(['%.' num2str(num_significant_digits) 'f(%d)'], rounded_mean, first_digit_error);
end
