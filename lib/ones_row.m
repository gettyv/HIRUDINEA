function X = ones_row(r,c, row)
%ONES_ROW Summary of this function goes here
%   Detailed explanation goes here
    X = zeros(r, c);
    X(row, :) = 1;
end

