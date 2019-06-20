%% Fill in the gaps by averaging the positions before and after the gap
%% Original coord needs to be a matrix with two columns. The first and the last
%% row needs to be numbers. Any rows with NaNs are changed with numbers that
%% connects the numbers at boundary.
function filled_coord = fill_gap(coord)

if size(coord, 2) ~=2
    display('The matrix should have two columns.');
    return
end
    
if any(any(isnan(coord(:,1:2))))
    nan_row_id = find(isnan(coord(:,1)));
    num_row_id = find(~isnan(coord(:,1)));
    if any(nan_row_id == 1) || any(nan_row_id == size(coord, 1))
        display('NaN at the first or last time point. Cannot determine the first position.');
        return
    end
    %%
    if numel(nan_row_id) ==1
        row_id = nan_row_id;
        coord(row_id, 1:2) =...
            (coord(row_id-1, 1:2) + coord(row_id +1, 1:2))/2;
    elseif numel(nan_row_id) > 1
        num_row_id_diff = diff(num_row_id) - 1;
        gap_id = find(num_row_id_diff > 0);
        for id = 1:numel(gap_id)
            gap_size = num_row_id_diff(gap_id(id));
            row_id = num_row_id(gap_id(id));
            coord((row_id+1):(row_id + gap_size), 1:2) =...
                ((gap_size:-1:1)'*coord(row_id, 1:2) + (1:gap_size)'*coord(row_id + 1+gap_size, 1:2))/(gap_size+1);
        end
    end
    %%
end
filled_coord = coord;
return
