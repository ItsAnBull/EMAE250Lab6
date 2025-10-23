function x = ZhouK_Lab6(A, b)

% display A and b
disp(A);
disp(b);

% calculate the number of rows of the matrix
[rows, ~] = size(A);

% check to see if the system is a first order system
if (rows == 1)

    % if the system is a first order system, solve for x
    x = [A \ b];

else

    % concatenate the b matrix to the end of the A matrix
    Astar = horzcat(A, b);

    % update the cols value accordingly
    cols = rows + 1;

    % identify the index of the row to row swap
    [~, pivot_target] = max(abs(Astar(:,1)));

    % perform partial pivoting
    temp = Astar(pivot_target,:);
    Astar(pivot_target,:) = Astar(1,:);
    Astar(1,:) = temp;

    % perform row ops for each row of the matrix
    for i = 2:rows

        % calculate the row factor
        fac = Astar(i,1) / Astar(1,1);

        % perform the elimination
        Astar(i,:) = Astar(i,:) - fac*Astar(1,:);

    end

    % define the child matrix
    child_Astar = Astar(2:rows,2:cols);
    
    % define the x_child vector as the solution vector x for the new A and 
    % b vectors to be passed into the child
    A_child = child_Astar(:,1:cols-2);
    b_child = child_Astar(:,cols-1);
    x_child = ZhouK_Lab6(A_child, b_child);

    % define the next x term to insert into the solution vector as this
    % math, which is essentially just solving the equation for the next
    % term
    x_calculate = (Astar(1,cols) - dot(Astar(1,1:cols-1),horzcat([0], transpose(x_child)))) / Astar(1,1);
    
    % define the solution vector
    x = [x_calculate; x_child];

end

