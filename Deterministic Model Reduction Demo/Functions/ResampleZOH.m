function Y = ResampleZOH(X, T, t_vector)
    Y = zeros(size(X,1),length(t_vector));
    j_Previous = 1;
    for i = 1 : length(t_vector)
        t_i = t_vector(i);
        for j = j_Previous : length(T)-1
            if (t_i >= T(j)) && (t_i < T(j+1))
                Y(:,i) = X(:,j);
                j_Previous = j;
                break;
            end
        end
    end
end

