function u = f_rec(t, u_0, r)
    if (t == 0)
        u = u_0;
    else
        u = r .* f_rec(t-1, u_0, r) .* (4 - f_rec(t-1, u_0, r).^2);
    end
end

