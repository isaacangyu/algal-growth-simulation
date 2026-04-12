function [l] = plot(mu, data)
    tt = 
    data = 
    delta_t = 0.416;
    B = @(t) B(t - delta_t) * exp(mu * delta_t);

    % calculate and return loss
    l = data - B;

end