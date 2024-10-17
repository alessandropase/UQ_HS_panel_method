function f = uq_VMCircular_invcdf(F, parameters)
    % Calcola l'inversa della CDF della distribuzione di von Mises
    % p: probabilità cumulativa (0 <= p <= 1)
    % mu: direzione media (in radianti)
    % kappa: parametro di concentrazione (controlla la diffusione della distribuzione)

    % Trova l'angolo corrispondente alla probabilità cumulativa p
    theta_range = linspace(0, 2*pi, 1000); % Intervallo di angoli
    cdf_values = uq_VMCircular_cdf(theta_range, parameters);
    
    % Rimuovi duplicati dai valori della CDF
    [unique_cdf, idx] = unique(cdf_values);
    unique_theta = theta_range(idx);

    % Interpolazione per ottenere un valore approssimato dell'angolo corrispondente a p
    f = interp1(unique_cdf, unique_theta, F, 'linear', 'extrap');
end
