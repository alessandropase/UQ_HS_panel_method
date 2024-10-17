function F = uq_VMCircular_cdf(X, parameters)
    % Calcola la CDF della distribuzione di von Mises
    % theta: vettore di angoli
    % mu: direzione media (in radianti)
    % kappa: parametro di concentrazione (controlla la diffusione della distribuzione)

    % Calcola la PDF della distribuzione di von Mises
    pdf_values = circ_vmpdf(X, parameters(1), parameters(2));

    % Calcola la CDF cumulativamente
    F = cumsum(pdf_values) * (2*pi) / numel(X); % Normalizzazione per la lunghezza dell'intervallo angolare
end
