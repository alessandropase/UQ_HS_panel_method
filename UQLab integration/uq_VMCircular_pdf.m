function F = uq_VMCircular_pdf (X, parameters)
% UQ_VMCIRCULAR_PDF(X, parameters) calculates the Probability Density Function
% values of samples X  that follow a Gaussian distribution with parameters 
% specified in the vector 'parameters'
F = circ_vmpdf(X, parameters(1), parameters(2));