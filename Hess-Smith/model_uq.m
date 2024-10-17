function L = model_uq(X)

U_vect = X(:, 1);
alpha_deg_vect = X(:, 2);

FileName = 'naca_23012_closed';
Cl = zeros(length(U_vect), 1);
L = zeros(length(U_vect), 1);


for ind = 1 : length(U_vect)
    U = U_vect(ind);
    alpha_deg = alpha_deg_vect(ind);

    u_inf = [U;
              0];                                                               % ASINTOTIC VELOCITY [m/s]
    
    
    alpha = alpha_deg;                                               % AoA [rad]
    
    
    
    nacapane_tab = importPanelisation(FileName);                                % IMPORT OF PANELISATION FROM FILE .txt
    array = table2array(nacapane_tab);                                          % CONVERTION OF THE TABLE TO AN ARRAY
    array(:, 1) = array(:, 1) - 0.25;
    array = ROT_alpha(array, alpha);
    
    array = flip(array);
    nacapane.x_pane = array(:, 1);                                              % DEFINING THE STRUCT OF POINTS - X
    nacapane.y_pane = array(:, 2);                                              % DEFINING THE STRUCT OF POINTS - Y
    
    N = length(nacapane.x_pane)-1;                                              % NUMBER OF PANELS
    
    
    %% MATRIX DEFINITION
    
    center  = zeros(N,2);                                                       % MEMORY ALLOCATION
    A_s     = zeros(N);                                                         % MEMORY ALLOCATION
    a_v_aux = zeros(N);                                                         % MEMORY ALLOCATION
    c_s_1   = zeros(1, N);                                                      % MEMORY ALLOCATION
    c_v_1   = zeros(1, N);                                                      % MEMORY ALLOCATION
    c_s_N   = zeros(1, N);                                                      % MEMORY ALLOCATION
    c_v_N   = zeros(1, N);                                                      % MEMORY ALLOCATION
    b_s     = zeros(N, 1);                                                      % MEMORY ALLOCATION
    a_v     = zeros(N, 1);                                                      % MEMORY ALLOCATION
    
    
    point_1 = [nacapane.x_pane(1); nacapane.y_pane(1)];                         % DEFINITION OF COORDINATES OF FIRST POINT OF FIRST PANEL 
    point_2 = [nacapane.x_pane(2); nacapane.y_pane(2)];                         % DEFINITION OF COORDINATES OF LAST POINT OF FIRST PANEL 
    point_N = [nacapane.x_pane(N); nacapane.y_pane(N)];                         % DEFINITION OF COORDINATES OF FIRST POINT OF LAST PANEL 
    point_NN = [nacapane.x_pane(N+1); nacapane.y_pane(N+1)];                    % DEFINITION OF COORDINATES OF LAST POINT OF LAST PANEL 
    
    [~, tau1] = versors(point_1, point_2);                                      % DEFINITION OF TANGENT VERSOR OF FIRST PANEL
    [~, tauN] = versors(point_N, point_NN);                                     % DEFINITION OF TANGENT VERSOR OF LAST PANEL
    
    for i = 1 : N
        point_i = [nacapane.x_pane(i); nacapane.y_pane(i)];                     % DEFINITION OF COORDINATES OF FIRST POINT OF PANEL i
        point_ii = [nacapane.x_pane(i+1); nacapane.y_pane(i+1)];                % DEFINITION OF COORDINATES OF LAST POINT OF PANEL i
        [n, ~] = versors (point_i, point_ii);                                   % DEFINITION OF NORMAL VERSOR OF PANEL i
        center_x = (point_i(1) + point_ii(1))*0.5;                              % DEFINITION OF X-COORDINATE OF CENTER POINT OF PANEL i
        center_y = (point_i(2) + point_ii(2))*0.5;                              % DEFINITION OF X-COORDINATE OF CENTER POINT OF PANEL i
        Centro = [center_x; center_y];                                          % DEFINITION OF COORDINATES OF CENTER POINT OF PANEL i
        center(i, :) = Centro;                                                  % STORING THE COORDINATE OF CENTER POINTS OF ALL PANELS
    
        for j = 1 : N
            Estremo_1 = [nacapane.x_pane(j); nacapane.y_pane(j)];               % DEFINITION OF COORDINATES OF FIRST POINT OF PANEL j               
            Estremo_2 = [nacapane.x_pane(j+1); nacapane.y_pane(j+1)];           % DEFINITION OF COORDINATES OF LAST POINT OF PANEL j
            beta = atan2((Estremo_2(2) - Estremo_1(2)) , ...
                         (Estremo_2(1) - Estremo_1(1)));                        % LOCAL SDR ANGLE [rad]
            Q = [ cos(beta), sin(beta);
                 -sin(beta), cos(beta)                                          % G2L MATRIX
                ];                                                              
    
            u_sorg = ViSorgente(Centro, Estremo_1, Estremo_2, Q', Q);           % COMPUTATION OF VELOCITY ON CENTER OF PANEL i DUE TO PANEL SOURCE j
            u_vort = ViVortice(Centro, Estremo_1, Estremo_2, Q', Q);            % COMPUTATION OF VELOCITY ON CENTER OF PANEL i DUE TO PANEL VORTEX j
    
            A_s(i,j) = u_sorg' * n;                                             % DEFINITION OF ELEMENT i,j OF A_S SUBMATRIX
            a_v(i) = a_v(i) + u_vort'*n;                                        % DEFINITION OF ELEMENT i OF a_v 
            
            if i == 1                                                           % IF-CYCLE TO DEFINE c_s AND c_v COEFFICIENT
              c_s_1(j) = u_sorg'*tau1;
              c_v_1(j) = u_vort'*tau1;
            end
            if i == N
              c_s_N(j) = u_sorg'*tauN;
              c_v_N(j) = u_vort'*tauN;
            end
        end
    
        b_s(i) = -u_inf' * n;                                                   % DEFINITION OF b_s VECTOR
    end
    b_v = - u_inf' * (tau1+tauN);                                               % DEFINITION OF b_v VARIABLE
    c_s = c_s_1+c_s_N;                                                          % FINAL CHARACTERIZATION OF c_s
    c_v = sum(c_v_1)+sum(c_v_N);                                                % FINAL CHARACTERIZATION OF c_s
    
    %% LINEAR SISTEM
    A = [A_s,  a_v;
         c_s, c_v];                                                             % ASSEMBLY OF LINEAR SYSTEM MATRIX
    b = [b_s; 
         b_v];                                                                  % ASSEMBLY OF CONSTANT TERM
    
    x = A\b;                                                                    % COMPUTATION OF LINEAR SISTEM SOLUTION
    q = x(1:(end-1));                                                           % INTENSITIES OF SURCES
    gamma = x(end);                                                             % INTENSITY OF VORTEX
    
    u_tot = @(x,y) TotalVelocity(q, gamma, u_inf, x, y, nacapane.x_pane, nacapane.y_pane); % COMPUTATION OF VELOCITY
    
    Cp = zeros(N, 1);                                                           % MEMORY ALLOCATION
    Circ = 0;                                                                   % VARIABLE INITIALIZATION
    Cl_new = 0;
    
    for i = 1 : N
    
        Estremo_1 = [nacapane.x_pane(i); nacapane.y_pane(i)];                   % DEFINITION OF COORDINATES OF FIRST POINT OF PANEL i           
        Estremo_2 = [nacapane.x_pane(i+1); nacapane.y_pane(i+1)];               % DEFINITION OF COORDINATES OF LAST POINT OF PANEL i           
    
        [n, tau] = versors (Estremo_1, Estremo_2);                              % DEFINITION OF TANGENT VERSOR OF PANEL i
        
        aux = u_tot(center(i,1),center(i,2))'*tau;
        Cp(i) = 1 - aux^2/(norm(u_inf))^2;                                      % PRESSURE COEFFICIENT COMPUTATION
        Cl_new = Cl_new +  Cp(i)*norm(Estremo_2-Estremo_1)* n(2);
        Circ = Circ + norm(Estremo_2 - Estremo_1)*gamma;                        % CIRCOLATION COMPUTATION
    end
    
    Cl(ind) = 2*Circ/norm(u_inf);                                                   % LIFT COEFFICIENT COMPUTATION
    L(ind) = 0.5*1.225*U^2*Cl(ind);
end
end