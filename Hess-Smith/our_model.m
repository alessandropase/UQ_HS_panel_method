function Cl = model_fun(U, alpha_deg)

FileName = 'naca_23012_closed';


u_inf = [U;
          0];                                                               % ASINTOTIC VELOCITY [m/s]

h_1 = 0.3;                                                                  % DISTANCE OF THE FOIL FROM THE GROUND


alpha = alpha_deg * pi / 180;                                               % AoA [rad]

nacapane_tab = importPanelisation(FileName);                                % IMPORT OF PANELISATION FROM FILE .txt
array = table2array(nacapane_tab);                                          % CONVERTION OF THE TABLE TO AN ARRAY

array = ROT_alpha(array, alpha);                                            % ROTATION OF THE FOIL OF AoA

array = flip(array);
nacapane.x_pane = array(:, 1);                                              % DEFINING THE STRUCT OF POINTS - X
nacapane.y_pane = array(:, 2) + h_1;                                        % DEFINING THE STRUCT OF POINTS - Y

N = length(nacapane.x_pane)-1;                                              % NUMBER OF PANELS

%% CREATING THE IMAGE OF THE FOIL 
nacapane.x_pane_m = array(:, 1);                                            % DEFINING THE STRUCT OF POINTS - X (image)
nacapane.y_pane_m = -(array(:, 2) + h_1);                                   % DEFINING THE STRUCT OF POINTS - Y (image)

%% MATRIX DEFINITION

center  = zeros(N,2);                                                       % MEMORY ALLOCATION
A_s     = zeros(N);                                                         % MEMORY ALLOCATION
a_v_aux = zeros(N);                                                         % MEMORY ALLOCATION
c_s_1   = zeros(1, N);                                                      % MEMORY ALLOCATION
c_v_1   = zeros(1, N);                                                      % MEMORY ALLOCATION
c_v_1_m = zeros(1, N);                                                      % MEMORY ALLOCATION
c_s_N   = zeros(1, N);                                                      % MEMORY ALLOCATION
c_v_N   = zeros(1, N);                                                      % MEMORY ALLOCATION
c_v_N_m = zeros(1, N);                                                      % MEMORY ALLOCATION
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
        Estremo_1_m = [nacapane.x_pane_m(j); nacapane.y_pane_m(j)];         % DEFINITION OF COORDINATES OF FIRST POINT OF PANEL j               
        Estremo_2_m = [nacapane.x_pane_m(j+1); nacapane.y_pane_m(j+1)];     % DEFINITION OF COORDINATES OF LAST POINT OF PANEL j
        beta_m = atan2((Estremo_2_m(2) - Estremo_1_m(2)) , ...
                     (Estremo_2_m(1) - Estremo_1_m(1)));                    % LOCAL SDR ANGLE [rad]
        Q_m = [ cos(beta_m), sin(beta_m);
             -sin(beta_m), cos(beta_m)                                      % G2L MATRIX
            ]; 


        u_sorg = ViSorgente(Centro, Estremo_1, Estremo_2, Q', Q);           % COMPUTATION OF VELOCITY ON CENTER OF PANEL i DUE TO PANEL SOURCE j
        u_sorg_m = ViSorgente(Centro, Estremo_1_m, Estremo_2_m, Q_m', Q_m); % COMPUTATION OF VELOCITY ON CENTER OF PANEL i DUE TO MIRROR PANEL SOURCE j        
        u_vort = ViVortice(Centro, Estremo_1, Estremo_2, Q', Q);            % COMPUTATION OF VELOCITY ON CENTER OF PANEL i DUE TO PANEL VORTEX j
        u_vort_m = ViVortice(Centro, Estremo_1_m, Estremo_2_m, Q_m', Q_m);  % COMPUTATION OF VELOCITY ON CENTER OF PANEL i DUE TO MIRROR PANEL VORTEX j

        A_s(i,j) = (u_sorg' + u_sorg_m') * n;                               % DEFINITION OF ELEMENT i,j OF A_S SUBMATRIX
        a_v(i) = a_v(i) + (u_vort' - u_vort_m') * n;                        % DEFINITION OF ELEMENT i OF a_v 
        
        if i == 1                                                           % IF-CYCLE TO DEFINE c_s AND c_v COEFFICIENT
          c_s_1(j) = (u_sorg'+u_sorg_m') * tau1;
          c_v_1(j) = u_vort'*tau1;
          c_v_1_m(j) = u_vort_m'*tau1;
        end
        if i == N
          c_s_N(j) = (u_sorg'+u_sorg_m') * tauN;
          c_v_N(j) = u_vort'*tauN;
          c_v_N_m(j) = u_vort_m'*tauN;
        end
    end

    b_s(i) = -u_inf' * n;                                                   % DEFINITION OF b_s VECTOR
    b_v = - u_inf' * (tau1+tauN);                                           % DEFINITION OF b_v VARIABLE

end
c_s = c_s_1+c_s_N;                                                          % FINAL CHARACTERIZATION OF c_s
c_v = sum(c_v_1)+sum(c_v_N) - (sum(c_v_1_m)+sum(c_v_N_m));                  % FINAL CHARACTERIZATION OF c_s

%% LINEAR SISTEM
A = [A_s,  a_v;
     c_s, c_v];                                                             % ASSEMBLY OF LINEAR SYSTEM MATRIX
b = [b_s; 
     b_v];                                                                  % ASSEMBLY OF CONSTANT TERM

x = A\b;                                                                    % COMPUTATION OF LINEAR SISTEM SOLUTION
q = x(1:(end-1));                                                           % INTENSITIES OF SURCES
gamma = x(end);                                                             % INTENSITY OF VORTEX

u_tot = @(x,y) TotalVelocity_ground(q, gamma, u_inf, x, y, nacapane.x_pane, nacapane.y_pane, nacapane.x_pane_m, nacapane.y_pane_m); % COMPUTATION OF VELOCITY

%% Coefficient computation

Cp = zeros(N, 1);                                                           % MEMORY ALLOCATION
Circ = 0;                                                                   % VARIABLE INITIALIZATION

for i = 1 : N

    Estremo_1 = [nacapane.x_pane(i); nacapane.y_pane(i)];                   % DEFINITION OF COORDINATES OF FIRST POINT OF PANEL i           
    Estremo_2 = [nacapane.x_pane(i+1); nacapane.y_pane(i+1)];               % DEFINITION OF COORDINATES OF LAST POINT OF PANEL i           

    [~, tau] = versors (Estremo_1, Estremo_2);                              % DEFINITION OF TANGENT VERSOR OF PANEL i
    
    aux = u_tot(center(i,1),center(i,2))'*tau;
    Cp(i) = 1 - aux^2/(norm(u_inf))^2;                                      % PRESSURE COEFFICIENT COMPUTATION

    Circ = Circ + norm(Estremo_2 - Estremo_1)*gamma;                        % CIRCOLATION COMPUTATION
end

Cl = -2*Circ/norm(u_inf);                                                   % LIFT COEFFICIENT COMPUTATION