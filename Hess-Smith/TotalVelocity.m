function u_tot = TotalVelocity (q, gamma, u_inf, x, y, x_pane, y_pane)

N = length(x_pane)-1;
u_tot = u_inf;

for i = 1 : N
    beta = atan2((y_pane(i+1) - y_pane(i)) , ...
                 (x_pane(i+1) - x_pane(i)));
    Q = [ cos(beta), sin(beta);
         -sin(beta), cos(beta)
        ];

    Estremo_1 = [x_pane(i); y_pane(i)];
    Estremo_2 = [x_pane(i+1); y_pane(i+1)];

    u_sorg = ViSorgente([x;y], Estremo_1, Estremo_2, Q', Q);
    u_vort = ViVortice([x;y], Estremo_1, Estremo_2, Q', Q);
    
    u_tot = u_tot + u_sorg*q(i) + gamma*u_vort;
end
