clear all;
clc;

mymodel= true;
if mymodel
    L = 130e-6;
    dsl = 0.625;
    dso = 0;
    st_turns= 4;
    rt_turns = 20;
    r = 0.0642;
    g0 = 1.6e-3;
    i=1;
    j=3;
else
    L = 0.076;
    dsl = 0.5;
    dso = -0.5;
    N= 52;
    r = 73.75e-3;
    g0 = 0.5e-3;
end


data_size = 100;
h_size = 100;
alpha = ((2*pi)/45);
gamma = ((11*pi)/180);
theta_r = linspace(0, 2*pi, data_size);
phi = linspace(0, 2*pi, data_size); 
z = linspace(0, L, data_size);
d = 0.3;
nri = ((alpha)/(2*pi));
mu_0 = 4* pi* 10e-7;

if mymodel
    Lxy_mat = cell(6, 6);
    turns_phases = cell(2, 3);
    winding_phases = cell(2, 3);
    for index = 1:3
        turns_phases{1, index} = my_stator_turns_function(h_size, phi, st_turns,index);
        turns_phases{2, index} = my_rotor_turns_function(h_size, data_size, phi, rt_turns, theta_r,index);
        [Z, Phi] = meshgrid(z, phi);
        d_z = (((dsl-dso)/L)*Z)+(dso);
        G = g0 * (1 - (d_z .* cos(Phi)));
        winding_phases{1, index} = winding_function(L, phi, z, G, turns_phases{1, index});
        winding_phases{2, index} = compute_rotorwinding_data(data_size, L, phi, z, G, turns_phases{2, index});
    end

    for i = 0:5
        if i > 2
            i_t = 2;
        else
            i_t = 1;
        end
        for j = 0:5
            j_t = mod(i, 3) + 1;
            j_w = mod(j, 3) + 1;
            if j > 2
                i_w = 2;
            else
                i_w = 1;
            end
            L_yx = zeros(data_size, data_size);
            t_f_temp = turns_phases{i_t, j_t};
            w_f_temp = winding_phases{i_w, j_w};
            for index_theta= 1:data_size
                for index_z = 1:data_size
                    integrand = zeros(data_size, data_size);
                    for index_phi = 1:data_size
                        if i_w == 1 && i_t ==1
                            temp = r * t_f_temp(index_theta, index_phi) * w_f_temp(index_theta,index_phi) / G(index_phi, index_z);
                        elseif i_w == 1 && i_t ==2
                            temp = r * t_f_temp(index_phi, index_theta) * w_f_temp(index_theta, index_phi) / G(index_phi, index_z);
                        elseif i_w == 2 && i_t ==1
                            temp = r * t_f_temp(index_theta, index_phi) * w_f_temp(index_phi, index_theta) / G(index_phi, index_z);
                        else
                            temp = r * t_f_temp(index_phi, index_theta) * w_f_temp(index_phi, index_theta) / G(index_phi, index_z);
                        end
                        for l = 1:data_size
                            integrand(index_phi, l) = temp;
                        end
                    end
                    temp = trapz(phi, trapz(z, integrand, 1));
                    L_yx(index_theta,index_z) = mu_0 * temp;
                end
            end

            [Z_i, Theta] = meshgrid(z, theta_r);
            % [L_yx, Z_i, Theta] = compute_inductance_old(data_size, winding_phases{m_w, n_w}, G, phi, z, theta_r, turns_phases{m_t, n_t}, r,  mu_0);
            Lxy_mat{i +1 ,j+1} = L_yx;
        end
    end
       
else
    n_si = turns_function(h_size, phi, N, i);
    
    nri_list = compute_nri_list(h_size, data_size, alpha, gamma, phi, j, theta_r);
    [Z, Phi] = meshgrid(z, phi);
    d_z = (((dsl-dso)/L)*Z)+(dso);
    G = g0 * (1 - (d_z .* cos(Phi)));
    statorwinding_data = winding_function(L, phi, z, G, n_si);
    rotorwinding_data = compute_rotorwinding_data(data_size, L, phi, z, G, nri_list);
    [L_yx, Z_i, Theta] = compute_inductance_old(data_size, statorwinding_data, G, phi, z, theta_r, nri_list, r,  mu_0);
end



%  plot_phi(phi, n_si, '\phi', 'n_{si}(\phi)');
% plot_phi(phi, nri_list(1, :), '\phi', 'n_{ri}(\phi)');
% plot_phi_for_each_z(phi, G, 'Phi (rad)', 'g (m)');
% plot_surf(Phi, Z, G, 'Phi (radians)', 'Z (meters)', 'G (meters)', 'G vs Phi vs Z');
% plot_phi(phi, statorwinding_data, '\phi', 'N_{s}(\phi)');
% plot_phi(phi, rotorwinding_data(20, :), '\phi', 'N_{r}(\phi)');
plot_surf(Theta, Z_i, Lxy_mat{1, 6}, '\theta (radians)', 'Z (meters)', 'L_{yx}', 'L_{yx} vs Z vs \theta');
% plot(Theta, L_yx,'red')
ylabel('L_{yx}');
xlabel('\theta (radians)');
title('Inductance at end of Stack Length')

plot_surf(Theta, Z_i, Lxy_mat{6, 1}, '\theta (radians)', 'Z (meters)', 'L_{yx}', 'L_{yx} vs Z vs \theta');
% plot(Theta, L_yx,'red')
ylabel('L_{yx}');
xlabel('\theta (radians)');
title('Inductance at end of Stack Length')


function res = double_integral(first_integral, second_integral, func)
first_result = trapz(first_integral, func, 2);
res = trapz(second_integral, first_result);
end

function n_si = turns_function(h_size, angle, N, i)
n_si = 3*N/2*ones(size(angle));
for h = 1:2:h_size
    term = (2*N/(h*pi))*sin(h*pi/3)*(1 + 2*cos(h*pi/9))*cos(2*h*(angle - (i-1)*pi/3));
    n_si = n_si + term;
end
end

function res = winding_function(L, phi, z, G , func)
int_n_over_g = double_integral(z, phi, func ./ G);
g_inv = 1./G;
g_avg = (1/(2*pi*L))*double_integral(z, phi, g_inv);
turns_func_avg = (1/(2*pi*L*g_avg))*int_n_over_g;
res = func - turns_func_avg;
end



function rotorwinding_data = compute_rotorwinding_data(data_size, L, phi, z, G, nri_list)
rotorwinding_data = zeros(data_size);
for index= 1: data_size
    rotorwinding_data(index, :) =  winding_function(L, phi, z, G, nri_list(index, :));
end
end



function n_si = my_stator_turns_function(h_size, phi, N,i)
n_si = ((38*N)/15)*ones(size(phi)); 

for h = 1:h_size 
    sum1 = (((N)/(pi*h)) * (cos(2*h*(phi-(((i-1)*pi)/3))))) *((2*sin((2*h*pi)/5) + (2*sin((2*h*4*pi)/15)) - (2*sin((2*h*11*pi)/15)) - (2*sin((2*h*4*pi)/5)) + (sin((2*h*pi)/3)) - (sin((2*2*h*pi)/3))));
    sum2 = (((N)/(pi*h)) * (sin(2*h*(phi-(((i-1)*pi)/3))))) * ((5 - (2*cos((2*h*pi)/5)) - (2*cos((2*h*4*pi)/15)) + (2*cos((2*h*11*pi)/15)) + (2*cos((2*h*4*pi)/5)) - (cos(2*h*pi/3)) + (cos((2*2*h*pi)/3)) - (5*cos(2*h*pi))));
    n_si = (n_si + sum1 + sum2)  ; 
end
n_si = repmat(n_si, numel(n_si), 1);
end

function nri_list = compute_nri_list(h_size, data_size, alpha, gamma, phi, j, theta_r)
    h = (1:h_size).'; % column vector
    nri_list = zeros(data_size, data_size); % preallocate result matrix
    
    for a = 1:data_size
        for b = 1:data_size
            sum1 = (4./(h.^2*pi*gamma)).*sin(h*alpha/2).*sin(h*gamma/2).*cos(h*(phi(b) - ((j-1)*alpha) - theta_r(a)));
            nri_list(a, b) = sum(sum1) + alpha/(2*pi);
        end
    end
end

function nri_list = my_rotor_turns_function(h_size, data_size, phi, N, theta_r,j)
nri_list = zeros(data_size);
for a = 1:data_size
    n_ri = 0;
    for h = 1:h_size 
        sum1 = (((N)/(pi*h)) * (cos(2*h*(phi- theta_r(a)-(((j-1)*pi)/3))))) *((sin((h*pi)/3) + (2*sin((h*pi)/2)) + (sin((2*h*pi)/3)) - (2*sin((h*3*pi)/2)) - (sin((4*h*pi)/3)) - (sin((5*h*pi)/3))));
        sum2 = (((N)/(pi*h)) * (sin(2*h*(phi-theta_r(a)-(((j-1)*pi)/3))))) * ((4 - (cos((h*pi)/3)) - (2*cos((h*pi)/2)) - (cos((2*h*pi)/3)) +(2*cos((3*h*pi)/2))+ (cos((h*4*pi)/3)) + (cos(5*h*pi/3)) - (4*cos(2*h*pi))));
        n_ri = (n_ri + sum1 + sum2);
    end

    nri_list(a, :) = n_ri+ (2*N);
end
end


function [L_yx, Z_i, Theta] = compute_inductance_old(data_size, statorwinding_data, G, phi, z, theta_r, nri_list, r, mu_0)
L_yx = zeros(data_size, data_size);
for index_theta= 1:data_size
    for index_z = 1:data_size
        integrand = zeros(data_size, data_size);

        for index_phi = 1:data_size
            temp = r * nri_list(index_phi, index_theta) * statorwinding_data(index_phi) / G(index_phi, index_z);
            for l = 1:data_size
                integrand(index_phi, l) = temp;
            end
        end
        temp = trapz(phi, trapz(z, integrand, 1));
        L_yx(index_theta,index_z) = mu_0 * temp;
    end
end

[Z_i, Theta] = meshgrid(z, theta_r);
end


function plot_phi(phi, data, xlabel_text, ylabel_text)
plot(phi, data);
xlabel(xlabel_text);
ylabel(ylabel_text);
xticks([0, pi, 2*pi]);
xticklabels({'0', '\pi', '2\pi'});
grid on;
end

function plot_phi_for_each_data_size(phi, nri_list, xlabel_text, ylabel_text)
figure;
hold on;
for t = 1:length(phi)
    plot(phi, nri_list(t, :));
end
hold off;
xlabel(xlabel_text);
ylabel(ylabel_text);
xticks([0, pi, 2*pi]);
xticklabels({'0', '\pi', '2\pi'});
end

function plot_phi_for_each_z(phi, G, xlabel_text, ylabel_text)
figure;
hold on;
for i = 1:length(phi)
    plot(phi, G(:,i));
end
hold off;
xlabel(xlabel_text);
ylabel(ylabel_text);
end

function plot_surf(x, y, z, xlabel_text, ylabel_text, zlabel_text, title_text)
figure;
surf(x, y, z);
xlabel(xlabel_text);
ylabel(ylabel_text);
zlabel(zlabel_text);
title(title_text);
colorbar;
end

