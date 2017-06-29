% Couette

close all;
clear all;
clc;

%% Options

path = {'/home/9-1-18-03_CouetteFlow_2D'};
% Case = {'couette_2D_x1L_dx2en5_viscM1997'};
Case = {'data_csv/'};
% File = {'profile.1.csv', ...
%         'profile.2.csv', ...
%         'profile.5.csv', ...
%         'profile.10.csv', ...
%         'profile.223.csv'};
File = {'current_simulation.csv'};

U_ref  = 1.25e-5;
L_ref  = 1.0e-3;
nu_ref = 1.0e-6;
%F      = 1.0e-4

% t      = [0.0225 0.045 0.1125 0.225 inf];
t      = [inf];
npy    = 50;

Opt    = {true,true};

Line     = {'.b'};
Line_ref = {'-k'};

XLabel = {'$y/L$'};
YLabel = {'$u/U$'};

XLim = {[0 1]};
YLim = {[0 1]};

XTick = {0:0.2:1};
YTick = {0:0.2:1};

Position = [0.35 0.3 0.3 0.4];

n_inf  = 1000;

%% Read

nt = length(t);
nc = length(Case);

u        = cell(nt,nc);
y        = cell(nt,nc);
M        = cell(nt,nc);
for i = 1:nt
    for j = 1:nc
        M{i,j}     = csvread(strcat(path{j},'/',Case{j},File{i}),1,0);
        u{i,j}     = M{i,j}(:,1)/U_ref;
        y{i,j}     = M{i,j}(:,19)/L_ref;
    end
end


%% Calculate

ny    = npy+2;
u_ref = cell(nt,1);
y_ref = cell(nt,1);
for it = 1:nt
    u_ref{it,1} = zeros(ny,1);
    y_ref{it,1} = zeros(ny,1);
    y_ref{it,1}(ny) = L_ref;
    for iy = 2:ny-1
        y_ref{it,1}(iy) = L_ref/npy*(iy-1.5);
    end
end

for it = 1:nt
	for iy = 1:ny
		u_ref{it,1}(iy) = U_ref/L_ref*y_ref{it,1}(iy);
		%u_ref{it,1}(iy) = 0.5*F*y_ref{it,1}(iy)*(y_ref{it,1}(iy)-L_ref)/nu_ref;
		for i = 1:n_inf
			%u_old    = u{it,1}(iy); %DEBUG
			u_ref{it,1}(iy) = u_ref{it,1}(iy) ...
                            + 2*U_ref/(i*pi)*(-1)^(i) ...
			                * sin(i*pi*y_ref{it,1}(iy)/L_ref) ...
                            * exp(-nu_ref*t(it)*(i*pi/L_ref)^2);
			%u_ref{it,1}(iy) = u_ref{it,1}(iy) ...
            %                + 4*F*L_ref^2/(nu_ref*(pi*(2*i+1))^3) ...
			%                * sin(pi*y_ref{it,1}(iy)*(2*i+1)/L_ref) ...
			%                * exp(-nu_ref*t(it)*((2*i+1)*pi/L_ref)^2);
        end
		%u(it,iy)/u_old %DEBUG
    end
    y_ref{it,1}(:) = y_ref{it,1}(:)/L_ref;
    u_ref{it,1}(:) = u_ref{it,1}(:)/U_ref;
end

%%Interpolation
xq = 0 : 0.01 : 1;
u_interp     = interp1(y{1},u{1},xq);
u_ref_interp = interp1(y_ref{1},u_ref{1},xq);

%% Get error
err=norm(u_interp-u_ref_interp);
limit = 0.1;

%% Print error
if(err > limit) 
    fprintf( 'Validation = !! PROBLEM !!, error too high compare to the reference data: %e', err);
else 
    fprintf( 'Validation = OK,  with an error compare to the reference data: %e', err);
end

h=figure('visible','off');
plot(xq, u_interp, xq, u_ref_interp);
legend('Current','Reference');
xlabel('Time');
ylabel('Velocity');
title('9-1-18-03 CouetteFlow 2D');
saveas(h,'/home/result/png/9-1-18-03_CouetteFlow_2D.png', 'png');

quit

%% Plot
% plot(y{1}, u{1}, y_ref{1}, u_ref{1});

% f = plotProfileRef( y(:,:), u(:,:), y_ref(:,:), u_ref(:,:), ...
%                     Line, Line_ref, Name, Name_ref, XLabel, YLabel, ...
%                     XLim, YLim, XTick, YTick, Legend, Opt, Position );

% print( f, 'png/couette', '-dpng', '-r300' );

