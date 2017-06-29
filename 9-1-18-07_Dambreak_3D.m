% Couette

close all;
clear all;
clc;

%% Options

path = {'/home/9-1-18-07_Dambreak_3D'};
% Case = {'couette_2D_x1L_dx2en5_viscM1997'};
Case = {'data_csv/'};
% File = {'profile.1.csv', ...
%         'profile.2.csv', ...
%         'profile.5.csv', ...
%         'profile.10.csv', ...
%         'profile.223.csv'};
File = {'current_simulation0.csv', ...
        'reference_simulation0.csv'};

U_ref  = 1.25e-5;
L_ref  = 0.5;
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

nt = length(File);
nc = length(Case);

u        = cell(nt,nc);
y        = cell(nt,nc);
M        = cell(nt,nc);
for i = 1:nt
    for j = 1:nc
        M{i,j}     = csvread(strcat(path{j},'/',Case{j},File{i}),1,0);
        u{i,j}     = M{i,j}(:,14)/0.3*-1;
        y{i,j}     = M{i,j}(:,64)*(9.81/0.3)^0.5;
    end
end

%%Interpolation
xq = 0 : 0.01 : 2.5;
u_interp     = interp1(y{1,1},u{1,1},xq);
u_ref_interp = interp1(y{2,1},u{2,1},xq);

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
xlabel('High');
ylabel('Velocity');
title('9-1-18-07 Dambreak 3D');
saveas(h,'/home/result/png/9-1-18-07_Dambreak_3D.png', 'png');

% plot(xq, u_interp, xq, u_ref_interp)

% quit
%% Plot
% plot(y{1}, u{1}, y_ref{1}, u_ref{1});

% f = plotProfileRef( y(:,:), u(:,:), y_ref(:,:), u_ref(:,:), ...
%                     Line, Line_ref, Name, Name_ref, XLabel, YLabel, ...
%                     XLim, YLim, XTick, YTick, Legend, Opt, Position );

% print( f, 'png/couette', '-dpng', '-r300' );

