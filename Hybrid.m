clear
clear global 
global topnodalvalue V_init count2  dt tt count3  count4 countstupid axon_r
global g_leak g_na g_kd V_T g_m t_max E_na E_k E_leak Cm cond g_na_alt g_kd_alt g_leak_alt
global count_head count_tail topnodalvalue_head topnodalvalue_tail

%%%%%%%%%%%% Simulation Variables %%%%%%%%%%%%

V_init = [-70.6878    0.0    1.0    0.0    0.5]; %mV
currentin = 2.0E-4;
currentinduration = 5; %ms
axon_r = 1.0; %um
swelling_r = 29.0; %um


simulating_time = 20; %ms
loging_interval_time = 0.1; %ms 

tot_one_d_length = 0.5; %cm
node_dist = 50; %um
dt=0.005; %ms
steps=ceil(simulating_time./dt);
minmesh = 3.0; %um


% Units used in computation are cm, mV
%%%%%%%%%%%% Axon System Constants (HH Type) %%%%%%%%%%%%
g_na = 50; % mS/cm^2
g_kd = 4.8; % mS/cm^2
V_T = -61.5; % mV
g_m = 0.13; % mS/cm^2  %not in classic H&H
t_max = 1123.5; % ms

g_leak = 0.1; % mS/cm^2
g_constant_leak = 1.0;
g_leak_alt = g_leak .* g_constant_leak;
g_constant_na = 1.0;
g_na_alt = g_na .* g_constant_na;
g_constant_kd = 1.0;
g_kd_alt = g_kd .* g_constant_kd;


E_na = 50; % mV
E_k = -90; % mV
E_leak = -70; % mV
% E_ca = 120; % mV % Calcium current

cond = 10.0; % mS/cm
Cm = 1; % uF/cm^-2



%%%%%%%%%%%% Simulation Variables %%%%%%%%%%%%


MonitoredNode = 1;


%dttt = 0.000001;
Cmin=-100;
Cmax=50;
count2=1;
count_head=1;
count_tail=1;

N = 1;
model = createpde(N);




%%%%%%%%%%%% Create Geometry %%%%%%%%%%%%
length_unit_ratio = 10000.0; %cm/um

node_dist = node_dist./length_unit_ratio; %cm
tot_n_node = ceil(tot_one_d_length./node_dist);

axon_r = axon_r./length_unit_ratio; %cm
axon_A = pi.*axon_r.*axon_r;

rec_upper_y = axon_r;
rec_right_x = 200.0./length_unit_ratio;
circle_center_x = 100.0./length_unit_ratio;
circle_center_y = 0.0./length_unit_ratio;
circle_radius = swelling_r./length_unit_ratio;




R1 = [3,4,0,rec_right_x,rec_right_x,0,rec_upper_y,rec_upper_y,0,0]';
R2 = [3,4,0,rec_right_x,rec_right_x,0,0,0,-2.*circle_radius,-2.*circle_radius]';

C1 = [1,circle_center_x,circle_center_y,circle_radius]';
C1 = [C1;zeros(length(R1) - length(C1),1)];
gm = [R1,C1,R2];
sf = '(R1+C1)-R2';

ns = char('R1','C1','R2');
ns = ns';
[g,bt] = decsg(gm,sf,ns);
[dl2,bt2] = csgdel(g,bt); 
geometryFromEdges(model,dl2);





%gdm = [1 
%    0
%    0
%    circle_size];
%g = decsg(gdm,'S1',('S1')');

% Create a geometry entity and append to the PDE Model
%geometryFromEdges(model,g);

figure;
pdegplot(model, 'EdgeLabels', 'on', 'FaceLabels', 'on');
title('Axon model (with edge labels)')
xlabel('X-coordinate, cm')
ylabel('Y-coordinate, cm')
axis equal

%%%%%%%%%%%%  2D PDE SETUP Mesh %%%%%%%%%%%%
minmesh = minmesh./length_unit_ratio;
                                           
generateMesh(model,'Hmax',minmesh,...
                         'GeometricOrder','quadratic',...
                         'MesherVersion','R2013a');
figure;
pdeplot(model);
xlabel(['Total mesh' num2str(size(model.Mesh.Nodes,2))]);
axis equal

% PDE SETUP

m = 0;
d = 0;
a = 0;
f = 0;
c = @(region,state) region.y;
specifyCoefficients(model, 'm', m,'d', d,'c', c, 'a', a, 'f', f);


% Boundary Condition


BCTop = applyBoundaryCondition(model,'neumann','Edge',[3,4,8,9],...
                                               'q',1.0E3,'g',@initBCimplicit_1);

BCtail = applyBoundaryCondition(model,'dirichlet','Edge', 1, 'u', @initBCimplicit_tail);
BChead = applyBoundaryCondition(model,'dirichlet','Edge', 2, 'u', @initBCimplicit_head);



R = solvepde(model);
u = R.NodalSolution;
%pdeplot(model,'XYData',u(:,end),'Contour','on','ColorMap','jet');


% update BC


%%%%%%%%%%%% Plottign parameter %%%%%%%%%%%%
logingiterval = ceil(loging_interval_time./dt);
plotinterval = logingiterval;
voltage_lower = -100;
voltage_upper = 50;


h = waitbar(0,'Calculating...');
set(h,'position',[100 750 300 60]); 

figure1 = figure(5);
set(figure1,'position',[100 50 1800 1100]); 



%axes1 = axes('Parent',figure1);
filename = 'testnew51.gif';


%%%%%%%%%%%% 1-D initialization %%%%%%%%%%%%



current_u_front = zeros(tot_n_node,5);
current_u_back = zeros(tot_n_node,5);
next_u_front = zeros(tot_n_node,1);
next_u_back = zeros(tot_n_node,1);
data_matrix_front{steps} = [];
data_matrix_back{steps} = [];

for ii=1:1:tot_n_node
current_u_front(ii,:) = V_init;
current_u_back(ii,:) = V_init;
end


const_1 = 2.0.*pi.*axon_r.*node_dist;
const_2 = 1.0./node_dist.*cond.*axon_A;
const_3 = 1.0./dt.*Cm.*2.0.*pi.*axon_r.*node_dist;

average_grad_head = 0.0;
average_grad_tail = 0.0;

for ii=1:1:tot_n_node
    
node_location(ii) = node_dist.*ii;    
    
end


%%%%%%%%%%%%  Initiate LOOP %%%%%%%%%%%%
count4=0;
countstupid = 0;

datalog{1}=zeros(ceil(steps./logingiterval),1);
for kk = 2:1:11
datalog{kk}=zeros(ceil(steps./logingiterval),size(topnodalvalue,1));
end

totalt=zeros(steps,1);
count3 = size(topnodalvalue,1);
influx = zeros(size(topnodalvalue,1),1);

q_c = @(region,state) region.y.*Cm./cond./dt;

overalltic=tic;
tot_toc_one_d = 0;
tot_toc_two_d = 0;
tot_toc_plot = 0;



%% %%%%%%%%%%%%   LOOP  %%%%%%%%%%%%



for tt = 1:1:steps
looptic=tic;
one_d_tic=tic;


%%%%%%%%%%%% 1D Head %%%%%%%%%%%%
alpha_value_front(:,1) = -0.32*(current_u_front(:,1) -V_T -13)./(exp(-(current_u_front(:,1) -V_T -13)./4) -1);
beta_value_front(:,1) = 0.28*(current_u_front(:,1) -V_T -40)./(exp((current_u_front(:,1) -V_T -40)./5) -1);
alpha_value_front(:,2) = 0.128*exp(-(current_u_front(:,1) -V_T -17)./18);
beta_value_front(:,2) = 4.0./(1.0+exp(-(current_u_front(:,1) -V_T -40)./5));
alpha_value_front(:,3) = -0.032*(current_u_front(:,1) -V_T -15)./(exp(-(current_u_front(:,1) -V_T -15)./5) -1);
beta_value_front(:,3) = 0.5*exp(-(current_u_front(:,1) -V_T -10)./40);
alpha_value_front(:,4) = 1.0./(1.0+exp(-(current_u_front(:,1) +35)./10));
beta_value_front(:,4) = t_max./(3.3.*exp((current_u_front(:,1) +35)./20) +exp(-(current_u_front(:,1)+35)./20));


partial_t_value_front(:,1:3) = alpha_value_front(:,1:3).*(1 -current_u_front(:,2:4)) -beta_value_front(:,1:3).*current_u_front(:,2:4);
partial_t_value_front(:,4) = (alpha_value_front(:,4) -current_u_front(:,5)) ./beta_value_front(:,4);


flux_value_front(:,1) = - g_leak.*(current_u_front(:,1)-E_leak) ;
flux_value_front(:,2) = - g_na.*current_u_front(:,2).^3.*current_u_front(:,3).*(current_u_front(:,1)-E_na);
flux_value_front(:,3) = - g_kd.*current_u_front(:,4).^4.*(current_u_front(:,1)-E_k) ;
flux_value_front(:,4) = - g_m.*current_u_front(:,5).*(current_u_front(:,1)-E_k);
net_flux_front = flux_value_front(:,1)+flux_value_front(:,2)+flux_value_front(:,3)+flux_value_front(:,4);
%net_flux = flux_value(:,1)+flux_value(:,2)+flux_value(:,3);

current_u_front(:,2:5) = current_u_front(:,2:5)+partial_t_value_front(:,1:4).*dt;    


if totalt(tt) < currentinduration
next_u_front(1)         = (currentin + net_flux_front(1).*const_1 ...
                    + (current_u_front(2,1)-current_u_front(1,1)).*const_2 ...
                    + current_u_front(1,1).*const_3)./const_3;
else
next_u_front(1)         = (net_flux_front(1).*const_1 ...
                    + (current_u_front(2,1)-current_u_front(1,1)).*const_2 ...
                    + current_u_front(1,1).*const_3)./const_3;    
end
                
                
                
for ii = 2:1:tot_n_node-1
    next_u_front(ii) = (net_flux_front(ii).*const_1 ...
                    - (current_u_front(ii,1)-current_u_front(ii-1,1)).*const_2 + (current_u_front(ii+1,1)-current_u_front(ii,1)).*const_2 ...
                   + current_u_front(ii,1).*const_3)./const_3;
end

next_u_front(tot_n_node)= (net_flux_front(tot_n_node).*const_1 - average_grad_head.*cond.*axon_A ... %% average_grad_head is on -x direction
                    - (current_u_front(tot_n_node,1)-current_u_front(tot_n_node-1,1)).*const_2 ...
                    + current_u_front(tot_n_node,1).*const_3)./const_3;

current_u_front(:,1) = next_u_front(:);



%%%%%%%%%%%% 1D Tail %%%%%%%%%%%%

alpha_value_back(:,1) = -0.32*(current_u_back(:,1) -V_T -13)./(exp(-(current_u_back(:,1) -V_T -13)./4) -1);
beta_value_back(:,1) = 0.28*(current_u_back(:,1) -V_T -40)./(exp((current_u_back(:,1) -V_T -40)./5) -1);
alpha_value_back(:,2) = 0.128*exp(-(current_u_back(:,1) -V_T -17)./18);
beta_value_back(:,2) = 4.0./(1.0+exp(-(current_u_back(:,1) -V_T -40)./5));
alpha_value_back(:,3) = -0.032*(current_u_back(:,1) -V_T -15)./(exp(-(current_u_back(:,1) -V_T -15)./5) -1);
beta_value_back(:,3) = 0.5*exp(-(current_u_back(:,1) -V_T -10)./40);
alpha_value_back(:,4) = 1.0./(1.0+exp(-(current_u_back(:,1) +35)./10));
beta_value_back(:,4) = t_max./(3.3.*exp((current_u_back(:,1) +35)./20) +exp(-(current_u_back(:,1)+35)./20));


partial_t_value_back(:,1:3) = alpha_value_back(:,1:3).*(1 -current_u_back(:,2:4)) -beta_value_back(:,1:3).*current_u_back(:,2:4);
partial_t_value_back(:,4) = (alpha_value_back(:,4) -current_u_back(:,5)) ./beta_value_back(:,4);


flux_value_back(:,1) = - g_leak.*(current_u_back(:,1)-E_leak) ;
flux_value_back(:,2) = - g_na.*current_u_back(:,2).^3.*current_u_back(:,3).*(current_u_back(:,1)-E_na);
flux_value_back(:,3) = - g_kd.*current_u_back(:,4).^4.*(current_u_back(:,1)-E_k) ;
flux_value_back(:,4) = - g_m.*current_u_back(:,5).*(current_u_back(:,1)-E_k);
net_flux_back = flux_value_back(:,1)+flux_value_back(:,2)+flux_value_back(:,3)+flux_value_back(:,4);
%net_flux = flux_value(:,1)+flux_value(:,2)+flux_value(:,3);

current_u_back(:,2:5) = current_u_back(:,2:5)+partial_t_value_back(:,1:4).*dt;    

next_u_back(1)         = (-average_grad_tail.*cond.*axon_A + net_flux_back(1).*const_1 ...
                    + (current_u_back(2,1)-current_u_back(1,1)).*const_2 ...
                    + current_u_back(1,1).*const_3)./const_3;
for ii = 2:1:tot_n_node-1
    next_u_back(ii) = (net_flux_back(ii).*const_1 ...
                    - (current_u_back(ii,1)-current_u_back(ii-1,1)).*const_2 + (current_u_back(ii+1,1)-current_u_back(ii,1)).*const_2 ...
                   + current_u_back(ii,1).*const_3)./const_3;
end

next_u_back(tot_n_node)= (net_flux_back(tot_n_node).*const_1 ...
                    - (current_u_back(tot_n_node,1)-current_u_back(tot_n_node-1,1)).*const_2 ...
                    + current_u_back(tot_n_node,1).*const_3)./const_3;

current_u_back(:,1) = next_u_back(:);


%%%%%%%%%%%% %%%%%%%%%%%%

tot_toc_one_d = tot_toc_one_d + toc(one_d_tic);
two_d_tic=tic;

%%%%%%%%%%%% 2D %%%%%%%%%%%%
count2=1;

setInitialConditions(model,R);

BCTop = applyBoundaryCondition(model,'neumann','Edge',[3,4,8,9],...
                                               'q',q_c,'g',@loopBC_1);
BCtail = applyBoundaryCondition(model,'dirichlet','Edge', 1, 'u', current_u_back(1,1));
BChead = applyBoundaryCondition(model,'dirichlet','Edge', 2, 'u', current_u_front(tot_n_node,1));
                                           
model.SolverOptions.MinStep = 0;
model.SolverOptions.MaxIterations = 50;
model.SolverOptions.ResidualTolerance = 1.0000e-10;
%model.SolverOptions.ResidualNorm = 2;
                                           
R = solvepde(model);

tot_toc_two_d = tot_toc_two_d + toc(two_d_tic);  %% TIME MARKER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

topnodalvalue(:,5) = interpolateSolution(R,topnodalvalue(:,1),topnodalvalue(:,2));

alpha_value(:,1) = -0.32*(topnodalvalue(:,5) -V_T -13)./(exp(-(topnodalvalue(:,5) -V_T -13)./4) -1);
beta_value(:,1) = 0.28*(topnodalvalue(:,5) -V_T -40)./(exp((topnodalvalue(:,5) -V_T -40)./5) -1);
alpha_value(:,2) = 0.128*exp(-(topnodalvalue(:,5) -V_T -17)./18);
beta_value(:,2) = 4.0./(1.0+exp(-(topnodalvalue(:,5) -V_T -40)./5));
alpha_value(:,3) = -0.032*(topnodalvalue(:,5) -V_T -15)./(exp(-(topnodalvalue(:,5) -V_T -15)./5) -1);
beta_value(:,3) = 0.5*exp(-(topnodalvalue(:,5) -V_T -10)./40);
alpha_value(:,4) = 1.0./(1.0+exp(-(topnodalvalue(:,5) +35)./10));
beta_value(:,4) = t_max./(3.3.*exp((topnodalvalue(:,5) +35)./20) +exp(-(topnodalvalue(:,5)+35)./20));

partial_t_value(:,2:4) = alpha_value(:,1:3).*(1 -topnodalvalue(:,6:8)) -beta_value(:,1:3).*topnodalvalue(:,6:8);
partial_t_value(:,5) = (alpha_value(:,4) -topnodalvalue(:,9)) ./beta_value(:,4);
topnodalvalue(:,6:9)=topnodalvalue(:,6:9)+partial_t_value(:,2:5).*dt;


    
[gradx_head,grady_head] = evaluateGradient(R,topnodalvalue_head(:,1),topnodalvalue_head(:,2));
normal_grad_head(:) = gradx_head.*topnodalvalue_head(:,3)+grady_head.*topnodalvalue_head(:,4);    
[gradx_tail,grady_tail] = evaluateGradient(R,topnodalvalue_tail(:,1),topnodalvalue_tail(:,2));
normal_grad_tail(:) = gradx_tail.*topnodalvalue_tail(:,3)+grady_tail.*topnodalvalue_tail(:,4);    
average_grad_head=mean(normal_grad_head);
average_grad_tail=mean(normal_grad_tail);

%%%%%%%%%%%%    %%%%%%%%%%%%
tot_toc_two_d = tot_toc_two_d + toc(two_d_tic);
plot_tic=tic;



%%%%%%%%%%%%  Logging  %%%%%%%%%%%%

    
if rem(tt,logingiterval) == 0

        [gradx,grady] = evaluateGradient(R,topnodalvalue(:,1),topnodalvalue(:,2));
    normal_grad(:) = gradx.*topnodalvalue(:,3)+grady.*topnodalvalue(:,4);

for ii=1:1:size(topnodalvalue,1)

    influx(ii) = -normal_grad(ii).*cond;

flux_value(:,2) = - g_leak.*(topnodalvalue(:,5)-E_leak) ;
flux_value(:,3) = - g_na.*topnodalvalue(:,6).^3.*topnodalvalue(:,7).*(topnodalvalue(:,5)-E_na);
flux_value(:,4) = - g_kd.*topnodalvalue(:,8).^4.*(topnodalvalue(:,5)-E_k) ;
flux_value(:,5) = - g_m.*topnodalvalue(:,9).*(topnodalvalue(:,5)-E_k);    
end

    
    
    
datalog{1}(ceil(tt./logingiterval),1) = totalt(tt);
datalog{2}(ceil(tt./logingiterval),:) = topnodalvalue(:,5);
datalog{3}(ceil(tt./logingiterval),:) = topnodalvalue(:,6);
datalog{4}(ceil(tt./logingiterval),:) = topnodalvalue(:,7);
datalog{5}(ceil(tt./logingiterval),:) = topnodalvalue(:,8);
datalog{6}(ceil(tt./logingiterval),:) = topnodalvalue(:,9);
datalog{7}(ceil(tt./logingiterval),:) = influx(:);
datalog{8}(ceil(tt./logingiterval),:) = flux_value(:,2);
datalog{9}(ceil(tt./logingiterval),:) = flux_value(:,3);
datalog{10}(ceil(tt./logingiterval),:) = flux_value(:,4);
datalog{11}(ceil(tt./logingiterval),:) = flux_value(:,5);

    data_matrix_back{tt./logingiterval}(:,:) = current_u_back(:,:);
    data_matrix_front{tt./logingiterval}(:,:) = current_u_front(:,:);


end
    
    
        totalt(tt+1)=dt+totalt(tt);


    
%%%%%%%%%%%% Plot %%%%%%%%%%%%    
    
    
if rem(tt,plotinterval) == 0


    
% Create axes
axes1 = subplot(3,5,[1 2 6 7]);

hold(axes1,'on');

pdeplot(model,'XYData',R.NodalSolution(:,end),'Contour','on','ColorMap','jet');
text(topnodalvalue(MonitoredNode,1),topnodalvalue(MonitoredNode,2),'\leftarrow Probe ','Color','red','FontSize',14)

if Cmin~=0 || Cmax~=0
set(axes1,'BoxStyle','full','CLim',[Cmin Cmax],'FontSize',16,'Layer',...
    'top','XGrid','on','YGrid','on');
end
axis equal
%drawnow

axes2 = subplot(3,5,3);
hold(axes2,'on');
plot2 = plot(datalog{1}(1:tt./plotinterval,1),datalog{2}(1:tt./plotinterval,MonitoredNode));
set(plot2(1),'Color','k');

xlabel('Time (ms)');
title('Membrane Voltage');
ylabel('(mV)');
axis([0 inf voltage_lower voltage_upper]);
box(axes2,'on');

axes3 = subplot(3,5,4);
hold(axes3,'on');

plot3 = plot(datalog{1}(1:tt./plotinterval,1),[datalog{3}(1:tt./plotinterval,MonitoredNode) datalog{4}(1:tt./plotinterval,MonitoredNode) datalog{5}(1:tt./plotinterval,MonitoredNode) datalog{6}(1:tt./plotinterval,MonitoredNode)],'Parent',axes3);
set(plot3(1),'DisplayName','m');
set(plot3(2),'DisplayName','h');
set(plot3(3),'DisplayName','n');
set(plot3(4),'DisplayName','p');
set(plot3(1),'Color','k');
set(plot3(2),'Color','r');
set(plot3(3),'Color','g');
set(plot3(4),'Color','b');
axis([0 inf -0.2 1.2]);
xlabel('Time (ms)');
title('Auxiliary Variables');
box(axes3,'on');
legend(axes3,'show');

axes4 = subplot(3,5,5);
hold(axes4,'on');
plot4 = plot(datalog{1}(1:tt./plotinterval,1),[datalog{7}(1:tt./plotinterval,MonitoredNode) datalog{8}(1:tt./plotinterval,MonitoredNode) datalog{9}(1:tt./plotinterval,MonitoredNode) datalog{10}(1:tt./plotinterval,MonitoredNode) datalog{11}(1:tt./plotinterval,MonitoredNode)],'Parent',axes4);
set(plot4(1),'DisplayName','Conduct');
set(plot4(2),'DisplayName','Leak');
set(plot4(3),'DisplayName','Na');
set(plot4(4),'DisplayName','Kd');
set(plot4(5),'DisplayName','K_slow');
set(plot4(1),'Color','k');
set(plot4(2),'Color','r');
set(plot4(3),'Color','g');
set(plot4(4),'Color','b');
set(plot4(5),'Color','y');

xlabel('Time (ms)');
title('Membrane Current');
ylabel('(uA/cm^2)');
box(axes4,'on');
legend1 = legend(axes4,'show');




axes6 = subplot(3,5,8);
plot6 = plot(node_location, data_matrix_front{1, tt./plotinterval}(:,1));
xlabel('coordinate (cm)');
title('Action Potential: Front');
ylabel('mV');
axis([0 inf voltage_lower voltage_upper]);
box(axes6,'on');

axes7 = subplot(3,5,9);
plot7 = plot(node_location, data_matrix_back{1, tt./plotinterval}(:,1));
xlabel('coordinate (cm)');
title('Action Potential: Back');
ylabel('mV');
axis([0 inf voltage_lower voltage_upper]);
box(axes7,'on');


axes8 = subplot(3,5,13);
plot8 = plot(node_location, data_matrix_front{1, tt./plotinterval}(:,2:5));
xlabel('coordinate (cm)');
title('Auxiliary Variables: Front');
set(plot8(1),'DisplayName','m');
set(plot8(2),'DisplayName','h');
set(plot8(3),'DisplayName','n');
set(plot8(4),'DisplayName','p');
set(plot8(1),'Color','k');
set(plot8(2),'Color','r');
set(plot8(3),'Color','g');
set(plot8(4),'Color','b');
axis([0 inf -0.2 1.2]);
box(axes8,'on');


axes9 = subplot(3,5,14);
plot6 = plot(node_location, data_matrix_back{1, tt./plotinterval}(:,2:5));
xlabel('coordinate (cm)');
title('Auxiliary Variables: Back');
set(plot6(1),'DisplayName','m');
set(plot6(2),'DisplayName','h');
set(plot6(3),'DisplayName','n');
set(plot6(4),'DisplayName','p');
set(plot6(1),'Color','k');
set(plot6(2),'Color','r');
set(plot6(3),'Color','g');
set(plot6(4),'Color','b');
axis([0 inf -0.2 1.2]);
box(axes9,'on');


axes10 = subplot(3,5,11);
plot10 = plot(node_location, flux_value_front(:,1:4));
xlabel('coordinate (cm)');
title('Membrane Current: Front');
ylabel('(uA/cm^2)');
set(plot10(1),'DisplayName','Leak');
set(plot10(2),'DisplayName','Na');
set(plot10(3),'DisplayName','Kd');
set(plot10(4),'DisplayName','K_slow');
set(plot10(1),'Color','r');
set(plot10(2),'Color','g');
set(plot10(3),'Color','b');
set(plot10(4),'Color','y');

box(axes10,'on');


axes11 = subplot(3,5,12);
plot11 = plot(node_location, flux_value_back(:,1:4));
xlabel('coordinate (cm)');
title('Membrane Current: Back');
ylabel('(uA/cm^2)');
set(plot11(1),'DisplayName','Leak');
set(plot11(2),'DisplayName','Na');
set(plot11(3),'DisplayName','Kd');
set(plot11(4),'DisplayName','K_slow');
set(plot11(1),'Color','r');
set(plot11(2),'Color','g');
set(plot11(3),'Color','b');
set(plot11(4),'Color','y');

box(axes11,'on');


axes5 = subplot(3,5,10);
cla(axes5);
set(axes5 ,'Visible','off')
text(0,1,['Current Simulated Steps:' num2str(tt)]);
text(0,0.9,['Current Simulated Time:' num2str(totalt(tt)) 'ms']);
text(0,0.8,['Wall Time of the Current Step:' num2str(toc(looptic)) 's']);
text(0,0.7,['Total Wall Time:' num2str(toc(overalltic)) 's']);
text(0,0.6,['Total 1D Time:' num2str(tot_toc_one_d) 's']);
text(0,0.5,['Total 2D Time:' num2str(tot_toc_two_d) 's']);
text(0,0.4,['Total Plot Time:' num2str(tot_toc_plot) 's']);
text(0,0.3,['Monitored Location: x' num2str(topnodalvalue(MonitoredNode,1)) ' y'  num2str(topnodalvalue(MonitoredNode,2))] );
text(0,0.2,['Total 2D Mesh Generated:' num2str(size(model.Mesh.Nodes,2))] );
text(0,0.1,['Max Mesh Size Allowed:' num2str(minmesh.*length_unit_ratio) 'um'] );
text(0,0.0,['1D mesh size:' num2str(node_dist) 'cm'] );
text(0,-0.1,['Current input:' num2str(currentin) 'uA'] );
text(0,-0.2,['Current input duration:' num2str(currentinduration) 'ms'] );





      frame = getframe(5);
      im = frame2im(frame);
      [imind,cm] = rgb2ind(im,256);
      if tt == plotinterval
          imwrite(imind,cm,filename,'gif', 'DelayTime',0.01,'Loopcount',inf);
      else
          imwrite(imind,cm,filename,'gif','DelayTime',0.01,'WriteMode','append');
      end
      
end

tot_toc_plot = tot_toc_plot + toc(plot_tic);

waitbar(tt / steps, h,sprintf('Current Step: %d, Simulated Time: %f',tt, totalt(tt)))    ;
end

close(h)
datalog{12} = topnodalvalue(:,1:4);
save data.mat datalog data_matrix_back data_matrix_front
save all.mat 
%createfigure_membrane_voltage(datalog{1}(:,1), datalog{2}(:,60))
%createfigure_aux(datalog{1}(:,1),[datalog{3}(:,60) datalog{4}(:,60) datalog{5}(:,60) datalog{6}(:,60)])
%createfigure_current(datalog{1}(:,1),[datalog{7}(:,60) datalog{8}(:,60) datalog{9}(:,60) datalog{10}(:,60) datalog{11}(:,60)])


