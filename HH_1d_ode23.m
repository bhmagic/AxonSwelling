
clear
clear global 

global dt node_dist tot_n_node axon_r axon_A const_1 const_2 currentin


V_init = [-70.6878    0.0    0.0    0.0    0.0]; %mV

length_unit_ratio = 10000.0; %cm/um


currentin=8E-4;
simulating_time = 5; %ms


total_length = 1; %cm
node_dist = 5; %um
node_dist = node_dist./length_unit_ratio; %cm

tot_n_node = ceil(total_length./node_dist);

axon_r(1:tot_n_node) = 2.0; %um
swelling_r = 50.0; %um
swelling_r = swelling_r./length_unit_ratio;
axon_r = axon_r./length_unit_ratio; %cm



for ii=1:1:tot_n_node
if abs(node_dist.*ii-total_length./2.0) < swelling_r 
    rrr = (swelling_r.*swelling_r-abs(node_dist.*ii-total_length./2.0).*abs(node_dist.*ii-total_length./2.0)).^0.5;
    if rrr > axon_r(ii)
    axon_r(ii) = rrr;
    end
end
end

axon_A = pi.*axon_r.*axon_r;
cond = 10.0; % mS/um +> mS/cm
const_1(1) = 2.0.*pi.*axon_r(1).*node_dist; %Total Suraface area
const_1(tot_n_node) = 2.0.*pi.*axon_r(tot_n_node).*node_dist; %Total Suraface area
for ii=2:1:tot_n_node-1
const_1(ii) = 2.0.*pi.*axon_r(ii).*(4.*node_dist.*node_dist+abs(axon_r(ii-1)-axon_r(ii+1)).^2.0).^0.5; %Total Suraface area
end
const_2 = 1.0./node_dist.*cond.*axon_A; %Conductivity 

initial_y(1:tot_n_node) = V_init(1);
initial_y(tot_n_node.*1+1:tot_n_node.*2) =  V_init(2);
initial_y(tot_n_node.*2+1:tot_n_node.*3) =  V_init(3);
initial_y(tot_n_node.*3+1:tot_n_node.*4) =  V_init(4);
initial_y(tot_n_node.*4+1:tot_n_node.*5) =  V_init(5);


for jj=0:1:ceil(simulating_time./2)-1
tic    
[t,y] = ode23(@HH_1d_ode,[0 2],initial_y);
for ii=1:1:200
y_slim(jj.*200+ii,:)=y(ii.*floor(size(t,1)./200),1:tot_n_node);
v1_slim(jj.*200+ii,:)=y(ii.*floor(size(t,1)./200),tot_n_node.*1+1:tot_n_node.*2);
v2_slim(jj.*200+ii,:)=y(ii.*floor(size(t,1)./200),tot_n_node.*2+1:tot_n_node.*3);
v3_slim(jj.*200+ii,:)=y(ii.*floor(size(t,1)./200),tot_n_node.*3+1:tot_n_node.*4);
v4_slim(jj.*200+ii,:)=y(ii.*floor(size(t,1)./200),tot_n_node.*4+1:tot_n_node.*5);
t_slim(jj.*200+ii)=t(ii.*floor(size(t,1)./200));

end
initial_y = y(size(y,1),:);

toc
end
%%




