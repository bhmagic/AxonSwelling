function cmatrix = loopBC_1(region,state)
global topnodalvalue count2 dt count3 count4 countstupid g_na_alt g_kd_alt g_leak_alt
global g_leak g_na g_kd V_T g_m t_max E_na E_k E_leak Cm  cond axon_r

nr = numel(region.x);
flux_value = zeros(nr,5);
alpha_value = zeros(nr,4);
beta_value = zeros(nr,4);
topnodalvalue_temp = zeros(nr,9);
if region.x ~= topnodalvalue(count2,1) || region.y ~= topnodalvalue(count2,2)
    error('mismatch')
else
    
    
    

alpha_value(:,1) = -0.32*(state.u -V_T -13)./(exp(-(state.u -V_T -13)./4) -1);
beta_value(:,1) = 0.28*(state.u -V_T -40)./(exp((state.u -V_T -40)./5) -1);
alpha_value(:,2) = 0.128*exp(-(state.u -V_T -17)./18);
beta_value(:,2) = 4.0./(1.0+exp(-(state.u -V_T -40)./5));
alpha_value(:,3) = -0.032*(state.u -V_T -15)./(exp(-(state.u -V_T -15)./5) -1);
beta_value(:,3) = 0.5*exp(-(state.u -V_T -10)./40);
alpha_value(:,4) = 1.0./(1.0+exp(-(state.u +35)./10));
beta_value(:,4) = t_max./(3.3.*exp((state.u +35)./20) +exp(-(state.u+35)./20));



partial_t_value(:,2:4) = alpha_value(:,1:3).*(1 -topnodalvalue(count2,6:8)) -beta_value(:,1:3).*topnodalvalue(count2,6:8);
partial_t_value(:,5) = (alpha_value(:,4) -topnodalvalue(count2,9)) ./beta_value(:,4);

topnodalvalue_temp(:,6:9) =topnodalvalue(count2,6:9)+partial_t_value(:,2:5).*dt;    


if region.y > axon_r
local_g_na = g_na_alt;
local_g_kd = g_kd_alt;
local_g_leak = g_leak_alt;
else
local_g_na = g_na;
local_g_kd = g_kd;
local_g_leak = g_leak;
end

flux_value(:,2) = - local_g_leak.*(state.u-E_leak) ;
flux_value(:,3) = - local_g_na.*topnodalvalue_temp(:,6).^3.*topnodalvalue_temp(:,7).*(state.u-E_na);
flux_value(:,4) = - local_g_kd.*topnodalvalue_temp(:,8).^4.*(state.u-E_k) ;
flux_value(:,5) = - g_m.*topnodalvalue_temp(:,9).*(state.u-E_k);

    
cmatrix = (flux_value(:,2)+flux_value(:,3)+flux_value(:,4)+flux_value(:,5))./cond + topnodalvalue(count2,5).*Cm./cond./dt;
cmatrix = cmatrix.*region.y;
    
        if state.u ~= 0 || isnan(state.u)
    else
        countstupid = countstupid + 1;
        end
    
end

if count2 == count3
    count2 = 1;
    count4=count4+1;

else
count2=count2+1;
end

end

