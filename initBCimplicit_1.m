function cmatrix = initBCimplicit_1(region,state)
global topnodalvalue
global count2
global V_init
nr = numel(region.x);
cmatrix = zeros(1,nr);
zmatrix1 = zeros(1,nr);
zmatrix2 = zeros(1,nr);
zmatrix3 = zeros(1,nr);
zmatrix4 = zeros(1,nr);

cmatrix(1,:) = V_init(1).*1.0E3;
zmatrix1(1,:) = V_init(1);
zmatrix2(1,:) = V_init(2);
zmatrix3(1,:) = V_init(3);
zmatrix4(1,:) = V_init(4);
zmatrix5(1,:) = V_init(5);
topnodalvalue(count2,:)=[region.x region.y region.nx region.ny zmatrix1 zmatrix2 zmatrix3 zmatrix4 zmatrix5];
count2=count2+1;


end

