function cmatrix = initBCimplicit_tail(region,state)
global topnodalvalue_tail
global count_tail
global V_init
nr = numel(region.x);
cmatrix = zeros(1,nr);
zmatrix1 = zeros(1,nr);
zmatrix2 = zeros(1,nr);
zmatrix3 = zeros(1,nr);
zmatrix4 = zeros(1,nr);

cmatrix(1,:) = V_init(1);
zmatrix1(1,:) = V_init(1);
zmatrix2(1,:) = V_init(2);
zmatrix3(1,:) = V_init(3);
zmatrix4(1,:) = V_init(4);
zmatrix5(1,:) = V_init(5);
topnodalvalue_tail(count_tail,:)=[region.x region.y region.nx region.ny zmatrix1 zmatrix2 zmatrix3 zmatrix4 zmatrix5];
count_tail=count_tail+1;


end

