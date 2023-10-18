


x0 = [12.1212121212121
12.1212121212121
4.84848484848485
12.1212121212121
12.1212121212121
2.42424242424242
9.69696969696970
2.42424242424242
9.69696969696970
4.84848484848485];

y0 = [36.3636363636364
25.4545454545455
22.0606060606061
6.54545454545455
19.6363636363636
16.9696969696970
38.5454545454546
34.1818181818182
30.5454545454545
26.4242424242424];

z0 = [5.09090909090909
4.12121212121212
38.3030303030303
39.0303030303030
32.2424242424242
36.8484848484849
26.4242424242424
37.5757575757576
29.8181818181818
7.03030303030303];

x = double((1:1:330));
y = double((1:1:165));
z = double((1:1:165));

x = x*(40/165); %from indices to real world coordinates in mm
y = y*(40/165);
z = z*(40/165);

x = transpose(x);
y = transpose(y);
z = transpose(z);

cd('C:\Users\iislam30\Documents\Bubble_DNS_Bruno\CODE')
load('p1.mat')
p1 = p_1;
cd('C:\Users\iislam30\Documents\Bubble_DNS_Bruno\CODE\Periodic Boundary Conditions')
[P] = GetPressure(p1,x0,y0,z0,x,y,z);




