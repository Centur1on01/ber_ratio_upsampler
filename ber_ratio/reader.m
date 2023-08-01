
signal = fopen('outup1.txt');
sig = textscan(signal,'%f\t%f');

x = real(sig{2});
y = imag(sig{2});

%sz = linspace(5, 5,1001);
%scatter(real(sig{1}),imag(sig{1}));
%figure;
%scatter(x,y); axis([-6 6 -6 6]);
%grid on;
%hold on;
%plot(x,-x);
%hold off;

berratio = fopen('out2.txt');
ber = textscan(berratio, '%f\t%f\t%f');
figure;
plot(ber{1}, ber{2});
hold on;
plot(ber{1}, ber{3});
hold off;

%real1 = fopen('real.txt');
%sig_re = textscan(real1,'%f\t%f\t\t\t%f');
%figure;
%plot(sig_re{1}, sig_re{2});
%hold on;
%plot(sig_re{1}, sig_re{3});
%hold off;

fclose(signal);
fclose(berratio);
%fclose(real1);


























