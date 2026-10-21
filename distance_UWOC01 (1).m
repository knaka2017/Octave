l = 0:1:1000;

%% Lambert-Beer's law
%%   T = exp(-A*L)
%%   T : transmittance  = Intensity of incident light / Intensity of output light
%%   A : absorption coefficient [1/m]
%%   L : distance [m]

a = [0.0105 0.023025 1 4];  % absorption coef., no scattering coef.
%% 1. 純水(450nm)
%% 2. 東京湾入り口（観音岬沖） 2000年1月(488nm)
%% 3. 東京湾入り口（観音岬沖） 1999年7月(488nm)
%% 4. 東京湾奥中央（幕張沖） 1999年7月(488nm)
%% 1: F. M. Sogandares, et al., "Absorption spectrum (340-640 nm) of pure water. I. Photothermal measurements.", Applied optics, Vol. 36, No. 33, pp. 8699-8709 (1997)
%% 2,3,4: 成田美穂 他：“東京湾における溶存態有機物と懸濁物による濁度の分布および各起因物質との関係”，東京海洋大学研究報告，Vol. 2，pp. 35-46 （2006）

rpr = -20; %% dB, transmittance
y01 = log(10^(0.1*rpr))./-a;
gain = 20; %% dB
g = log(10^(0.1* gain))./a;
y10 = g+y01;

%figure(1);
%semilogy(1:length(a), g, 'o;;');
%xticks([1 2 3 4]);

%figure(2);
%semilogy(1:length(a), y01, 'o;;');
%xticks([1 2 3 4]);

figure(3);
%semilogy(1:length(a), y01, sprintf("+;%ddB of transmission;",rpr), 1:length(a), y10, sprintf("o;w/ %ddB gain of PMT;",gain));
semilogy(1:length(a), y01, sprintf("+;%ddB of transmission;",rpr));
ylabel('distance [m]');
xticks([1 2 3 4]);
grid on;
text((1:length(a))+0.05, y01, {sprintf("%6.2f m",y01(1)) sprintf("%6.2f m",y01(2)) sprintf("%6.2f m",y01(3)) sprintf("%6.2f m",y01(4))});
%text((1:length(a))+0.05, y10, {sprintf("%6.2f m",y10(1)) sprintf("%6.2f m",y10(2)) sprintf("%6.2f m",y10(3)) sprintf("%6.2f m",y10(4))});
%text((1:length(a))+0.08, y01+g/2.6, {sprintf("%+6.2f m",g(1)) sprintf("%+6.2f m",g(2)) sprintf("%+6.2f m",g(3)) sprintf("%+6.2f m",g(4))});
