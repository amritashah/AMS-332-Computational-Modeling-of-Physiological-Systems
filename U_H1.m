w = 0.5;
u = 0.1;
Xp = 1;
Xr = 1;
Xprot = 0:20;
Xrna1 = (Xp/w).*Xprot;
Xrna2 = ones(size(Xprot));
Xrna2 = Xrna2.*(Xr/u);

plot(Xprot, Xrna1)
hold on
plot(Xprot, Xrna2)

xlabel('[X_{prot}]')
ylabel('[X_{rna}]')
legend('[X_{rna}] = (X_{prot} / w) [X_{prot}]','[X_{rna}] = X_{rna} / u', 'Location', 'Northwest')
hold off

%%
clear
w = 0.5;
u = 5;
Xp = 1;
Xr = 1;
K = 1;
h = 1;

Xprot = 0:0.1:10;
Xrna1 = (Xp/w).*Xprot;
Xrna2 = (u/Xr).*((Xprot.^h) ./ ((K^h) + (Xprot.^h)));

plot(Xprot, Xrna1)
hold on
plot(Xprot, Xrna2)

xlabel('[X_{prot}]')
ylabel('[X_{rna}]')
hold off