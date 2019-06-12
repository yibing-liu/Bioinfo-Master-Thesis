embryoData = csvread(filename,1,2);

D = embryoData;

D(isnan(D)) = 0;
[m,n] = size(D);

% normalize
for i = 1:n
    l = round(m/4);
    d = medfilt1(D(:,i),l);
    D(:,i) = D(:,i) - d;
end

% estimate sigma
d = D(abs(D)<2);
sig = 1.48*median(abs(d(:)-median(d(:))));

% Run PLA
alpha1 = 0.5*sqrt(m)*sig;
alpha2 = 0.01*alpha1; %0.01*alpha1;
beta = 2*sig;
tol = 1e-6;

tic;
[B,E] = RPLA(D,alpha1,alpha2,beta,tol);
t_pla = toc;


lrc_filename = ['lrc_',filename];
csvwrite(lrc_filename, B);

sc_filename = ['sc_',filename];
csvwrite(sc_filename, E);


