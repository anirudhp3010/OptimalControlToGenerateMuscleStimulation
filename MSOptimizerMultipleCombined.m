%% Problem setup
toms t

Vrel_pro = out.Vrel_pro.data;
lce_pro = out.lce_pro.data;
Te_pro = out.Te_pro.data;
Ve_pro = out.Ve_pro.data;
STIM_pro = out.STIM_pro.data;
J_pro = trapz(out.tout, STIM_pro.^2);

%% Plot reference trajectory
figure()
plot(out.tout, STIM_pro);
ylabel('STIM')
xlabel('t (s)')
title('STIM reference trajectory');

figure()
subplot(4,1,1)
plot(out.tout, Vrel_pro);
hold on
ylabel('Vrel');
xlabel('t (s)');
title('Vrel reference trajectory');

subplot(4,1,2)
plot(out.tout, lce_pro);
hold on
ylabel('lce (m)');
xlabel('t (s)');
title('lce reference trajectory');

subplot(4,1,3)
plot(out.tout, Te_pro*180/pi);
hold on
ylabel('Te (deg)');
xlabel('t (s)');
title('Te reference trajectory');

subplot(4,1,4)
plot(out.tout, Ve_pro);
hold on
ylabel('Ve (rad/s)');
xlabel('t (s)');
title('Ve reference trajectory');

%% Solve the problem 1
nvec = [50];
for i=1:length(nvec);
    n = nvec(i);

    p=tomPhase('p',t,0,T_fin,n);
    setPhase(p);
    tomStates Vrel lce Te Ve
    tomControls STIM

    % Initial guess
    if i==1
        Vrelg = interp1(1:length(Vrel_pro), Vrel_pro, linspace(1,length(Vrel_pro),n));
        lceg = interp1(1:length(lce_pro), lce_pro, linspace(1,length(lce_pro),n));
        Teg = interp1(1:length(Te_pro), Te_pro, linspace(1,length(Te_pro),n));
        Veg = interp1(1:length(Ve_pro), Ve_pro, linspace(1,length(Ve_pro),n));
        STIMg = interp1(1:length(STIM_pro), STIM_pro, linspace(1,length(STIM_pro),n));
        xgs =  icollocate({ Vrel == Vrelg; lce == lceg; Te == Teg; Ve == Veg});
        xgu = collocate(STIM == STIMg);

        x0 = {xgs
            xgu};
    else
        x0 = {icollocate({Vrel == Vrel_opt; lce == lce_opt
            Te == Te_opt; Ve == Ve_opt})
            collocate(STIM == STIM_opt)};
    end

    % Box constraints
    cbox = {icollocate({0 <= Vrel <= 1 ...
            0 <= lce <= 0.1140 ...
            0 <= Te <= 90/180*pi ...
            0 <= Ve <= 2}) ...
            0 <= collocate(STIM) <= 1};

    % Boundary constraints
    cbnd = {initial({
            Vrel == 0;
            lce == 0.1140;
            Te == 0;
            Ve == 0})
            final({Te == ThetaFin})
            final({Ve == 0})
            };

    % ODEs and path constraints
    ceq = collocate({
            dot(Vrel) == m*(STIM - Vrel);
            dot(lce) == (brel*lce_opt*((625*(a2e*Te^2 + a1e*Te + a0 - lce - lse_0)^2)/lse_0^2 + ((q0 + (Vrel^3*c^3*lce^3*u^3*(k - 1)^3)/(lce_opt^3*(k - lce/lce_opt)^3))*(1/width^2 - (2*lce)/(lce_opt*width^2) + lce^2/(lce_opt^2*width^2) - 1))/((Vrel^3*c^3*lce^3*u^3*(k - 1)^3)/(lce_opt^3*(k - lce/lce_opt)^3) + 1)))/((625*(a2e*Te^2 + a1e*Te + a0 - lce - lse_0)^2)/lse_0^2 - (arel*(1/(exp(100 - (100*lce)/lce_opt) + 1) - 1/((exp((100*lce)/lce_opt - 100) + 1)*(1/width^2 - (2*lce)/(lce_opt*width^2) + lce^2/(lce_opt^2*width^2) - 1)))*(q0 + (Vrel^3*c^3*lce^3*u^3*(k - 1)^3)/(lce_opt^3*(k - lce/lce_opt)^3))*(1/width^2 - (2*lce)/(lce_opt*width^2) + lce^2/(lce_opt^2*width^2) - 1))/((Vrel^3*c^3*lce^3*u^3*(k - 1)^3)/(lce_opt^3*(k - lce/lce_opt)^3) + 1));
            dot(Te) == Ve;
            dot(Ve) == (3*(Cload*Ve + Kload*Te - (625*Fmax*a1e*(a2e*Te^2 + a1e*Te + a0 - lce - lse_0)^2)/lse_0^2 + (Llarm*Marm*g*sin(Te))/2 + Llarm*Mload*g*sin(Te)))/(Llarm^2*Marm);});

    % Objective
    J = STIM.^2;
    %J = (Te-ThetaFin).^2 + STIM.^2;
    %J = 2*(Te-ThetaFin).^2 + STIM.^2;
    %J = (Te-ThetaFin).^2 + 2*STIM.^2;
    objective = integrate(J);

    %% Solve the problem
    options = struct;
    options.name = 'Musculoskeletal Optimizer';
    if i==1
        options.solver = 'multiMin';
        options.xInit = 100;
    end
    %options.scale = 'auto'
    solution = ezsolve(objective, {cbox, cbnd, ceq}, x0, options);

    % Optimal v and more to use as starting guess
    Vrel_opt = subs(Vrel, solution);
    lce_opt = subs(lce, solution);
    Te_opt = subs(Te, solution);
    Ve_opt = subs(Ve, solution);
    STIM_opt = subs(STIM, solution);
end
t1 = subs(collocate(t),solution);
Vrel1 = subs(collocate(Vrel_opt),solution);
lce1 = subs(collocate(lce_opt),solution);
Te1 = subs(collocate(Te_opt),solution);
Ve1 = subs(collocate(Ve_opt),solution);
STIM1 = subs(collocate(STIM_opt),solution);
J_opt1 = trapz(t1, STIM1.^2);

%% Solve the problem 2
toms t
nvec = [50];
for i=1:length(nvec);
    n = nvec(i);

    p=tomPhase('p',t,0,T_fin,n);
    setPhase(p);
    tomStates Vrel lce Te Ve
    tomControls STIM

    % Initial guess
    if i==1
        Vrelg = interp1(1:length(Vrel_pro), Vrel_pro, linspace(1,length(Vrel_pro),n));
        lceg = interp1(1:length(lce_pro), lce_pro, linspace(1,length(lce_pro),n));
        Teg = interp1(1:length(Te_pro), Te_pro, linspace(1,length(Te_pro),n));
        Veg = interp1(1:length(Ve_pro), Ve_pro, linspace(1,length(Ve_pro),n));
        STIMg = interp1(1:length(STIM_pro), STIM_pro, linspace(1,length(STIM_pro),n));
        xgs =  icollocate({ Vrel == Vrelg; lce == lceg; Te == Teg; Ve == Veg});
        xgu = collocate(STIM == STIMg);

        x0 = {xgs
            xgu};
    else
        x0 = {icollocate({Vrel == Vrel_opt; lce == lce_opt
            Te == Te_opt; Ve == Ve_opt})
            collocate(STIM == STIM_opt)};
    end

    % Box constraints
    cbox = {icollocate({0 <= Vrel <= 1 ...
            0 <= lce <= 0.1140...
            0 <= Te <= 90/180*pi ...
            0 <= Ve <= 2}) ...
            0 <= collocate(STIM) <= 1};

    % Boundary constraints
    cbnd = {initial({
            Vrel == 0;
            lce == 0.1140;
            Te == 0;
            Ve == 0})
            final({Te == ThetaFin})
            final({Ve == 0})
            };

    % ODEs and path constraints
    ceq = collocate({
            dot(Vrel) == m*(STIM - Vrel);
            dot(lce) == (brel*lce_opt*((625*(a2e*Te^2 + a1e*Te + a0 - lce - lse_0)^2)/lse_0^2 + ((q0 + (Vrel^3*c^3*lce^3*u^3*(k - 1)^3)/(lce_opt^3*(k - lce/lce_opt)^3))*(1/width^2 - (2*lce)/(lce_opt*width^2) + lce^2/(lce_opt^2*width^2) - 1))/((Vrel^3*c^3*lce^3*u^3*(k - 1)^3)/(lce_opt^3*(k - lce/lce_opt)^3) + 1)))/((625*(a2e*Te^2 + a1e*Te + a0 - lce - lse_0)^2)/lse_0^2 - (arel*(1/(exp(100 - (100*lce)/lce_opt) + 1) - 1/((exp((100*lce)/lce_opt - 100) + 1)*(1/width^2 - (2*lce)/(lce_opt*width^2) + lce^2/(lce_opt^2*width^2) - 1)))*(q0 + (Vrel^3*c^3*lce^3*u^3*(k - 1)^3)/(lce_opt^3*(k - lce/lce_opt)^3))*(1/width^2 - (2*lce)/(lce_opt*width^2) + lce^2/(lce_opt^2*width^2) - 1))/((Vrel^3*c^3*lce^3*u^3*(k - 1)^3)/(lce_opt^3*(k - lce/lce_opt)^3) + 1));
            dot(Te) == Ve;
            dot(Ve) == (3*(Cload*Ve + Kload*Te - (625*Fmax*a1e*(a2e*Te^2 + a1e*Te + a0 - lce - lse_0)^2)/lse_0^2 + (Llarm*Marm*g*sin(Te))/2 + Llarm*Mload*g*sin(Te)))/(Llarm^2*Marm);});

    % Objective
    %J = STIM.^2;
    J = (Te-ThetaFin).^2 + STIM.^2;
    %J = 2*(Te-ThetaFin).^2 + STIM.^2;
    %J = (Te-ThetaFin).^2 + 2*STIM.^2;
    objective = integrate(J);

    %% Solve the problem
    options = struct;
    options.name = 'Musculoskeletal Optimizer';
    if i==1
        options.solver = 'multiMin';
        options.xInit = 100;
    end
    %options.scale = 'auto'
    solution = ezsolve(objective, {cbox, cbnd, ceq}, x0, options);

    % Optimal v and more to use as starting guess
    Vrel_opt = subs(Vrel, solution);
    lce_opt = subs(lce, solution);
    Te_opt = subs(Te, solution);
    Ve_opt = subs(Ve, solution);
    STIM_opt = subs(STIM, solution);
end
t2 = subs(collocate(t),solution);
Vrel2 = subs(collocate(Vrel_opt),solution);
lce2 = subs(collocate(lce_opt),solution);
Te2 = subs(collocate(Te_opt),solution);
Ve2 = subs(collocate(Ve_opt),solution);
STIM2 = subs(collocate(STIM_opt),solution);
J_opt2 = trapz(t2, (Te2-ThetaFin).^2 + STIM2.^2);

%% Solve the problem 3
toms t
nvec = [50];
for i=1:length(nvec);
    n = nvec(i);

    p=tomPhase('p',t,0,T_fin,n);
    setPhase(p);
    tomStates Vrel lce Te Ve
    tomControls STIM

    % Initial guess
    if i==1
        Vrelg = interp1(1:length(Vrel_pro), Vrel_pro, linspace(1,length(Vrel_pro),n));
        lceg = interp1(1:length(lce_pro), lce_pro, linspace(1,length(lce_pro),n));
        Teg = interp1(1:length(Te_pro), Te_pro, linspace(1,length(Te_pro),n));
        Veg = interp1(1:length(Ve_pro), Ve_pro, linspace(1,length(Ve_pro),n));
        STIMg = interp1(1:length(STIM_pro), STIM_pro, linspace(1,length(STIM_pro),n));
        xgs =  icollocate({ Vrel == Vrelg; lce == lceg; Te == Teg; Ve == Veg});
        xgu = collocate(STIM == STIMg);

        x0 = {xgs
            xgu};
    else
        x0 = {icollocate({Vrel == Vrel_opt; lce == lce_opt
            Te == Te_opt; Ve == Ve_opt})
            collocate(STIM == STIM_opt)};
    end

    % Box constraints
    cbox = {icollocate({0 <= Vrel <= 1 ...
            0 <= lce <= 0.1140...
            0 <= Te <= 90/180*pi ...
            0 <= Ve <= 2}) ...
            0 <= collocate(STIM) <= 1};

    % Boundary constraints
    cbnd = {initial({
            Vrel == 0;
            lce == 0.1140;
            Te == 0;
            Ve == 0})
            final({Te == ThetaFin})
            final({Ve == 0})
            };

    % ODEs and path constraints
    ceq = collocate({
            dot(Vrel) == m*(STIM - Vrel);
            dot(lce) == (brel*lce_opt*((625*(a2e*Te^2 + a1e*Te + a0 - lce - lse_0)^2)/lse_0^2 + ((q0 + (Vrel^3*c^3*lce^3*u^3*(k - 1)^3)/(lce_opt^3*(k - lce/lce_opt)^3))*(1/width^2 - (2*lce)/(lce_opt*width^2) + lce^2/(lce_opt^2*width^2) - 1))/((Vrel^3*c^3*lce^3*u^3*(k - 1)^3)/(lce_opt^3*(k - lce/lce_opt)^3) + 1)))/((625*(a2e*Te^2 + a1e*Te + a0 - lce - lse_0)^2)/lse_0^2 - (arel*(1/(exp(100 - (100*lce)/lce_opt) + 1) - 1/((exp((100*lce)/lce_opt - 100) + 1)*(1/width^2 - (2*lce)/(lce_opt*width^2) + lce^2/(lce_opt^2*width^2) - 1)))*(q0 + (Vrel^3*c^3*lce^3*u^3*(k - 1)^3)/(lce_opt^3*(k - lce/lce_opt)^3))*(1/width^2 - (2*lce)/(lce_opt*width^2) + lce^2/(lce_opt^2*width^2) - 1))/((Vrel^3*c^3*lce^3*u^3*(k - 1)^3)/(lce_opt^3*(k - lce/lce_opt)^3) + 1));
            dot(Te) == Ve;
            dot(Ve) == (3*(Cload*Ve + Kload*Te - (625*Fmax*a1e*(a2e*Te^2 + a1e*Te + a0 - lce - lse_0)^2)/lse_0^2 + (Llarm*Marm*g*sin(Te))/2 + Llarm*Mload*g*sin(Te)))/(Llarm^2*Marm);});

    % Objective
    %J = STIM.^2;
    %J = (Te-ThetaFin).^2 + STIM.^2;
    J = 2*(Te-ThetaFin).^2 + STIM.^2;
    %J = (Te-ThetaFin).^2 + 2*STIM.^2; 
    objective = integrate(J);

    %% Solve the problem
    options = struct;
    options.name = 'Musculoskeletal Optimizer';
    if i==1
        options.solver = 'multiMin';
        options.xInit = 100;
    end
    %options.scale = 'auto'
    solution = ezsolve(objective, {cbox, cbnd, ceq}, x0, options);

    % Optimal v and more to use as starting guess
    Vrel_opt = subs(Vrel, solution);
    lce_opt = subs(lce, solution);
    Te_opt = subs(Te, solution);
    Ve_opt = subs(Ve, solution);
    STIM_opt = subs(STIM, solution);
end
t3 = subs(collocate(t),solution);
Vrel3 = subs(collocate(Vrel_opt),solution);
lce3 = subs(collocate(lce_opt),solution);
Te3 = subs(collocate(Te_opt),solution);
Ve3 = subs(collocate(Ve_opt),solution);
STIM3 = subs(collocate(STIM_opt),solution);
J_opt3 = trapz(t3, 2*(Te3-ThetaFin).^2 + STIM3.^2);

%% Solve the problem 4
toms t
nvec = [50];
for i=1:length(nvec);
    n = nvec(i);

    p=tomPhase('p',t,0,T_fin,n);
    setPhase(p);
    tomStates Vrel lce Te Ve
    tomControls STIM

    % Initial guess
    if i==1
        Vrelg = interp1(1:length(Vrel_pro), Vrel_pro, linspace(1,length(Vrel_pro),n));
        lceg = interp1(1:length(lce_pro), lce_pro, linspace(1,length(lce_pro),n));
        Teg = interp1(1:length(Te_pro), Te_pro, linspace(1,length(Te_pro),n));
        Veg = interp1(1:length(Ve_pro), Ve_pro, linspace(1,length(Ve_pro),n));
        STIMg = interp1(1:length(STIM_pro), STIM_pro, linspace(1,length(STIM_pro),n));
        xgs =  icollocate({ Vrel == Vrelg; lce == lceg; Te == Teg; Ve == Veg});
        xgu = collocate(STIM == STIMg);

        x0 = {xgs
            xgu};
    else
        x0 = {icollocate({Vrel == Vrel_opt; lce == lce_opt
            Te == Te_opt; Ve == Ve_opt})
            collocate(STIM == STIM_opt)};
    end

    % Box constraints
    cbox = {icollocate({0 <= Vrel <= 1 ...
            0 <= lce <= 0.1140...
            0 <= Te <= 90/180*pi ...
            0 <= Ve <= 2}) ...
            0 <= collocate(STIM) <= 1};

    % Boundary constraints
    cbnd = {initial({
            Vrel == 0;
            lce == 0.1140;
            Te == 0;
            Ve == 0})
            final({Te == ThetaFin})
            final({Ve == 0})
            };

    % ODEs and path constraints
    ceq = collocate({
            dot(Vrel) == m*(STIM - Vrel);
            dot(lce) == (brel*lce_opt*((625*(a2e*Te^2 + a1e*Te + a0 - lce - lse_0)^2)/lse_0^2 + ((q0 + (Vrel^3*c^3*lce^3*u^3*(k - 1)^3)/(lce_opt^3*(k - lce/lce_opt)^3))*(1/width^2 - (2*lce)/(lce_opt*width^2) + lce^2/(lce_opt^2*width^2) - 1))/((Vrel^3*c^3*lce^3*u^3*(k - 1)^3)/(lce_opt^3*(k - lce/lce_opt)^3) + 1)))/((625*(a2e*Te^2 + a1e*Te + a0 - lce - lse_0)^2)/lse_0^2 - (arel*(1/(exp(100 - (100*lce)/lce_opt) + 1) - 1/((exp((100*lce)/lce_opt - 100) + 1)*(1/width^2 - (2*lce)/(lce_opt*width^2) + lce^2/(lce_opt^2*width^2) - 1)))*(q0 + (Vrel^3*c^3*lce^3*u^3*(k - 1)^3)/(lce_opt^3*(k - lce/lce_opt)^3))*(1/width^2 - (2*lce)/(lce_opt*width^2) + lce^2/(lce_opt^2*width^2) - 1))/((Vrel^3*c^3*lce^3*u^3*(k - 1)^3)/(lce_opt^3*(k - lce/lce_opt)^3) + 1));
            dot(Te) == Ve;
            dot(Ve) == (3*(Cload*Ve + Kload*Te - (625*Fmax*a1e*(a2e*Te^2 + a1e*Te + a0 - lce - lse_0)^2)/lse_0^2 + (Llarm*Marm*g*sin(Te))/2 + Llarm*Mload*g*sin(Te)))/(Llarm^2*Marm);});

    % Objective
    %J = STIM.^2;
    %J = (Te-ThetaFin).^2 + STIM.^2;
    %J = 2*(Te-ThetaFin).^2 + STIM.^2;
    J = (Te-ThetaFin).^2 + 2*STIM.^2;
    objective = integrate(J);

    %% Solve the problem
    options = struct;
    options.name = 'Musculoskeletal Optimizer';
    if i==1
        options.solver = 'multiMin';
        options.xInit = 100;
    end
    %options.scale = 'auto'
    solution = ezsolve(objective, {cbox, cbnd, ceq}, x0, options);

    % Optimal v and more to use as starting guess
    Vrel_opt = subs(Vrel, solution);
    lce_opt = subs(lce, solution);
    Te_opt = subs(Te, solution);
    Ve_opt = subs(Ve, solution);
    STIM_opt = subs(STIM, solution);
end
t4 = subs(collocate(t),solution);
Vrel4 = subs(collocate(Vrel_opt),solution);
lce4 = subs(collocate(lce_opt),solution);
Te4 = subs(collocate(Te_opt),solution);
Ve4 = subs(collocate(Ve_opt),solution);
STIM4 = subs(collocate(STIM_opt),solution);
J_opt4 = trapz(t4, (Te4-ThetaFin).^2 + 2*STIM4.^2);

%% Plot outputs
figure()
plot(t1, STIM1);
hold on
plot(t2, STIM2);
plot(t3, STIM3);
plot(t4, STIM4);
ylabel('STIM')
xlabel('t (s)')
title('STIM optimal control');
legend('Q0R1', 'Q1R1', 'Q2R1', 'Q1R2')

figure()
subplot(4,1,1)
plot(t1, Vrel1);
hold on
plot(t2, Vrel2);
plot(t3, Vrel3);
plot(t4, Vrel4);
ylabel('Vrel');
xlabel('t (s)');
title('Vrel optimal trajectory');
legend('Q0R1', 'Q1R1', 'Q2R1', 'Q1R2')

subplot(4,1,2)
plot(t1, lce1);
hold on
plot(t2, lce2);
plot(t3, lce3);
plot(t4, lce4);
ylabel('lce (m)');
xlabel('t (s)');
title('lce optimal trajectory');

subplot(4,1,3)
plot(t1, Te1*180/pi);
hold on
plot(t2, Te2*180/pi);
plot(t3, Te3*180/pi);
plot(t4, Te4*180/pi);
ylabel('Te (deg)');
xlabel('t (s)');
title('Te optimal trajectory');

subplot(4,1,4)
plot(t1, Ve1);
hold on
plot(t2, Ve2);
plot(t3, Ve3);
plot(t4, Ve4);
ylabel('Ve (rad/s)');
xlabel('t (s)');
title('Ve optimal trajectory');

%% Export STIM output

clear simin;
simin.time = t';
simin.signals.values = STIM;

a0 = 0.286;
a1e = -0.014;
a2e = -3.96e-3;

Llarm = 0.5;
Marm = 2.5;
Iarm = Marm*(Llarm^2)/3;
Mload = 5;
Cload = -10;
Kload = -10;
g = -9.8;
c = 1.373e-4;
u = 52700;
k = 2.90;
q0 = 0.005;
m = 11.25;

width = 0.66;
arel = 0.41;
brel = 5.2;
qcrit = 0.03;
a = 1/(width^2);
brel = 5.2;

lse_0 = 0.172; % Tendon natural length (m)
lce_opt = 0.092; % Optimal CE length (m)
lpe_0 = 1.4*lce_opt; % PE natural length (m)
%lpe_0 = 1*lce_opt; % PE natural length (m)
Fmax = 1420; % MVC value (N)

kse = Fmax/((0.04*lse_0)^2); % SE Stiffness (N/m)
kpe = (0.5*Fmax)/(1 + width - (lpe_0/lce_opt)); % PE Stiffness (N/m)

lmtc_init = a0 + (a1e*ThetaInit) + (a2e*(ThetaInit^2));
lce_init = lmtc_init-lse_0;
lce_rel_init = lce_init/lce_opt;