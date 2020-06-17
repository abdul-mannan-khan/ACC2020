global r n_flips t0 max_w
h1 = create_axis([1 1 1;2 2 2],15);
h2 = create_axis([1 1 1],15,'BottomMargin',0.2);
h3 = create_axis([1 1 1],15,'BottomMargin',0.2);
for s = 5e1%logspace(1,3,3)
    [t,j,xi] = noise(s);
    norm_e = zeros(numel(t),2);
    e = zeros(6,numel(t),2);
    inner_r_rd = zeros(numel(t),2);
    e_omega = zeros(numel(t),2);
    T = zeros(numel(t),2);
    for I = 1:numel(t)
        [pd,vd,ddot_pd,d3dot_pd] = reference(t(I));
        e(:,I,1) = xi(I,1:6)'-[pd;vd];
        norm_e(I,1) = sqrt(sum(e(1:3,I,1).^2));
        R1 = reshape(xi(I,7:15)',[3 3]);
        inner_r_rd(I,1) = 1-r'*R1'*rd(ddot_pd,e(:,I,1));
        T(I,1) = T_1(ddot_pd,e(:,I,1),R1);
        R2 = reshape(xi(I,27:35)',[3 3]);
        e(:,I,2) = xi(I,21:26)'-[pd;vd];
        T(I,2) = T_2(ddot_pd,e(:,I,2));
        norm_e(I,2) = sqrt(sum(e(1:3,I,2).^2));
        q1 = xi(I,16);
        q2 = xi(I,36:38)';
        inner_r_rd(I,2) = 1-r'*R2'*rd(ddot_pd,e(:,I,2));%V2(R2'*rd(ddot_pd,e(:,I,2)),q2);
        [~,omegad1] = F1_torque(xi(I,1:20)',xi(I,43),xi(I,44));
        [~,omegad2] = F2_torque(xi(I,21:42)',xi(I,45),xi(I,46));
        e_omega(I,1) = norm(xi(I,17:19)');
        e_omega(I,2) = norm(xi(I,39:41)');
    end
    %
    axes(h3)
    t1 = t0+n_flips*2*pi/max_w;
    ylim = enlarge([min(min(T,[],1)) max(max(T,[],1))],1.1);
    patch([t0 t1 t1 t0],[ylim(1) ylim(1) ylim(2) ylim(2)],0.9*ones(1,3),'edgecolor','none')
    hold on
    leg0 = plot(t,T,'linewidth',1);
    set(leg0(1),'linewidth',2)
    grid on
    ylabel('$Thrust$ [N]')
    xlabel('$t [s]$')
    set(gca,'ylim',ylim)
    %
    axes(h2)
    t1 = t0+n_flips*2*pi/max_w;
    ylim = enlarge([min(min(norm_e,[],1)) max(max(norm_e,[],1))],1.1);
    patch([t0 t1 t1 t0],[ylim(1) ylim(1) ylim(2) ylim(2)],0.9*ones(1,3),'edgecolor','none')
    hold on
    leg1 = plot(t,norm_e,'linewidth',1);
    set(leg1(1),'linewidth',2)
    grid on
    ylabel('$\norm{\e[p]}$ [m]')
    xlabel('$t [s]$')
    set(gca,'ylim',ylim)
    %
    axes(h1(1))
    ylim = [-0.1 2.1];
    patch([t0 t1 t1 t0],[ylim(1) ylim(1) ylim(2) ylim(2)],0.9*ones(1,3),'edgecolor','none')
    hold on
    leg2 = plot(t,inner_r_rd,'linewidth',1);
    leg2(1).LineWidth = 2;
    grid on
    set(gca,'xlim',[0.8 2],'xticklabel','','ylim',ylim)
    ylabel('$1-r\tp\R\rd(\ddot\pd,\e)$')
    %
    axes(h1(2))
    ylim = enlarge([min(min(e_omega,[],1)) max(max(e_omega,[],1))],1.1);
    patch([t0 t1 t1 t0],[ylim(1) ylim(1) ylim(2) ylim(2)],0.9*ones(1,3),'edgecolor','none')
    hold on
    leg3 = plot(t,e_omega,'linewidth',1);
    leg3(1).LineWidth = 2;
    grid on
    set(gca,'xlim',[0.8 2],'ylim',ylim)
    ylabel('$\norm{\omega}$ [rad/s]')
    xlabel('$t [s]$')
end
[~,~,~,txt] = legend(leg3,{'Section V-A','Section V-B'},...
    'position', [0.6208    0.3357    0.3603    0.1599]);
axes(h2)
[~,~,~,txt] = legend(leg1,{'Section V-A','Section V-B'},...
    'position',[0.5944    0.6572    0.3889    0.3074]);
[~,~,~,txt] = legend(leg0,{'Section V-A','Section V-B'},...
    'position',[0.5873    0.4248    0.3889    0.3074]);