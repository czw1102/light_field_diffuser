function [vk, f] = LF_ADMM_solver(psf,b,solverSettings)
% ADMM solver to compute 4D light fields
% H: Impulse stack 4D array
% b: measurement image from camera
% solverSettings: user defined params.
%
% The original version was developed for 3D imaging by Nick Antipa et al.
% The current version was modified for 4D light-field imaging.
% Modified by Zewei Cai - 2020-07-14   
%

assert(size(psf,1) == size(b,1) || size(psf,2) == size(b,2),'image and impulse have different dimensions');

if ~isfield(solverSettings,'print_interval')
    solverSettings.print_interval = 1;
end

if ~isfield(solverSettings,'disp_figs')
    solverSettings.disp_figs = 1;
end

if ~isfield(solverSettings,'autotune')
    solverSettings.autotune = 1;
end

if ~isfield(solverSettings,'maxIter')
    solverSettings.maxIter = 200;
end

if ~isfield(solverSettings,'regularizer')
    solverSettings.regularizer = 'TV';
end

if ~isfield(solverSettings,'display_func')
    solverSettings.display_func = @(x)x;
end

if ~isfield(solverSettings,'fig_num')
    solverSettings.fig_num = 1;
end

if solverSettings.disp_figs ~= 0
    solverSettings.fighandle = figure(solverSettings.fig_num);
    clf
end


mu1 = solverSettings.mu1;   %Set initial ADMM parameters
mu2 = solverSettings.mu2;
mu3 = solverSettings.mu3;

[Nt, Ns, Nv, Nu] = size(psf);   %Get problem size

% Setup convolutional forward op
vec = @(X)reshape(X,numel(X),1);
crop4d = @(x)x(:,:,1,1);
pad4d = @(x)padarray(x,[0 0 Nv-1 Nu-1],'post');
psf = circshift(flip(flip(psf,3),4),[0,0,ceil(Nv/2)+1,ceil(Nu/2)+1])/norm(psf(:));  %Shift impulse stack
Hs = fftn(ifftshift(psf));
Hs_conj = conj(Hs);
clear psf
Hfor = @(x)real(ifftshift(ifftn(Hs.*fftn(ifftshift(x)))));
Hadj = @(x)real(ifftshift(ifftn(Hs_conj.*fftn(ifftshift(x)))));
HtH = abs(Hs.*Hs_conj);

vk = zeros(Nt,Ns,Nv,Nu);   %Initialize variables. vk is the primal (this is the image you want to find)
xi = zeros(size(vk));  % Dual associated with Mv = nu (boundary condition variables)
rho = zeros(size(vk));  % Dual associated with v = w   (nonnegativity)
Dtb = pad4d(b);

switch lower(solverSettings.regularizer)
    case('tv')
        PsiTPsi = generate_laplacian(Nt, Ns, Nv, Nu);
        eta_1 = zeros(Nt-1,Ns,Nv,Nu);  %Duals associatd with Psi v = u (TV sparsity)
        eta_2 = zeros(Nt,Ns-1,Nv,Nu);
        eta_3 = zeros(Nt,Ns,Nv-1,Nu);
        eta_4 = zeros(Nt,Ns,Nv,Nu-1);
        PsiT = @(P1,P2,P3,P4)cat(1,P1(1,:,:,:),diff(P1,1,1),-P1(end,:,:,:)) + ...
            cat(2,P2(:,1,:,:),diff(P2,1,2),-P2(:,end,:,:)) + ...
            cat(3,P3(:,:,1,:),diff(P3,1,3),-P3(:,:,end,:)) +...
            cat(4,P4(:,:,:,1),diff(P4,1,4),-P4(:,:,:,end));
       
        % Sparsifying map
        Psi = @(x)deal(-diff(x,1,1),-diff(x,1,2),-diff(x,1,3),-diff(x,1,4));
        [uk1, uk2, uk3, uk4] = Psi(zeros(Nt, Ns,Nv,Nu));
        Lvk1 = uk1;
        Lvk2 = uk2;
        Lvk3 = uk3;
        Lvk4 = uk4;
        
    case('tv_native')
        PsiTPsi = generate_laplacian(Nt, Ns, Nv, Nu);
        PsiT = @(P1,P2,P3,P4,P5)cat(1,P1(1,:,:,:),diff(P1,1,1),-P1(end,:,:,:)) + ...
            cat(2,P2(:,1,:,:),diff(P2,1,2),-P2(:,end,:,:)) + ...
            cat(3,P3(:,:,1,:),diff(P3,1,3),-P3(:,:,end,:)) + ...
            cat(4,P4(:,:,:,1),diff(P4,1,4),-P4(:,:,:,end)) + ...
            solverSettings.tau_n*P5;
                
        %Sparsifying with gradient and l1
        Psi = @(x)deal(-diff(x,1,1),-diff(x,1,2),-diff(x,1,3),-diff(x,1,4),solverSettings.tau_n*x);
        [uk1, uk2, uk3, uk4, uk5] = Psi(zeros(Nt,Ns,Nv,Nu));
        Lvk1 = uk1;
        Lvk2 = uk2;
        Lvk3 = uk3;
        Lvk4 = uk4;
        Lvk5 = uk5;
        eta_1 = zeros(Nt-1,Ns,Nv,Nu);  %Duals associatd with Psi v = u (TV sparsity)
        eta_2 = zeros(Nt,Ns-1,Nv,Nu);
        eta_3 = zeros(Nt,Ns,Nv-1,Nu);
        eta_4 = zeros(Nt,Ns,Nv,Nu-1);
        eta_5 = zeros(Nt,Ns,Nv,Nu);
        PsiTPsi = PsiTPsi + solverSettings.tau_n^2;
        
    case('native')        
        PsiTPsi = 1;
        PsiT = @(x)x;
        Psi = @(x)x;   %Identity operator for native sparsity
        uk = zeros(size(vk));
        Lvk = uk;
        eta = uk;
end

v_mult = 1./(mu1*HtH + mu2*PsiTPsi + mu3);  %Denominator of v update (in 4D frequency space)

DtD = pad4d(ones(size(b)));   %D'D(b)
nu_mult = 1./(DtD + mu1);   %denominator of nu update

n = 0;  %Initialize number of steps to 0

% Store solver parameters in structure, f
% Initialize residuals with NaNs
f.dual_resid_s = zeros(1,solverSettings.maxIter)./0;   
f.primal_resid_s = zeros(1,solverSettings.maxIter)./0;
f.dual_resid_u = f.dual_resid_s;
f.primal_resid_u = f.dual_resid_u;
f.dual_resid_w = f.dual_resid_s;
f.primal_resid_w = f.dual_resid_s;
f.objective = f.primal_resid_u;   
f.data_fidelity = f.primal_resid_u;
f.regularizer_penalty = f.primal_resid_u;
f.run_time = f.dual_resid_s;

Hvkp = zeros(size(vk));

while n<solverSettings.maxIter
    
    tic;    
    n = n+1;
    Hvk = Hvkp;
    nukp = nu_mult.*(mu1*(xi/mu1 + Hvk) + Dtb);
    wkp = max(rho/mu3 + vk,0);
    switch lower(solverSettings.regularizer)
        case('tv')
            [uk1, uk2, uk3, uk4] = soft_threshold_4D(Lvk1+eta_1/mu2, Lvk2+eta_2/mu2, Lvk3+eta_3/mu2, Lvk4+eta_4/mu2,solverSettings.tau/mu2);
            vkp_numerator = mu3*(wkp-rho/mu3) + ...
                mu2*PsiT(uk1 - eta_1/mu2,uk2 - eta_2/mu2, uk3 - eta_3/mu2, uk4 - eta_4/mu2) + ...
                mu1*Hadj(nukp - xi/mu1);
            
        case('tv_native')
            [uk1, uk2, uk3, uk4, uk5] = soft_threshold_4D(Lvk1 + eta_1/mu2, Lvk2 + eta_2/mu2, ...
                Lvk3 + eta_3/mu2, Lvk4 + eta_4/mu2, solverSettings.tau/mu2, Lvk5 + eta_5/mu2);
            vkp_numerator = mu3*(wkp-rho/mu3) + ...
                mu2*PsiT(uk1 - eta_1/mu2,uk2 - eta_2/mu2, uk3 - eta_3/mu2, uk4 - eta_4/mu2, uk5 - eta_5/mu2) + ...
                mu1*Hadj(nukp - xi/mu1);
            
        case('native')
            uk = soft_threshold_4D([],[],[],[],solverSettings.tau_n/mu2,Lvk + eta/mu2);
            vkp_numerator = mu3*(wkp-rho/mu3) + mu2*PsiT(uk - eta/mu2) + mu1*Hadj(nukp - xi/mu1);
    end
    
    
    vkp = real(ifftshift(ifftn(v_mult .* fftn(ifftshift(vkp_numerator)))));
    
    %Update dual and parameter for Hs=v constraint
    Hvkp = Hfor(vkp);
    r_sv = Hvkp-nukp;
    xi = xi + mu1*r_sv;
    f.dual_resid_s(n) = mu1*norm(vec(Hvk - Hvkp));
    f.primal_resid_s(n) = norm(vec(r_sv));
    [mu1, mu1_update] = update_param(mu1,solverSettings.resid_tol,solverSettings.mu_inc,solverSettings.mu_dec,f.primal_resid_s(n),f.dual_resid_s(n));
    
    % Update dual and parameter for Ls=v
    f.data_fidelity(n) = .5*norm(crop4d(Hvkp)-b,'fro')^2;
    switch lower(solverSettings.regularizer)
        case('tv')
            Lvk1_ = Lvk1;
            Lvk2_ = Lvk2;
            Lvk3_ = Lvk3;
            Lvk4_ = Lvk4;
            [Lvk1, Lvk2, Lvk3, Lvk4] = Psi(vkp);
            r_su_1 = Lvk1 - uk1;
            r_su_2 = Lvk2 - uk2;
            r_su_3 = Lvk3 - uk3;
            r_su_4 = Lvk4 - uk4;
            eta_1 = eta_1 + mu2*r_su_1;
            eta_2 = eta_2 + mu2*r_su_2;
            eta_3 = eta_3 + mu2*r_su_3;
            eta_4 = eta_4 + mu2*r_su_4;
            f.dual_resid_u(n) = mu2*sqrt(norm(vec(Lvk1_ - Lvk1))^2 + norm(vec(Lvk2_ - Lvk2))^2 + norm(vec(Lvk3_ - Lvk3))^2 + norm(vec(Lvk4_ - Lvk4))^2);
            f.primal_resid_u(n) = sqrt(norm(vec(r_su_1))^2 + norm(vec(r_su_2))^2 + norm(vec(r_su_3))^2 + norm(vec(r_su_4))^2);
            f.regularizer_penalty(n) = solverSettings.tau*(sum(vec(abs(Lvk1))) + sum(vec(abs(Lvk2))) + sum(vec(abs(Lvk3))) + sum(vec(abs(Lvk4))));
            
        case('tv_native')
            Lvk1_ = Lvk1;
            Lvk2_ = Lvk2;
            Lvk3_ = Lvk3;
            Lvk4_ = Lvk4;
            Lvk5_ = Lvk5;
            [Lvk1, Lvk2, Lvk3, Lvk4, Lvk5] = Psi(vkp);
            r_su_1 = Lvk1 - uk1;
            r_su_2 = Lvk2 - uk2;
            r_su_3 = Lvk3 - uk3;
            r_su_4 = Lvk4 - uk4;
            r_su_5 = Lvk5 - uk5;
            eta_1 = eta_1 + mu2*r_su_1;
            eta_2 = eta_2 + mu2*r_su_2;
            eta_3 = eta_3 + mu2*r_su_3;
            eta_4 = eta_4 + mu2*r_su_4;
            eta_5 = eta_5 + mu2*r_su_5;
            f.dual_resid_u(n) = mu2*sqrt(norm(vec(Lvk1_ - Lvk1))^2 + norm(vec(Lvk2_ - Lvk2))^2 + ...
                norm(vec(Lvk3_ - Lvk3))^2 + norm(vec(Lvk4_ - Lvk4))^2 + norm(vec(Lvk5_ - Lvk5))^2);
            f.primal_resid_u(n) = sqrt(norm(vec(r_su_1))^2 + norm(vec(r_su_2))^2 + ...
                norm(vec(r_su_3))^2 + norm(vec(r_su_4))^2 + norm(vec(r_su_5))^2);
            f.regularizer_penalty(n) = solverSettings.tau*(sum(vec(abs(Lvk1))) + sum(vec(abs(Lvk2))) + ...
               sum(vec(abs(Lvk3))) + sum(vec(abs(Lvk4)))) + solverSettings.tau_n*sum(vec(abs(Lvk5)));
           
        case('native')
            Lvk_ = Lvk;
            Lvk = Psi(vkp);
            r_su = Lvk - uk;
            eta = eta + mu2*r_su;
            f.dual_resid_u(n) = mu2*norm(vec(Lvk_ - Lvk));
            f.primal_resid_u(n) = norm(vec(r_su));
            f.regularizer_penalty(n) = solverSettings.tau_n*(sum(vec(abs(Lvk))));
    end
    f.objective(n) = f.data_fidelity(n) + f.regularizer_penalty(n);
    
    
    [mu2, mu2_update] = update_param(mu2,solverSettings.resid_tol,...
        solverSettings.mu_inc,solverSettings.mu_dec,...
        f.primal_resid_u(n),f.dual_resid_u(n));
    
    % Update nonnegativity dual and parameter (s=w)
    r_sw = vkp-wkp;
    rho = rho + mu3*r_sw;
    f.dual_resid_w(n) = mu3*norm(vec(vk - vkp));
    f.primal_resid_w(n) = norm(vec(r_sw));
    [mu3, mu3_update] = update_param(mu3,solverSettings.resid_tol,solverSettings.mu_inc,solverSettings.mu_dec,f.primal_resid_w(n),f.dual_resid_w(n));
    
    %Update filters
    if mu1_update || mu2_update || mu3_update
        mu_update = 1;
    else
        mu_update = 0;
    end
    if mu_update
        v_mult = 1./(mu1*HtH + mu2*PsiTPsi + mu3);  %This is the frequency space division fo S update
        nu_mult = 1./(DtD + mu1);
    end
    
    
    vk = vkp;
    
    if mod(n,solverSettings.disp_figs) == 0
        draw_figures(vk,solverSettings)
    end
    
    f.run_time(n) = toc;
    
    if mod(n,solverSettings.print_interval) == 0
        time_iter = sum(f.run_time(n-solverSettings.print_interval+1:n));
         fprintf('iter: %i \t t: %.2g \t cost: %.2g \t data_fidelity: %.2g \t norm: %.2g \t Primal v: %.2g \t Dual v: %.2g \t Primal u: %.2g \t Dual u: %.2g \t Primal w: %.2g \t Dual w: %.2g \t mu1: %.2g \t mu2: %.2g \t mu3: %.2g \n',...
            n,time_iter,f.objective(n),f.data_fidelity(n),f.regularizer_penalty(n),f.primal_resid_s(n), f.dual_resid_s(n),f.primal_resid_u(n), f.dual_resid_u(n),f.primal_resid_w(n), f.dual_resid_w(n),mu1,mu2,mu3)
        
    end
end

end


% Private function to display figures
function draw_figures(xk, solverSettings)

set(0,'CurrentFigure',solverSettings.fighandle)
img = sum(sum(xk, 3), 4);
imagesc(solverSettings.display_func(img));
axis image
colormap gray
colorbar
caxis([0 prctile(img(:),solverSettings.disp_percentile)])
set(gca,'fontSize',10)
axis off
drawnow

end


function PsiTPsi = generate_laplacian(Nt,Ns,Nv,Nu)

lapl = zeros(Nt,Ns,Nv,Nu);    %Compute laplacian in closed form. This is the kernal to compute Psi'Psi
lapl(1) = 8;    
lapl(2,1,1,1) = -1;
lapl(end,1,1,1) = -1;
lapl(1,2,1,1) = -1;
lapl(1,end,1,1) = -1;
lapl(1,1,2,1) = -1;
lapl(1,1,end,1) = -1;
lapl(1,1,1,2) = -1;
lapl(1,1,1,end) = -1;
PsiTPsi = abs(fftn(lapl));   %Compute power spectrum of laplacian

end


function [mu_out, mu_update] = update_param(mu,resid_tol,mu_inc,mu_dec,r,s)

if r > resid_tol*s
    mu_out = mu*mu_inc;
    mu_update = 1;
elseif r*resid_tol < s
    mu_out = mu/mu_dec;
    mu_update = -1;
else
    mu_out = mu;
    mu_update = 0;
end

end

    
function [varargout] = soft_threshold_4D(v,h,d,k,tau,varargin)

%Perform isotropic soft thresholding on LF differences, v, h, d, and k
%using parameter tau. If a 5th input it added, assume it's the original
%LF and soft threshold that as well (for TV + sparsity).

if size(v,1) ~= 0  %If no v, h, d, or k passed in, skip gradient thresholding
    mag = sqrt(cat(1,v,zeros(1,size(v,2),size(v,3),size(v,4))).^2 + ...
        cat(2,h,zeros(size(h,1),1,size(h,3),size(h,4))).^2 + ...
        cat(3,d,zeros(size(d,1),size(d,2),1,size(d,4))).^2 +...
        cat(4,k,zeros(size(k,1),size(k,2),size(k,3),1)).^2);
    magt = soft_threshold(mag,tau);
    mmult = magt./mag;
    mmult(mag==0) = 0;
    varargout{1} = v.*mmult(1:end-1,:,:,:);
    varargout{2} = h.*mmult(:,1:end-1,:,:);
    varargout{3} = d.*mmult(:,:,1:end-1,:);
    varargout{4} = k.*mmult(:,:,:,1:end-1);
    if ~isempty(varargin)  %5th argument is native sparsity
        varargout{5} = soft_threshold(varargin{1},tau);
    end
else
    varargout{1} = soft_threshold(varargin{1},tau);
end

end


function threshed = soft_threshold(x,tau)

threshed = max(abs(x)-tau,0);
threshed = threshed.*sign(x);

end