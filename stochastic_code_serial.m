% Parallelized microscale model code to simulate simple stochastic
%  dynamics of actin filaments including polymerization, depolymerization,
%  and branching from a single nucleation site.
% Options that can be turned on/off include:
% 1) capping of filaments
% 2) limited Arp23
% 3) limited actin monomers
%
% This code does "run_number" independent simulations of the micro model
% and then computes the average actin density across all those simulations
% Authors: Calina Copos, Minghao Rostami, Brittany Bannish, Kelsey Gasior, Adriana Dawes, Rebecca Pinals
% Last updated: 4/23/23

% Please contact Minghao W. Rostami for any problem you may encounter when
% running the code.

% When using this code, please cite the following reference:
%
%   MW Rostami, BE Bannish, K Gasior, RL Pinals, C Copos, and A Dawes.
%   Inferring local molecular dynamics from the global actin network
%   structure: a case study of 2D synthetic branching actin networks.

close all
clear all


run_number = 100; %number of independent simulations
box_size = 0.2; 
array_size=(10.0-box_size)/box_size;

Factin_array_T = zeros(array_size*array_size,run_number);

%variables to set to 1 to choose what to include in the model. a
%variable with a value of 1 will include the corresponding feature
%(capping, limited Arp2/3, etc.) in the model. Note: if
%Arp23yes and monomeryes are both set to 0, we have the unlimited
%resources case
capyes     = 1;       % 1 if we want to include capping, 0 if we don't
Arp23yes   = 0;       % 1 if we want to include limited Arp2/3, 0 if we don't
monomeryes = 1;       % 1 if we want to include limited actin monomer, 0 if we don't
plotnetworkyes = 1;      % 1 if we want to make a movie of the actin network, 0 if we don't


% Set discretization
%
dt     = 0.005;                  % temporal discretization
dt0    = 0.005;                  % (engineered) temporal discretization
Tend   = 10.0;                   % end time (a.u.)
Nt     = Tend/dt;                % number of timesteps
dx     = 0.0027;                 % size of actin monomer (um)
tplot  = 10;                     % visualize every 10 time units


% Set actin filament parameters
%
ppoly0   = 0.324*dt/dx;  % probability to polymerize
pdepoly = 0.0134*dt/dx;  % probability to depolymerize
Lbranch = 0.2;           % critical length (um)
mu0 = 2.0;               % controls mean of Gaussian distribution for branching probability
sigma0 = 1.0;            % controls width of Gaussian distribution for branching probability
arp23angle = 70*pi/180; % 70 degree branching angle by Arp2/3
pcap    = 0.05*dt*capyes; % probability of capping
NArp0 = 24;             % initial number of Arp2/3 (in limited case)
Nmonomers0 = 10000;       % number of G-actin monomers in the system

% Pre-allocation of vectors for speed (optional)
%
Nmax = 100000;



tic
for stats=1:run_number
    rng(stats,'twister')
    
    % Set discretization
    %
    
    b = zeros(Nmax,2);
    t = zeros(Nmax,2);
    branchloc = zeros(Nmax,2);
    L = zeros(Nmax,1);
    Lmin = zeros(Nmax,1);
    th = zeros(Nmax,1);
    lengthActin = zeros(Nt,Nmax);
    captip = ones(Nmax,1);
    tlim=zeros(Nt,Nmax);
    
    thetamat = zeros(Nt,Nmax);
    
    % Initialize nucleation site
    %
    N = 1;                                      % # filaments
    xnucl = 0.0; ynucl = 0.0;                   % nucleation site
    th(1) = rand(1,1)*2*pi;                     % angle of filament growth
    b(1,1) = xnucl; b(1,2) = ynucl;             % base location
    t(1,1) = b(1,1); t(1,2) = b(1,2);           % tip location
    branchloc(1,1) = 0.0; branchloc(1,2) = 0.0; % last branch location
    lengthActin(1,1) = sqrt( (t(1,1)-b(1,1))^2 + (t(1,2)-b(1,2))^2 );
    
    if stats==1 %only make a movie for the first run
        % following 4 pieces are for saving network movie
        w = VideoWriter('network_movie.avi');
        w.FrameRate = 10;
        w.Quality = 100;
        open(w);
        % end variables for network movie
    end
    
    % Plot domain
    %
    xmin = -5; xmax = 5;
    ymin = -5; ymax = 5;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                    Run microscale simulation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %tic
    
    NArp = NArp0;
    Nmonomers = Nmonomers0;
    for tau=1:Nt
        
        %generate random numbers to be used this time step
        rvec = rand(N,4);
                    
        % For each filament determine if it gets capped by capping protein
        captip(1:N) = (1-(rvec(:,4)<pcap)).*captip(1:N);
        
        if monomeryes==1
            % For each filament i grow or decay the ActinLength of the filament
            % (No preallocation -- growing list of filaments)
            ppoly = ppoly0*exp((Nmonomers-Nmonomers0)/Nmonomers0)*Nmonomers/Nmonomers0; %multiply by ratio so that probability is smaller when there's fewer monomers
            
            L(1:N) =  sqrt( (t(1:N,1)-b(1:N,1)).^2 + (t(1:N,2)-b(1:N,2)).^2 );
            Lmin(1:N) = sqrt( (t(1:N,1)-branchloc(1:N,1)).^2 + (t(1:N,2)-branchloc(1:N,2)).^2 );
            remove_check = sum(rvec(:,1)<ppoly);
            testvec = ((rvec(:,1)<ppoly).*(Nmonomers>=remove_check) - ((L(1:N)>dx)&(Lmin(1:N)>dx)).*(rvec(:,2)>=(1-pdepoly)))*dx.*captip(1:N); 
            Nmonomers = Nmonomers - remove_check*(Nmonomers>=remove_check) + sum(((L(1:N)>dx)&(Lmin(1:N)>dx)).*(rvec(:,2)>=(1-pdepoly)));
            t(1:N,1) = t(1:N,1) + cos(th(1:N)).*testvec;
            t(1:N,2) = t(1:N,2) + sin(th(1:N)).*testvec;
        else
            % For each filament i grow or decay the ActinLength of the filament
            % (No preallocation -- growing list of filaments)
            L(1:N) =  sqrt( (t(1:N,1)-b(1:N,1)).^2 + (t(1:N,2)-b(1:N,2)).^2 );
            Lmin(1:N) = sqrt( (t(1:N,1)-branchloc(1:N,1)).^2 + (t(1:N,2)-branchloc(1:N,2)).^2 );
            testvec = ((rvec(:,1)<ppoly0) - ((L(1:N)>dx)&(Lmin(1:N)>dx)).*(rvec(:,2)>=(1-pdepoly)))*dx.*captip(1:N);
            t(1:N,1) = t(1:N,1) + cos(th(1:N)).*testvec;
            t(1:N,2) = t(1:N,2) + sin(th(1:N)).*testvec;
        end %end if monomeryes==1
        
        
        if mod(tau,dt0/dt)==0
            % For each filament i decide branching event
            % (No preallocation -- growing list of filaments)
            
            if Arp23yes==1
                mu = mu0*(2*NArp0-NArp)/NArp0;
                sigma = sigma0;
                pbranchvec = normcdf(Lmin(1:N),mu,sigma).*(Lmin(1:N)>Lbranch)*NArp/NArp0; %multiply by ratio so that branching probability is lower when there's less free Arp2/3
                BranchSites = find(rvec(:,3)<pbranchvec);
                add_check = length(BranchSites)*(NArp>=length(BranchSites));
                NArp = NArp - add_check;
            else
                pbranchvec = normcdf(Lmin(1:N),mu0,sigma0).*(Lmin(1:N)>Lbranch);
                BranchSites = find(rvec(:,3)<pbranchvec);
            end %end if Arp23yes==1
            
            if ~isempty(BranchSites)
                Nadd = length(BranchSites);
                b(N+1:N+Nadd,1) = t(BranchSites,1);
                b(N+1:N+Nadd,2) = t(BranchSites,2);
                branchloc(BranchSites,1) = t(BranchSites,1);
                branchloc(BranchSites,2) = t(BranchSites,2);
                branchrvec = rand(Nadd,1);
                th(N+1:N+Nadd) = (branchrvec<0.5).*(th(BranchSites)+arp23angle) + (branchrvec>=0.5).*(th(BranchSites)-arp23angle);
                t(N+1:N+Nadd,1) = b(N+1:N+Nadd,1) + dx*cos(th(N+1:N+Nadd));
                t(N+1:N+Nadd,2) = b(N+1:N+Nadd,2) + dx*sin(th(N+1:N+Nadd));
                branchloc(N+1:N+Nadd,1) = t(N+1:N+Nadd,1);
                branchloc(N+1:N+Nadd,2) = t(N+1:N+Nadd,2);
                Lmin(N+1:N+Nadd) = sqrt( (t(N+1:N+Nadd,1)-b(N+1:N+Nadd,1)).^2 + (t(N+1:N+Nadd,2)-b(N+1:N+Nadd,2)).^2 );
                N = N + Nadd;
            end
        end
        
        
        if stats==1
            if plotnetworkyes==1
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % For plotting the network
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                figure(2);
                if mod(tau,tplot)==0
                    for i=1:N
                        scatter(b(i,1),b(i,2),50,'ro','filled'); hold on;
                        plot([b(i,1),t(i,1)],[b(i,2),t(i,2)],'-k','linewidth',2);
                    end
                    scatter(b(1,1),b(1,2),200,'bs','filled'); hold on;
                    xlim([-3 3]); ylim([-3 3]);
                    title(['time = ',num2str(tau*dt)]);
                    set(gca,'fontname','times','fontsize',30); box on;
                    drawnow;
                    hold off;
                    frame = getframe(gcf);
                    writeVideo(w,frame);
                end
            end %end if plotnetworkyes==1
        end %end if stats=1
        
    end %end for tau
    
    Factin_array=Factin_density(N,b,t,tau,tplot,box_size);
    Factin_array_T(:,stats) = reshape(Factin_array,array_size*array_size,1);
    
end %end stats

close(w) %for network movie
toc

%Make a gray scale image of the density data from the jth independent run
%
j=1; %run number you want to plot
densityMat = reshape(Factin_array_T(:,j),49,49); % obtain density matrix
densityGray = mat2gray(densityMat); % map density to grayscale between [0,1]
imshow(densityGray) % show the grayscale image

%change the name of the file below for each set of runs:
save Factin_array_cap_monomer Factin_array_T

