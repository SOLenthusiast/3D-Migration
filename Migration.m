
v_mig=2050;
nparams=5;params=nan*ones(1,nparams);
n_longueur=3;
% for ik=1:3
for ik=1:3

  if ik==1  
    str_file=(['data' num2str(ik) '.mat']);
    file_data=load(str_file);
    plotimage(file_data.seis,file_data.t,file_data.x);
    saveas(gcf,(['zero-phase.png']));

 %%Simple migration
     [arymig,tmig,xmig]=kirk(file_data.seis,v_mig,file_data.t,file_data.x,params);
     plotimage(arymig,file_data.t,file_data.x);   
     saveas(gcf,(['kirch_v_mig' num2str(params(1)) '.png']));
    
 %%Vrms calculation
     nt=size(file_data.seis,1);
     t1=file_data.tv(1);t2=file_data.tv(end);
     tout=linspace(t1,t2,nt);
     vrms=vint2vrms(file_data.v,file_data.tv,tout);

 %%Kirchoff migration    
     ratio=1;
     elapsedTime=zeros(1,n_longueur);
     for j=1:n_longueur
       
         mig_aper=(abs(max(file_data.x)-min(file_data.x)))/ratio;
         ratio=ratio*2; 
         params(1)=mig_aper;
         tic
         [arymig_rms,tmig_rms_tot,xmig_rms_tot]=kirk(file_data.seis,vrms,file_data.t,file_data.x,params);
         elapsedTime(j) = toc;
         plotimage(arymig_rms,file_data.t,file_data.x);   
         saveas(gcf,(['kirchoff_' num2str(params(1)) '.png']));  
% %         Elapsed time is 10.057669 seconds.
% %         Elapsed time is 5.830910 seconds.
% %         Elapsed time is 5.830910 seconds.     
     end
     elapsedTime

    
%%FD migration
         dtau=[1e-2 3e-2 2e-3];
         for i=1:length(dtau);
             [arymig_rms_df,tmig_rms_df,xmig_rms_df]=fd15mig(file_data.seis,vrms,file_data.t,file_data.x,dtau(i));
             plotimage(arymig_rms_df,file_data.t,file_data.x);
             saveas(gcf,(['fd_' num2str(dtau(i)) '.png']));             
         end

         %% Gazdag
         params_gaz=params(1:4);
         params_gaz(3)=200;
         ratio=1;      
         for kj=1:2
             params_gaz(4)=90/ratio;
             [out,cputime]=ps_migt(file_data.seis,file_data.t,file_data.x,vrms,params_gaz);
             plotimage(out,file_data.t,file_data.x);
             saveas(gcf,(['gazdag_nopad_' num2str(params_gaz(4)) '.png']));
             ratio=ratio*6;
         end
         %% en limitant � 15 deg (gasdag.png) les reflecteurs profonds (probablement � fort ...
         %% pendage) sont tronqu�s.        
         clf
         params_gaz=params(1:4);
         params_gaz(4)=90;
         params_gaz(2)=60;
         params_gaz(3)=200;
 
         for ij=1:2
             if ij==2
                     params_gaz(1)=0.1;
             end
             [out,cputime]=ps_migt(file_data.seis,file_data.t,file_data.x,vrms,params_gaz);
             plotimage(out,file_data.t,file_data.x);
             saveas(gcf,(['gazdag_' num2str(params_gaz(2)) '_' num2str(params_gaz(1)) '.png']));
         end
         %% Pas de diff�rence observable entre la 0 et 0.1 tpad, mais changement de phase de l'ondelette
         %% d'une forme antisymm�trique � une forme sym�trique
        
           %% Stolt migation
           params_fk=[ params params params(1:3) ];
           [seismig,tmig,xmig]=fkmig(file_data.seis,file_data.t,file_data.x,v_mig,params_fk);
           plotimage(seismig,file_data.t,file_data.x);
           saveas(gcf,'Stolt.png');
 
             %% En prof ni stolt ni gazdag ne retrouvent les reflecteurs, gazdag marche mieux en surface,
             %% stolt semble de moindre qualit� en gle.

  else if ik==2  
     str_file=(['data' num2str(ik) '.mat']);
     file_data=load(str_file);    
     
     %%Vrms calculation
     nt=size(file_data.seis,1);
     t1=file_data.tv(1);t2=file_data.tv(end);
     tout=linspace(t1,t2,nt);
     vrms=vint2vrms(file_data.v,file_data.tv,tout);
     
%         %% PS
         params_gaz=params(1:4);
         params_gaz(4)=90;
         params_gaz(2)=60;
         params_gaz(3)=200;
 
         [out,cputime]=ps_migt(file_data.seis,file_data.t,file_data.x,vrms,params_gaz);
         plotimage(out,file_data.t,file_data.x);
         saveas(gcf,(['ps_dat2_' num2str(params_gaz(2)) '_' num2str(params_gaz(1)) '.png']));
 
         %% PSPI    
         val_init=tout(2);
         tout(2)=2e-3;
         [zosmig,exzos]=pspi_stack_tmig(file_data.seis,file_data.t,file_data.x,vrms,0,tout,[0 inf],[]);
         plotimage(zosmig,file_data.t,file_data.x);
         saveas(gcf,(['PSPI_dat2.png']));
         tout(2)=val_init;
         
         %% pspi retrouve les reflecteurs � la mi-profondeur sans laisser de train��e blanche,...
         %% gazdag ne retrouve que les 2 premiers points les plus superficiels + aug du ray. de courb. des...
         %% hyperboles avec la profondeur
     
      else
           str_file=(['data' num2str(ik) '.mat']);
           file_data=load(str_file);
           
           %% Vrms calculation
           nt=size(file_data.seis,1);
           t1=file_data.tv(1);t2=file_data.tv(end);
           tout=linspace(t1,t2,nt);
           vrms=vint2vrms(file_data.v,file_data.tv,tout);
           
           %% PSPI
           val_init=tout(2);
          tout(2)=file_data.t(2);
           [zosmig,exzos]=pspi_stack_tmig(file_data.seis,file_data.t,file_data.x,vrms,0,tout,[0 inf],[]);
           plotimage(zosmig,file_data.t,file_data.x);
           saveas(gcf,(['PSPI_dat3.png']));
           tout(2)=val_init;
            
           %% Kirchoff migration
           mig_aper=abs(max(file_data.x)-min(file_data.x));
           params(1)=mig_aper;
           [arymig_rms,tmig_rms_tot,xmig_rms_tot]=kirk(file_data.seis,vrms,file_data.t,file_data.x,params);
           plotimage(arymig_rms,file_data.t,file_data.x);
           saveas(gcf,(['Kirch_dat3' num2str(params(1)) '.png']));
     
           %% FD migration
           %dtau=[1e-2 2e-3 3e-2];
           dtau=2e-3;          
           [arymig_rms_df,tmig_rms_df,xmig_rms_df]=fd15mig(file_data.seis,vrms,file_data.t,file_data.x,dtau);
           plotimage(arymig_rms_df,file_data.t,file_data.x);
           saveas(gcf,(['FD_dat3' num2str(dtau) '.png']));
           %% FD ne marche pas, PSPI retrouve toutes les interfaces, Kirch resout moins bien les interf ...
           %% superficielles du � la longueur de fen�tre un peu grande qui fait que le bruit en profond ...
           %% affecte les interface pr�coces en temps, avec cette train�e blanche verticale,...
           %% Pour le reste de la section le bruit est remplac� par des smiles
      end
  end
end


% figure(2);plot(file_data.v,file_data.tv,'ro');
% set(gca,'YDir','reverse')
