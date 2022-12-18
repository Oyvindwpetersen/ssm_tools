% %%
% 
% clc
% close all
% 
% P01=eye(nx)*1e-3;
% 
% close all
% 
% x0=zeros(nx,1);
% [x,y]=ssmod_forward_stoch(A,[],G,[],[],x0,[],w,v);
% 
% 
% 
% 
% [x_filt p_filt P_ss Pp_ss M_ss L_ss] = JIS_trunc_ss(A,B,G,J,y,x0,R,Q,S,P01);
% [x_filtb p_filtb P_ssb Pp_ssb M_ssb L_ssb] = JIS_trunc_ss(A,B,G,J,y,x0,R,Q,S*0,P01);
% 
% 
% [x_filt2 p_filt2 P_ss2 Pp_ss2 M_ss L_ss] = JIS_trunc_ss2(A,B,G,J,y,x0,R,Q,S,P01);
% [x_filt2b p_filt2b P_ss2b Pp_ss2b M_ssb L_ssb] = JIS_trunc_ss2(A,B,G,J,y,x0,R,Q,S*0,P01);
% 
% 
% close all
% plotTime(t,x,x_filt,x_filtb,x_filt2,x_filt2b);xlimall(gcf,[0 1000])
% plotTime(t,p,p_filt,p_filtb,p_filt2,p_filt2b);xlimall(gcf,[0 1000])
% 
% 
% delta_p=sum((p_filt-p_filt2).^2,2)
% delta_pb=sum((p_filtb-p_filt2b).^2,2)

%%

% Note 2020-05-15:

% JIS with S appears to be the same for Maes (2016) and Yan Niu (2011). 

% KF with S appears to be wrong for Maes (2016) (missing effect of S on Kalman gain) and correct for sources in books etc.

