function [...
    Hyp,Hyw,Hyp_e,...
    delta_Txp_jis,delta_Txw_jis,delta_Txv_jis,...
    delta_Tpp_jis,delta_Tpw_jis,delta_Tpv_jis...
    ]=jis_tf_err( ...
    omega_axis_plot,dt,...
    Ad,Bd,Gd,Jd,Q,R,S,...
    Ad_e,Bd_e,Gd_e,Jd_e,Q_e,R_e,S_e ...
    )
 
%% Transfer function from p,w,v to errors in state and input estimates
%
% Inputs:
% omega_axis: frequency axis
% dt: time step
% Ad: state matrix
% Bd: input matrix
% Gd: output matrix
% Jd: direct transmission matrix
% Q: state noise covariance
% R: output noise covariance
% S: mixed noise covariance
% Ad_e: state matrix (erroneous)
% Bd_e: input matrix (erroneous)
% Gd_e: output matrix (erroneous)
% Jd_e: direct transmission matrix (erroneous)
% Q_e: state noise covariance (erroneous)
% R_e: output noise covariance (erroneous)
% S_e: mixed noise covariance (erroneous)
%
% Outputs:
% Hyp: transfer from force to output
% Hyw: transfer from process noise to output
% Hyp_e: transfer from force to output (erroneous)
% delta_Txp_jis: transfer from force to state estimate
% delta_Txw_jis: transfer from process noise to state estimate
% delta_Txv_jis: transfer from output noise to state estimate
% delta_Tpp_jis: transfer from force to input estimate
% delta_Tpw_jis: transfer from process input to state estimate
% delta_Tpv_jis: transfer from output noise to input estimate
%

%% Parse inputs

%% Forward transfer functions

% Hxp=ssmod_tf(Ad,Bd,Gd,Jd,omega_axis_plot,dt,'type','is'); % not used
Hyp=ssmod_tf(Ad,Bd,Gd,Jd,omega_axis_plot,dt,'type','io');

Hyw=ssmod_tf(Ad,eye(size(Ad)),Gd,zeros(size(Gd,1),size(Ad,2)),omega_axis_plot,dt,'type','io');

Hyp_e=ssmod_tf(Ad_e,Bd_e,Gd_e,Jd_e,omega_axis_plot,dt,'type','io');

%% JIS

[T_py_jis T_x0y_jis T_x1y_jis]=JIS_tf(Ad,Bd,Gd,Jd,Q,R,S,dt,omega_axis_plot);
[T_py_jis_e T_x0y_jis_e T_x1y_jis_e]=JIS_tf(Ad_e,Bd_e,Gd_e,Jd_e,Q_e,R_e,S_e,dt,omega_axis_plot);

% Exception: erroneous model contains more states (e.g. to compensate for bias)
% Reduce to the same states as the true model
if size(T_x0y_jis_e,1)>size(T_x0y_jis,1)
    nx=size(T_x0y_jis);
    T_x0y_jis_e=T_x0y_jis_e(1:nx,:,:);
    T_x1y_jis_e=T_x1y_jis_e(1:nx,:,:);
end

% Transfer function from p, w, and v to state and input estimate
delta_Txp_jis=mtimes3(T_x0y_jis-T_x0y_jis_e,Hyp);
delta_Txw_jis=mtimes3(T_x0y_jis-T_x0y_jis_e,Hyw);
delta_Txv_jis=T_x0y_jis-T_x0y_jis_e;

delta_Tpp_jis=mtimes3(T_py_jis-T_py_jis_e,Hyp);
delta_Tpw_jis=mtimes3(T_py_jis-T_py_jis_e,Hyw);
delta_Tpv_jis=T_py_jis-T_py_jis_e;

% Check
% plottf(omega_axis_plot,abs(mtimes3(T_py_jis,Hyp)),'xlim',[0 15]);
% plottf(omega_axis_plot,abs(mtimes3(T_py_jis_e,Hyp_e)),'xlim',[0 15]);


