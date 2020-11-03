close all
clear
clc
%% Reading Mesh file

meshfile_name = sprintf('lid_driven_unstruct.msh');     % reading mesh file name

readMeshData    % reading mesh data

preProcess      % pre processing the mesh data

%%  Initializing parameters

Re = 100 ;              % Reynolds number
dt = 0.01 ;             % time step

u_curr = zeros(num_cells+sum(num_bcFaces),1) ; u_pred = zeros(num_cells+sum(num_bcFaces),1) ; u_Next = zeros(num_cells+sum(num_bcFaces),1) ;
v_curr = zeros(num_cells+sum(num_bcFaces),1) ; v_pred = zeros(num_cells+sum(num_bcFaces),1) ; v_Next = zeros(num_cells+sum(num_bcFaces),1) ;
p_Next = zeros(num_cells+sum(num_bcFaces),1) ;

uc_pred = zeros(number_nodes,1) ; vc_pred = zeros(number_nodes,1) ;
pc_Next = zeros(number_nodes,1) ;
uc_Next = zeros(number_nodes,1) ; vc_Next = zeros(number_nodes,1) ;

flux = zeros(num_cells,4) ;
uu = zeros(4,1) ; vv = zeros(4,1) ;

pdiff = zeros(num_cells,4) ;
udiff = zeros(num_cells,4) ;    vdiff = zeros(num_cells,4) ;
PSx   = zeros(num_cells,4) ;    PSy   = zeros(num_cells,4) ;
flux_noPressureCorrection = zeros(num_cells,4) ;

iter    = 0;
error   = 1;
fprintf('  Starting time loop \n\n  ');

while error > 1e-3      % Main error loop for Checking steady state or Transient time
    
    iter = iter + 1 ;   % Number of iterations
    fprintf(' Time = %f \n', iter * dt );
    
    %% Predicted Velocity
    
    residual_error = 1 ;
    
    while ( residual_error > 1e-06 )
        
        uPred_residual = 0;
        vPred_residual = 0;
        
        for  i = 1 : num_cells
            convection_CentralCoeff = 0 ;
            uPredConvection = 0 ;       vPredConvection = 0 ;
            % Convection Term
            for k = 1:4
                if ( flux(i,k) >= 0)         % % Upwind for Convection
                    uu(k) = u_pred(i);
                    vv(k) = v_pred(i);
                    convection_CentralCoeff = convection_CentralCoeff + flux(i,k);
                else
                    uu(k) = u_pred(fout(i,k));
                    vv(k) = v_pred(fout(i,k));
                end
                % Total Convection for each cell for U and V
                uPredConvection = uPredConvection + flux(i,k) * uu(k)  ;
                vPredConvection = vPredConvection + flux(i,k) * vv(k)  ;
                
                uc_pred(lnc(i,k)) = 0 ; vc_pred(lnc(i,k)) = 0 ;
                for z = 1:length(lcn(1,:))
                    if ( lcn(lnc(i,k),z) > 0 )
                        % Volume Interpolation for Nodal Values of U and V
                        uc_pred(lnc(i,k)) = uc_pred(lnc(i,k)) + pointWt(lnc(i,k),z)*u_pred(lcn(lnc(i,k),z)) ;
                        vc_pred(lnc(i,k)) = vc_pred(lnc(i,k)) + pointWt(lnc(i,k),z)*v_pred(lcn(lnc(i,k),z)) ;
                    end
                end
            end
            
            uPredDiff_CentralCoeff = 0 ;    uTdiff = 0 ;
            vPredDiff_CentralCoeff = 0 ;    vTdiff = 0 ;
            
            % Diffusion Term
            for k = 1:4
                udiff(i,k) = normCoff(i,k) * ( u_pred(fout(i,k)) - u_pred(i) ) ...
                    + crossCoff(i,k) *snSign(i,k)* ( uc_pred(faces(lcf(i,k),2)) - uc_pred(faces(lcf(i,k),1)) );
                vdiff(i,k) = normCoff(i,k) * ( v_pred(fout(i,k)) - v_pred(i) ) ...
                    + crossCoff(i,k) *snSign(i,k)* ( vc_pred(faces(lcf(i,k),2)) - vc_pred(faces(lcf(i,k),1)) ) ;
                % Central Coefficient for U and V Predicted Velocity diffusion term
                uPredDiff_CentralCoeff = uPredDiff_CentralCoeff + normCoff(i,k) ;
                vPredDiff_CentralCoeff = vPredDiff_CentralCoeff + normCoff(i,k) ;
                % Total Diffsion for U and V Predicted Velocity
                uTdiff  =  uTdiff + udiff(i,k)  ;
                vTdiff  =  vTdiff + vdiff(i,k)  ;
            end
            
            % Total Central Coefficient for U and V Predicted Velocity
            uPredTotal_CentralCoeff	= 1 + dt/Vp(i) * convection_CentralCoeff + dt/Re/Vp(i) * uPredDiff_CentralCoeff ;
            vPredTotal_CentralCoeff = 1 + dt/Vp(i) * convection_CentralCoeff + dt/Re/Vp(i) * vPredDiff_CentralCoeff ;
            
            % Calculating Error or Residue and Correcting the value using Gauss Seidal iterative method
            uPred_error   	= u_curr(i) - u_pred(i) - dt/Vp(i) * uPredConvection + dt/Re * uTdiff/Vp(i) ;                              % residue for predicetd X velocity
            uPred_residual 	= uPred_residual + uPred_error * uPred_error;
            u_pred(i)       = uPred_error/uPredTotal_CentralCoeff + u_pred(i) ;
            
            vPred_error  	= v_curr(i) - v_pred(i) - dt/Vp(i) * vPredConvection + dt/Re * vTdiff/Vp(i) ;                              % residue for predicetd Y velocity
            vPred_residual	= vPred_residual + vPred_error * vPred_error;
            v_pred(i)    	= vPred_error/vPredTotal_CentralCoeff + v_pred(i) ;
            
        end
        residual_error = sqrt( ( uPred_residual + vPred_residual )/(num_cells)) ;                                                       % RMS for predicted velocity
    end
    % BC for Predicted Velocity ( Must implement according to order in the mesh file )
    for j = interior_faces +1 : interior_faces+num_bcFaces(1)
        u_pred(faces(j,4)) = 1 ;
        v_pred(faces(j,4)) = 0 ;
    end
    for j = interior_faces+num_bcFaces(1) +1 : interior_faces+sum(num_bcFaces)
        u_pred(faces(j,4)) = 0 ;
        v_pred(faces(j,4)) = 0 ;
    end
    
    %% Pressure Poisson
    
    residual_error = 1 ;
    
    while ( residual_error > 1e-06 )
        
        pressure_residual = 0;
        for  i = 1 : num_cells
            % Pressure Nodal Values Calculation
            for k = 1:4
                pc_Next(lnc(i,k)) = 0 ;
                for z = 1:length(lcn(1,:))
                    if ( lcn(lnc(i,k),z) > 0 )
                        % Volume Interpolation for Nodal Values of P
                        pc_Next(lnc(i,k)) = pc_Next(lnc(i,k)) + pointWt(lnc(i,k),z)*p_Next(lcn(lnc(i,k),z)) ;
                    end
                end
            end
            
            pNext_CentralCoeff = 0 ;    pTdiff = 0 ;
            % Diffusion Term
            for k = 1:4
                pdiff(i,k) = normCoff(i,k) * ( p_Next(fout(i,k)) - p_Next(i) ) ...
                    + crossCoff(i,k) *snSign(i,k)* ( pc_Next(faces(lcf(i,k),2)) - pc_Next(faces(lcf(i,k),1)) ) ;
                % Central Coefficient for Pressure Diffusion
                pNext_CentralCoeff = pNext_CentralCoeff + normCoff(i,k) ;
                % Total Diffusion of presuure
                pTdiff  =  pTdiff + pdiff(i,k)  ;
            end
            
            totalflux_noPressureCorrection = 0 ;
            for k = 1 : 4
                % Calculating Flux without pressure correction for each face
                flux_noPressureCorrection(i,k) =  ( u_pred(i)*faceWt(i,k) + u_pred(fout(i,k))*(1-faceWt(i,k)) ) * sn1(i,k) ...
                    +  ( v_pred(i)*faceWt(i,k) + v_pred(fout(i,k))*(1-faceWt(i,k)) ) * sn2(i,k)  ;
                
                totalflux_noPressureCorrection = totalflux_noPressureCorrection + flux_noPressureCorrection(i,k) ;
            end
            % Calculating Error or Residue and Correcting the value using Gauss Seidal iterative method
            p_Next_error        =  totalflux_noPressureCorrection - dt * pTdiff ;
            pressure_residual   =  pressure_residual + p_Next_error * p_Next_error;
            p_Next(i)           = -p_Next_error/( dt*pNext_CentralCoeff ) + p_Next(i);
        end
        residual_error = sqrt(pressure_residual/(num_cells));
    end
    
    % BC for Pressure
    for j = interior_faces +1 : interior_faces+num_bcFaces(1)
        p_Next(faces(j,4)) = p_Next(faces(j,3)) ;
    end
    for j = interior_faces+num_bcFaces(1) +1 : interior_faces+sum(num_bcFaces)
        p_Next(faces(j,4)) = p_Next(faces(j,3)) ;
    end
    
    %% Corrected Velocity
    
    residual_error  = 1 ;
    
    while ( residual_error > 1e-06 )
        uNext_residual = 0;
        vNext_residual = 0;
        
        for  i = 1 : num_cells
            convection_CentralCoeff = 0 ;
            uConvection = 0 ;       vConvection = 0 ;
            
            % Convection Term
            for k = 1:4 % flux Calculation
                flux(i,k) =  ( u_pred(i)*faceWt(i,k) + u_pred(fout(i,k))*(1-faceWt(i,k)) ) * sn1(i,k)...
                    +  ( v_pred(i)*faceWt(i,k) + v_pred(fout(i,k))*(1-faceWt(i,k)) ) * sn2(i,k) -dt * pdiff(i,k) ;
                
                if ( flux(i,k) >= 0)         % Upwind for Convection
                    uu(k) = u_Next(i);
                    vv(k) = v_Next(i);
                    convection_CentralCoeff = convection_CentralCoeff + flux(i,k);
                else
                    uu(k) = u_Next(fout(i,k));
                    vv(k) = v_Next(fout(i,k));
                end
                
                % Total Convection for each cell for U and V
                uConvection = uConvection + flux(i,k) * uu(k)  ;
                vConvection = vConvection + flux(i,k) * vv(k)  ;
                
                uc_Next(lnc(i,k)) = 0 ; vc_Next(lnc(i,k)) = 0 ;
                for z = 1:length(lcn(1,:))
                    if ( lcn(lnc(i,k),z) > 0 )
                        % Volume Interpolation for Nodal Values of U and V
                        uc_Next(lnc(i,k)) = uc_Next(lnc(i,k)) + pointWt(lnc(i,k),z)*u_Next(lcn(lnc(i,k),z)) ;
                        vc_Next(lnc(i,k)) = vc_Next(lnc(i,k)) + pointWt(lnc(i,k),z)*v_Next(lcn(lnc(i,k),z)) ;
                    end
                end
                
            end
            
            uNextDiff_CentralCoeff = 0 ; vNextDiff_CentralCoeff = 0 ;
            uTdiff = 0 ;    vTdiff = 0 ;
            
            % Diffusion Term in Momentum for U and V
            for k = 1:4
                udiff(i,k) = normCoff(i,k) * ( u_Next(fout(i,k)) - u_Next(i) ) ...
                    + crossCoff(i,k) *snSign(i,k)* ( uc_Next(faces(lcf(i,k),2)) - uc_Next(faces(lcf(i,k),1)) );
                
                vdiff(i,k) = normCoff(i,k) * ( v_Next(fout(i,k)) - v_Next(i) ) ...
                    + crossCoff(i,k) *snSign(i,k)* ( vc_Next(faces(lcf(i,k),2)) - vc_Next(faces(lcf(i,k),1)) );
                % Central Coefficient for U and V Corrected Velocity diffusion term
                uNextDiff_CentralCoeff = uNextDiff_CentralCoeff + normCoff(i,k) ;
                vNextDiff_CentralCoeff = vNextDiff_CentralCoeff + normCoff(i,k) ;
                
                % Total Diffsion for U and V Corrected Velocity
                uTdiff  =  uTdiff + udiff(i,k)  ;
                vTdiff  =  vTdiff + vdiff(i,k)  ;
            end
            
            % Using gauss divergence theorem integral dp/dn*dV = Sigma(pfsfn)
            % where f is face and n is direction(x,y,z)
            SigmaPfSfx = 0 ; SigmaPfSfy = 0 ;
            for k =1:4
                % Pressure for U Velocity
                PSx(i,k) = ( p_Next(i)*faceWt(i,k) + p_Next(fout(i,k))*(1-faceWt(i,k)) ) * sn1(i,k) ;
                SigmaPfSfx = SigmaPfSfx + PSx(i,k) ;
                % Pressure for V Velocity
                PSy(i,k) = ( p_Next(i)*faceWt(i,k) + p_Next(fout(i,k))*(1-faceWt(i,k)) ) * sn2(i,k) ;
                SigmaPfSfy = SigmaPfSfy + PSy(i,k) ;
            end
            
            % Total Central Coefficient for U and V Corrected Velocity
            uNextTotal_CentralCoeff = 1 + dt/Vp(i) * convection_CentralCoeff + dt/Re/Vp(i) * uNextDiff_CentralCoeff ;
            vNextTotal_CentralCoeff = 1 + dt/Vp(i) * convection_CentralCoeff + dt/Re/Vp(i) * vNextDiff_CentralCoeff ;
            
            % Calculating Error or Residue and Correcting the value using Gauss Seidal iterative method
            uNext_error  	= u_curr(i) - u_Next(i) - dt/Vp(i) * uConvection + dt/Re * uTdiff/Vp(i) - ( dt/Vp(i) * SigmaPfSfx );   % residue for predicetd X velocity
            uNext_residual	= uNext_residual + uNext_error * uNext_error;
            u_Next(i)   	= uNext_error/uNextTotal_CentralCoeff + u_Next(i) ;
            
            vNext_error  	= v_curr(i) - v_Next(i) - dt/Vp(i) * vConvection + dt/Re * vTdiff/Vp(i) - ( dt/Vp(i) * SigmaPfSfy );   % residue for predicetd Y velocity
            vNext_residual	= vNext_residual + vNext_error * vNext_error;
            v_Next(i)   	= vNext_error/vNextTotal_CentralCoeff + v_Next(i) ;
            
        end
        residual_error = sqrt( (uNext_residual+vNext_residual)/(num_cells));
        
    end
    % BC for Corrected velocities
    for j = interior_faces +1 : interior_faces+num_bcFaces(1)
        u_Next(faces(j,4)) = 1 ;
        v_Next(faces(j,4)) = 0 ;
    end
    for j = interior_faces+num_bcFaces(1) +1 : interior_faces+sum(num_bcFaces)
        u_Next(faces(j,4)) = 0 ;
        v_Next(faces(j,4)) = 0 ;
    end
    
    
    error = 0;
    
    for  i = 1 : num_cells + sum(num_bcFaces)
        error = error + (u_Next(i) - u_curr(i))^2 + (v_Next(i) - v_curr(i))^2 ;     % Residue for difflux_erence in current and next time step values
        u_curr(i)  = u_Next(i);                                                              % updating values for next time step
        v_curr(i)  = v_Next(i);
    end
    
    error = sqrt(error/(num_cells));
    error = error/dt;
    
end

%% Exporting to .plt for Tecplot or Paraview
fid=fopen('lid_driven_uns.plt','w');             % file name
fprintf(fid,'VARIABLES= X,Y,U,V,P \n');    % parameters to export
fprintf(fid,'ZONE T="U", F=FEPOINT, N=');
fprintf(fid,'%d, ',number_nodes);
fprintf(fid,'E=%d, ET=QUADRILATERAL \n',num_cells);
for i=1:number_nodes
    fprintf(fid,'%f %f %f %f %f \n',vertices(i,1),vertices(i,2),uc_Next(i),vc_Next(i),pc_Next(i) );
end
fprintf(fid,'\n\n');
for i=1:num_cells
    fprintf(fid,'%f %f %f %f \n',cellOrder(i,1),cellOrder(i,2),cellOrder(i,3),cellOrder(i,4));
end
fclose(fid);