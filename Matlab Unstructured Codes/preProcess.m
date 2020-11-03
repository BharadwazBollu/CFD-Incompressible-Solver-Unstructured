
faces = [ int_faces ; TOP ; BOTTOM ; LEFT ; RIGHT ] ;   % Storing all faces

Vp = zeros(num_cells+sum(num_bcFaces),1) ;          % Volume of each cell
sn1     = zeros(num_cells,4) ;                      % Surface Normals
sn2     = zeros(num_cells,4) ;
normCoff= zeros(num_cells,4) ;                      % Normal Coefficent
snSign  = zeros(num_cells,4) ;                      % Surface Normal Sign
crossCoff = zeros(num_cells,4) ;                    % Cross Coefficent
faceWt  = zeros(num_cells,4) ;                      % Face Weight
pointWt = zeros(number_nodes,4) ;                   % Point or Node Weight
fout    = zeros(num_cells,4) ;
pV  = zeros(number_nodes,4) ;
lcn     = zeros(number_nodes,4) ;                   % linking cell to node
pointNode = zeros(number_nodes,8) ;
lcf = zeros(num_cells,4) ;                          % linking cell to face
lnc = zeros(num_cells,4) ;                          % linking node to cell
xc = zeros(num_cells,1) ;   yc = zeros(num_cells,1) ;   % cell centers
cellOrder = zeros(num_cells,4) ;                    % Cell Ordering for Post Process to Tecplot/ParaView
%% Computing cell Centers

for j = interior_faces+1:length(faces)
    faces(j,4)          = num_cells + k ;
    xc(num_cells + k)   = 0.5*(vertices(faces(j,1),1) + vertices(faces(j,2),1)) ;
    yc(num_cells + k)   = 0.5*(vertices(faces(j,1),2) + vertices(faces(j,2),2)) ;
    k = k + 1 ;
end

for i = 1:num_cells
    k = 1 ;
    for j = 1:length(faces)
        if      ( faces(j,3) == i  )
            lcf(i,k) = j ;              % linking cell to face
            xc(i) = xc(i) + 0.5 * 0.25 * ( vertices(faces(j,1),1) + vertices(faces(j,2),1) ) ;
            yc(i) = yc(i) + 0.5 * 0.25 * ( vertices(faces(j,1),2) + vertices(faces(j,2),2) ) ;
            k = k + 1 ;
        elseif  ( faces(j,4) == i   )
            lcf(i,k) = j ;
            xc(i) = xc(i) + 0.5 * 0.25 * ( vertices(faces(j,1),1) + vertices(faces(j,2),1) ) ;
            yc(i) = yc(i) + 0.5 * 0.25 * ( vertices(faces(j,1),2) + vertices(faces(j,2),2) ) ;
            k = k + 1 ;
        end
    end
end

%% Computing Surface Normal, Normal and Cross Coefficients
for i = 1:num_cells
    Vp(i) = 0 ;
    for k = 1:4
        
        sn1(i,k) = vertices(faces(lcf(i,k),2),2) - vertices(faces(lcf(i,k),1),2) ;
        sn2(i,k) = vertices(faces(lcf(i,k),1),1) - vertices(faces(lcf(i,k),2),1) ;
        
        dx  = 0.5 * (vertices(faces(lcf(i,k),1),1) + vertices(faces(lcf(i,k),2),1) ) - xc(i) ;
        dy  = 0.5 * (vertices(faces(lcf(i,k),1),2) + vertices(faces(lcf(i,k),2),2) ) - yc(i) ;
        
        dot = sn1(i,k) * dx + sn2(i,k) * dy ;
        
        if ( dot < 0 )
            snSign(i,k) = -1 ;
            sn1(i,k)    = -sn1(i,k);
            sn2(i,k)    = -sn2(i,k);
            fout(i,k)   = faces(lcf(i,k),3) ;
        else
            fout(i,k)   = faces(lcf(i,k),4) ;
            snSign(i,k) =  1 ;
        end
        
        Vp(i) = Vp(i) + sn1(i,k) * 0.5*(vertices(faces(lcf(i,k),1),1) + vertices(faces(lcf(i,k),2),1)) ;
        
        area_sf = sqrt( sn1(i,k)^2 + sn2(i,k)^2 ) ;
        
        dist = sqrt( (xc(faces(lcf(i,k),4)) - xc(faces(lcf(i,k),3)))^2 + (yc(faces(lcf(i,k),4)) - yc(faces(lcf(i,k),3)))^2 ) ;
        ex  = (xc(faces(lcf(i,k),4)) - xc(faces(lcf(i,k),3)))/dist ;
        ey  = (yc(faces(lcf(i,k),4)) - yc(faces(lcf(i,k),3)))/dist ;
        
        normX = (sn1(i,k)^2 + sn2(i,k)^2 ) * ex/(sn1(i,k)*ex + sn2(i,k)*ey )  ;
        normY = (sn1(i,k)^2 + sn2(i,k)^2 ) * ey/(sn1(i,k)*ex + sn2(i,k)*ey )  ;
        
        normCoff(i,k)   = sqrt( normX^2 + normY^2 )/dist ;
        
        tx = sn1(i,k) - normX ;
        ty = sn2(i,k) - normY ;
        
        crossCoff(i,k)  = sqrt( tx^2 + ty^2 )/area_sf ;
        
    end
end
%% Computing Face Weight
for i = 1:num_cells
    for k = 1:4
        if ( snSign(i,k) > 0)
            xx = 3 ; yy = 4 ;
        else
            xx = 4 ; yy = 3 ;
        end
        faceWt(i,k) = Vp(faces(lcf(i,k),yy))/( Vp(faces(lcf(i,k),3)) + Vp(faces(lcf(i,k),4)) ) ;
    end
end
%% Computing Node Weight
for j = 1:number_nodes
    zz = 1 ;
    for i = 1:num_cells
        for k = 1:4
            if ( faces(lcf(i,k),1) == j  )
                pointNode(j,zz) = i ;
                zz = zz + 1 ;
            elseif ( faces(lcf(i,k),2) == j  )
                pointNode(j,zz) = i ;
                zz = zz + 1 ;
            end
        end
    end
    k = 1 ; temp = 0 ;
    for zz = 1:length(pointNode(j,:))-1
        if ( lcn(j,k) ~= temp )
            temp = pointNode(j,zz) ;
            k = k + 1 ;
        end
        lcn(j,k) = pointNode(j,zz) ;
    end
    
    temp = 0 ;
    
    for k = 1:length(lcn(1,:))
        if ( lcn(j,k) > 0 )
            pV(j,k) = Vp(lcn(j,k)) ;
        end
        if ( pV(j,k) > 0 )
            temp = temp + 1/pV(j,k) ;
        end
    end
    
    for k = 1:length(lcn(1,:))
        if ( pV(j,k) > 0 )
            pointWt(j,k) = (1/pV(j,k))/temp ;
        end
    end
end

for i = 1:num_cells
    zz = 1 ;
    for j = 1:number_nodes
        for k = 1:length(lcn(1,:))
            if ( lcn(j,k) == i )
                lnc(i,zz) = j ;     % linking node to cells
                zz = zz + 1;
            end
        end
    end
end

for i=1:num_cells
    angle = atan2(vertices(lnc(i,:),2) - yc(i), vertices(lnc(i,:),1) - xc(i));
    [~, order] = sort(angle);
    cellOrder(i,:) = lnc(i,order) ;
end