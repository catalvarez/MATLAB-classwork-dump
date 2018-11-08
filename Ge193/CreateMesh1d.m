function  [gcoord,connectivity,Nnodes,Nelements]=CreateMesh1d(x1,x2,dx,nod)

% create connectivity
    switch nod
        case 2
            Nnodes=(x2-x1)/dx+1; % total number of nodes
            gcoord=[x1:dx:x2]';
            nodes=[1:length(gcoord)]';
            Nelements=length(nodes)-1;
            connectivity=zeros(Nelements,nod);
            for I=1:Nelements
                connectivity(I,1:3)=[I I+1];
            end
        case 3
            Nelements=round((x2-x1)/dx);
            Nnodes=2*Nelements+1; % total number of nodes
            gcoord=linspace(x1,x2,Nnodes)';
            nodes=[1:length(gcoord)]';
            connectivity=zeros(Nelements,nod);
            for I=1:Nelements
                ii=2*(I-1)+1;
                connectivity(I,:)=[ii ii+1 ii+2 ];
            end
        otherwise
            disp(' not yet implemented ')
    end
    
end
    
    