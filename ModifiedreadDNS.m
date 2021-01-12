cd 'F:\Thesis\DNS data\K6_k1';
fclose all;
f=fopen('Save_00029600.ess','rb');

plt=1; % nx needed only when plotting
nx=256;
ny=nx/2;
nz=nx/2;

while 1
    varName = fread(f,20,'*char');
    varType = fread(f,10,'*char');
    l1=strfind(varName',char(0));
    l2=strfind(varType',char(0));
    varName=varName(1:l1(1)-1)';
    varType=varType(1:l2(1)-1)';
    disp([varType '  ' varName]);   % var type
    if strcmp(varType,'int');       type='int';
    elseif strcmp(varType,'FLOAT'); type='double';
    else type=''; error('unknown type');
    end
    n    = fread(f,1,'int'); % number of data
    disp(['  Number of data: ' num2str(n)]);
    if n>1
        str = input(['n=' num2str(n) ', continue? (y/n) '],'s');
        if ~strcmp(str,'y')
            break;
        end
    end
    data = fread(f,n,type);
    disp(['  Value: ' num2str(data(1))] );
    
    % Plot a plane
    if n>1 && plt
        data3D=reshape(data,nz,ny,nx);
        planeData=zeros(nx,nz);
        for i=1:nx
        for j=1:nz
            planeData(i,j)=data3D(j,100,i);
        end
        end

        figure()
        imagesc(planeData')
        set(gca,'Ydir','Normal');
        colormap(jet);
        title(varName)
    end
    
    % Vol plot
    if 0
        vol3d('cdata',data3D);
        axis equal off
        set(gcf, 'color', 'w');
        colormap(wb)
        alphamap(linspace(0., 0.008, 255));
    end

    save(varName,'data')

   
end