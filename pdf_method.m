function [edf_liu, J]=pdf_method(data,mask,N)

data=mask.*data;
data(isnan(data))=0;
P=data(:);

%% Dipole kernel using the functional form in kspace

C=dipole2D(size(data));
c=conj(C);

Mn=1-mask;
W=ones(size(mask));

%% Initialize variables and set iterations

x=zeros(size(data));
r=W.*mask.*data; 

F=fftshift(fftn(ifftshift(W.*r))); 
d=Mn.*fftshift((ifftn(ifftshift(c.*F)))); 
z0_ip=sum(reshape(conj(d).*d, [numel(data) 1]));

TMP=[];
for i=1:N;
    Ad=mask.*W.*fftshift(ifftn(ifftshift(C.*fftshift(fftn(ifftshift(Mn.*d))))));
    Ad_ip=sum(reshape(conj(Ad).*Ad,[numel(data) 1]));
    alpha=(z0_ip)/(Ad_ip);
    x=x+alpha*d;
    r=r-alpha*Ad;
    s=Mn.*fftshift(ifftn(ifftshift(c.*fftshift(fftn(ifftshift(.*W.*r)))))); 
    
    z0_new_ip=sum(reshape(conj(s).*s,[numel(data) 1]));
    beta=(z0_new_ip)/(z0_ip);
    z0_ip=z0_new_ip;
    d=s+beta*d;
    
    tmp=norm(r(:));
    
    disp(['||r|| = ' num2str(tmp)])
    TMP=[TMP;tmp];
    
end

field_in= fftshift(ifftn(ifftshift(C.*fftshift(fftn(ifftshift(Mn.*x))))));
edf_liu=real(mask.*(data - field_in));
J= real(mask.* field_in);


