
function [imgEnderezada, numValCarta, imgAreaBusqueda,posIMGNUM, posIMGSim,posIMREY,areainfo] = procesar(img)
%% debe ingresar carta recortada ( solo 1)
posIMREY=0;
[M,N,P]= size(img);
mask=zeros(M,N);

ip=find(img(:,:,1)>40 | img(:,:,2)>100 | img(:,:,3)>100);  %Buscar objetos que no sean el fondo negro
mask(ip)=255;                                      %Se crea una matriz de 255s en donde se encontraron los objetos

areainfo(1)=0;
areainfo(2)=0;
areainfo(3)=0;
areainfo(4)=0;


% izquierda
[y_xmin]=find(mask(:,1)==255);  
y_xmin = round(sum(y_xmin)/length(y_xmin));

% derecha
[y_xmax]=find(mask(:,end)==255);  
y_xmax = round(sum(y_xmax)/length(y_xmax));

% arriba
[x_ymin]=find(mask(1,:)==255);  
x_ymin = round(sum(x_ymin)/length(x_ymin));

% abajo
[x_ymax]=find(mask(end,:)==255);  
x_ymax = round(sum(x_ymax)/length(x_ymax));




% determinar el lado mas largo para enderezar
lado_1 = [abs(1 - x_ymin)  abs(y_xmin - 1 )];  %izquierda - arriba
lado_2 = [abs(1 - x_ymax)  abs(y_xmin - M)];  %izquierda - abajo

if (abs( lado_1(1)^2 + lado_1(2)^2)  )^(1/2) > (abs( lado_2(1)^2 + lado_2(2)^2 ) )^(1/2)
    punto_1 = [1 x_ymin];  %arriba
    punto_2 = [y_xmin 1];  %izquierda
    cw =  1;
else
    punto_1 = [y_xmin   1]; %izquierda
    punto_2 = [M   x_ymax]; %abajo
    cw = -1;
end




% aproximar a recta, ya que se tomaron puntos de esquinas es un primer lugar

% si la horizontal es mas grande que la vertical
if abs(punto_2(1) - punto_1(1)) <= abs(punto_2(2) - punto_1(2))
    
    if cw == -1
        r_x1 = punto_1(2)+10;
        r_y1 = find(mask(:,r_x1)==255, 1, 'last' );
        
        r_x2 = punto_2(2)-10;
        r_y2 = find(mask(:,r_x2)==255, 1, 'last' );
    else
        r_x1 = punto_1(2)-10;
        r_y1 = find(mask(:,r_x1)==255, 1 );
        
        r_x2 = punto_2(2)+10;
        r_y2 = find(mask(:,r_x2)==255, 1 );
        
        
    end
else
    
    if cw == 1
        r_y1 = punto_1(1)+10;
        r_x1 = find(mask(r_y1,:)==255, 1 );
        
        r_y2 = punto_2(1)-10;
        r_x2 = find(mask(r_y2,:)==255, 1 );
        
    else
        r_y1 = punto_1(1)+10;
        r_x1 = find(mask(r_y1,:)==255, 1 );
        
        r_y2 = punto_2(1)-10;
        r_x2 = find(mask(r_y2,:)==255, 1 );
    
    end
end



angle = atan(  (-r_x1 + r_x2) / (-r_y1+r_y2) );



 

[rowsi,colsi,z]= size(img); 

rads = -angle;



rowsf=ceil(rowsi*abs(cos(rads))+colsi*abs(sin(rads)));                      
colsf=ceil(rowsi*abs(sin(rads))+colsi*abs(cos(rads)));                     


C=uint8(zeros([rowsf colsf 3 ]));


xo=ceil(rowsi/2);                                                            
yo=ceil(colsi/2);

midx=ceil((size(C,1))/2);
midy=ceil((size(C,2))/2);


for i=1:size(C,1)
    for j=1:size(C,2)                                                       

         x= (i-midx)*cos(rads)+(j-midy)*sin(rads);                                       
         y= -(i-midx)*sin(rads)+(j-midy)*cos(rads);                             
         x=round(x)+xo;
         y=round(y)+yo;

         if (x>=1 && y>=1 && x<=size(img,1) &&  y<=size(img,2) ) 
              C(i,j,:)=img(x,y,:);  
         end

    end
end



[r, c, n]= size(C);
if (r<c)
%     for i=1:1:M
%         for j=1:1:N
%             newC(i,j)=C(j,N-i+1);
%         end
%     end
    newC=imrotate(C,90);
    
    C=newC;
end



imgEnderezada=C;


C_GRAY=rgb2gray(C);
pout = imadjust(C_GRAY);
C_BN=255*imbinarize(pout);
[Iobjeto2,numobjeto2] = bwlabel(C_BN);
posC_BN=cell(1,numobjeto2);
posC_CL=cell(1,numobjeto2);
[M,N]= size(C_BN);
if numobjeto2>0 
    for num=1:numobjeto2
        
        [xTemp,yTemp]=find(Iobjeto2==num);
        posUbi{num}=[xTemp,yTemp];
        
        
       If=zeros(M,N); 
       for x=1:M
           for y=1:N
               if Iobjeto2(x,y)==num
                  If(x,y)=255;    
               end   
           end    
       end

    end
end


tem = zeros(1,numobjeto2);
% Detección de regiones indeseadas por tamaño
if numobjeto2 > 0
    for s = 1:numobjeto2
        tem(s) = length(Iobjeto2(Iobjeto2==s));  % Almacena el tamaño de cada region
        if tem(s) < 1500 % Menores a este valor son regiones indeseadas
            % Quitamos las regiones indeseadas
            for i = 1:M
                for j = 1:N
                    if Iobjeto2(i,j) == s
                        C_BN(i,j) = 0; % Ponemos a 0 las regiones indeseadas
                    end
                end
            end
        end
    end
end

[Iobjeto3,numobjeto3] = bwlabel(C_BN);

if numobjeto3>0 
    for num=1:numobjeto3
        
        [xTemp,yTemp]=find(Iobjeto3==num);
        posUbi{num}=[xTemp,yTemp];
        
        
       If=zeros(M,N); 
       for x=1:M
           for y=1:N
               if Iobjeto3(x,y)==num
                  If(x,y)=255;    
               end   
           end    
       end
        posC_BN{num}=If;


    end
end
      

for i=1:numobjeto3
    x=posUbi{i}(:,1);
    y=posUbi{i}(:,2);
    im_prueba=C(min(x):max(x),min(y):max(y),:,:,:);
    posC_CL{i}=im_prueba;
    
end





%imagen recortada para contar numero carta.
posIMBN=posC_CL{1}(20:280,35:170,:,:,:);

%imagen recortada del numero.
posIMGNUM=posC_CL{1}(5:50,5:40,:,:,:);

%imagen recortada del simbolo.
posIMGSim=posC_CL{1}(50:80,5:40,:,:,:);

%%verificacion de color
%Color Rojo
detectColorRojo=find(posIMGSim(:,:,1)>=100 & posIMGSim(:,:,2)<=40 &posIMGSim(:,:,3)<=50);
numValCarta(1)=length(detectColorRojo);
%Color Amarillo
detectColorAmarillo=find(posIMBN(:,:,1)>=140 & posIMBN(:,:,1)<=190 & posIMBN(:,:,2)>=120 & posIMBN(:,:,2)<=180 & posIMBN(:,:,3)<=90);
esRey=length(detectColorAmarillo);

%Color Azul
detectColorAzul=find(posIMBN(:,:,1)>=50 & posIMBN(:,:,1)<=120 & posIMBN(:,:,2)>=70 & posIMBN(:,:,2)<=110 & posIMBN(:,:,3)>=70 & posIMBN(:,:,3)<=180);
estaVolteada=length(detectColorAzul);



%deteccion de simbolo 
%vector de cimagen con areaDeBusqueda
imgAreaBusqueda=posIMBN; propsGray=rgb2gray(posIMGSim); pout_imadjust = imadjust(propsGray); BW = 255*imbinarize(pout_imadjust); test50=regionprops(BW,'All'); test1=[test50.Area]; test2=[test50.Perimeter];
circularities =  test2.^2 ./ (4*pi*test1);


numValCarta(3)=max(circularities);
areainfo(1)  = esRey;

areainfo(2) = estaVolteada;

if (esRey<100)
    numValCarta(4)=0;
   if (estaVolteada>8000)
       
       numValCarta(2)=0;
       
   else
       
    
    
    
    %deteccion de numero 1 al 10
    textybn=rgb2gray(posIMBN);  textybn1=255-255*imbinarize(textybn);
    [M,N]= size(textybn1);
    [Iobjeto,numobjeto] = bwlabel(textybn1);

    tem = zeros(1,numobjeto);
    % Detección de regiones indeseadas por tamaño
    if numobjeto > 0
        for s = 1:numobjeto
            tem(s) = length(Iobjeto(Iobjeto==s));  % Almacena el tamaño de cada region
            if tem(s) < 400 % Menores a este valor son regiones indeseadas
                % Quitamos las regiones indeseadas
                for i = 1:M
                    for j = 1:N
                        if Iobjeto(i,j) == s
                            textybn1(i,j) = 0; % Ponemos a 0 las regiones indeseadas
                        end
                    end
                end
            end
        end

    end
    
    [Iobjetonew,numobjetonew] = bwlabel(textybn1);


    numValCarta(2)=numobjetonew;
    
   end

else
    numValCarta(4)=1;
    
    %cortamos areaSombrero
     posIMJSom=posC_CL{1}(27:46,75:146,:,:,:);
     
     %Detectamos area roja sombrero
     detectColorRojoSombrero=find(posIMJSom(:,:,1)>=110 & posIMJSom(:,:,1)<=180 & posIMJSom(:,:,2)<=50 & posIMJSom(:,:,3)<=80);
     %Detectamos area amarillo sombrero
     detectColorAmarilloSombrero=find(posIMJSom(:,:,1)>=130 & posIMJSom(:,:,1)<=190 & posIMJSom(:,:,2)>=130 & posIMJSom(:,:,2)<=180 & posIMJSom(:,:,3)>=20 & posIMJSom(:,:,3)<=90);
         
     posIMREY=posIMJSom;
     
     tamSomRojo=length(detectColorRojoSombrero);
     tamSomAmarillo=length(detectColorAmarilloSombrero);
     areainfo(3) = tamSomRojo;
     areainfo(4) = tamSomAmarillo;
     
     if (tamSomRojo > 200 && tamSomAmarillo < 130)
         numValCarta(2)=11;
     elseif (tamSomRojo > 10 && tamSomRojo < 85  && tamSomAmarillo < 170)
         numValCarta(2)=12;    
     elseif (tamSomRojo < 170  && tamSomAmarillo > 178)    
         numValCarta(2)=13;  
     else
         numValCarta(2)=randi([12,13]);
     end
     %numValCarta(1)=length(detectColorRojo);
    
end



end