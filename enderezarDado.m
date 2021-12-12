function [dadoEnderezado] = enderezarDado(img)

[M,N,P]= size(img);
mask=zeros(M,N);

ip=find(img(:,:,1)>40 | img(:,:,2)>100 | img(:,:,3)>100);  %Buscar objetos que no sean el fondo negro
mask(ip)=255;                                      %Se crea una matriz de 255s en donde se encontraron los objetos




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

Daddogray=rgb2gray(C);
BNDado=255*imbinarize(Daddogray);
[Iobjeto3,numobjeto3] = bwlabel(BNDado);

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
    im_prueba=C(min(x)+5:max(x)-5,min(y)+5:max(y)-5,:,:,:);
    dadoEnderezado=im_prueba;
    
end

end