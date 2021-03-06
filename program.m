
function varargout = program(varargin)

% PROGRAM MATLAB code for program.fig
%      PROGRAM, by itself, creates a new PROGRAM or raises the existing
%      singleton*.
%
%      H = PROGRAM returns the handle to a new PROGRAM or the handle to
%      the existing singleton*.
%
%      PROGRAM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PROGRAM.M with the given input arguments.
%
%      PROGRAM('Property','Value',...) creates a new PROGRAM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before program_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to program_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help program

% Last Modified by GUIDE v2.5 12-Nov-2021 13:59:56

% Begin initialization code - DO NOT EDIT

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @program_OpeningFcn, ...
                   'gui_OutputFcn',  @program_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before program is made visible.
function program_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to program (see VARARGIN)

% Choose default command line output for program
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes program wait for user response (see UIRESUME)
% uiwait(handles.figure1);

global WBC_ID WBC_FORMAT indiceDiamante indiceEspada defaultPanel defaultFigure 
defaultPanel = get(handles.uipanel1,'OuterPosition');
defaultFigure = get(handles.figure1,'OuterPosition');

info=imaqhwinfo('winvideo');
WBC={'Seleccionar Video'};
WBC_ID = "0";
WBC_FORMAT=["0"; 'RGB24_1280x960';'YUY2_640x480'];
indiceDiamante = 1.456;
indiceEspada = 1.52;

upclogo=imread('UPC.jpg','jpg');
axes(handles.axes14)
imshow(upclogo)

for i=1:length(info.DeviceIDs)
    temp = imaqhwinfo('winvideo',i);
    WBC = [WBC; temp.DeviceName];
end
for i=1:length(info.DeviceIDs)
    temp = imaqhwinfo('winvideo',i);
    WBC_ID = [WBC_ID; temp.DeviceID];
end
set(handles.popupmenu1,'String',WBC)

% --- Outputs from this function are returned to the command line.
function varargout = program_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global WBC_ID WBC_FORMAT vid
dispoEntrada = get(handles.popupmenu1,"value");
handles.output = hObject;
axes(handles.axes1);
sizevid={[0 0];[1280 960]; [640 480]};
if dispoEntrada ~= 1
    closepreview;
    vid=videoinput('winvideo',WBC_ID(dispoEntrada),WBC_FORMAT(dispoEntrada));
    
    hImage=image(zeros(sizevid{dispoEntrada}(2),sizevid{dispoEntrada}(1),3));
    preview(vid,hImage);
    guidata(hObject, handles);


end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global vid I
%cd 'data'
    
posicion = get(handles.uipanel2,'Position');
set(handles.uipanel2,'Position',[posicion(1) posicion(2) posicion(3) 33.7 ])

I = getsnapshot(vid);
closepreview();
axes(handles.axes1)
imshow(I)



% --- Executes on button press in UploadFile.
function UploadFile_Callback(hObject, eventdata, handles)
% hObject    handle to UploadFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global I
I=0;
[filename, pathname] = uigetfile( {'*.jpg;*.png', 'Archivo de imagen (*.jpg,*.png'}, 'Seleccione una imagen');
direccion=fullfile(pathname,filename);


if isequal(filename,0)
   disp('Cancelado por el usuario');
else
   disp(['Se ha abierto el archivo con exito: ', fullfile(pathname,filename)]);

    posicion = get(handles.uipanel2,'Position');
    set(handles.uipanel2,'Position',[posicion(1) posicion(2) posicion(3) 33.7 ])

   I=imread(direccion,'jpg');
   closepreview();
   axes(handles.axes1);
   imshow(I)
   
end


% --- Executes on button press in BotonProcesar.
function BotonProcesar_Callback(hObject, eventdata, handles)
% hObject    handle to BotonProcesar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global I imgEnderezada numValCarta imgAreaBusqueda posIMGNUM posIMGSim index numobjeto indiceDiamante indiceEspada posIMGcolor posIMREY areainfo

set(handles.uipanel2,'Visible','off')

Igray=rgb2gray(I);

[M,N]= size(Igray); % Determina el tama?o de la imagen.
Ibin = zeros(M,N); % Imagen modificada
for x = 1:M
    for y = 1:N
        if Igray(x,y) <= 70 % Intensidades. Estaba  170.
            Ibin(x,y) = 0;
        else
            Ibin(x,y) = 255;
        end
    end
end




% Primera matriz de etiquetado (Para detectar regiones peque?as indeseadas)
[ob,nob] = bwlabel(Ibin);

tem = zeros(1,nob);
% Detecci?n de regiones indeseadas por tama?o
if nob > 0
    for s = 1:nob
        tem(s) = length(ob(ob==s));  % Almacena el tama?o de cada region
        if tem(s) < 1500 % Menores a este valor son regiones indeseadas
            % Quitamos las regiones indeseadas
            for i = 1:M
                for j = 1:N
                    if ob(i,j) == s
                        Ibin(i,j) = 0; % Ponemos a 0 las regiones indeseadas
                    end
                end
            end
        end
    end
end


% Segunda matriz de etiquetado (Para detectar todas las fichas)
[Iobjeto,numobjeto] = bwlabel(Ibin);

pos=cell(1,numobjeto); %imagen individual en tama?o de captura
posIMGcolor=cell(1,numobjeto); %imagen de carta cortada ( solo la carta )
posUbi=cell(1,numobjeto); %posicion de carta detectada
imgEnderezada=cell(1,numobjeto); %imagen enderezada
numValCarta=cell(1,numobjeto);
imgAreaBusqueda=cell(1,numobjeto); %imagen de busqueda de numero
posIMGNUM=cell(1,numobjeto); %Imagen recortada numero
posIMGSim=cell(1,numobjeto); %Imagen recortada simbolo
posIMREY=cell(1,numobjeto); %Imagen recortada cabeza rey
areainfo=cell(1,numobjeto); %area info

if numobjeto>0 
    for num=1:numobjeto
        
        [xTemp,yTemp]=find(Iobjeto==num);
        % 
        posUbi{num}=[xTemp,yTemp];
        
        
       If=zeros(M,N); 
       for x=1:M
           for y=1:N
               if Iobjeto(x,y)==num
                  If(x,y)=255;    
               end   
           end    
       end
        pos{num}=If;

    end


for i=1:numobjeto
    x=posUbi{i}(:,1);
    y=posUbi{i}(:,2);
    im_prueba=I(min(x):max(x),min(y):max(y),:,:,:);
    posIMGcolor{i}=im_prueba;
%       figure()
%       imshow(im_prueba)
%       impixelinfo;
end


for i=1:numobjeto
    [imgEnderezada{i}, numValCarta{i}, imgAreaBusqueda{i}, posIMGNUM{i}, posIMGSim{i}, posIMREY{i},areainfo{i} ]=procesar(posIMGcolor{i});
end




index = 1;
axes(handles.axes2);
imshow(imgEnderezada{index})
   
axes(handles.axes3);
imshow(imgAreaBusqueda{index})

axes(handles.axes4);
imshow(posIMGNUM{index})
   
axes(handles.axes5);
imshow(posIMGSim{index})   

textoNum=numValCarta{index}(2);

if (numValCarta{index}(1)>10)
     textoSim = "Rojo";
else
     textoSim ="Negro";
end       

textoProp=numValCarta{index}(3);

if (textoSim=="Rojo")
    if (textoProp<indiceDiamante)
        
        textFormSim="Diamente";
    else
        textFormSim="Corazon";
    end
else
    if (textoProp<indiceEspada)
        textFormSim="Espada";
    else
        textFormSim="Trebol";
    end
end

    if (numValCarta{index}(4) == 1)
        if textoNum == 12
            esRey="Es reyna";
        else
            esRey="Es rey";
        end
    else
        esRey="No es rey";
    end

if (textoNum <= 10 && textoNum > 0)
            numeroCarta=num2str(textoNum);
    elseif (textoNum == 11)
           numeroCarta="J";
    elseif (textoNum == 12)
           numeroCarta="Q";
    elseif (textoNum == 13)
           numeroCarta="K";
else
    numeroCarta="Volteado";
    textoSim = "Azul";
    textFormSim = "No hay";
end




set(handles.text2,'String',strcat('El valor de la carta es : ',numeroCarta));   
set(handles.text3,'String',strcat('Color : ',textoSim));   
set(handles.text4,'String',strcat('Simbolo : ',textFormSim));   
set(handles.text5,'String',esRey); 
    set(handles.pushbutton5,'Enable','off');
    set(handles.pushbutton6,'Enable','on');

end


% --- Executes on button press in pushbutton5.
%CARTA ANTERIOR
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global index imgEnderezada imgAreaBusqueda posIMGNUM posIMGSim numValCarta indiceDiamante indiceEspada

index = index-1;

pushbutton12_Callback(hObject, eventdata, handles)
if (index>0)

        axes(handles.axes2);
        imshow(imgEnderezada{index})

        axes(handles.axes3);
        imshow(imgAreaBusqueda{index})

        axes(handles.axes4);
        imshow(posIMGNUM{index})

        axes(handles.axes5);
        imshow(posIMGSim{index})   

        textoNum=numValCarta{index}(2);
       

        if (numValCarta{index}(1)>10)
            textoSim = "Rojo";
        else
            textoSim ="Negro";
        end       
        
        textoProp=numValCarta{index}(3);
        if (textoSim=="Rojo")
            if (textoProp<indiceDiamante)

                textFormSim="Diamente";
            else
                textFormSim="Corazon";
            end
        else
            if (textoProp<indiceEspada)
                textFormSim="Espada";
            else
                textFormSim="Trebol";
            end
        end
        if (numValCarta{index}(4) == 1)
            if textoNum == 12
                esRey="Es reyna";
            else
                esRey="Es rey";
            end
        else
            esRey="No es rey";
        end
        
        if (textoNum <= 10 && textoNum > 0)
            numeroCarta=num2str(textoNum);
            elseif (textoNum == 11)
                   numeroCarta="J";
            elseif (textoNum == 12)
                   numeroCarta="Q";
            elseif (textoNum == 13)
                   numeroCarta="K";
        else
            numeroCarta="Volteado";
            textoSim = "Azul";
            textFormSim = "No hay";
        end
        
        set(handles.text2,'String',strcat('El valor de la carta es : ',numeroCarta));   
        set(handles.text3,'String',strcat('Color : ',textoSim));  
        set(handles.text4,'String',strcat('Simbolo : ',textFormSim)); 
        set(handles.text5,'String',esRey);
        if (get(handles.pushbutton6,'Enable')== "off")
            set(handles.pushbutton6,'Enable','on');
        end
        
        if ( index-1 <=0)
            set(handles.pushbutton5,'Enable','off');
        end

else

        set(handles.pushbutton5,'Enable','off');
        index = index+1;

end






% --- Executes on button press in pushbutton6.
%CARTA SIGUIENTE
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global index imgEnderezada imgAreaBusqueda posIMGNUM posIMGSim numValCarta numobjeto indiceDiamante indiceEspada
index = index+1;

pushbutton12_Callback(hObject, eventdata, handles)
if (index<=numobjeto)

        axes(handles.axes2);
        imshow(imgEnderezada{index})

        axes(handles.axes3);
        imshow(imgAreaBusqueda{index})

        axes(handles.axes4);
        imshow(posIMGNUM{index})

        axes(handles.axes5);
        imshow(posIMGSim{index})   

        textoNum=numValCarta{index}(2);
        
        if (numValCarta{index}(1)>10)
            textoSim = "Rojo";
        else
            textoSim ="Negro";
        end   
        
        textoProp=numValCarta{index}(3);
        
        if (textoSim=="Rojo")
            if (textoProp<indiceDiamante)

                textFormSim="Diamente";
            else
                textFormSim="Corazon";
            end
        else
            if (textoProp<indiceEspada)
                textFormSim="Espada";
            else
                textFormSim="Trebol";
            end
        end
        
        if (numValCarta{index}(4) == 1)
            if textoNum == 12
                esRey="Es reyna";
            else
                esRey="Es rey";
            end
        else
            esRey="No es rey";
        end
        
        if (textoNum <= 10 && textoNum > 0)
            numeroCarta=num2str(textoNum);
            elseif (textoNum == 11)
                   numeroCarta="J";
            elseif (textoNum == 12)
                   numeroCarta="Q";
            elseif (textoNum == 13)
                   numeroCarta="K";
        else
            numeroCarta="Volteado";
            textoSim = "Azul";
            textFormSim = "No hay";
        end
        
        set(handles.text2,'String',strcat('El valor de la carta es : ',numeroCarta));   
        set(handles.text3,'String',strcat('Color : ',textoSim));  
        set(handles.text4,'String',strcat('Simbolo : ',textFormSim));
        set(handles.text5,'String',esRey);
        temp=get(handles.pushbutton5,'Enable');
        if (temp == "off")
            set(handles.pushbutton5,'Enable','on');
        end
        if ( index+1 >numobjeto)
            set(handles.pushbutton6,'Enable','off');
        end
else

        set(handles.pushbutton6,'Enable','off');
        index = index-1;

end


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global defaultFigure
set(handles.figure1,'OuterPosition',[defaultFigure(1) defaultFigure(2) 240 defaultFigure(4)]);
set(handles.uipanel4,'visible','off');
%set(handles.uipanel1,'OuterPosition',[2.8 0.923 205 59]);

% --- Executes on button press in radiobutton3.
function radiobutton3_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton3
if ( get(handles.radiobutton3,'Value')==1)
       set(handles.UploadFile,'Visible','off');
       set(handles.popupmenu1,'Visible','on');
       set(handles.pushbutton1,'Visible','on');
       set(handles.pushbutton2,'Visible','on');
       set(handles.BotonProcesar,'Position',[10.8 0.95 25 2.5]);
   end


% --- Executes on button press in radiobutton4.
function radiobutton4_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton4
   if ( get(handles.radiobutton4,'Value')==1)
       set(handles.UploadFile,'Visible','on');
       set(handles.popupmenu1,'Visible','off');
       set(handles.pushbutton1,'Visible','off');
       set(handles.pushbutton2,'Visible','off');
       set(handles.BotonProcesar,'Position',[10.8 3.46 25 2.5]);
   end


% --- Executes during object creation, after setting all properties.
function uibuttongroup2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uibuttongroup2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global I turnoJuego rondaJuego 
set(handles.uipanel2,'Visible','off')
set(handles.pushbutton5,'Enable','off')
set(handles.pushbutton6,'Enable','off')
Igray=rgb2gray(I);

%pout_imadjust = imadjust(Igray);

BW = 255*imbinarize(Igray,0.3);
[M,N]= size(BW);
[ob,nob] = bwlabel(BW);

 
axes(handles.axes2);
imshow(BW)
 if nob>0 && 1<=nob
    for num=1:nob
        if(length(find(ob==num))>800)
      [xTemp,yTemp]=find(ob==num);
        % 
      posUbiDado{num}=[xTemp,yTemp];
        
        
     If=zeros(M,N); 
     for x=1:M
         for y=1:N
             if ob(x,y)==num
                  If(x,y)=255;    
             end   
         end    
     end
        end
    end
[ob,nob] = bwlabel(If);
 
 for i=1:nob
    x=posUbiDado{i}(:,1);
    y=posUbiDado{i}(:,2);
    im_prueba=I(min(x):max(x),min(y):max(y),:,:,:);
    posIMGcolorDADO=im_prueba;

end
axes(handles.axes3);
imshow(posIMGcolorDADO)


[dadoEnderezado] = enderezarDado(posIMGcolorDADO);

axes(handles.axes4);
imshow(dadoEnderezado)




Daddogray=rgb2gray(dadoEnderezado);
BNDado=255- 255*imbinarize(Daddogray);

[obdado,nobdado] = bwlabel(BNDado);
[M,N]= size(BNDado);


 axes(handles.axes5);
imshow(BNDado)

set(handles.text12,'String',strcat('Tu dado : ',num2str(nobdado)));

%dado de la pc
pause(1)
dadoPC=randi([1,6]);
set(handles.text11,'String',strcat('Dado de PC : ',num2str(dadoPC)));

if( nobdado>dadoPC || nobdado==dadoPC)
    set(handles.text13,'String','Usuario comienza');
    set(handles.text10,'String','TURNO USUARIO');
    turnoJuego=1;   
    rondaJuego=1;
else
    set(handles.text13,'String','PC comienza');
    set(handles.text10,'String','TURNO PC');
    turnoJuego=0;
    rondaJuego=1;
end
set(handles.pushbutton8,'Enable', 'off');
set(handles.pushbutton9,'Enable','on');

end

% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global I turnoJuego rondaJuego indiceDiamante indiceEspada cartaUsuario cartaPC

Igray=rgb2gray(I);

[M,N]= size(Igray); % Determina el tama?o de la imagen.
Ibin = zeros(M,N); % Imagen modificada
for x = 1:M
    for y = 1:N
        if Igray(x,y) <= 70 % Intensidades. Estaba  170.
            Ibin(x,y) = 0;
        else
            Ibin(x,y) = 255;
        end
    end
end




% Primera matriz de etiquetado (Para detectar regiones peque?as indeseadas)
[ob,nob] = bwlabel(Ibin);

tem = zeros(1,nob);
% Detecci?n de regiones indeseadas por tama?o
if nob > 0
    for s = 1:nob
        tem(s) = length(ob(ob==s));  % Almacena el tama?o de cada region
        if tem(s) < 1500 % Menores a este valor son regiones indeseadas
            % Quitamos las regiones indeseadas
            for i = 1:M
                for j = 1:N
                    if ob(i,j) == s
                        Ibin(i,j) = 0; % Ponemos a 0 las regiones indeseadas
                    end
                end
            end
        end
    end
end





% Segunda matriz de etiquetado (Para detectar todas las fichas)
[Iobjeto,numobjeto] = bwlabel(Ibin);

pos=cell(1,numobjeto); %imagen individual en tama?o de captura
posIMGcolor=cell(1,numobjeto); %imagen de carta cortada ( solo la carta )
posUbi=cell(1,numobjeto); %posicion de carta detectada
imgEnderezada=cell(1,numobjeto); %imagen enderezada
numValCarta=cell(1,numobjeto);
imgAreaBusqueda=cell(1,numobjeto); %imagen de busqueda de numero
posIMGNUM=cell(1,numobjeto); %Imagen recortada numero
posIMGSim=cell(1,numobjeto); %Imagen recortada simbolo


% Segunda matriz de etiquetado (Para detectar todas las fichas)
[Iobjeto,numobjeto] = bwlabel(Ibin);

pos=cell(1,numobjeto); %imagen individual en tama?o de captura
posIMGcolor=cell(1,numobjeto); %imagen de carta cortada ( solo la carta )
posUbi=cell(1,numobjeto); %posicion de carta detectada
imgEnderezada=cell(1,numobjeto); %imagen enderezada
numValCarta=cell(1,numobjeto);
imgAreaBusqueda=cell(1,numobjeto); %imagen de busqueda de numero
posIMGNUM=cell(1,numobjeto); %Imagen recortada numero
posIMGSim=cell(1,numobjeto); %Imagen recortada simbolo


if numobjeto>0 
    for num=1:numobjeto
        
        [xTemp,yTemp]=find(Iobjeto==num);
        % 
        posUbi{num}=[xTemp,yTemp];
        
        
       If=zeros(M,N); 
       for x=1:M
           for y=1:N
               if Iobjeto(x,y)==num
                  If(x,y)=255;    
               end   
           end    
       end
        pos{num}=If;

    end


for i=1:numobjeto
    x=posUbi{i}(:,1);
    y=posUbi{i}(:,2);
    im_prueba=I(min(x):max(x),min(y):max(y),:,:,:);
    posIMGcolor{i}=im_prueba;
%       figure()
%       imshow(im_prueba)
%       impixelinfo;
end


for i=1:numobjeto
    [imgEnderezada{i}, numValCarta{i}, imgAreaBusqueda{i}, posIMGNUM{i}, posIMGSim{i}]=procesar(posIMGcolor{i});
end




index = 1;
axes(handles.axes2);
imshow(imgEnderezada{index})
   
axes(handles.axes3);
imshow(imgAreaBusqueda{index})

axes(handles.axes4);
imshow(posIMGNUM{index})
   
axes(handles.axes5);
imshow(posIMGSim{index})   


if (turnoJuego==1)
    set(handles.axes6,'Visible','on');
    axes(handles.axes6);
    imshow(posIMGNUM{index}) 
    rondaJuego=rondaJuego+1;
    turnoJuego=0;
    set(handles.text10,'String', 'Turno PC');
    cartaUsuario=numValCarta{index}(2);
    if (rondaJuego>=3)
       if (cartaUsuario > cartaPC)
           set(handles.text10,'String', '!GANA USUARIO!');
           set(handles.pushbutton8,'Enable', 'on');
           set(handles.pushbutton9,'Enable', 'off');
       end
    end
    
elseif (turnoJuego==0)
    set(handles.axes7,'Visible','on');
    axes(handles.axes7);
    imshow(posIMGNUM{index}) 
    rondaJuego=rondaJuego+1;
    turnoJuego=0;
    cartaPC=numValCarta{index}(2);
   set(handles.text10,'String', 'Turno Usuario');
    if (rondaJuego>=3)
       if (cartaUsuario < cartaPC)
           set(handles.text10,'String', '!GANA PC!');
           set(handles.pushbutton8,'Enable', 'on');
           set(handles.pushbutton9,'Enable', 'off');
       end
    end
end

textoNum=numValCarta{index}(2);

if (numValCarta{index}(1)>10)
     textoSim = "Rojo";
else
     textoSim ="Negro";
end       

textoProp=numValCarta{index}(3);

if (textoSim=="Rojo")
    if (textoProp<indiceDiamante)
        
        textFormSim="Diamente";
    else
        textFormSim="Corazon";
    end
else
    if (textoProp<indiceEspada)
        textFormSim="Espada";
    else
        textFormSim="Trebol";
    end
end

    if (numValCarta{index}(4) == 1)
        if textoNum == 12
            esRey="Es reyna";
        else
            esRey="Es rey";
        end
    else
        esRey="No es rey";
    end

if (textoNum <= 10 && textoNum > 0)
            numeroCarta=num2str(textoNum);
    elseif (textoNum == 11)
           numeroCarta="J";
    elseif (textoNum == 12)
           numeroCarta="Q";
    elseif (textoNum == 13)
           numeroCarta="K";
else
    numeroCarta="Volteado";
    textoSim = "Azul";
    textFormSim = "No hay";
end




set(handles.text2,'String',strcat('El valor de la carta es : ',numeroCarta));   
set(handles.text3,'String',strcat('Color : ',textoSim));   
set(handles.text4,'String',strcat('Simbolo : ',textFormSim));   
set(handles.text5,'String',esRey); 
    set(handles.pushbutton5,'Enable','off');
    set(handles.pushbutton6,'Enable','off');



end



% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.axes6,'visible','off');
set(handles.axes7,'visible','off');
cla(handles.axes6);
cla(handles.axes7);
set(handles.pushbutton8,'Enable','on');
set(handles.pushbutton9,'Enable','off');
set(handles.text10,'String', '');
set(handles.text11,'String', '');
set(handles.text12,'String', '');
set(handles.text13,'String', '');



% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global defaultPanel defaultFigure

set(handles.figure1,'OuterPosition',defaultFigure);
set(handles.axes6,'visible','off');
set(handles.axes7,'visible','off');
cla(handles.axes6);
cla(handles.axes7);
set(handles.pushbutton8,'Enable','on');
set(handles.pushbutton9,'Enable','off');
set(handles.text10,'String', '');
set(handles.text11,'String', '');
set(handles.text12,'String', '');
set(handles.text13,'String', '');
%set(handles.uipanel1,'OuterPosition',defaultPanel);


% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global index imgEnderezada imgAreaBusqueda posIMGNUM posIMGSim numValCarta posIMREY posIMGcolor defaultFigure areainfo

set(handles.figure1,'OuterPosition',[defaultFigure(1) defaultFigure(2) 240 defaultFigure(4)]);
set(handles.uipanel4,'visible','on');
length(posIMREY{index})
axes(handles.axes8);
imshow(posIMGcolor{index})
axes(handles.axes9);
imshow(imgEnderezada{index})
axes(handles.axes10);
imshow(imgAreaBusqueda{index})
axes(handles.axes11);
imshow(posIMGSim{index})
axes(handles.axes12);
imshow(posIMGNUM{index})
set(handles.text20,'String',strcat('Area Amarilla : ',num2str(areainfo{index}(1))));
set(handles.text21,'String',strcat('Area Azul : ',num2str(areainfo{index}(2))));
if(length(posIMREY{index})> 1)
    set(handles.axes13,'visible','on');
    axes(handles.axes13);
    imshow(posIMREY{index})
    set(handles.text22,'String',strcat('Area Amarilla : ',num2str(areainfo{index}(4))));
    set(handles.text23,'String',strcat('Area Rojo : ',num2str(areainfo{index}(3))));
    set(handles.text19,'visible','on');
else
    cla(handles.axes13)
    set(handles.axes13,'visible','off');
    set(handles.text22,'String','');
    set(handles.text23,'String','');
    set(handles.text19,'visible','off');
    
end
