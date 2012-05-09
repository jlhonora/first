% Used to write in png and eps the different images obtained with FIRST and
% IDEAL

plotting = 1;
saving = 1;
%acq = 1; % 1 = thigh, 2 = brain
acq = 2; fileprefix = './exportimgs/thigh/thighR205';
%acq = 1; fileprefix = './exportimgs/brain/brainR2';
% For file formats refer to
% http://www.mathworks.com/help/techdoc/ref/print.html
fileformat = '-depsc'; postproc = 0; % 0 for eps
%fileformat = '-dpng'; filesuffix = '.png'; postproc = 1; % 1 for png

% Bounding box for png images 
% For 128x128 png
% Box = [1347 531;    %X,Y upper left
%        3450 2634;]; %X,Y lower right
% For 256x256 png
Box = [1053 317;    %X,Y upper left
       3747 3011;]; %X,Y lower right

% Boundig box for images with colorbar
% For 256x256 png and Colorbar   
BoxColorbar = [937, 566;    %X,Y upper left
               3552 2764;]; %X,Y lower right   
   
ndpi = '-r600';

% Position of the identifying letter
letterposition = [8 235];%[5 115];
fontsize = 18;
withtext = 0;

witharrows = 0;

% Image Limits
imgbox = [1 256;  % Rows
          1 256]; % Cols
      
% Intensity Ranges
% For brain
if(acq==1)
fatabs = [0 6];
waterabs = [0 6];
sumrange = [0 8];
end

% For thigh
if(acq==2)
fatabs = [0 1];
waterabs = [0 1];
sumrange = [0 1.4];
end

phirange = [-100 100];
R2range = [0 350;];

% thres = 0.1*mean(mean(mean(m_total)));
% fatmask = ((m_total(imgbox(1):imgbox(3),imgbox(2):imgbox(4),2)>thres)...
%          + (m_ideal(imgbox(1):imgbox(3),imgbox(2):imgbox(4),2)>thres))>=1;

fatmask = ones(size(phi_total));
Ncols = size(phi_total,2);
Nrows = size(phi_total,1);

if(0)
    figure, % Object
    setwhitebackground
    subplot(2,2,1)
    imshow(abs(m_total(imgbox(1):imgbox(3),imgbox(2):imgbox(4),1)),[0,12]) % Water
    title('Water')
    ylabel('Optim')
    subplot(2,2,2)
    imshow(abs(m_total(imgbox(1):imgbox(3),imgbox(2):imgbox(4),2)),[0,12]) % Fat
    title('Fat')
    colorbar

    subplot(2,2,3)
    imshow(abs(m_ideal(imgbox(1):imgbox(3),imgbox(2):imgbox(4),1)),[0,12]) % Water
    ylabel('IDEAL')
    subplot(2,2,4)
    imshow(abs(m_ideal(imgbox(1):imgbox(3),imgbox(2):imgbox(4),2)),[0,12]) % Fat
    colorbar

    figure, % Phi
    setwhitebackground
    subplot(1,2,1)
    imshow((phi_total(imgbox(1):imgbox(3),imgbox(2):imgbox(4))),[-130 130])
    title('B0 - Optim')
    colorbar
    subplot(1,2,2)
    imshow((phi_ideal(imgbox(1):imgbox(3),imgbox(2):imgbox(4))),[-130 130])
    title('B0 - IDEAL')
    colorbar

    figure, % R2*
    setwhitebackground
    subplot(1,2,1)
    fatmask = ((m_total(imgbox(1):imgbox(3),imgbox(2):imgbox(4),2)>0.1)...
             + (m_ideal(imgbox(1):imgbox(3),imgbox(2):imgbox(4),2)>0.1))>=1;
    imshow((R2_total(imgbox(1):imgbox(3),imgbox(2):imgbox(4))).*fatmask,[0 250])
    colorbar
    title('R2* - Optim')
    subplot(1,2,2)
    imshow((R2_map(imgbox(1):imgbox(3),imgbox(2):imgbox(4))).*fatmask,[0 250])
    title('R2* - IDEAL')
    colorbar

    figure, % Species sum
    setwhitebackground
    subplot(1,2,1)
    imshow(abs(m_total(imgbox(1):imgbox(3),imgbox(2):imgbox(4),1))...
         + abs(m_total(imgbox(1):imgbox(3),imgbox(2):imgbox(4),2)),[0 12])
    colorbar
    title('Species sum - Optim')
    subplot(1,2,2)
    imshow(abs(m_ideal(imgbox(1):imgbox(3),imgbox(2):imgbox(4),1))...
         + abs(m_ideal(imgbox(1):imgbox(3),imgbox(2):imgbox(4),2)),[0 12])
    title('Species sum - IDEAL')
    colorbar
end

% Saves in files
if(saving)
    if(gcf==1)
        figure(gcf), % Object - FIRST
    else
        figure(gcf+1)
    end
    setwhitebackground
    imshow(abs(m_total(imgbox(1):imgbox(3),imgbox(2):imgbox(4),1)),waterabs) % Water
    if(withtext)
        text(letterposition(1),letterposition(2),'\bfa','FontSize',fontsize,'Color',[1 1 1], 'EdgeColor', [0 0 0])
    end
    if(witharrows)
        if(acq==1)
            x = [28 38]; y = 256-[156 146];
            [x,y] = dsxy2figxy(x,y);
           annotation(gcf,'arrow','X',x,'Y',y, 'HeadStyle', 'plain', 'HeadWidth', 4, 'HeadLength',4,...
               'LineWidth', 1.5, 'Color', [1 1 1]) 
            x = [217 201]; y = 256-[101 99];
            [x,y] = dsxy2figxy(x,y);
           annotation(gcf,'arrow','X',x,'Y',y, 'HeadStyle', 'plain', 'HeadWidth', 4, 'HeadLength',4,...
               'LineWidth', 1.5, 'Color', [1 1 1]) 
            x = [208 221]; y = 256-[201 211];
            [x,y] = dsxy2figxy(x,y);
           annotation(gcf,'arrow','X',x,'Y',y, 'HeadStyle', 'plain', 'HeadWidth', 4, 'HeadLength',4,...
               'LineWidth', 1.5, 'Color', [1 1 1]) 
        end
    end
    filename = [fileprefix 'waterabsfirst'];
    print('-f',ndpi,fileformat,filename)
    if(postproc)
        I = imread([filename filesuffix]);
        imwrite(I(Box(3):Box(4),Box(1):Box(2)),[filename filesuffix]);
    end
    figure(gcf+plotting)
    setwhitebackground
    imshow(abs(m_total(imgbox(1):imgbox(3),imgbox(2):imgbox(4),2)),fatabs) % Fat
    if(withtext)
        text(letterposition(1),letterposition(2),'\bfb','FontSize',fontsize,'Color',[1 1 1], 'EdgeColor', [0 0 0])
    end
    if(witharrows)
%         x = [150 150]/Ncols; y = [160, 160]/Nrows;
%        annotation('arrow',x,y, 'HeadStyle', 'vback1', 'LineWidth', 1.2, 'Color', [1 1 1]) 
    end    
    filename = [fileprefix 'fatabsfirst'];
    print('-f',ndpi,fileformat,filename)
    if(postproc)
        I = imread([filename filesuffix]);
        imwrite(I(Box(3):Box(4),Box(1):Box(2)),[filename filesuffix]);
    end


    figure(gcf+plotting), % Object - IDEAL
    setwhitebackground
    imshow(abs(m_ideal(imgbox(1):imgbox(3),imgbox(2):imgbox(4),1)),waterabs) % Water
    if(withtext)
        text(letterposition(1),letterposition(2),'\bfc','FontSize',fontsize,'Color',[1 1 1]);%, 'BackgroundColor', [0 0 0])
    end
    if(witharrows)
       if(acq==1)
            x = [28 38]; y = 256-[156 146];
            [x,y] = dsxy2figxy(x,y);
           annotation(gcf,'arrow','X',x,'Y',y, 'HeadStyle', 'plain', 'HeadWidth', 4, 'HeadLength',4,...
               'LineWidth', 1.5, 'Color', [1 1 1]) 
            x = [217 201]; y = 256-[101 99];
            [x,y] = dsxy2figxy(x,y);
           annotation(gcf,'arrow','X',x,'Y',y, 'HeadStyle', 'plain', 'HeadWidth', 4, 'HeadLength',4,...
               'LineWidth', 1.5, 'Color', [1 1 1]) 
            x = [208 221]; y = 256-[201 211];
            [x,y] = dsxy2figxy(x,y);
           annotation(gcf,'arrow','X',x,'Y',y, 'HeadStyle', 'plain', 'HeadWidth', 4, 'HeadLength',4,...
               'LineWidth', 1.5, 'Color', [1 1 1]) 
        end
    end
    filename = [fileprefix 'waterabsideal'];
    print('-f',ndpi,fileformat,filename)
    if(postproc)
        I = imread([filename filesuffix]);
        imwrite(I(Box(3):Box(4),Box(1):Box(2)),[filename filesuffix]);
    end
    figure(gcf+plotting)
    setwhitebackground
    imshow(abs(m_ideal(imgbox(1):imgbox(3),imgbox(2):imgbox(4),2)),fatabs) % Fat
    if(withtext)
        text(letterposition(1),letterposition(2),'\bfd','FontSize',fontsize,'Color',[1 1 1]);%, 'BackgroundColor', [0 0 0])
    end
    if(witharrows)
       x = [150 150]/Ncols; y = [160, 160]/Nrows;
       annotation('arrow',x,y, 'HeadStyle', 'vback1', 'LineWidth', 1.2, 'Color', [1 1 1]) 
    end
    filename = [fileprefix 'fatabsideal'];
    print('-f',ndpi,fileformat,filename)
    if(postproc)
        I = imread([filename filesuffix]);
        imwrite(I(Box(3):Box(4),Box(1):Box(2)),[filename filesuffix]);
    end

    figure(gcf+plotting) % Phi - FIRST
    setwhitebackground    
    imshow((phi_total(imgbox(1):imgbox(3),imgbox(2):imgbox(4))),phirange)    
    colorbar, mBox = BoxColorbar;
    if(withtext)
        text(letterposition(1),letterposition(2),'\bfa','FontSize',fontsize,'Color',[1 1 1]);%, 'BackgroundColor', [0 0 0])
    end
    filename = [fileprefix 'phifirst'];
    print('-f',ndpi,fileformat,filename)
    if(postproc)
        I = imread([filename filesuffix]);
        imwrite(I(mBox(3):mBox(4),mBox(1):mBox(2)),[filename filesuffix]);
    end
    
    figure(gcf+plotting) % Phi - IDEAL
    setwhitebackground 
    imshow((phi_ideal(imgbox(1):imgbox(3),imgbox(2):imgbox(4))),phirange)
    colorbar, mBox = BoxColorbar;
    if(withtext)
        text(letterposition(1),letterposition(2),'\bfb','FontSize',fontsize,'Color',[1 1 1], 'BackgroundColor', [0 0 0])
    end
    filename = [fileprefix 'phiideal'];
    print('-f',ndpi,fileformat,filename)
    if(postproc)
        I = imread([filename filesuffix]);
        imwrite(I(mBox(3):mBox(4),mBox(1):mBox(2)),[filename filesuffix]);
    end
    
    figure(gcf+plotting) % R2* - FIRST
    setwhitebackground    
    imshow((R2_total(imgbox(1):imgbox(3),imgbox(2):imgbox(4))).*fatmask,R2range)   
    colorbar, mBox = BoxColorbar;
    if(withtext)
        text(letterposition(1),letterposition(2),'\bfa','FontSize',fontsize,'Color',[1 1 1]);%, 'BackgroundColor', [0 0 0])
    end
    if(witharrows)
       x = [150 150]/Ncols; y = [160, 160]/Nrows;
       annotation('textarrow',x,y,'String','R2* = 1',...
           'HeadStyle', 'vback1', 'LineWidth', 1.2, 'Color', [1 1 1], 'FontSize', 18, 'TextColor', [1 1 1])
       %annotation('rectangle', [100,100,10,20]) % [x y w h]
    end
    filename = [fileprefix 'R2first'];
    print('-f',ndpi,fileformat,filename)
    if(postproc)
        I = imread([filename filesuffix]);
        imwrite(I(mBox(3):mBox(4),mBox(1):mBox(2)),[filename filesuffix]);
    end
    
    figure(gcf+plotting) % R2* - IDEAL
    setwhitebackground 
    imshow((R2_map(imgbox(1):imgbox(3),imgbox(2):imgbox(4))).*fatmask,R2range)
    colorbar, mBox = BoxColorbar;
    if(witharrows)
       x = [150 150]/Ncols; y = [160, 160]/Nrows;
       annotation('textarrow',x,y,'HeadStyle', 'vback1', ...
           'LineWidth', 1.2, 'Color', [1 1 1], 'FontSize', 18, 'TextColor', [1 1 1])
       %annotation('rectangle', [100,100,10,20]) % [x y w h]
    end
    if(withtext)
        text(letterposition(1),letterposition(2),'\bfb','FontSize',fontsize,'Color',[1 1 1], 'BackgroundColor', [0 0 0])
    end
    filename = [fileprefix 'R2ideal'];
    print('-f',ndpi,fileformat,filename)
    if(postproc)
        I = imread([filename filesuffix]);
        imwrite(I(mBox(3):mBox(4),mBox(1):mBox(2)),[filename filesuffix]);
    end
    
    figure(gcf+plotting) % R2* - Fitted
    setwhitebackground 
    imshow((R2_fitted(imgbox(1):imgbox(3),imgbox(2):imgbox(4))).*fatmask,R2range)
    colorbar, mBox = BoxColorbar;
    if(witharrows)
       x = [150 150]/Ncols; y = [160, 160]/Nrows;
       annotation('textarrow',x,y,'HeadStyle', 'vback1', ...
           'LineWidth', 1.2, 'Color', [1 1 1], 'FontSize', 18, 'TextColor', [1 1 1])
       %annotation('rectangle', [100,100,10,20]) % [x y w h]
    end
    if(withtext)
        text(letterposition(1),letterposition(2),'\bfc','FontSize',fontsize,'Color',[1 1 1], 'BackgroundColor', [0 0 0])
    end
    filename = [fileprefix 'R2fitted'];
    print('-f',ndpi,fileformat,filename)
    if(postproc)
        I = imread([filename filesuffix]);
        imwrite(I(mBox(3):mBox(4),mBox(1):mBox(2)),[filename filesuffix]);
    end
    
    figure(gcf+plotting), % Species sum - FIRST
    setwhitebackground
    imshow(abs(m_total(imgbox(1):imgbox(3),imgbox(2):imgbox(4),1))...
         + abs(m_total(imgbox(1):imgbox(3),imgbox(2):imgbox(4),2)),sumrange)
    if(withtext)
        text(letterposition(1),letterposition(2),'\bfa','FontSize',fontsize,'Color',[1 1 1], 'BackgroundColor', [0 0 0])
    end
    if(witharrows)
       x = [150 150]/Ncols; y = [160, 160]/Nrows;
       annotation('textarrow',x,y, 'HeadStyle', 'vback1', 'LineWidth', 1.2, 'Color', [1 1 1]) 
    end
    filename = [fileprefix 'sumfirst'];
    print('-f',ndpi,fileformat,filename)
    if(postproc)
        I = imread([filename filesuffix]);
        imwrite(I(Box(3):Box(4),Box(1):Box(2)),[filename filesuffix]);
    end

    figure(gcf+plotting), % Species sum - IDEAL
    setwhitebackground
    imshow(abs(m_ideal(imgbox(1):imgbox(3),imgbox(2):imgbox(4),1))...
         + abs(m_ideal(imgbox(1):imgbox(3),imgbox(2):imgbox(4),2)),sumrange)
    if(witharrows)
       x = [150 150]/Ncols; y = [160, 160]/Nrows;
       annotation('textarrow',x,y, 'HeadStyle', 'vback1', 'LineWidth', 1.2, 'Color', [1 1 1]) 
    end
    if(withtext)
        text(letterposition(1),letterposition(2),'\bfb','FontSize',fontsize,'Color',[1 1 1], 'BackgroundColor', [0 0 0])
    end
    filename = [fileprefix 'sumideal'];
    print('-f',ndpi,fileformat,filename)
    if(postproc)
        I = imread([filename filesuffix]);
        imwrite(I(Box(3):Box(4),Box(1):Box(2)),[filename filesuffix]);
    end
    
    if(plotting==0)
        close gcf
    end
end


%%% ANOTHER VERSION
% 
% plotting = 0;
% saving = 1;
% fileprefix = './exportimgs/thighs';
% % For file formats refer to
% % http://www.mathworks.com/help/techdoc/ref/print.html
% fileformat = '-depsc';
% fileformat = '-dpng'; filesuffix = '.png';
% postproc = 1;
% Box = [1347 531;    %X,Y upper left
%        3450 2634;]; %X,Y lower right
% ndpi = '-r600';
% 
% letterposition = [5 115];
% fontsize = 14;
% 
% % Image Limits
% imgbox = [1 128;  % Rows
%           1 128]; % Cols
% 
% if(0)
%     figure, % Object
%     setwhitebackground
%     subplot(2,2,1)
%     imshow(abs(m_total(imgbox(1):imgbox(3),imgbox(2):imgbox(4),1)),[0,12]) % Water
%     title('Water')
%     ylabel('Optim')
%     subplot(2,2,2)
%     imshow(abs(m_total(imgbox(1):imgbox(3),imgbox(2):imgbox(4),2)),[0,12]) % Fat
%     title('Fat')
%     colorbar
% 
%     subplot(2,2,3)
%     imshow(abs(m_ideal(imgbox(1):imgbox(3),imgbox(2):imgbox(4),1)),[0,12]) % Water
%     ylabel('IDEAL')
%     subplot(2,2,4)
%     imshow(abs(m_ideal(imgbox(1):imgbox(3),imgbox(2):imgbox(4),2)),[0,12]) % Fat
%     colorbar
% 
%     figure, % Phi
%     setwhitebackground
%     subplot(1,2,1)
%     imshow((phi_total(imgbox(1):imgbox(3),imgbox(2):imgbox(4))),[-130 130])
%     title('B0 - Optim')
%     colorbar
%     subplot(1,2,2)
%     imshow((phi_ideal(imgbox(1):imgbox(3),imgbox(2):imgbox(4))),[-130 130])
%     title('B0 - IDEAL')
%     colorbar
% 
%     figure, % R2*
%     setwhitebackground
%     subplot(1,2,1)
%     fatmask = ((m_total(imgbox(1):imgbox(3),imgbox(2):imgbox(4),2)>0.1)...
%              + (m_ideal(imgbox(1):imgbox(3),imgbox(2):imgbox(4),2)>0.1))>=1;
%     imshow((R2_total(imgbox(1):imgbox(3),imgbox(2):imgbox(4))).*fatmask,[0 250])
%     colorbar
%     title('R2* - Optim')
%     subplot(1,2,2)
%     imshow((R2_map(imgbox(1):imgbox(3),imgbox(2):imgbox(4))).*fatmask,[0 250])
%     title('R2* - IDEAL')
%     colorbar
% 
%     figure, % Species sum
%     setwhitebackground
%     subplot(1,2,1)
%     imshow(abs(m_total(imgbox(1):imgbox(3),imgbox(2):imgbox(4),1))...
%          + abs(m_total(imgbox(1):imgbox(3),imgbox(2):imgbox(4),2)),[0 12])
%     colorbar
%     title('Species sum - Optim')
%     subplot(1,2,2)
%     imshow(abs(m_ideal(imgbox(1):imgbox(3),imgbox(2):imgbox(4),1))...
%          + abs(m_ideal(imgbox(1):imgbox(3),imgbox(2):imgbox(4),2)),[0 12])
%     title('Species sum - IDEAL')
%     colorbar
% end
% 
% % Saves in files
% if(saving)
%     if(gcf==1)
%         figure(gcf), % Object - FIRST
%     else
%         figure(gcf+1)
%     end
%     object = abs(m_total(:,:,1));
%     range = [0 12];
%     filename = 'waterabsfirst';
%     letter = 'a';
%     lettercolor = [1 1 1];
%     ShowAndWrite(object, range, imgbox, ndpi, letter, letterposition, ...
%              lettercolor, fontsize, filename, fileformat, fileprefix, ...
%              filesuffix, postproc, Box, plotting);
%          
%     object = abs(m_total(:,:,2));
%     range = [0 12];
%     filename = 'fatabsfirst';
%     letter = 'b';
%     lettercolor = [1 1 1];
%     ShowAndWrite(object, range, imgbox, ndpi, letter, letterposition, ...
%              lettercolor, fontsize, filename, fileformat, fileprefix, ...
%              filesuffix, postproc, Box, plotting);
% 
% 
% 	% Object - IDEAL
%     object = abs(m_ideal(:,:,1));
%     range = [0 12];
%     filename = 'waterabsideal';
%     letter = 'c';
%     lettercolor = [1 1 1];
%     ShowAndWrite(object, range, imgbox, ndpi, letter, letterposition, ...
%              lettercolor, fontsize, filename, fileformat, fileprefix, ...
%              filesuffix, postproc, Box, plotting);    
% 
%     object = abs(m_ideal(:,:,2));
%     range = [0 12];
%     filename = 'fatabsideal';
%     letter = 'd';
%     lettercolor = [1 1 1];
%     ShowAndWrite(object, range, imgbox, ndpi, letter, letterposition, ...
%              lettercolor, fontsize, filename, fileformat, fileprefix, ...
%              filesuffix, postproc, Box, plotting);  
% 
%     % Phi - FIRST
%     object = phi_total;
%     range = [-130 130];
%     filename = 'phifirst';
%     letter = 'a';
%     lettercolor = [1 1 1];
%     ShowAndWrite(object, range, imgbox, ndpi, letter, letterposition, ...
%              lettercolor, fontsize, filename, fileformat, fileprefix, ...
%              filesuffix, postproc, Box, plotting);  
%     
%     % Phi - IDEAL
%     object = phi_ideal;
%     range = [-130 130];
%     filename = 'phiideal';
%     letter = 'b';
%     lettercolor = [1 1 1];
%     ShowAndWrite(object, range, imgbox, ndpi, letter, letterposition, ...
%              lettercolor, fontsize, filename, fileformat, fileprefix, ...
%              filesuffix, postproc, Box, plotting); 
%     
%     fatmask = ((m_total(:,:,2)>0.1)...
%              + (m_ideal(:,:,2)>0.1))>=1;
%     % R2* - FIRST
%     object = R2_total.*fatmask;
%     range = [0 250];
%     filename = 'R2first';
%     letter = 'a';
%     lettercolor = [1 1 1];
%     ShowAndWrite(object, range, imgbox, ndpi, letter, letterposition, ...
%              lettercolor, fontsize, filename, fileformat, fileprefix, ...
%              filesuffix, postproc, Box, plotting); 
%     
%     % R2* - IDEAL
%     object = R2_map.*fatmask;
%     range = [0 250];
%     filename = 'R2ideal';
%     letter = 'b';
%     lettercolor = [1 1 1];
%     ShowAndWrite(object, range, imgbox, ndpi, letter, letterposition, ...
%              lettercolor, fontsize, filename, fileformat, fileprefix, ...
%              filesuffix, postproc, Box, plotting); 
%     
%     % Species sum - FIRST
%     object = abs(m_total(:,:,1))...
%            + abs(m_total(:,:,2));
%     range = [0 12];
%     filename = 'sumfirst';
%     letter = 'a';
%     lettercolor = [1 1 1];
%     ShowAndWrite(object, range, imgbox, ndpi, letter, letterposition, ...
%              lettercolor, fontsize, filename, fileformat, fileprefix, ...
%              filesuffix, postproc, Box, plotting); 
% 
%     figure(gcf+plotting), % Species sum - IDEAL
%     object = abs(m_ideal(:,:,1))...
%            + abs(m_ideal(:,:,2));
%     range = [0 12];
%     filename = 'sumideal';
%     letter = 'b';
%     lettercolor = [1 1 1];
%     ShowAndWrite(object, range, imgbox, ndpi, letter, letterposition, ...
%              lettercolor, fontsize, filename, fileformat, fileprefix, ...
%              filesuffix, postproc, Box, plotting); 
%     
%     if(plotting==0)
%         close gcf
%     end
% end

% function [] = ShowAndWrite(object, range, imgbox, ndpi, letter, letterposition, lettercolor, fontsize, filename, fileformat, fileprefix, filesuffix, postproc, Box, plotting)
% % Auxiliar function - not used
% figure(gcf+plotting),
% setwhitebackground
% imshow(object(imgbox(1):imgbox(3),imgbox(2):imgbox(4)),range)
% text(letterposition(1),letterposition(2),['\bf' letter],'FontSize',fontsize,'Color',lettercolor)
% filename = [fileprefix filename];
% print('-f',ndpi,fileformat,filename)
% if(postproc)
%     I = imread([filename filesuffix]);
%     imwrite(I(Box(3):Box(4),Box(1):Box(2)),[filename filesuffix]);
% end
    
    