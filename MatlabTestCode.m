%% PIV CODE - Lab 3
% Brigham Ostergaard
% 9-15-2022

%% Section to change inputs
clc
clear all
close all

t = cputime;

% picture1 = "PIVlab_gen_2_01.tif";
% picture2 = "PIVlab_gen_2_02.tif";

picture1 = "A001_1.bmp";
picture2 = "A001_2.bmp";

M = 50;  %pixels for M/N of the search window size
N = M;
info = imfinfo(picture1);

height = info.Width;
width = info.Height;

overlap = 0.25;     %overlap percentage     
overlap = 1-overlap;

method = "MQD";      %three types: CC, MQD, FFT 

%Plots the first interrogation window, phi vs [m,n];
PlotFirst = true;

% The search window is set at a 10x10, which is m,n inside the function.

%Could just do pixels/frame, not pixels/second...
fprintf("Images Loaded. Beginning Processing...\n")

run_Code(width, height, M, N, overlap, method, picture1, picture2, PlotFirst, t)
%% This is some code I think will help
%This is the function we will turn in, along with the related functions.
function run_Code(width, height, M, N, overlap, method, picture1, picture2, PlotFirst, t)     
    fprintf("Spacing Interrogation Windows...\n")

    % We need to check to make sure n_x_windows (and y) don't allow us to 
    %   go outside of the pictures dimensions ever.


    n_x_windows = 0;
    current_pixels = 0;  
    overlap_pixels_M = floor(overlap*M);  %round down
    while (current_pixels + overlap_pixels_M + 1 <= width)
        if(current_pixels + overlap_pixels_M + 1 <= width)
            current_pixels = current_pixels + (overlap_pixels_M);
            n_x_windows = n_x_windows + 1;
        end
    end
    fprintf("Finished X, moving to Y...\n")

    n_y_windows = 0;
    current_pixels = 0;
    overlap_pixels_N = floor(overlap*N);  %round down
    while(current_pixels + overlap_pixels_N + 1 <= height)
        if(current_pixels + overlap_pixels_N + 1 <= height)
            current_pixels = current_pixels + (overlap_pixels_N);
            n_y_windows = n_y_windows + 1;
        end
    end
    fprintf("Picture Matrix Established...\n")
    %Now, we have the # of interrogation windows along the x and y axis. We can
    %do a loop now. This should account for the non-complete boxes on the right
    %and bottom edge, as we will ignore them (due to the <= in the setup above.
    
    m = 10; % x size of search window
    n = m; % y size of search window   
    
    fprintf("Running Coorelation Functions...\n")
    Picture_Matrix1 = double(imread(picture1)); %First picture
    Picture_Matrix2 = double(imread(picture2)); %second picture

    showGridPlot = false;   %Shows a second plot of the interrogation windows all on the picture
    %Now we will start a loop to go through each interrogation windows.
    if(showGridPlot == true)
        figure('Name', 'Test Search Window')
        imshow(picture1)
        hold on
    end
    x = zeros(n_x_windows,n_y_windows);
    y = zeros(n_x_windows,n_y_windows); 
    u = zeros(n_x_windows,n_y_windows); 
    v = zeros(n_x_windows,n_y_windows); 

    for i = 1:n_y_windows 
        starty = ((i-1)*overlap_pixels_N) + 1;
        endy = starty + N;
        for j = 1:n_x_windows
            %In this loop, we set g1 and g2 inside the interrogation window,
            %and run the algorithms to get the correlation values.

            startx = ((j-1)*overlap_pixels_M) + 1;
            endx = startx + M;

            if(showGridPlot == true)
                plot([starty, starty, endy, endy, starty], [startx,endx,endx,startx,startx], 'Color', 'g');
            end
            if(endx < width && endy < height)
                image1 = Picture_Matrix1(startx:endx, starty:endy);
                image2 = Picture_Matrix2(startx:endx, starty:endy);
    
                returnables = Search_Interrogation_Window(image1,image2, M, N, m,n, method, PlotFirst);
                PlotFirst = false;    
    
                x(i,j) = returnables(1);
                y(i,j) = returnables(2);
                u(i,j) = returnables(3);
                v(i,j) = returnables(4);
            %Now, we did set x and y, but we don't have them in the right spot, so we 
            % shift them, as they are referenced to the middle of the first window.
                x(i,j) = x(i,j) + (i-1)*(overlap_pixels_M); 
                y(i,j) = y(i,j) + (j-1)*(overlap_pixels_N);
            end
        end
    end
    if(showGridPlot == true)
        hold off
    end
    
    %We add the image and plot the quiver over top here.
    Name = "Vector Field - " + method;
    figure('Name', Name)
    imshow(picture1)
    hold on
    quiver(x,y,u,v, 'Color', 'Red', 'LineWidth', 2);

    hold off

    %Computational time of the program. Outputted to the screen
    total_time = cputime-t


end
%%%%%%%%%%%%%%%%%%%%--Functions--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function returnables = Search_Interrogation_Window(image1, image2, M, N, m,n, method, PlotFirst)
    
    %Search here, store the values of the middle in x and y. Return in Loc.
    cx = floor(M/2)+1; %if M were 100, it gives 51, if it were 101, it gives 51.
    cy = floor(N/2)+1;
    
    %g1 is the center of the matrix. Where is it?
    firstx = cx-floor(m/2);
    firsty = cy-floor(m/2);
    g1 = image1(firstx:firstx+m, firsty:firsty+n);

    Int_Value_Matrix = zeros(M-m,N-n);  %May be incorrect size. (m-1) and (n-1)?
    for i = 1:(M-m)
        % the M-m takes into account the width of the interrogation window
        for j = 1:(N-n)
%             fprintf("m = %f\tn = %f\n", m, n)
            g2 = image2(i:(i+m),j:(j+n));
            %The line above is a problem. It goes outside the range.
            % I can't subtract 1, because then g2 is smaller than g1.
            if(method == "CC")
                phi = CrossCorrelationFunc(g2,g1);
            elseif (method == "MQD")
                phi = MQDFunc(g2,g1);
            elseif (method == "FFT")
                Int_Value_Matrix = FFTFunc(g2,g1,PlotFirst);
            end

            if(method ~= "FFT")
                if(PlotFirst)
                    figure('Name', 'First Search Window')
                    surf(phi);
                    xlabel('m');
                    ylabel('n');
                    zlabel('phi');
                end

                if(method == "CC") 
                    phi = sum(sum(phi));
                elseif(method == "MQD") 
                    phi = 1/sum(sum(phi));
                end

                Int_Value_Matrix(i,j) = phi;

            end
            %Take the phi value and put in matrix
            PlotFirst = false;
        end
    end
    %Testing the max function to see how it works
    [V,i] = max(Int_Value_Matrix);
    [M,j] = max(max(Int_Value_Matrix));

    %cx/cy are the coordinates of the center of the int. window.
    dx = (i(1)-cx);
    dy = (j-cy);
    x = cx;
    y = cy;
    u = dx;
    v = dy;
    returnables = [x,y,u,v];
    %Returns the location, inside a matrix 
end
%% Correlation Functions
function phi = CrossCorrelationFunc(g2,g1)
    phi = (g1.*g2);
end
function phi = MQDFunc(g2,g1)
    phi = (g1-g2).^2; 
end
function phi = FFTFunc(g2, g1,plotSurf)
    g1p = fft2(g1);
    g2p = fft2(g2);
    phi = fftshift(ifft2(g1p.*conj(g2p)));
    if(plotSurf == true)
        figure('Name', 'First Search Window')
        surf(phi);
        xlabel('m');
        ylabel('n');
        zlabel('phi');
    end
    phi = (fftshift(ifft2(g1p.*conj(g2p))));
end