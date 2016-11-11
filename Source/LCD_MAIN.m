%% import data-------------------------------------------------------------
clear;
clc;
% user inputs the number of the respective distributions
num_neg = input('Enter the number of negative distributions. 0-10. ');
num_pos = input('Enter the number of positive distributions. 0-10. ');

% user inputs whether the flash is IC or not
IC = input('Enter 1 if the flash is intracloud. ');

%get which field mills to use for analysis
useMills = input('Enter the field mills to use. ', 's');
mills = importdata('mills.txt'); %file containing field mill coordinates
useMills = str2num(useMills);

% put the coordinates of the field mills in matrix coords
for k = useMills(:)
    coords = [mills(k,1), mills(k,2), mills(k,3)];
end 
[numMills, ~] = size(coords);

%import E-field data for field mills
eData = importdata('field.txt');

% create cell structure to hold the negative and positive LDAR/LMA data
% for each distribution.
negatives = cell(1,num_neg);
positives = cell(1,num_pos);

% import negative and positive distributions
% each distribution's coordinate data is held in .txt files denoted:
% neg1.txt, neg2.txt, pos1.txt, pos2.txt
for k = 1:num_neg
    filename = sprintf('neg%d.txt',k);
    negatives{k} = importdata(filename);
    negatives{k}(:,3) = negatives{k}(:,3); 
end

for k = 1:num_pos
    filename = sprintf('pos%d.txt',k);
    positives{k} = importdata(filename);
    positives{k}(:,3) = positives{k}(:,3); 
end


%% compute image charges---------------------------------------------------

% create cell to hold the image charges for each positive/negative
% distribution
neg_images = cell(1,num_neg);
pos_images = cell(1,num_pos);

% for each image charge, just flip each coordinate point about the xy-plane
for k = 1:num_neg
    p = negatives{k};
    p = [p(:,1), p(:,2), p(:,3) .* -1];
    neg_images{k} = p;
end

for k = 1:num_pos
    p = positives{k};
    p = [p(:,1), p(:,2), p(:,3) .* -1];
    pos_images{k} = p;
end

%% run PCA on the points and the images.-----------------------------------

% the layout for the pcaPoints and pcaImages cells is:
% {k,1}     {k,2}   {k,3}   {k,4}   {k,5}
%  mean    dirVect    t    endpts    len
%  k corresponds to the particular charge distribution, 1-10.

% create cells for the pca results for the distribution and their image
% charges
pcaNegatives = cell(num_neg,5);
pcaNegImages = cell(num_neg,5);
pcaPositives = cell(num_pos,5);
pcaPosImages = cell(num_pos,5);

%pca for negatives
for k = 1:num_neg
    [a, b, c, d, e] = princom(negatives{k});
    pcaNegatives{k,1} = a;
    pcaNegatives{k,2} = b;
    pcaNegatives{k,3} = c;
    pcaNegatives{k,4} = d;
    pcaNegatives{k,5} = e;
end

%pca for negative images
for k = 1:num_neg
    [a, b, c, d, e] = princom(neg_images{k});
    pcaNegImages{k,1} = a;
    pcaNegImages{k,2} = b;
    pcaNegImages{k,3} = c;
    pcaNegImages{k,4} = d;
    pcaNegImages{k,5} = e;
end

%pca for positives
for k = 1:num_pos
    [a, b, c, d, e] = princom(positives{k});
    pcaPositives{k,1} = a;
    pcaPositives{k,2} = b;
    pcaPositives{k,3} = c;
    pcaPositives{k,4} = d;
    pcaPositives{k,5} = e;
end

%pca for positive images
for k = 1:num_pos
    [a, b, c, d, e] = princom(pos_images{k});
    pcaPosImages{k,1} = a;
    pcaPosImages{k,2} = b;
    pcaPosImages{k,3} = c;
    pcaPosImages{k,4} = d;
    pcaPosImages{k,5} = e;
end

%%  integrate over PCA lines-----------------------------------------------

% initialize matrix to hold the results from the integration
coeff = zeros(numMills, num_pos+num_neg);

%integrate over negatives
for k = 1:numMills
    for j = 1:num_neg
        A = intergrl(j, coords(k,:), pcaNegatives);
        B = intergrl(j, coords(k,:), pcaNegImages);
        coeff(k,j) = -A + B;
    end
end

%integrate over positives
for k = 1:numMills
    for j = 1:num_pos
        A = intergrl(j, coords(k,:), pcaPositives);
        B = intergrl(j, coords(k,:), pcaPosImages);
        coeff(k,j+num_neg) = A - B;
    end
end

%% solve for the charge on each line---------------------------------------
if IC == 1
    % if the flash is IC, perform weighted solution
    
    % create weight vector w:
    w = zeros(1, numMills + 1);
    w(1:numMills) = 0.0001; % weight of each data point
    w(numMills + 1) = 1000000; % weight of conserved charge
    
    % create a row on the bottom of coeff matrix for the constraint for 
    % charge conservation. the negative charges plus the positive charges 
    % must equal zero:
    % -Qneg1 - Qneg2 - ... - Q(negs) + Qpos1 + Qpos2 + ... + Q(pos) = 0
    for count = 1:num_neg
        coeff(numMills + 1, count) = -1;
    end
    for count = num_neg+1:num_neg+num_pos
        coeff(numMills + 1, count) = 1;
    end
    eData(numMills + 1) = 0;
    
    % solve the weighted system:
    Q = lscov(coeff, eData, w');
else
    % for CG flashes, perform linear least squares solution with bounding
    % constraints: a negative charge must live inside (-inf, 0) and a
    % positive charge must live inside (0, inf):
    
    % the lower bounds for each q in Q:
    lb = [ -inf .* ones(1,num_neg), zeros(1,num_pos) ];
    ub = [ zeros(1,num_neg), inf .* ones(1,num_pos) ];
    Q = lsqlin(coeff, eData, [], [], [], [], lb, ub);
end

%% output results----------------------------------------------------------

% output data to out.txt
file = fopen('out.txt', 'w');
if IC ~= 1
    fprintf(file, 'CG Flash.\n\n');
else 
    fprintf(file, 'IC Flash. \n\n');
end

fprintf(file, 'Negatives:\n');
for i = 1:num_neg
   fprintf(file, 'Q%d = -%f C\n', i, Q(i)); 
end
fprintf(file, '\nPositives:\n');
for i = num_neg+1:num_neg+num_pos
    fprintf(file, 'Q%d = %f C\n', i, Q(i));
end

if IC ~= 1
    fprintf(file, '\nNet Charge:\n');
    fprintf(file, '%f C\n\n', sum(Q));
end 
fprintf(file, '\nNormalized E Fields:\n');
for j = 1:numMills
    fprintf(file, '%f V/m\n', coeff(j,:) * Q);
end

fprintf(file, '\nResiduals:\n');
for j = 1:numMills
   fprintf(file, '%f V/m\n', abs( (coeff(j,:) * Q) - eData(j))); 
end





