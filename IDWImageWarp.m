function im2 = IDWImageWarp(im, psrc, pdst)
%% basic image manipulations
% get image (matrix) size
[h,w,dim] = size(im);
im2 = im;
[X,Y] = meshgrid(1:w,1:h);

%% Compute sigma and initialze some variables
n = size(psrc,1);
distance_psrc = squareform(pdist(psrc));
distance_psrc_square = distance_psrc.^2 + diag(ones(1,n));
processed_ones = ones(n) - diag(diag(ones(n)));
sigma_ij = processed_ones ./ distance_psrc_square;

coef_TxTy = zeros(2,2);
b_Tx = zeros(2,1);
b_Ty = zeros(2,1);
t1 = zeros(n,1);
t2 = zeros(n,1);
t3 = zeros(n,1);
t4 = zeros(n,1);

%% Compute Tx and Ty
for i=1:n
    psrc_i = psrc - repmat(psrc(i,:),n,1);
    pdst_i = pdst - repmat(pdst(i,:),n,1);
    coef_TxTy(1,1) = sum(sigma_ij(i,:)'.* psrc_i(:,1).* psrc_i(:,1));
    coef_TxTy(1,2) = sum(sigma_ij(i,:)'.* psrc_i(:,1).* psrc_i(:,2));
    coef_TxTy(2,2) = sum(sigma_ij(i,:)'.* psrc_i(:,2).* psrc_i(:,2));
    coef_TxTy(2,1) = coef_TxTy(1,2);
    b_Tx(1,1) = sum(sigma_ij(i,:)'.* pdst_i(:,1).* psrc_i(:,1));
    b_Tx(2,1) = sum(sigma_ij(i,:)'.* pdst_i(:,1).* psrc_i(:,2));
    b_Ty(1,1) = sum(sigma_ij(i,:)'.* pdst_i(:,2).* psrc_i(:,1));
    b_Ty(2,1) = sum(sigma_ij(i,:)'.* pdst_i(:,2).* psrc_i(:,2));
    Tx = coef_TxTy \ b_Tx;
    Ty = coef_TxTy \ b_Ty;
    t1(i) = Tx(1);
    t2(i) = Tx(2);  
    t3(i) = Ty(1);
    t4(i) = Ty(2);
end
all_Tx = [t1 t2];
all_Ty = [t3 t4];

clear psrc_i pdst_i coef_TxTy b_Tx Tx Ty t1 t2 t3 t4

%% Compute warpped image(IDW)
for ii=1:h
    for jj=1:w
        p_minus_psrc = repmat([ii jj],n,1) - psrc;
        p_test = logical(sum(logical(p_minus_psrc),2));
        
        if sum(p_test) < n
            index_zero_row = find(p_test == 0);
            new_p_x = round(pdst(index_zero_row,1));
            new_p_y = round(pdst(index_zero_row,2));
            im3(new_p_x,new_p_y,:) = im2(ii,jj,:);
            continue
            
        else
            f_ix = sum(p_minus_psrc.* all_Tx,2) + pdst(:,1);
            f_iy = sum(p_minus_psrc.* all_Ty,2) + pdst(:,2);
            rep_p = repmat([ii jj],n,1);
            distance_psrc_p_square = sum((psrc - rep_p).^2,2);
            sigma_p = ones(n,1)./ distance_psrc_p_square;
            w_i = sigma_p / sum(sigma_p);
            new_p_x = round(sum(w_i.* f_ix));
            new_p_y = round(sum(w_i.* f_iy));
            
            if  new_p_x > h
                 new_p_x = h;
            end
            if new_p_x < 1
                new_p_x = 1;
            end
            if new_p_y > w
                new_p_y = w;
            end
            if new_p_y < 1
                new_p_y = 1;
            end
            im3(new_p_x,new_p_y,:) = im2(ii,jj,:);
            X(ii,jj) = new_p_x;
            Y(ii,jj) = new_p_y;
        
        end
    end
end

figure(2)
axis off
Z = zeros(h,w);
S = surface(X,Y,Z,im2,'FaceColor','texturemap','EdgeColor','none');
im2 = im3;

im2(:,:,1) = medfilt2(im3(:,:,1),[5,5]);
im2(:,:,2) = medfilt2(im3(:,:,2),[5,5]);
im2(:,:,3) = medfilt2(im3(:,:,3),[5,5]);
