function [mixedsig, mixedfilters, CovEvals, covtrace, movm, ...
    movtm] = Cellsort3dPCAmat(mov, nPCs)
% Cellsort 3d‚ðtiff‚Å‚È‚­”z—ñ(XYZT)‚É‘Î‚µ‚ÄŠ|‚¯‚é
% [mixedsig, mixedfilters, CovEvals, covtrace, movm, movtm] = CellsortPCA(fn, flims, nPCs, dsamp, outputdir, badframes)
%
% CELLSORT
% Read TIFF movie data and perform singular-value decomposition (SVD)
% dimensional reduction.
%
% Inputs:
%   fn - movie file name. Must be in TIFF format.
%   nPCs - number of principal components to be returned
%   pixz - number of z slices
%   outputdir - directory in which to store output .mat files
%
% Outputs:
%   mixedsig - N x T matrix of N temporal signal mixtures sampled at T
%   points.
%   mixedfilters - N x X x Y x Z array of N spatial signal mixtures sampled at
%   X x Y x Z spatial points.
%   CovEvals - largest eigenvalues of the covariance matrix
%   covtrace - trace of covariance matrix, corresponding to the sum of all
%   eigenvalues (not just the largest few)
%   movm - average of all movie time frames at each pixel
%   movtm - average of all movie pixels at each time frame, after
%   normalizing each pixel deltaF/F

tic
fprintf('run CellsortPCA')

%-----------------------
[pixw, pixh, pixz, nt] = size(mov);
nf = nt*pixz;
npix = pixw*pixh*pixz;% 3d image

fprintf('   %d pixels x %d time frames;', npix, nt);

% Create covariance matrix
if nt < npix
    fprintf(' using temporal covariance matrix.\n')
    [covmat, mov, movm, movtm] = create_tcov(mov, pixw, pixh, pixz, nt);
else
    fprintf(' using spatial covariance matrix.\n')
    [covmat, mov, movm, movtm] = create_xcov(mov, pixw, pixh, pixz, nt);
end

covtrace = trace(covmat) / npix;
movm = reshape(movm, pixw, pixh, pixz);

if nt < npix
    % Perform SVD on temporal covariance
    [mixedsig, CovEvals] = cellsort_svd(covmat, nPCs, nt, npix);
    
    % Load the other set of principal components
    [mixedfilters] = reload_moviedata(pixw*pixh*pixz, mov, mixedsig, CovEvals);
else
    % Perform SVD on spatial components
    [mixedfilters, CovEvals] = cellsort_svd(covmat, nPCs, nt, npix);
    mixedfilters = mixedfilters' * npix;
    
    % Load the other set of principal components
    [mixedsig] = reload_moviedata(nt, mov', mixedfilters', CovEvals);
    mixedsig = mixedsig' / npix^2;
end

mixedfilters = reshape(mixedfilters, pixw, pixh, pixz, nPCs);

%------------
% Save the output data
toc

    function [covmat, mov, movm, movtm] = create_xcov(mov,pixw,pixh,pixz,nt)
        %-----------------------
        % Load movie data to compute the spatial covariance matrix
        
        npix1 = pixw*pixh*pixz; % 3d
        
        nframes = nt*pixz;
        
        mov = reshape(mov, npix1, nt); % allpixels x time
        
        % DFoF normalization of each pixel
        movm = mean(mov,2); % Average over time
        movmzero = (movm==0);
        movm(movmzero) = 1;
        mov = mov ./ (movm * ones(1,nt)) - 1; % Compute Delta F/F
        mov(movmzero, :) = 0;
        
        movtm = mean(mov,1); % Average over space
        mov = mov - ones(size(mov,1),1)*movtm;
        
        covmat = (mov*mov')/size(mov,2);
        covmat = covmat * size(mov,2)/size(mov,1); % Rescale to gree with create_tcov
        toc
    end

    function [covmat, mov, movm, movtm] = create_tcov(mov,pixw,pixh,pixz,nt)
        %-----------------------
        % Load movie data to compute the temporal covariance matrix
        npix1 = pixw*pixh*pixz;
        
        nframes = nt*pixz;
        mov = reshape(mov, npix1, nt);
        
        % DFoF normalization of each pixel
        movm = mean(mov,2); % Average over time
        movmzero = (movm==0); % Avoid dividing by zero
        movm(movmzero) = 1;
        mov = mov ./ (movm * ones(1,nt)) - 1;
        mov(movmzero, :) = 0;
        
        c1 = (mov'*mov)/npix1;
        movtm = mean(mov,1); % Average over space
        covmat = c1 - movtm'*movtm;
        clear c1
    end

    function [mixedsig, CovEvals, percentvar] = cellsort_svd(covmat, nPCs, nt, npix1)
        %-----------------------
        % Perform SVD
        
        covtrace1 = trace(covmat) / npix1;
        
        opts.disp = 0;
        opts.issym = 'true';
        if nPCs<size(covmat,1)
            [mixedsig, CovEvals] = eigs(covmat, nPCs, 'LM', opts);  % pca_mixedsig are the temporal signals, mixedsig
        else
            [mixedsig, CovEvals] = eig(covmat);
            CovEvals = diag( sort(diag(CovEvals), 'descend'));
            nPCs = size(CovEvals,1);
        end
        CovEvals = diag(CovEvals);
        if nnz(CovEvals<=0)
            nPCs = nPCs - nnz(CovEvals<=0);
            fprintf(['Throwing out ',num2str(nnz(CovEvals<0)),' negative eigenvalues; new # of PCs = ',num2str(nPCs),'. \n']);
            mixedsig = mixedsig(:,CovEvals>0);
            CovEvals = CovEvals(CovEvals>0);
        end
        
        mixedsig = mixedsig' * nt;
        CovEvals = CovEvals / npix1;
        
        percentvar = 100*sum(CovEvals)/covtrace1;
        fprintf([' First ',num2str(nPCs),' PCs contain ',num2str(percentvar,3),'%% of the variance.\n'])
    end

    function [mixedfilters] = reload_moviedata(npix1, mov, mixedsig, CovEvals)
        %-----------------------
        % Re-load movie data
        nPCs1 = size(mixedsig,1);
        
        Sinv = inv(diag(CovEvals.^(1/2)));

        movtm1 = mean(mov,1); % Average over space
        movuse = mov - ones(npix1,1) * movtm1;
        mixedfilters = reshape(movuse * mixedsig' * Sinv, npix1, nPCs1);
    end

    function j = tiff_frames(fn)
        %
        % n = tiff_frames(filename)
        %
        % Returns the number of slices in a TIFF stack.
        %
        % Modified April 9, 2013 for compatibility with MATLAB 2012b
        
        j = length(imfinfo(fn));            
    end
end