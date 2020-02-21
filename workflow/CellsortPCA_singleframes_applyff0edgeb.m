function [mixedsig, mixedfilters, CovEvals, covtrace, movm, ...
    movtm] = CellsortPCA_singleframes_applyff0edge(fn, flims, nPCs, dsamp, outputdir, badframes,forceload,edgecrop)

% [mixedsig, mixedfilters, CovEvals, covtrace, movm, movtm] = CellsortPCA(fn, flims, nPCs, dsamp, outputdir, badframes)
%
% CELLSORT
% Read TIFF movie data and perform singular-value decomposition (SVD)
% dimensional reduction.
% 
%   NOTE: This is the nimmerjahn code except modified to work with single
%   tiff frames as opposed to tiff stacks, fn in this case corresponds to
%   the name of the first frame i.e. img_00000_00.tif
%
% Inputs:
%   fn - movie file name. Must be in TIFF format.
%   flims - 2-element vector specifying the endpoints of the range of
%   frames to be analyzed. If empty, default is to analyze all movie
%   frames.
%   nPCs - number of principal components to be returned
%   dsamp - optional downsampling factor. If scalar, specifies temporal
%   downsampling factor. If two-element vector, entries specify temporal
%   and spatial downsampling, respectively.
%   outputdir - directory in which to store output .mat files
%   badframes - optional list of indices of movie frames to be excluded
%   from analysis
%
% Outputs:
%   mixedsig - N x T matrix of N temporal signal mixtures sampled at T
%   points.
%   mixedfilters - N x X x Y array of N spatial signal mixtures sampled at
%   X x Y spatial points.
%   CovEvals - largest eigenvalues of the covariance matrix
%   covtrace - trace of covariance matrix, corresponding to the sum of all
%   eigenvalues (not just the largest few)
%   movm - average of all movie time frames at each pixel
%   movtm - average of all movie pixels at each time frame, after
%   normalizing each pixel deltaF/F
%
% Eran Mukamel, Axel Nimmerjahn and Mark Schnitzer, 2009
% Email: eran@post.harvard.edu, mschnitz@stanford.edu
%

tic
fprintf('-------------- CellsortPCA %s: %s -------------- \n', date, fn)

m = dir;

%-----------------------

% Note, I am calling this function regardless because then if I define a
% flims outside of function, I can just refer to the vals/dir_use vectors
% as indices
%[dir_use frame_nt nt_full] = tiff_frames(fn);

dir_use = dir('*.tif');
nt_full = length(dir_use);

    if (nargin<2)||(isempty(flims))
        flims = [1,nt_full];
    end


if nargin<4 || isempty(dsamp)
    dsamp = [1,1]; 
end
    
if size(dsamp,2) == 1
        dsamp_time = dsamp(1);
        dsamp_space = 1;
        
elseif size(dsamp,2) == 2
        dsamp_time = dsamp(1);
        dsamp_space = dsamp(2);        
end

useframes = setdiff((flims(1):flims(2)), badframes);
nt = length(flims(1):flims(2));

if nargin<3 || isempty(nPCs)
    nPCs = min(150, nt);
end

if nargin<5 || isempty(outputdir)
    outputdir = [pwd,'/cellsort_preprocessed_data/'];
end
if nargin<6
    badframes = [];
end
if isempty(dir(outputdir))
    mkdir(pwd, '/cellsort_preprocessed_data/')
end
if outputdir(end)~='/';
    outputdir = [outputdir, '/'];
end



[fpath, fname] = fileparts(fn);
    if isempty(badframes)
        fnmat = [outputdir, fname, '_',num2str(flims(1)),',',num2str(flims(2)), '_', date,'.mat'];
    else
        fnmat = [outputdir, fname, '_',num2str(flims(1)),',',num2str(flims(2)),'_selframes_', date,'.mat'];
    end
    
    if ~isempty(dir(fnmat))
        fprintf('CELLSORT: Movie %s already processed;', ...
            fn)
        %forceload = input(' Re-load data? [0-no/1-yes] ');
        % make choice when calling function.
        if isempty(forceload) || forceload==0
            load(fnmat)
            return
        end
    end

fncovmat = [outputdir, fname, '_cov_', num2str(flims(1)), ',', num2str(flims(2)), '_', date,'.mat'];

[pixw,pixh] = size(imread(fn,1));
npix = pixw*pixh;

fprintf('   %d pixels x %d time frames;', npix, nt*dsamp_time)
    if nt<npix
        fprintf(' using temporal covariance matrix.\n')
    else
        fprintf(' using spatial covariance matrix.\n')
    end


%% note amending th code to always use tcov


% Create covariance matrix
    if nt < npix
        [covmat, mov, movm, movtm] = create_tcov(fn, pixw, pixh, useframes, nt, dsamp,edgecrop);
    else
    %Note: This should be xcov and not tcov    
        [covmat, mov, movm, movtm] = create_tcov(fn, pixw, pixh, useframes, nt, dsamp,edgecrop);
    end
    
if edgecrop == 0
    covtrace = trace(covmat) / npix;
    movm = reshape(movm, pixw/dsamp_space, pixh/dsamp_space);
    %pass to SVD if edgecrop = 0
    [pixw_edge, pixh_edge] = size(movm);
    npix_edge = pixw_edge * pixh_edge;
else
    tmp = imread(fn,1);
    tmp = tmp(edgecrop:end-edgecrop,edgecrop:end-edgecrop);
    tmp = imresize(tmp, 1/dsamp_space, 'bilinear' );
    
    [pixw_edge, pixh_edge] = size(tmp);
    npix_edge = pixw_edge * pixh_edge;  
    covtrace = trace(covmat) / npix;
    movm = reshape(movm, pixw_edge, pixh_edge);
end

% if nt < npix

    % Perform SVD on temporal covariance
    [mixedsig, CovEvals, percentvar] = cellsort_svd(covmat, nPCs, nt, npix_edge);

    % Load the other set of principal components
	[mixedfilters] = reload_moviedata(npix_edge, mov, mixedsig, CovEvals);
      mixedfilters = reshape(mixedfilters, pixw_edge,pixh_edge,nPCs);
           

           if (dsamp_space>1) 
            mixedfilters = imresize(mixedfilters,dsamp_space);
            movm = imresize(movm,dsamp_space);
           end


    firstframe_full = imread(fn,1);
    firstframe = firstframe_full;
    if dsamp_space>1
        firstframe = imresize(firstframe, size(mov(:,:,1)),'bilinear');
    end

%------------
% Save the output data
save(fnmat,'mixedfilters','CovEvals','mixedsig', ...
    'movm','movtm','covtrace', 'fn', 'flims', 'nPCs', 'dsamp', 'outputdir', 'badframes')
fprintf(' CellsortPCA: saving data and exiting; ')
toc


function [covmat, mov, movm, movtm] = create_tcov(fn, pixw, pixh, useframes, nt, dsamp,edgecrop)
        %-----------------------
       % Load movie data to compute the temporal covariance matrix

        % Downsampling

            dsamp_time = dsamp(1);
            dsamp_space = dsamp(2); % Spatial downsample
        
        %NO SPATIAL DOWNSAMPLING
        if (dsamp_space==1)
            if edgecrop == 0 %No edgecroppping
                [pixw_dsamp,pixh_dsamp] = size(imresize( imread(fn,1), 1/dsamp_space, 'bilinear' ));
                npix = pixw_dsamp * pixh_dsamp;
                mov = ones(pixw_dsamp, pixh_dsamp, nt);
            else
            %changing size of mov for edge removal
                tmp = imread(fn,1);
                tmp = tmp(edgecrop:end-edgecrop,edgecrop:end-edgecrop);
                tmp = imresize(tmp, 1/dsamp_space, 'bilinear' );
            
                [pixw_dsamp pixh_dsamp] = size(tmp);
                npix = pixw_dsamp * pixh_dsamp;
                mov = ones(pixw_dsamp, pixh_dsamp, nt);
            end
        
            kk = 1;
            for jjind=flims(1):flims(2)
               temp_fname = dir_use(jjind).name;

               if edgecrop == 0
                mov(:,:,kk) = imresize(imread(temp_fname), 1/dsamp_space, 'bilinear' );
               else
                tmp = imread(temp_fname,1);
                tmp = tmp(edgecrop:end-edgecrop,edgecrop:end-edgecrop);
                tmp = imresize(tmp, 1/dsamp_space, 'bilinear' );                  
  
                mov(:,:,kk) = tmp;
               end  
                   
                kk = kk+1;
                if mod(jjind,500)==1
                    fprintf(' Read frame %4.0f out of %4.0f; ', jjind, nt*dsamp_time)
                    toc
                end

            end
         
         
            
        else
        %YES SPATIAL DOWNSAMPLING
            if edgecrop == 0
                [pixw_dsamp,pixh_dsamp] = size(imresize( imread(fn,1), 1/dsamp_space, 'bilinear' ));
                npix = pixw_dsamp * pixh_dsamp;
                mov = ones(pixw_dsamp, pixh_dsamp, nt);
                
                for jjind=1:length(useframes)
                 jj = useframes(jjind);
                 temp_fname = dir_use(jj).name;
                 mov(:,:,jjind) = imresize( imread(temp_fname), 1/dsamp_space, 'bilinear' );
                    if mod(jjind,500)==1
                     fprintf(' Read frame %4.0f out of %4.0f; ', jjind, nt*dsamp_time)
                     toc
                    end
                end
            else
            %changing size of mov for edge removal
                
                tmp = imread(fn,1);
                tmp = tmp(edgecrop:end-edgecrop,edgecrop:end-edgecrop);
                tmp = imresize(tmp, 1/dsamp_space, 'bilinear' );

                [pixw_dsamp, pixh_dsamp] = size(tmp);
                npix = pixw_dsamp * pixh_dsamp;
                mov = ones(pixw_dsamp, pixh_dsamp, nt);
            
                for jjind=1:length(useframes)
                 jj = useframes(jjind);
                 temp_fname = dir_use(jj).name;
                	tmp = imread(temp_fname,1);
                    tmp = tmp(edgecrop:end-edgecrop,edgecrop:end-edgecrop);
                    tmp = imresize(tmp, 1/dsamp_space, 'bilinear' );
                    mov(:,:,jjind) = tmp;
                    
                    if mod(jjind,500)==1
                     fprintf(' Read frame %4.0f out of %4.0f; ', jjind, nt*dsamp_time)
                     toc
                    end
                end            
            end            
         end

        fprintf(' Read frame %4.0f out of %4.0f; ', jjind, nt*dsamp_time)
        toc
        mov = reshape(mov, npix, nt);
      
         % DFoF normalization of each pixel
         movm = mean(mov,2); % Average over time
         movmzero = (movm==0); % Avoid dividing by zero
         movm(movmzero) = 1;
         mov = mov ./ (movm * ones(1,nt)) - 1;
         mov(movmzero, :) = 0;

         if dsamp_time>1
             mov = filter(ones(dsamp_time,1)/dsamp_time, 1, mov, [], 2);
             mov = downsample(mov', dsamp_time)';
         end


        c1 = (mov'*mov)/npix;
        movtm = mean(mov,1); % Average over space

        covmat = c1 - movtm'*movtm;
        clear c1
end

%%
    function [mixedsig, CovEvals, percentvar] = cellsort_svd(covmat, nPCs, nt, npix)
        %-----------------------
        % Perform SVD

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
      
   %   mixedsig = interp(mixedsig,dsamp_time^2);
      
        CovEvals = CovEvals / npix;

        percentvar = 100*sum(CovEvals)/covtrace;
        fprintf([' First ',num2str(nPCs),' PCs contain ',num2str(percentvar,3),'%% of the variance.\n'])
    end

    function [mixedfilters] = reload_moviedata(npix, mov, mixedsig, CovEvals)
        %-----------------------
        % Re-load movie data
        nPCs = size(mixedsig,1);

        Sinv = inv(diag(CovEvals.^(1/2)));
        
        movuse = mov - ones(npix,1) * movtm;
         
        mixedfilters = reshape(movuse * mixedsig' * Sinv, npix, nPCs);

        
    end

 end