%initialize distributed computing toolbox
%matlabpool


%%

%fill in the directories here

%social/novel, analyzed as seperate files.
dir_data_social_individual = ...
   {
   '/sohal1/sohal/Frost/inscopix/socialnovel2017-0915/2017-0915-1-homecage/recording_20170915_161115/dstif'
   '/sohal1/sohal/Frost/inscopix/socialnovel2017-0915/2017-0915-1-object/recording_20170915_162536/dstif'
   '/sohal1/sohal/Frost/inscopix/socialnovel2017-0915/2017-0915-1-social/recording_20170915_163825/dstif'
   '/sohal1/sohal/Frost/inscopix/socialnovel2017-0915/2017-0915-2-homecage/recording_20170915_174410/dstif'
   '/sohal1/sohal/Frost/inscopix/socialnovel2017-0915/2017-0915-2-object/recording_20170915_175608/dstif'
   '/sohal1/sohal/Frost/inscopix/socialnovel2017-0915/2017-0915-2-social/recording_20170915_180756/dstif'
   '/sohal1/sohal/Frost/inscopix/socialnovel2017-0915/2017-0915-3-homecage/recording_20170915_182545/dstif'
   '/sohal1/sohal/Frost/inscopix/socialnovel2017-0915/2017-0915-3-object/recording_20170915_183734/dstif'
   '/sohal1/sohal/Frost/inscopix/socialnovel2017-0915/2017-0915-3-social/recording_20170915_184921/dstif'
   '/sohal1/sohal/Frost/inscopix/socialnovel2017-0915/2017-0915-4-homecage/recording_20170915_192150/dstif'
   '/sohal1/sohal/Frost/inscopix/socialnovel2017-0915/2017-0915-4-object/recording_20170915_193400/dstif'
   '/sohal1/sohal/Frost/inscopix/socialnovel2017-0915/2017-0915-4-social/recording_20170915_195415/dstif'
   '/sohal1/sohal/Frost/inscopix/socialnovel2017-0915/2017-0915-5-homecageredo/recording_20170915_211359/dstif'
   '/sohal1/sohal/Frost/inscopix/socialnovel2017-0915/2017-0915-5-novelredo/recording_20170915_210216/dstif'
   '/sohal1/sohal/Frost/inscopix/socialnovel2017-0915/2017-0915-5-social/recording_20170915_204958/dstif'
   '/sohal1/sohal/Frost/inscopix/socialnovel2017-0915/2017-0915-6-homecage/recording_20170915_170004/dstif'
   '/sohal1/sohal/Frost/inscopix/socialnovel2017-0915/2017-0915-6-object/recording_20170915_171200/dstif'
   '/sohal1/sohal/Frost/inscopix/socialnovel2017-0915/2017-0915-6-social/recording_20170915_172430/dstif'

}  

%Social-Novel, standard protocol(10/10/10min) concatenated
dir_data = {
    '/sohal1/sohal/Frost/inscopix/socialnovel/mouse1H-O-S'
    '/sohal1/sohal/Frost/inscopix/socialnovel/mouse2H-O-S'
    '/sohal1/sohal/Frost/inscopix/socialnovel/mouse3H-O-S'
    '/sohal1/sohal/Frost/inscopix/socialnovel/mouse4H-O-S'
    '/sohal1/sohal/Frost/inscopix/socialnovel/mouse5H-O-S'
    '/sohal1/sohal/Frost/inscopix/socialnovel/mouse6H-O-S'
    '/sohal1/sohal/Frost/inscopix/socialnovel/mouse9H-O-S'
    '/sohal1/sohal/Frost/inscopix/socialnovel/mouse10H-O-S'
    '/sohal1/sohal/Frost/inscopix/socialnovel/mouse124H-O-S'
}

%Homecage-Openfield-EPM(15/20/20min)

dir_data = ...
{
 
'/sohal1/sohal/Frost/inscopix/Context-vs-Ensembles-Day0/mouse1-HC-OF-EPM'
'/sohal1/sohal/Frost/inscopix/Context-vs-Ensembles-Day0/mouse2-HC-OF-EPM'
'/sohal1/sohal/Frost/inscopix/Context-vs-Ensembles-Day0/mouse3-HC-OF-EPM'
'/sohal1/sohal/Frost/inscopix/Context-vs-Ensembles-Day0/mouse4-HC-OF-EPM'
'/sohal1/sohal/Frost/inscopix/Context-vs-Ensembles-Day0/mouse5-HC-OF-EPM'
'/sohal1/sohal/Frost/inscopix/Context-vs-Ensembles-Day0/mouse6-HC-OF-EPM'
'/sohal1/sohal/Frost/inscopix/Context-vs-Ensembles-Day0/mouse9-HC-OF-EPM'
'/sohal1/sohal/Frost/inscopix/Context-vs-Ensembles-Day0/mouse10-HC-OF-EPM'
'/sohal1/sohal/Frost/inscopix/Context-vs-Ensembles-Day0/mouse11-HC-OF-EPM'
'/sohal1/sohal/Frost/inscopix/Context-vs-Ensembles-Day0/mouse124-HC-OF-EPM'

}

%Day 1 home cage vs homecage + MMN (0.5 Hz)
dir_data = {...
'/sohal1/sohal/Frost/inscopix/Context-vs-Ensembles-Day1/mouse1-d0-HC-HCMMN'
'/sohal1/sohal/Frost/inscopix/Context-vs-Ensembles-Day1/mouse2-d0-HC-HCMMN'
'/sohal1/sohal/Frost/inscopix/Context-vs-Ensembles-Day1/mouse3-d0-HC-HCMMN'
'/sohal1/sohal/Frost/inscopix/Context-vs-Ensembles-Day1/mouse4-d0-HC-HCMMN'
'/sohal1/sohal/Frost/inscopix/Context-vs-Ensembles-Day1/mouse5-d0-HC-HCMMN'
'/sohal1/sohal/Frost/inscopix/Context-vs-Ensembles-Day1/mouse6-d0-HC-HCMMN'
'/sohal1/sohal/Frost/inscopix/Context-vs-Ensembles-Day1/mouse9-d0-HC-HCMMN'
'/sohal1/sohal/Frost/inscopix/Context-vs-Ensembles-Day1/mouse10-d0-HC-HCMMN'
'/sohal1/sohal/Frost/inscopix/Context-vs-Ensembles-Day1/mouse124-d0-HC-HCMMN'
}  

%Day 2: home cage vs context A vs Context B
dir_data = { ...
'/sohal1/sohal/Frost/inscopix/Context-vs-Ensembles-Day2/Mouse1-Home-A-B'
'/sohal1/sohal/Frost/inscopix/Context-vs-Ensembles-Day2/Mouse2-Home-A-B'
'/sohal1/sohal/Frost/inscopix/Context-vs-Ensembles-Day2/Mouse3-Home-A-B'
'/sohal1/sohal/Frost/inscopix/Context-vs-Ensembles-Day2/Mouse4-Home-A-B'
'/sohal1/sohal/Frost/inscopix/Context-vs-Ensembles-Day2/Mouse5-Home-A-B'
'/sohal1/sohal/Frost/inscopix/Context-vs-Ensembles-Day2/Mouse6-Home-A-B'
'/sohal1/sohal/Frost/inscopix/Context-vs-Ensembles-Day2/Mouse9-Home-A-B'
'/sohal1/sohal/Frost/inscopix/Context-vs-Ensembles-Day2/Mouse10-Home-A-B'
'/sohal1/sohal/Frost/inscopix/Context-vs-Ensembles-Day2/Mouse124-Home-A-B'
}

%day3: home cage vs context A vs Context C
dir_data = { ...
'/sohal1/sohal/Frost/inscopix/Context-vs-Ensembles-Day3/Mouse1-Home-A-C-nohomecage'
'/sohal1/sohal/Frost/inscopix/Context-vs-Ensembles-Day3/Mouse2-Home-A-C'
'/sohal1/sohal/Frost/inscopix/Context-vs-Ensembles-Day3/Mouse4-Home-A-C-Amissing'
'/sohal1/sohal/Frost/inscopix/Context-vs-Ensembles-Day3/Mouse5-Home-A-C'
'/sohal1/sohal/Frost/inscopix/Context-vs-Ensembles-Day3/Mouse6-Home-A-C'
'/sohal1/sohal/Frost/inscopix/Context-vs-Ensembles-Day3/Mouse9-Home-A-C'
'/sohal1/sohal/Frost/inscopix/Context-vs-Ensembles-Day3/Mouse10-Home-A-C'
'/sohal1/sohal/Frost/inscopix/Context-vs-Ensembles-Day3/Mouse124-Home-A-C'
}


%Day4 - Homecage vs Context B vs Context B + Juvenile
dir_data = {...
'/sohal1/sohal/Frost/inscopix/Context-vs-Ensembles-Day4/Mouse1-Home-B-BJuv'
'/sohal1/sohal/Frost/inscopix/Context-vs-Ensembles-Day4/Mouse2-Home-B-BJuv'
'/sohal1/sohal/Frost/inscopix/Context-vs-Ensembles-Day4/Mouse3-Home-B-BJuv'
'/sohal1/sohal/Frost/inscopix/Context-vs-Ensembles-Day4/Mouse4-Home-B-BJuv'
'/sohal1/sohal/Frost/inscopix/Context-vs-Ensembles-Day4/Mouse5-Home-B-BJuv'
'/sohal1/sohal/Frost/inscopix/Context-vs-Ensembles-Day4/Mouse6-Home-B-BJuv'
'/sohal1/sohal/Frost/inscopix/Context-vs-Ensembles-Day4/Mouse9-Home-B-BJuv'
'/sohal1/sohal/Frost/inscopix/Context-vs-Ensembles-Day4/Mouse10-Home-B-BJuv'
'/sohal1/sohal/Frost/inscopix/Context-vs-Ensembles-Day4/Mouse124-Home-B-BJuv'
                 
}

%Day 5 - Homecage vs context C vs COntext C


%Kara RS Dataset
dir_data = {...
'/sohal1/sohal/Frost/inscopix/2017-1109-mouse124RS/recording_20171109_154618/dstif'
}

%anthony files
dir_data = {...
'/sohal1/sohal/Frost/aleeinscopix/20170413_VIPhalomch_SG6_m297_OF/recording_20170413_144255/dstif'
'/sohal1/sohal/Frost/aleeinscopix/20170413_VIPhalomch_SG6_m325_OF/recording_20170413_141141/dstif'
'/sohal1/sohal/Frost/aleeinscopix/20170414_VIPhalomch_SG6_m297_RTEPM2/recording_20170414_145945/dstif'
'/sohal1/sohal/Frost/aleeinscopix/20170414_VIPhalomch_SG6_m325_RTEPM/recording_20170414_131009/dstif'
'/sohal1/sohal/Frost/aleeinscopix/20170418_VIPhalomch_SG6_m297_RTEPM3/recording_20170418_094524/dstif'
'/sohal1/sohal/Frost/aleeinscopix/20170602_VIPhalomch_SG6_m149_OF_Partial/recording_20170602_155516/dstif'
'/sohal1/sohal/Frost/aleeinscopix/20170603_VIPhalomch_SG6_m148_RTEPM/recording_20170603_162221/dstif'
'/sohal1/sohal/Frost/aleeinscopix/20170603_VIPhalomch_SG6_m149_RTEPM/recording_20170603_152929/dstif'
'/sohal1/sohal/Frost/aleeinscopix/20170603_VIPhalomch_SG6_m297_round2_RTEPM/recording_20170603_155457/dstif'
'/sohal1/sohal/Frost/aleeinscopix/201704014_VIPhalomch_SG6_m297_RTEPM/recording_20170414_124552/dstif'

'/sohal1/sohal/Frost/aleeinscopix/20160513_VIPhaloSG6m_RTEPM_m2049/trimmedFiles/singleframes/dstif'
'/sohal1/sohal/Frost/aleeinscopix/20160513_VIPhaloSG6m_RTEPM_m2225/trimmedFiles/singleframes/dstif'
'/sohal1/sohal/Frost/aleeinscopix/20160513_VIPhaloSG6m_RTEPM_m2383/trimmedFiles/singleframes/dstif'


}

%anthony controls
dir_data = {...
    '/sohal1/sohal/Frost/aleeinscopix/controls/20160117_m2083_VIPhalo_SynG6m_EPM_inscopix/unprocessedSingles/dstif'
    '/sohal1/sohal/Frost/aleeinscopix/controls/20170215_VIPeYFP_SG6m_RTEPM_m6011/recording_20170215_164630/dstif'
    '/sohal1/sohal/Frost/aleeinscopix/controls/20170215_VIPeYFP_SG6m_RTEPM_m6196/recording_20170215_162024/dstif'
    '/sohal1/sohal/Frost/aleeinscopix/controls/recording_20171109_114711_4320/dstif'
    '/sohal1/sohal/Frost/aleeinscopix/controls/recording_20171109_121359_4513/dstif'
    '/sohal1/sohal/Frost/aleeinscopix/controls/recording_20171109_130616_4321/dstif'
    '/sohal1/sohal/Frost/aleeinscopix/controls/recording_20171201_112318_4520EPM/dstif'
    '/sohal1/sohal/Frost/aleeinscopix/controls/recording_20171201_114910_4518EPM/dstif'
    '/sohal1/sohal/Frost/aleeinscopix/controls/recording_20171201_122826_4339EPM/dstif'
    '/sohal1/sohal/Frost/aleeinscopix/controls/recording_20171201_130003_4521EPM/dstif'
    '/sohal1/sohal/Frost/aleeinscopix/controls/recording_20171201_132807_4519EPM/dstif'
    '/sohal1/sohal/Frost/aleeinscopix/controls/recording_20171201_135630_4517EPM/dstif'
    '/sohal1/sohal/Frost/aleeinscopix/controls/recording_20171201_145209_EGFPEPM/dstif'
}

%Mismatch + SSA 0.5 Hz
dir_data = {...
    '/sohal1/sohal/Frost/inscopix/mismatch vs SSA/2017-0902-2-MMN-SSA'
    '/sohal1/sohal/Frost/inscopix/mismatch vs SSA/2017-0902-4-MMN-SSA'    
    '/sohal1/sohal/Frost/inscopix/mismatch vs SSA/2017-0902-5-MMN-SSA'
    '/sohal1/sohal/Frost/inscopix/mismatch vs SSA/2017-0902-6-MMN-SSA'
    '/sohal1/sohal/Frost/inscopix/mismatch vs SSA/2017-1027-124-MMN-SSA'

}

%Social replay
dir_data = {...
    '/sohal1/sohal/Frost/inscopix/socialreplay/2017-1026-2-socialreplay'
    '/sohal1/sohal/Frost/inscopix/socialreplay/2017-1026-4-socialreplay'
    '/sohal1/sohal/Frost/inscopix/socialreplay/2017-1026-5-socialreplay'
    '/sohal1/sohal/Frost/inscopix/socialreplay/2017-1026-6-socialreplay'
    '/sohal1/sohal/Frost/inscopix/socialreplay/2017-1026-124-Socialpreplay'
    '/sohal1/sohal/Frost/inscopix/socialreplay/2018-0104-4513-Socialpreplay'
    '/sohal1/sohal/Frost/inscopix/socialreplay/2018-0104-4517-Socialpreplay'
    '/sohal1/sohal/Frost/inscopix/socialreplay/2018-0104-4518-Socialpreplay'
    '/sohal1/sohal/Frost/inscopix/socialreplay/2018-0104-4520-Socialpreplay'    
}


%BLA ANXIETY
dir_data = {...
    '/sohal1/sohal/Frost/inscopix/BLA ANXIETY/2018-0222-4451-HC_OF_EPM'
    '/sohal1/sohal/Frost/inscopix/BLA ANXIETY/2018-0222-4452-HC_OF_EPM'
    '/sohal1/sohal/Frost/inscopix/BLA ANXIETY/2018-0222-4453-HC_OF_EPM'
    '/sohal1/sohal/Frost/inscopix/BLA ANXIETY/2018-0222-4454-HC_OF_EPM'  
}


%social defeat days 1 and 2
dir_data = {...
    '/sohal1/sohal/Frost/inscopix/Social Defeat/social defeat day 1/2018-0130-124-Day1'
    '/sohal1/sohal/Frost/inscopix/Social Defeat/social defeat day 1/2018-0130-4513-Day1'
    '/sohal1/sohal/Frost/inscopix/Social Defeat/social defeat day 1/2018-0130-4517-Day1'
    '/sohal1/sohal/Frost/inscopix/Social Defeat/social defeat day 1/2018-0130-4518-Day1'
    '/sohal1/sohal/Frost/inscopix/Social Defeat/social defeat day 2/2018-0131-124-Day2'
    '/sohal1/sohal/Frost/inscopix/Social Defeat/social defeat day 2/2018-0131-4513-Day2'
    '/sohal1/sohal/Frost/inscopix/Social Defeat/social defeat day 2/2018-0131-4517-Day2'
    '/sohal1/sohal/Frost/inscopix/Social Defeat/social defeat day 2/2018-0131-4518-Day2'
   
}

%second round of HC-OF-EPM
dir_data = {...
    '/sohal1/sohal/Frost/inscopix/Homecage-OF-EPM/Segmented/2018-0625-1764full'
    '/sohal1/sohal/Frost/inscopix/Homecage-OF-EPM/Segmented/2018-0624-1764a-hc+of'
    '/sohal1/sohal/Frost/inscopix/Homecage-OF-EPM/Segmented/2018-0624-1833'
    '/sohal1/sohal/Frost/inscopix/Homecage-OF-EPM/Segmented/2018-0624-1962'
    '/sohal1/sohal/Frost/inscopix/Homecage-OF-EPM/Segmented/2018-0625-1768'
    '/sohal1/sohal/Frost/inscopix/Homecage-OF-EPM/Segmented/2018-0824-4546'
    '/sohal1/sohal/Frost/inscopix/Homecage-OF-EPM/Segmented/2018-1115-1896wsync'
    '/sohal1/sohal/Frost/inscopix/Homecage-OF-EPM/Segmented/2018-1115-1900wsync'
    '/sohal1/sohal/Frost/inscopix/Homecage-OF-EPM/Segmented/2018-1024-vpa2'
    '/sohal1/sohal/Frost/inscopix/Homecage-OF-EPM/Segmented/2018-1024-vpa4-secondOFbout'
    '/sohal1/sohal/Frost/inscopix/Homecage-OF-EPM/Segmented/2018-1024-vpa5'
    '/sohal1/sohal/Frost/inscopix/Homecage-OF-EPM/Segmented/2018-1024-vpa7'
    '/sohal1/sohal/Frost/inscopix/Homecage-OF-EPM/Segmented/2018-1028-vpa17'
    '/sohal1/sohal/Frost/inscopix/Homecage-OF-EPM/Segmented/2018-1028-vpa18'
    }

%second round of HC-O-S
dir_data = {...
    '/sohal1/sohal/Frost/inscopix/socialnovel/segmented/2018-1101-1764'
    '/sohal1/sohal/Frost/inscopix/socialnovel/segmented/2018-1101-1833'
    '/sohal1/sohal/Frost/inscopix/socialnovel/segmented/2018-1101-4546'
    '/sohal1/sohal/Frost/inscopix/socialnovel/segmented/2018-1101-vpa2'
    '/sohal1/sohal/Frost/inscopix/socialnovel/segmented/2018-1101-vpa5'
    '/sohal1/sohal/Frost/inscopix/socialnovel/segmented/2018-1101-vpa7'
    '/sohal1/sohal/Frost/inscopix/socialnovel/segmented/2018-1102-vpa17juvrepeated'
    '/sohal1/sohal/Frost/inscopix/socialnovel/segmented/2018-1102-vpa18'
    
}



%second round of social replay
dir_data = {...
    '/sohal1/sohal/Frost/inscopix/socialreplay/segmented/2018-0718-1764-socialreplay'
    '/sohal1/sohal/Frost/inscopix/socialreplay/segmented/2018-0718-1833-socialreplay'
    '/sohal1/sohal/Frost/inscopix/socialreplay/segmented/2018-0718-1962-socialreplay'
    '/sohal1/sohal/Frost/inscopix/socialreplay/segmented/2018-0828-4546-socialreplay'
    '/sohal1/sohal/Frost/inscopix/socialreplay/segmented/2018-1116-1896-socialreplay'
    '/sohal1/sohal/Frost/inscopix/socialreplay/segmented/2018-1116-1900-socialreplay'
    '/sohal1/sohal/Frost/inscopix/socialreplay/segmented/2018-1030-vpa02-socialreplay'
    '/sohal1/sohal/Frost/inscopix/socialreplay/segmented/2018-1030-vpa02-socialreplay-usingfullHC2'
    '/sohal1/sohal/Frost/inscopix/socialreplay/segmented/2018-1030-vpa04-socialreplay'
    '/sohal1/sohal/Frost/inscopix/socialreplay/segmented/2018-1030-vpa05-socialreplay'
    '/sohal1/sohal/Frost/inscopix/socialreplay/segmented/2018-1030-vpa07-socialreplay'
    '/sohal1/sohal/Frost/inscopix/socialreplay/segmented/2018-1030-vpa17-socialreplay'
    '/sohal1/sohal/Frost/inscopix/socialreplay/segmented/2018-1030-vpa18-socialreplay'  
    
 }   
    

%%   

fn_data = sprintf('20190816-mm.mat',date);
save(['/sohal1/sohal/Frost/inscopixanalysis/' fn_data],'dir_data')

%%

for kk =1:size(dir_data,1)

        kk
        tempdir = dir_data{kk};
        cd(tempdir);
            try load Numbrframes; end
            try load movies; end

        m = dir;
        mm = dir('*.tif');
        for i = 1:10
            if strfind(m(i).name, '.tif')
                fn = m(i).name
                break
            else
                i = i+1;
            end
        end
        
    
        framerate = 20; 
        if size(mm) > 75000; flims = [1 75000]; else flims = []; end
        nPCs = 125;
        %consider increase to 130 from 100 to see if we get more PC/ICs
        
        dsamp = [1 1]; 
        badframes = []; 
        outputdir = []; 
        % Note dsamp vector first element is the time vector down sample whereas
        % the second element is the spatial downsample
        
        edgecrop = 0; %will insert border of pixel value one around edge during PCA
        forceload = 1; %if 1 will force reperforming PCA;otherwise will simply load saved data
        %will only use movie before and after treatment for all cases
        %dsamp - first term is temporal and second term is spatial
        %downsampling factor; given that pixels are already binned (4x4) during
        %acquisition will use no spatial downsampling.
        %[mixedsig, mixedfilters, CovEvals, covtrace, movm, movtm] = CellsortPCA_singleframes(fn, flims, nPCs, dsamp, outputdir ,badframes,forceload);
        [mixedsig, mixedfilters, CovEvals, covtrace, movm, movtm] = CellsortPCA_singleframes_applyff0edgeb(fn, flims, nPCs, dsamp, outputdir ,badframes,forceload,edgecrop);

        if edgecrop > 1
            mixedfilters_cropped = mixedfilters;
            mixedfilters = zeros(size(movm,1)+2*edgecrop-1, size(movm,2)+2*edgecrop-1, size(mixedfilters,3));
            
            movm_cropped = movm;
            movm = zeros(size(mixedfilters,1), size(mixedfilters,2));    
            movm(edgecrop:end-edgecrop,edgecrop:end-edgecrop) = movm_cropped;
            for nm = 1:size(mixedfilters,3)
                mixedfilters(edgecrop:end-edgecrop,edgecrop:end-edgecrop, nm) = mixedfilters_cropped(:,:,nm);
            end
        end
        
        PCuse  = []; mu = 0.15; nIC = 125; ica_A_guess = []; termtol = [10^-6];  maxrounds = [1000];
        [ica_sig, ica_filters, ica_A, numiter] = CellsortICA(mixedsig, mixedfilters, CovEvals, PCuse, mu, nIC, ica_A_guess, termtol, maxrounds);
        
        %pure spatial ICA (mu = 0) always superior to pure temporal ICA (mu
        %= 1). spatialtemporal ica with mu < 0.5 usually superior
        %'With real data we found by exploration that spatio-temporal ICA (� ? 0.1�0.2) extracted the most components resembling Purkinje cells, so we
        %habitually used this approach(mukamel et Schnitzer)'
        
        
        smwidth = 1.5;thresh = 7.5;arealims = [10 100];plotting = 0; %nearly all ICs < 250 (ie ~ 15x15)
        [ica_segments_all, segmentlabel_all, segcentroid_all, seg_majoraxis_all] = CellsortSegmentation(ica_filters, smwidth, thresh, arealims, plotting);

        max_length =100; [ica_segments, seg_majoraxis, segcentroid] = majoraxisfilter(ica_segments_all, seg_majoraxis_all, max_length, segcentroid_all); %francisco used 40
        edge_distance = 10;  [ ica_segments, segcentroid, exclude_edge ] = edgefilter( ica_segments, segcentroid, edge_distance, movm );
        
        min_distance = 6;[ica_segments, segcentroid] = distancefilter( ica_segments, segcentroid, min_distance );

        %note:on olympus set up (slice rig) 7 pixels corresponds to ~2.5*7 um = 16.5 um; 8 pixels is ~ 20
        %um; there is substantial crosstalk at < 7 pix when cells are
        %dense/slice is active
        
        %similarly there is an increase in correlations for cells located
        %less than ~ 12 pixels in inscopix data using PCA/ICA method
        %looking at either subtracted fluorescence or raster analysis
        
        %however simply removing all cells spaced < 30 um appears too
        %conservative as often results in near 50% of cells being removed
        %January 2019 edit to remove cells spaced less than 6 pixels (15
        %um) or cells < 12 pixels(30 um that have correlation coefficient
        %(from raster) > 0.3.
        
        %backgroundmask - will use for surround subtraction
        backgroundmask = zeros(size(ica_segments,2),size(ica_segments,3));
        for i = 1:size(ica_segments_all,1); 
        backgroundmask = backgroundmask + squeeze(ica_segments_all(i,:,:) > 1); 
        end

        %subtractmin = 0;flims = []; cell_sig = CellsortApplyFilter_singleframes_subtractmin0(fn, ica_segments, flims, movm, subtractmin);
       	subtractsurround = 1;flims = []; [cell_sig,cell_sig_subtracted,cell_sig_surround] = CellsortApplyFilter_singleframes_subtractsurroundmask(fn, ica_segments, flims, movm, subtractsurround,backgroundmask);

            %for some reason keep getting nans for cells that are close to
            %edge and near other cells
            killsig = find(isnan(cell_sig_subtracted(:,1)))
                ica_segments(killsig,:) = [];
                segcentroid(killsig,:) = [];
                seg_majoraxis(:,killsig) = [];
                cell_sig(killsig,:) = [];
                cell_sig_subtracted(killsig,:) = [];
                cell_sig_surround(killsig,:) = [];

        %using subtract surround will subtract median value from pixels
        %bordering each IC to remove non-cell intrinsic signals. Somewhat
        %computationally intensive but useful if activity rather dense or
        %if there are fluctuations in background fluorescence
        
        % **********Note(**************: Remember to initialize workers if using the
        % distributed computing toolbox (i.e. parfor), just type 'matlabpool' into
        % the command prompt 

        %Code to convert edge subtracted signal to F/F0 and low pass filter
        fff = designfilt('lowpassfir', 'PassbandFrequency', 0.5, 'StopbandFrequency', .65, 'PassbandRipple', 1, 'StopbandAttenuation', 25);
        sig = [];
        for i = 1:size(cell_sig_subtracted,1)
        sig(i,:) = (cell_sig_subtracted(i,:)-prctile(cell_sig_subtracted(i,:)', .5))/mean(cell_sig_surround(i,:)');
        sig(i,:) =filtfilt(fff, sig(i,:));
        end

        temp_ppc = csaps([1:size(sig,2)],sig,0.05);
        C = fnval(temp_ppc,[1:length(sig)]);             

        thresh = 4;thresh2 = 15; thresh3 =250; thresh_abs = 0.0125; spikewin = 40;mindur = 1; baseline_var = 0.5;
        spiketimes = [];
        cellid = [];
        spt_beg = [];
        spt_end = [];
        s = [];
        b = [];
        e = [];
        new_spiketimes = [];
        new_cellid = [];

                    parfor i = 1:size(sig,1);    
                        [s b e] = detectspikes12(C(i,:)', thresh, thresh2,thresh3, thresh_abs, spikewin, mindur, baseline_var);
                        spiketimes = [spiketimes s];
                        spt_beg = [spt_beg b];
                        spt_end = [spt_end e];
                        cellid = [cellid, i*ones(1,length(s))];
                        i
                    end

    for i = 1:numel(spiketimes)
        new_spiketimes = [spt_beg(i):spiketimes(i) new_spiketimes];
        new_cellid = [cellid(i)*ones(1,spiketimes(i) - spt_beg(i) + 1) new_cellid];
        i;
    end

% Code for saving data file
temp = pwd; temp = strrep(temp, ' ', '_'); temp = strrep(temp, '/', '_'); 
raster = sparse(new_cellid, new_spiketimes, 1, size(sig,1),size(sig,2));
raster = full(raster > 0); 

%remove closely spaced cells with high correlation coefficients
jnk = clean_raster(segcentroid, raster,C, 12, 0.25);
if isempty (jnk) == 0
                cell_sig(jnk,:) = [];
                ica_segments(jnk,:) = [];
                cell_sig_subtracted(jnk,:) = [];
                cell_sig_surround(jnk,:) = [];
                segcentroid(jnk,:) = [];
                sig(jnk,:) = [];
                C(jnk,:) = [];
                
                %Remove spikes found in duplicate cells
                tmp = []
                for i = 1:size(jnk,1)
                    a = find(new_cellid(:) == jnk(i));
                    tmp = [tmp;a];
                end
                new_cellid(tmp) = [];
                new_spiketimes(tmp) = [];
                    
                    jnk = sort(jnk,'descend');
                    for j = 1: size(jnk,1)
                        [a ~] = find(new_cellid(:) >= jnk(j));
                        new_cellid(a) = new_cellid(a)-1;
                    end
end            

                
               
inscopix{kk} = raster;

final_raster = full(raster); final_raster(manual_badcells,:) = [];
final_C = C; final_C(manual_badcells,:) = [];


save('results_20190806.mat')
save(['/sohal1/sohal/Frost/inscopixanalysis/' fn_data], 'inscopix' , '-append')

clearvars -except dir_data inscopix fn_data

end

