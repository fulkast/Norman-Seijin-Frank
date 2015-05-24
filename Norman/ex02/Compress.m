function I_comp = Compress(I, method)

% Default settings that will be used by the submission machine.
if (nargin<2)
    method = 'byChannelSVD';
end

switch method
    case 'allChannels'
        % Results on submission system: 
        % Run-time:             1.166679s
        % Mean squared error:   0.0084671
        % Compression rate:     0.18311
        disp('Method used: CompressAllChannels');
        I_comp = CompressAllChannels(I);
    case 'byChannel'
        % Results on submission system: 
        % 
        % For colorSpace 'native' (rgb, gray)
        % Run-time:             1.444166s
        % Mean squared error:   0.012878
        % Compression rate:     0.012878
        disp('Method used: CompressByChannel');
        I_comp = CompressByChannel(I);
    case 'byChannelSVD'
        disp('Method used: CompressByChannel');
        I_comp = CompressByChannelSVD(I);
    case 'wavelet'
        disp('Method used: wavlet');
        [I_comp.orig.h, I_comp.orig.w, dummy] = size(I)
        % Round to next lower power of two:
        h = pow2(floor(log2(I_comp.orig.h)))
        w = pow2(floor(log2(I_comp.orig.w)))
        I = im2uint8(imresize(I, [h, w]));
        [CR,BPP] = wcompress('c',I,'compressed.wtc','spiht','maxloop',12,'wname','bior4.4');
        I_comp.cr = CR;
        I_comp.bpp = BPP;
        disp(['Compression ratio: ', num2str(CR), '%'])
    otherwise
        assert(false, 'This method is not supported')
end
