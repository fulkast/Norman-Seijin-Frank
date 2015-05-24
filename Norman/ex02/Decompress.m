function I_rec = Decompress(I_comp, method)

% Default settings that will be used by the submission machine.
if (nargin<2)
    method = 'byChannelSVD';
end

switch method
    case 'allChannels'
        I_rec = DecompressAllChannels(I_comp);
    case 'byChannel'
        I_rec = DecompressByChannel(I_comp);
    case 'byChannelSVD'
        I_rec = DecompressByChannelSVD(I_comp);
    case 'wavelet'
        I_rec = wcompress('u','compressed.wtc');
        I_rec = imresize(I_rec, [I_comp.orig.h, I_comp.orig.w]);
        c = size(I_rec, 3);
        if c == 1
            I_rec = double(I_rec)/max(I_rec(:));
        else
            I_rec = im2double(I_rec);
        end
    otherwise
        assert(false, 'This method is not supported')
end