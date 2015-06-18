function [upintrude,uphost,downintrude,downhost] = neighbourhood_layout_up_down(blocksize,intrusion);


    topsection = 1:blocksize;
    topsection = topsection-1;
    topsection = topsection * blocksize;
    topsection = repmat(topsection,(blocksize-intrusion),1);
    topsection = topsection + repmat([1:blocksize-intrusion]',1,blocksize);
    topbase = reshape(topsection,numel(topsection),1);
    topsection = 1:blocksize;
    topsection = topsection-1;
    topsection = topsection * blocksize;
    topsection = repmat(topsection,(intrusion),1);
    topsection = topsection + repmat([1:intrusion]',1,blocksize);
    topsection = reshape(topsection,numel(topsection),1);
    topintrusion = topsection;
%     topintrusion = reshape(repmat(blocksize*(0:intrusion-1)',1,blocksize-2*intrusion)+repmat(intrusion+1:blocksize-intrusion,intrusion,1),(blocksize-2*intrusion)*intrusion,1);
    bottomintrusion = reshape(blocksize*blocksize - (repmat(fliplr(1:blocksize)-1,intrusion,1)*blocksize + repmat(fliplr(0:intrusion-1)',1,blocksize)),blocksize*intrusion,1);
%     bottomintrusion = reshape(blocksize*blocksize - (repmat(fliplr(intrusion+1:blocksize-intrusion)-1,intrusion,1)*blocksize + repmat(fliplr(0:intrusion-1)',1,blocksize-2*intrusion)),(blocksize-2*intrusion)*intrusion,1);
    bottombase = blocksize*blocksize - reshape((repmat(fliplr(1:blocksize)-1,blocksize-intrusion,1)*blocksize + repmat(fliplr(0:blocksize-intrusion-1)',1,blocksize)),blocksize*(blocksize-intrusion),1);

    
    downhost = bottombase;
    downintrude = bottomintrusion;
    upintrude = topintrusion;
    uphost = topbase;
    
    