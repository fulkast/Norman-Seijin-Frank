function mask = random_mask(size,missing)
% function mask = random_mask(size,missing)
% Returns a binary mask of requested size with missing per cent of pixels

nelm = prod(size);
sel = randperm(nelm);
mask = true(size);
mask(sel(1:round(missing*end))) = false;
