vome_unstag=0*vomecrty;
vome_unstag(:,:,1:$y.size-1,:)=0.5*(vomecrty(:,:,1:$y.size-1,:) + vomecrty(:,:,0:$y.size-2,:));
