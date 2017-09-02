function result = interleave(m, dim)
if dim==2
    result = reshape([m(:) m(:)]',2*size(m,1), []);
else
    result=interleave(interleave(m,dim/2), 2); 
end
end



