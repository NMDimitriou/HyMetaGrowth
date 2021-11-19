function sim=cosine_sim(P,Q)
%  sim = cosine_sim(P,Q) cosine similarity measure of two discrete probability
%  distributions 
% P =  1 x nbins
% Q =  1 x nbins 
% sim = 1x1
if size(P,2)~=size(Q,2)
    error('the number of columns in P and Q should be the same');
end
if sum(~isfinite(P(:))) + sum(~isfinite(Q(:)))
   error('the inputs contain non-finite values!')
end

    dP = sqrt(sum(P.^2,2));
    dQ = sqrt(sum(Q.^2,2));
    
    dotPQ = sum(P.*Q,2);
    
    sim = dotPQ/(dP*dQ);
end