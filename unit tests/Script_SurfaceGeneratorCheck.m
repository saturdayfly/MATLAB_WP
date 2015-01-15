% compute z profile directly from ciro's data
% run after part 1bis of main program

% parameters 

    nq = 5;
    lx = 1;
    ly = 1;
    amp = 0.1/2/nq;   % check this before runing

% LOAD ah, dh_iq, dh_jq and fh from Ciro's data

    count = 1;

    
xvec = X(:);
yvec = Y(:);

clear lh;
for q=1:1:length(xvec)

        xx = xvec(q);
        yy = yvec(q);
        lh(count) = 0;
        
        for idx = 1:1:(2*nq+1)^2
                if dh_iq(idx)~= 0 && dh_jq(idx) ~= 0
                    lh(count) = lh(count) + amp*ah(idx)*sin(2*pi*((dh_iq(idx)*xx/lx)+(dh_jq(idx)*yy/ly)+fh(idx)));
                else
                    % fprintf('q vector is zero for idx = %d\n',idx);
                end
        end
        count = count+1;

end

hold on;
plot3(xvec,yvec,lh,'o'); axis equal;

dlmwrite('surface.dat',[xvec yvec lh'],'delimiter','\t');